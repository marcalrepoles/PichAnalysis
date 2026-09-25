"""Offline comparison of frozen, independent experimental datasets.

Presence sets are descriptive. Differential labels use only previously
persisted classifications and never recalculate statistics.
"""
from __future__ import annotations

import hashlib
import json
import re
import shutil
import uuid
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from .column_mapping import suggest_columns
from .entity_identity import tokens
from .identifier_detection import IdentifierType
from .importer import read_table


class ComparisonError(ValueError):
    pass


TYPE_ALIASES = {"entrez": "ncbi_gene", "refseq_protein": "refseq"}
IDENTITY_COLUMNS = {"uniprot": "uniprot_accession", "gene_symbol": "gene_symbol",
    "ncbi_gene": "ncbi_gene_id", "ensembl_gene": "ensembl_gene_id",
    "ensembl_protein": "ensembl_protein_id", "refseq": "refseq_id"}
SUPPORTED_TYPES = {item.value for item in IdentifierType} - {"unknown"} | {"original"}
MASTER_COLUMNS = ("comparison_entity", "comparison_entity_type", "A_present", "B_present",
    "A_original_identifier", "B_original_identifier", "A_source_rows", "B_source_rows",
    "A_mapping_status", "B_mapping_status", "match_basis", "mapping_provenance")
DIFF_COLUMNS = ("comparison_entity", "A_direction", "A_effect", "A_log2FC", "A_P.Value",
    "A_adj.P.Val", "A_combined_significant", "B_direction", "B_effect", "B_log2FC",
    "B_P.Value", "B_adj.P.Val", "B_combined_significant", "direction_pattern",
    "significance_pattern", "comparison_class")


@dataclass(frozen=True)
class TableSpec:
    path: str
    identifier_column: str
    identifier_type: str
    sheet: str | None = None
    alias: str = ""
    source_kind: str = "external"


@dataclass(frozen=True)
class ComparisonSpec:
    a: TableSpec
    b: TableSpec
    comparison_type: str
    mapping_path: str | None = None
    differential_a: str | None = None
    differential_b: str | None = None
    comparable_contrasts: bool = False
    differential_match_a: bool = False
    differential_match_b: bool = False


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def suggestions(frame: pd.DataFrame) -> dict:
    """Suggestions only: callers must request explicit user confirmation."""
    return suggest_columns(frame)


def canonical_type(value: str) -> str:
    return TYPE_ALIASES.get(value, value)


def _normal(value: str, kind: str) -> str:
    value = str(value).strip()
    return value.upper() if kind != "original" else value


def _blank(value: object) -> bool:
    return pd.isna(value) or str(value).strip().casefold() in {"", "na", "nan", "none", "null"}


def _read_source(path: Path, sheet: str | None):
    try:
        return read_table(path, sheet)
    except Exception as error:
        # The shared importer deliberately rejects delimiter-free CSV. A valid
        # one-column experimental table is still sufficient for presence sets.
        if path.suffix.lower() not in {".csv", ".tsv"} or "delimiter" not in str(error).lower():
            raise
        frame = pd.read_csv(path, sep="\t" if path.suffix.lower() == ".tsv" else ",",
            encoding="utf-8-sig")
        if len(frame.columns) != 1 or frame.empty or frame.columns[0].startswith("Unnamed:"):
            raise
        return frame, None

def _mapping_index(frame: pd.DataFrame, source_type: str, target_type: str) -> dict[str, set[str]]:
    source, target = IDENTITY_COLUMNS.get(source_type), IDENTITY_COLUMNS.get(target_type)
    if not source or not target or source not in frame or target not in frame:
        return {}
    index: dict[str, set[str]] = {}
    for _, row in frame.iterrows():
        for left in tokens(row[source]):
            for right in tokens(row[target]):
                index.setdefault(_normal(left, source_type), set()).add(_normal(right, target_type))
    return index


def _resolve(frame: pd.DataFrame, spec: TableSpec, target_type: str,
             mapping: pd.DataFrame | None, provenance: str):
    if spec.identifier_column not in frame:
        raise ComparisonError(f"Identifier column not found: {spec.identifier_column}")
    source_type = canonical_type(spec.identifier_type)
    lookup = _mapping_index(mapping, source_type, target_type) if mapping is not None else {}
    entities: dict[str, dict] = {}
    ambiguous, unmapped = [], []
    for row_number, raw in enumerate(frame[spec.identifier_column].tolist(), 1):
        if _blank(raw):
            unmapped.append({"source_row": row_number, "original_identifier": "", "status": "Unmapped"})
            continue
        original = str(raw).strip()
        # A protein group remains one experimental row. Members contribute only
        # to qualitative set membership; intensities are never duplicated.
        members = tokens(original) if source_type == "uniprot" and target_type != "original" else (original,)
        for member in members:
            if source_type == target_type or target_type == "original":
                targets, status = {_normal(member, source_type)}, "Exact"
            elif mapping is None:
                unmapped.append({"source_row": row_number, "original_identifier": original,
                    "status": "Mapping required"})
                continue
            else:
                targets, status = lookup.get(_normal(member, source_type), set()), "Mapped"
            if len(targets) > 1:
                ambiguous.append({"source_row": row_number, "original_identifier": original,
                    "source_member": member, "candidate_targets": ";".join(sorted(targets)),
                    "status": "Ambiguous", "mapping_provenance": provenance})
                continue
            if not targets:
                unmapped.append({"source_row": row_number, "original_identifier": original,
                    "status": "Unmapped", "mapping_provenance": provenance})
                continue
            entity = next(iter(targets))
            item = entities.setdefault(entity, {"rows": set(), "original": set(), "status": set(),
                "basis": set(), "provenance": set()})
            item["rows"].add(row_number)
            item["original"].add(original)
            item["status"].add(status)
            item["basis"].add("group membership" if len(members) > 1 else
                "exact identifier" if status == "Exact" else "persisted mapping")
            if status == "Mapped":
                item["provenance"].add(provenance)
    original_ids = [str(value).strip() for value in frame[spec.identifier_column] if not _blank(value)]
    duplicate_rows = len(original_ids) - len(set(original_ids))
    return entities, ambiguous, unmapped, duplicate_rows


def _joined(item: dict | None, field: str) -> str:
    return ";".join(map(str, sorted(item[field], key=str))) if item else ""


def _master(a: dict, b: dict, kind: str) -> pd.DataFrame:
    rows = []
    for entity in sorted(a.keys() | b.keys()):
        left, right = a.get(entity), b.get(entity)
        rows.append({"comparison_entity": entity, "comparison_entity_type": kind,
            "A_present": bool(left), "B_present": bool(right),
            "A_original_identifier": _joined(left, "original"),
            "B_original_identifier": _joined(right, "original"),
            "A_source_rows": _joined(left, "rows"), "B_source_rows": _joined(right, "rows"),
            "A_mapping_status": _joined(left, "status") if left else "Not present",
            "B_mapping_status": _joined(right, "status") if right else "Not present",
            "match_basis": ";".join(filter(None, (_joined(left, "basis"), _joined(right, "basis")))),
            "mapping_provenance": ";".join(filter(None, (_joined(left, "provenance"),
                _joined(right, "provenance"))))})
    return pd.DataFrame(rows, columns=MASTER_COLUMNS)


def _differential_index(frame: pd.DataFrame, entities: set[str], kind: str):
    eligible = ("original_identifier", "original_id", "feature_id") if kind == "original" else (
        IDENTITY_COLUMNS.get(kind, ""),)
    columns = [name for name in eligible if name in frame]
    index, collisions = {}, set()
    for _, row in frame.iterrows():
        matches = {part.upper() for name in columns for part in tokens(row[name])
            if part.upper() in entities}
        for entity in matches:
            if entity in index:
                collisions.add(entity)
            else:
                index[entity] = row.to_dict()
    for entity in collisions:
        index.pop(entity, None)
    return index, collisions


def _sig(row: dict | None) -> bool | None:
    if row is None:
        return None
    value = str(row.get("combined_significant", "")).casefold()
    return True if value in {"true", "1", "yes"} else False if value in {"false", "0", "no"} else None


def _direction(row: dict | None) -> str:
    if row is None:
        return "Not tested"
    raw = str(row.get("effect_direction", "")).casefold()
    return "Up" if raw in {"up", "positive", "increase", "higher_in_condition_a"} else "Down" if raw in {
        "down", "negative", "decrease", "higher_in_condition_b"} else "Unknown"


def compare_differential(frame_a: pd.DataFrame, frame_b: pd.DataFrame,
                         entities: set[str], comparable: bool, identity_type: str = "uniprot") -> tuple[pd.DataFrame, dict]:
    a, ambiguous_a = _differential_index(frame_a, entities, identity_type)
    b, ambiguous_b = _differential_index(frame_b, entities, identity_type)
    rows = []
    for entity in sorted(entities):
        left, right = a.get(entity), b.get(entity)
        da, db, sa, sb = _direction(left), _direction(right), _sig(left), _sig(right)
        if not comparable or entity in ambiguous_a | ambiguous_b:
            pattern, cls = "Not comparable / not tested", "Not comparable"
        elif sa and sb and da in {"Up", "Down"} and db in {"Up", "Down"}:
            pattern = f"{da} / {db}"
            cls = "Significant opposite direction" if da != db else f"Significant {da.lower()} in both"
        elif sa:
            pattern, cls = "Changed only in A", "Significant only in A"
        elif sb:
            pattern, cls = "Changed only in B", "Significant only in B"
        elif left is None or right is None:
            pattern, cls = "Not comparable / not tested", "Not comparable"
        else:
            pattern, cls = "No significant change in either", "Not significant in either"
        result = {"comparison_entity": entity, "A_direction": da, "B_direction": db,
            "direction_pattern": pattern, "significance_pattern":
            f"{sa if sa is not None else 'NA'} / {sb if sb is not None else 'NA'}",
            "comparison_class": cls}
        for prefix, row in (("A", left), ("B", right)):
            for output, names in (("effect", ("effect",)), ("log2FC", ("log2FC",)),
                    ("P.Value", ("P.Value", "p_value")), ("adj.P.Val", ("adj.P.Val", "fdr")),
                    ("combined_significant", ("combined_significant",))):
                result[f"{prefix}_{output}"] = next((row[name] for name in names if row is not None
                    and name in row and str(row[name]).strip()), pd.NA)
        rows.append(result)
    return pd.DataFrame(rows, columns=DIFF_COLUMNS), {"ambiguous_A": sorted(ambiguous_a),
        "ambiguous_B": sorted(ambiguous_b)}


def _workbook(tables: dict[str, pd.DataFrame], summary: dict, path: Path) -> None:
    labels = {"master_comparison": "Master", "a_only": "A only", "b_only": "B only",
        "shared": "Shared", "ambiguous": "Ambiguous", "unmapped_a": "Unmapped A",
        "unmapped_b": "Unmapped B", "differential_comparison": "Differential",
        "significant_up_both": "Up both", "significant_down_both": "Down both",
        "opposite_direction": "Opposite", "significant_only_a": "Only A changed",
        "significant_only_b": "Only B changed",
        "qualitative_a": "Qualitative A", "qualitative_b": "Qualitative B"}
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        pd.DataFrame(summary.items(), columns=["metric", "value"]).to_excel(writer,
            sheet_name="Summary", index=False)
        for name, frame in tables.items():
            frame.to_excel(writer, sheet_name=labels[name], index=False)


def _bar_graph(title: str, values: dict[str, int], path: Path) -> None:
    from PySide6.QtCore import Qt
    from PySide6.QtGui import QColor, QFont, QImage, QPainter
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])
    image = QImage(900, max(360, 105 + len(values) * 60), QImage.Format_RGB32)
    image.fill(Qt.white)
    painter = QPainter(image)
    painter.setFont(QFont("Arial", 15))
    painter.drawText(25, 35, title)
    maximum = max(values.values(), default=1) or 1
    for position, (label, count) in enumerate(values.items()):
        y = 75 + position * 60
        painter.setBrush(QColor("#4682b4"))
        painter.drawRect(240, y, int(530 * count / maximum), 35)
        painter.drawText(20, y + 25, label)
        painter.drawText(795, y + 25, str(count))
    painter.end()
    image.save(str(path), "PNG")


def compare(project, spec: ComparisonSpec, *, run_id: str | None = None) -> Path:
    kind = canonical_type(spec.comparison_type)
    if kind not in set(IDENTITY_COLUMNS) | {"original"}:
        raise ComparisonError("Select a supported comparison identity.")
    if spec.a.identifier_type not in SUPPORTED_TYPES or spec.b.identifier_type not in SUPPORTED_TYPES:
        raise ComparisonError("Confirm a supported identifier type for each dataset.")
    if kind == "original" and canonical_type(spec.a.identifier_type) != canonical_type(spec.b.identifier_type):
        raise ComparisonError("Original identifiers require the same type on both sides.")
    if bool(spec.differential_a) != bool(spec.differential_b):
        raise ComparisonError("Select two Differential runs or neither.")
    if spec.differential_a and spec.differential_b:
        for side, differential_run_id, confirmed in ((spec.a, spec.differential_a, spec.differential_match_a),
                                        (spec.b, spec.differential_b, spec.differential_match_b)):
            if side.source_kind.startswith("differential:") and side.source_kind != f"differential:{differential_run_id}":
                raise ComparisonError("Selected Differential run does not match the frozen table source.")
            if not confirmed and side.source_kind != f"differential:{differential_run_id}":
                raise ComparisonError("Confirm that each Differential run corresponds to its dataset.")
    frame_a, sheet_a = _read_source(Path(spec.a.path), spec.a.sheet)
    frame_b, sheet_b = _read_source(Path(spec.b.path), spec.b.sheet)
    mapping = pd.read_csv(spec.mapping_path, dtype=str, keep_default_na=False) if spec.mapping_path else None
    provenance = str(spec.mapping_path or "")
    a, ambiguous_a, unmapped_a, duplicates_a = _resolve(frame_a, spec.a, kind, mapping, provenance)
    b, ambiguous_b, unmapped_b, duplicates_b = _resolve(frame_b, spec.b, kind, mapping, provenance)
    master = _master(a, b, kind)
    tables = {"master_comparison": master,
        "a_only": master.loc[master.A_present & ~master.B_present].copy(),
        "b_only": master.loc[~master.A_present & master.B_present].copy(),
        "shared": master.loc[master.A_present & master.B_present].copy(),
        "ambiguous": pd.DataFrame([{"side": side, **row} for side, records in
            (("A", ambiguous_a), ("B", ambiguous_b)) for row in records]),
        "unmapped_a": pd.DataFrame(unmapped_a), "unmapped_b": pd.DataFrame(unmapped_b)}
    differential_meta = None
    if spec.differential_a and spec.differential_b:
        from .differential_analysis import load_run
        output_a, output_b = load_run(project, spec.differential_a), load_run(project, spec.differential_b)
        diff, conflicts = compare_differential(output_a["tables"]["all_results"],
            output_b["tables"]["all_results"], set(master.comparison_entity),
            spec.comparable_contrasts, kind)
        tables["differential_comparison"] = diff
        tables["qualitative_a"] = output_a["tables"]["qualitative_candidates"].copy()
        tables["qualitative_b"] = output_b["tables"]["qualitative_candidates"].copy()
        differential_meta = {"A": output_a["metadata"], "B": output_b["metadata"],
            "A_qualitative_rows": len(output_a["tables"]["qualitative_candidates"]),
            "B_qualitative_rows": len(output_b["tables"]["qualitative_candidates"]), **conflicts}
        if spec.comparable_contrasts:
            for name, cls in (("significant_up_both", "Significant up in both"),
                    ("significant_down_both", "Significant down in both"),
                    ("opposite_direction", "Significant opposite direction"),
                    ("significant_only_a", "Significant only in A"),
                    ("significant_only_b", "Significant only in B")):
                tables[name] = diff[diff.comparison_class == cls].copy()
    summary = {"Table A original rows": len(frame_a), "Table B original rows": len(frame_b),
        "A unique comparison entities": len(a), "B unique comparison entities": len(b),
        "A duplicate identifier rows": duplicates_a, "B duplicate identifier rows": duplicates_b,
        "A only": len(tables["a_only"]), "B only": len(tables["b_only"]),
        "Shared": len(tables["shared"]), "Ambiguous": len(tables["ambiguous"]),
        "Unmapped A": len(unmapped_a), "Unmapped B": len(unmapped_b)}
    now = datetime.now(timezone.utc)
    identifier = run_id or f"{now.strftime('%Y%m%dT%H%M%SZ')}-{uuid.uuid4().hex[:8]}"
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]{0,100}", identifier):
        raise ComparisonError("Invalid comparison run ID.")
    root = Path(project.root) / "analyses/experiment_comparison/runs" / identifier
    if root.exists():
        raise ComparisonError("Comparison run ID already exists.")
    for folder in ("inputs/dataset_a", "inputs/dataset_b", "tables", "graphs"):
        (root / folder).mkdir(parents=True)
    try:
        for side, item in (("a", spec.a), ("b", spec.b)):
            shutil.copy2(item.path, root / f"inputs/dataset_{side}" / Path(item.path).name)
        if spec.mapping_path:
            (root / "inputs/mapping").mkdir()
            shutil.copy2(spec.mapping_path, root / "inputs/mapping" / Path(spec.mapping_path).name)
        for name, frame in tables.items():
            frame.to_csv(root / "tables" / f"{name}.csv", index=False)
        config = asdict(spec)
        for side, sheet in (("a", sheet_a), ("b", sheet_b)):
            config[side]["sheet"] = sheet
            config[side]["alias"] = config[side]["alias"] or f"Dataset {side.upper()}"
        (root / "comparison_config.json").write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
        manifest = {"run_id": identifier, "created_at": now.isoformat(), "status": "Ready",
            "input_hashes": {side: sha256(root / f"inputs/dataset_{side}" / Path(item.path).name)
                for side, item in (("a", spec.a), ("b", spec.b))},
            "input_lineages": {"A": sha256(Path(spec.a.path)), "B": sha256(Path(spec.b.path))},
            "summary": summary, "differential": differential_meta,
            "mapping_sha256": sha256(Path(spec.mapping_path)) if spec.mapping_path else None}
        (root / "comparison_manifest.json").write_text(json.dumps(manifest, indent=2,
            default=str) + "\n", encoding="utf-8")
        _workbook(tables, summary, root / "experiment_comparison.xlsx")
        _bar_graph(f"{config['a']['alias']} vs {config['b']['alias']} - membership counts", {name: len(tables[key]) for name, key in
            (("A only", "a_only"), ("Shared", "shared"), ("B only", "b_only"))},
            root / "graphs/membership.png")
        if "differential_comparison" in tables:
            _bar_graph(f"{config['a']['alias']} vs {config['b']['alias']} - Differential classes", tables["differential_comparison"]
                .comparison_class.value_counts().to_dict(), root / "graphs/differential_patterns.png")
    except Exception:
        expected_parent = (Path(project.root) / "analyses/experiment_comparison/runs").resolve()
        if root.resolve().parent == expected_parent and root.is_dir():
            shutil.rmtree(root)
        raise
    return root


def list_runs(project) -> list[str]:
    folder = Path(project.root) / "analyses/experiment_comparison/runs"
    return sorted((path.name for path in folder.iterdir() if path.is_dir()) if folder.is_dir()
        else [], reverse=True)


def load_run(project, run_id: str) -> dict:
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]{0,100}", str(run_id)):
        raise ComparisonError("Invalid comparison run ID.")
    root = Path(project.root) / "analyses/experiment_comparison/runs" / run_id
    config = json.loads((root / "comparison_config.json").read_text(encoding="utf-8"))
    manifest = json.loads((root / "comparison_manifest.json").read_text(encoding="utf-8"))
    if manifest.get("status") != "Ready":
        raise ComparisonError("Comparison run is not Ready.")
    if config.get("mapping_path"):
        mapping = root / "inputs/mapping" / Path(config["mapping_path"]).name
        if sha256(mapping) != manifest.get("mapping_sha256"):
            raise ComparisonError("Frozen Mapping catalog integrity mismatch.")
    for side in ("a", "b"):
        frozen = root / f"inputs/dataset_{side}" / Path(config[side]["path"]).name
        if sha256(frozen) != manifest["input_hashes"][side]:
            raise ComparisonError(f"Frozen dataset {side.upper()} integrity mismatch.")
    tables = {}
    for path in (root / "tables").glob("*.csv"):
        try:
            tables[path.stem] = pd.read_csv(path, dtype=str, keep_default_na=False)
        except pd.errors.EmptyDataError:
            tables[path.stem] = pd.DataFrame()
    return {"run_root": root, "config": config, "manifest": manifest, "tables": tables}
