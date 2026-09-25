"""Descriptive comparison of frozen functional-analysis results.

No scientific module is executed here. Stable database IDs, run-local tables,
snapshots and classifications are retained separately for A and B.
"""
from __future__ import annotations

import json
import re
import uuid
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from .experiment_comparison import load_run as load_base


class FunctionalComparisonError(ValueError):
    pass


MODULES = ("go", "kegg", "reactome", "mitocarta", "mtdna", "domains", "complexes", "string")
UNITS = {
    "go": ("BP", "MF", "CC"), "kegg": ("Pathways",),
    "reactome": ("Pathways",), "mitocarta": ("Membership", "Subcompartments", "MitoPathways"),
    "mtdna": ("Categories", "Sources"),
    "domains": ("InterPro", "Pfam", "InterPro architectures", "Pfam architectures"),
    "complexes": ("Complexes",), "string": ("Nodes", "Edges"),
}
PARAMETER_KEYS = ("evidence_filter", "minimum_overlap", "min_count", "fdr_cutoff",
    "p_cutoff", "mapping_policy", "network_type", "combined_score_threshold",
    "threshold_preset", "max_hop", "degree1_selection_mode", "annotation_set_id")


@dataclass(frozen=True)
class FunctionalSpec:
    comparison_run_id: str
    module: str
    run_a: str
    run_b: str


def _safe(value: str) -> str:
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]{0,100}", str(value)):
        raise FunctionalComparisonError("Invalid run ID.")
    return str(value)


def list_source_runs(project, module: str) -> list[str]:
    if module not in MODULES:
        raise FunctionalComparisonError("Unknown functional module.")
    from . import (go_analysis, kegg_analysis, reactome_analysis, mitocarta_analysis,
        mtdna_analysis, interpro_pfam_analysis, complex_analysis, string_analysis)
    listing = {"go": go_analysis.list_go_runs, "kegg": kegg_analysis.list_kegg_runs,
        "reactome": reactome_analysis.list_reactome_runs,
        "mitocarta": mitocarta_analysis.list_mitocarta_runs,
        "mtdna": mtdna_analysis.list_mtdna_runs,
        "domains": interpro_pfam_analysis.list_interpro_pfam_runs,
        "complexes": complex_analysis.list_complex_runs,
        "string": string_analysis.list_string_runs}[module]
    return sorted(listing(project), reverse=True)


def _read(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FunctionalComparisonError(f"Incomplete run: {path}")
    try:
        return pd.read_csv(path, dtype=str, keep_default_na=False)
    except (OSError, ValueError, pd.errors.EmptyDataError) as error:
        raise FunctionalComparisonError(f"Incomplete run: {path}: {error}") from error


def _go(project, run_id: str):
    root = Path(project.root) / "analyses/GO/runs" / _safe(run_id)
    metadata = json.loads((root / "metadata.json").read_text(encoding="utf-8"))
    annotations = _read(root / "annotation/go_annotations.csv")
    tables = {}
    for ontology in ("BP", "MF", "CC"):
        base = ontology.lower()
        frequency = root / "frequency" / f"go_{base}_frequency.csv"
        if frequency.is_file():
            tables[ontology] = {"frequency": _read(frequency),
                "enrichment": _read(root / "enrichment" / f"go_{base}_all.csv"),
                "significant": _read(root / "enrichment" / f"go_{base}_significant.csv"),
                "membership": annotations[annotations.ontology.astype(str) == ontology].copy()}
    if not tables:
        raise FunctionalComparisonError("Incomplete GO run: no ontology results.")
    return metadata, tables


def _source(project, module: str, run_id: str):
    _safe(run_id)
    if module == "go":
        return _go(project, run_id)
    if module == "kegg":
        from .kegg_analysis import read_kegg_outputs
        out = read_kegg_outputs(project, run_id)
        return out.metadata, {"Pathways": {"frequency": out.frequency,
            "enrichment": out.enrichment,
            "significant": _read(out.root / "enrichment/kegg_pathway_significant.csv"),
            "membership": out.mapping}}
    if module == "reactome":
        from .reactome_analysis import read_reactome_outputs
        out = read_reactome_outputs(project, run_id)
        return out.metadata, {"Pathways": {"frequency": out.frequency,
            "enrichment": out.enrichment, "significant": out.significant,
            "membership": out.membership, "hierarchy": out.hierarchy}}
    if module == "mitocarta":
        from .mitocarta_analysis import read_mitocarta_outputs
        out = read_mitocarta_outputs(project, run_id)
        return out.metadata, {
            "Membership": {"frequency": out.membership[
                out.membership.is_mitocarta.astype(str).str.lower().isin({"true", "1", "t"})].copy()},
            "Subcompartments": {"frequency": out.subcompartment_frequency,
                "enrichment": out.subcompartment_enrichment,
                "significant": out.subcompartment_significant},
            "MitoPathways": {"frequency": out.pathway_frequency,
                "enrichment": out.pathway_enrichment,
                "significant": out.pathway_significant}}
    if module == "mtdna":
        from .mtdna_analysis import read_mtdna_outputs
        out = read_mtdna_outputs(project, run_id)
        t = out["tables"]
        return out["metadata"], {"Categories": {"frequency": t["category_frequency"],
            "enrichment": t["enrichment_all"], "significant": t["enrichment_significant"],
            "membership": t["category_to_entities"]},
            "Sources": {"frequency": t["source_frequency"],
                "membership": t["source_to_entities"]}}
    if module == "domains":
        from .interpro_pfam_analysis import read_interpro_pfam_outputs
        out = read_interpro_pfam_outputs(project, run_id)
        return out.metadata, {"InterPro": {"frequency": out.interpro_frequency,
            "enrichment": out.interpro_enrichment, "significant": out.interpro_significant},
            "Pfam": {"frequency": out.pfam_frequency,
                "enrichment": out.pfam_enrichment, "significant": out.pfam_significant},
            "InterPro architectures": {"frequency": out.interpro_architecture_frequency},
            "Pfam architectures": {"frequency": out.pfam_architecture_frequency}}
    if module == "complexes":
        from .complex_analysis import read_complex_outputs
        out = read_complex_outputs(project, run_id)
        t = out["tables"]
        return out["metadata"], {"Complexes": {"frequency": t["complex_frequency"],
            "enrichment": t["complex_enrichment_all"],
            "significant": t["complex_enrichment_significant"],
            "coverage": t["complex_coverage"],
            "membership": t["complex_to_proteins"],
            "alternative_groups": t["alternative_component_groups"]}}
    if module == "string":
        from .string_analysis import read_string_outputs
        out = read_string_outputs(project, run_id)
        t = out["tables"]
        return out["metadata"], {"Nodes": {"frequency": t["expanded_nodes"]},
            "Edges": {"frequency": t["expanded_edges"]}}
    raise FunctionalComparisonError("Unknown functional module.")


def _first(row, names: tuple[str, ...]):
    for name in names:
        if name in row and str(row[name]).strip() not in {"", "nan", "NA"}:
            return str(row[name]).strip()
    return ""


def _columns(module: str, unit: str):
    if module == "go":
        return (("GO_ID",), ("Description", "GO_term"), ("Protein_count", "Count"),
            ("Protein_fraction",), ("entity_id",))
    if module == "kegg":
        return (("pathway_id",), ("pathway_name", "name", "description"),
            ("gene_count", "target_gene_count"), ("gene_fraction",), ("genes", "gene_id"))
    if module == "reactome":
        return (("Reactome_ID",), ("Pathway",), ("Protein_count", "Target_count"),
            ("Protein_fraction",), ("reactome_entity_key", "UniProt"))
    if module == "mitocarta":
        if unit == "Membership":
            return ("canonical_gene_key",), ("gene_symbol",), (), (), ("canonical_gene_key",)
        if unit == "Subcompartments":
            return (("Subcompartment", "subcompartment"), ("Subcompartment", "subcompartment"),
                ("Gene_count", "Target_count"), ("Fraction_of_MitoCarta_target",), ("Genes",))
        return (("MitoPathway", "Pathway", "pathway"), ("MitoPathway", "Pathway", "pathway"),
            ("Gene_count", "Target_count"), ("Fraction_of_MitoCarta_target",), ("Genes",))
    if module == "mtdna":
        key = ("category", "Category", "evidence_category") if unit == "Categories" else (
            "source", "Source", "evidence_source")
        return (key, ("display_name",) + key, ("target_entity_count", "entity_count", "Target_count"), ("fraction_of_target",),
            ("target_entities", "entities", "entity_key"))
    if module == "domains":
        if "architectures" in unit:
            return (("Architecture",), ("Architecture",), ("Protein_count",),
                ("Fraction_of_target",), ("Proteins",))
        key = ("InterPro_ID", "ID") if unit == "InterPro" else ("Pfam_ID", "ID")
        return (key, ("Name", "Description"), ("Protein_count", "Target_count"),
            ("Fraction_of_target",), ("Proteins",))
    if module == "complexes":
        return (("complex_id",), ("complex_name",), ("target_member_count", "Target_count"),
            ("fraction_of_target",), ("target_member_proteins", "canonical_uniprot"))
    if module == "string":
        if unit == "Nodes":
            return ("string_protein_id",), ("preferred_name", "gene_symbol"), (), (), ()
        return ("protein_a", "protein_b"), (), (), (), ()
    raise FunctionalComparisonError("Unsupported comparison unit.")


def _record_index(module: str, unit: str, tables: dict):
    ids, names, counts, fractions, member_columns = _columns(module, unit)
    frequency = tables.get("frequency", pd.DataFrame())
    if frequency.empty:
        return {}
    if any(name not in frequency for name in ids) and not any(name in frequency for name in ids):
        raise FunctionalComparisonError(f"No stable {unit} identifier in persisted run.")
    enrichment = tables.get("enrichment", pd.DataFrame())
    significant = tables.get("significant", pd.DataFrame())
    membership = tables.get("membership", pd.DataFrame())
    coverage = tables.get("coverage", pd.DataFrame())
    hierarchy = tables.get("hierarchy", pd.DataFrame())
    def key(row):
        if module == "string" and unit == "Edges":
            pair = (_first(row, ("protein_a",)), _first(row, ("protein_b",)))
            return "|".join(sorted(pair)) if all(pair) else ""
        return _first(row, ids)
    enrichment_by = {key(row): row for _, row in enrichment.iterrows()} if not enrichment.empty else {}
    significant_ids = {key(row) for _, row in significant.iterrows()} if not significant.empty else set()
    records = {}
    for _, row in frequency.iterrows():
        entity = key(row)
        if not entity:
            continue
        raw_count = _first(row, counts)
        if counts and raw_count:
            try:
                if float(raw_count) <= 0:
                    continue
            except ValueError:
                pass
        enrichment_row = enrichment_by.get(entity, {})
        members = set()
        for field in member_columns:
            if field in row:
                members.update(part for part in re.split(r"[;,|]", str(row[field])) if part.strip())
        if not membership.empty:
            id_field = next((field for field in ids if field in membership), None)
            if id_field:
                matching = membership[membership[id_field].astype(str) == entity]
                for field in member_columns:
                    if field in matching:
                        members.update(value for cell in matching[field] for value in
                            re.split(r"[;,|]", str(cell)) if value.strip())
        members = {value.strip() for value in members if value.strip() and value.strip() != "nan"}
        extra = {}
        if module == "complexes" and not coverage.empty and "complex_id" in coverage:
            hit = coverage[coverage.complex_id.astype(str) == entity]
            if not hit.empty:
                extra = {field: hit.iloc[0][field] for field in (
                    "coverage_class", "protein_component_coverage", "covered_component_groups",
                    "total_protein_component_groups") if field in hit}
        if module == "string":
            extra = {field: row[field] for field in ("combined_score", "hop_level", "edge_type", "node_type") if field in row}
        if module == "reactome" and not hierarchy.empty and "Reactome_ID" in hierarchy:
            hit = hierarchy[hierarchy.Reactome_ID.astype(str) == entity]
            if not hit.empty:
                extra["hierarchy"] = hit.to_json(orient="records")
        records[entity] = {"name": _first(row, names), "count": raw_count,
            "fraction": _first(row, fractions), "members": members,
            "p_value": _first(enrichment_row, ("p_value", "pvalue", "P.Value")),
            "FDR": _first(enrichment_row, ("FDR", "p_adjust", "p.adjust", "adj.P.Val")),
            "overlap": _first(enrichment_row, ("Target_count", "Count", "target_gene_count")),
            "tested": entity in enrichment_by, "enriched": entity in significant_ids if entity in enrichment_by else None,
            "extra": extra, "raw": row.to_dict()}
    return records


def _snapshot(metadata):
    return _first(metadata, ("snapshot_id", "database_snapshot_id", "annotation_set_id",
        "go_snapshot_id", "reactome_snapshot_id", "string_snapshot_id"))


def _lineage(metadata):
    return _first(metadata, ("frozen_input_hash", "input_hash", "input_sha256"))


def _compare_unit(module, unit, left, right, snap_a, snap_b):
    a, b = _record_index(module, unit, left), _record_index(module, unit, right)
    rows, members = [], []
    for entity in sorted(a.keys() | b.keys()):
        x, y = a.get(entity), b.get(entity)
        present_class = "Shared" if x and y else "A only" if x else "B only"
        ex, ey = x.get("enriched") if x else None, y.get("enriched") if y else None
        if x and y and (x["tested"] or y["tested"]):
            if not x["tested"]:
                enrichment_class = "Not tested in A"
            elif not y["tested"]:
                enrichment_class = "Not tested in B"
            elif ex and ey:
                enrichment_class = "Enriched in both"
            elif ex:
                enrichment_class = "Enriched only in A"
            elif ey:
                enrichment_class = "Enriched only in B"
            else:
                enrichment_class = "Not enriched in either"
        else:
            enrichment_class = "Not applicable"
        row = {"entity_id": entity, "name_A": x["name"] if x else "",
            "name_B": y["name"] if y else "", "present_A": bool(x), "present_B": bool(y),
            "A_count": x["count"] if x else "", "B_count": y["count"] if y else "",
            "A_fraction": x["fraction"] if x else "", "B_fraction": y["fraction"] if y else "",
            "A_overlap": x["overlap"] if x else "", "B_overlap": y["overlap"] if y else "",
            "A_p_value": x["p_value"] if x else "", "B_p_value": y["p_value"] if y else "",
            "A_FDR": x["FDR"] if x else "", "B_FDR": y["FDR"] if y else "",
            "A_enriched": ex, "B_enriched": ey, "presence_class": present_class,
            "comparison_class": enrichment_class if enrichment_class != "Not applicable" else present_class,
            "snapshot_A": snap_a, "snapshot_B": snap_b,
            "A_score": x["extra"].get("combined_score", "") if x else "",
            "B_score": y["extra"].get("combined_score", "") if y else "",
            "A_details": json.dumps(x["extra"], default=str) if x else "",
            "B_details": json.dumps(y["extra"], default=str) if y else ""}
        rows.append(row)
        left_members, right_members = x["members"] if x else set(), y["members"] if y else set()
        for value in sorted(left_members | right_members):
            members.append({"entity_id": entity, "member_id": value,
                "A_member": value in left_members, "B_member": value in right_members,
                "member_class": "Shared" if value in left_members & right_members else
                    "A only" if value in left_members else "B only"})
    master_columns = ("entity_id", "name_A", "name_B", "present_A", "present_B",
        "A_count", "B_count", "A_fraction", "B_fraction", "A_overlap", "B_overlap",
        "A_p_value", "B_p_value", "A_FDR", "B_FDR", "A_enriched", "B_enriched",
        "presence_class", "comparison_class", "snapshot_A", "snapshot_B",
        "A_score", "B_score", "A_details", "B_details")
    return pd.DataFrame(rows, columns=master_columns), pd.DataFrame(members, columns=(
        "entity_id", "member_id", "A_member", "B_member", "member_class"))


def _bar(path: Path, counts: dict):
    from PySide6.QtGui import QImage, QPainter, QColor
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])
    image = QImage(700, 290, QImage.Format_RGB32)
    image.fill(0xFFFFFF)
    painter = QPainter(image)
    maximum = max(counts.values(), default=1) or 1
    for index, (label, count) in enumerate(counts.items()):
        y = 45 + index * 75
        painter.drawText(20, y + 20, label)
        painter.fillRect(170, y, int(420 * count / maximum), 30, QColor("#4682b4"))
        painter.drawText(610, y + 20, str(count))
    painter.end()
    image.save(str(path), "PNG")


def compare(project, spec: FunctionalSpec) -> Path:
    if spec.module not in MODULES:
        raise FunctionalComparisonError("Unsupported functional module.")
    _safe(spec.comparison_run_id)
    _safe(spec.run_a)
    _safe(spec.run_b)
    base = load_base(project, spec.comparison_run_id)
    meta_a, tables_a = _source(project, spec.module, spec.run_a)
    meta_b, tables_b = _source(project, spec.module, spec.run_b)
    if spec.module == "string" and _first(meta_a, ("network_type",)) != _first(meta_b, ("network_type",)):
        raise FunctionalComparisonError("Not directly comparable: STRING network types differ.")
    snap_a, snap_b = _snapshot(meta_a), _snapshot(meta_b)
    warnings = []
    if snap_a and snap_b and snap_a != snap_b:
        warnings.append("These runs were generated using different database snapshots.")
    if spec.module == "string" and _first(meta_a, ("combined_score_threshold",)) != _first(meta_b, ("combined_score_threshold",)):
        warnings.append("STRING thresholds differ; node and edge counts are descriptive.")
    parameters = {side: {key: meta.get(key) for key in PARAMETER_KEYS if key in meta}
        for side, meta in (("A", meta_a), ("B", meta_b))}
    for key in PARAMETER_KEYS:
        if key in parameters["A"] and key in parameters["B"] and parameters["A"][key] != parameters["B"][key]:
            warnings.append(f"Parameter differs: {key}.")
    lineages = {"A": _lineage(meta_a), "B": _lineage(meta_b)}
    base_hashes = base["manifest"].get("input_hashes", {})
    for side, hash_key in (("A", "a"), ("B", "b")):
        if lineages[side] and lineages[side] != base_hashes.get(hash_key):
            warnings.append(f"Run {side} input lineage differs from or cannot be verified against dataset {side}.")
        elif not lineages[side]:
            warnings.append(f"Run {side} input lineage is not recorded; dataset correspondence is unverified.")
    now = datetime.now(timezone.utc)
    identifier = f"{now.strftime('%Y%m%dT%H%M%SZ')}-{uuid.uuid4().hex[:8]}"
    root = base["run_root"] / "functional" / spec.module / identifier
    root.mkdir(parents=True, exist_ok=False)
    counts = {}
    for unit in UNITS[spec.module]:
        if unit not in tables_a or unit not in tables_b:
            continue
        master, members = _compare_unit(spec.module, unit, tables_a[unit], tables_b[unit], snap_a, snap_b)
        folder = root / re.sub(r"[^a-z0-9]+", "_", unit.lower()).strip("_")
        folder.mkdir()
        enriched = master[master.comparison_class.isin(("Enriched only in A", "Enriched only in B", "Enriched in both"))] if not master.empty else master
        for name, frame in (("enrichment", enriched), ("master", master), ("a_only", master[master.presence_class == "A only"] if not master.empty else master),
                ("b_only", master[master.presence_class == "B only"] if not master.empty else master),
                ("shared", master[master.presence_class == "Shared"] if not master.empty else master),
                ("members", members)):
            frame.to_csv(folder / f"{name}.csv", index=False)
        counts[unit] = {label: int((master.presence_class == label).sum()) for label in
            ("A only", "Shared", "B only")}
        _bar(folder / "counts.png", counts[unit])
        enrichment_counts = {label: int((master.comparison_class == label).sum()) for label in ("Enriched only in A", "Enriched in both", "Enriched only in B")}
        if any(enrichment_counts.values()):
            _bar(folder / "enrichment_counts.png", enrichment_counts)
    if not counts:
        raise FunctionalComparisonError("Incomplete run: no comparable persisted units.")
    config = asdict(spec) | {"snapshot_A": snap_a, "snapshot_B": snap_b,
        "parameters": parameters, "input_lineages": lineages,
        "source_metadata_A": meta_a, "source_metadata_B": meta_b,
        "comparison_identity": base["config"].get("comparison_type"), "warnings": warnings}
    (root / "functional_config.json").write_text(json.dumps(config, indent=2, default=str) + "\n", encoding="utf-8")
    with pd.ExcelWriter(root / "functional_comparison.xlsx", engine="openpyxl") as writer:
        overview = pd.DataFrame([{"unit": unit, **values} for unit, values in counts.items()])
        overview.to_excel(writer, sheet_name="Overview", index=False)
        for unit in counts:
            folder = root / re.sub(r"[^a-z0-9]+", "_", unit.lower()).strip("_")
            for label in ("master", "a_only", "b_only", "shared", "members", "enrichment"):
                pd.read_csv(folder / f"{label}.csv").to_excel(writer,
                    sheet_name=f"{unit[:18]} {label[:10]}"[:31], index=False)
    (root / "functional_manifest.json").write_text(json.dumps({"status": "Ready",
        "created_at": now.isoformat(), "units": list(counts), "counts": counts,
        "source_comparison_run": spec.comparison_run_id}, indent=2) + "\n", encoding="utf-8")
    return root


def list_functional_runs(project, comparison_run_id: str, module: str) -> list[str]:
    base = load_base(project, _safe(comparison_run_id))["run_root"]
    folder = base / "functional" / module
    return sorted((path.name for path in folder.iterdir() if path.is_dir()) if folder.is_dir()
        else [], reverse=True)


def load_functional(project, comparison_run_id: str, module: str, functional_run_id: str):
    base = load_base(project, _safe(comparison_run_id))["run_root"]
    root = base / "functional" / module / _safe(functional_run_id)
    config = json.loads((root / "functional_config.json").read_text(encoding="utf-8"))
    manifest = json.loads((root / "functional_manifest.json").read_text(encoding="utf-8"))
    if manifest.get("status") != "Ready":
        raise FunctionalComparisonError("Incomplete functional comparison.")
    tables = {}
    for unit in manifest["units"]:
        folder = root / re.sub(r"[^a-z0-9]+", "_", unit.lower()).strip("_")
        tables[unit] = {name: _read(folder / f"{name}.csv") for name in
            ("master", "a_only", "b_only", "shared", "members", "enrichment")}
    return {"root": root, "config": config, "manifest": manifest, "tables": tables}
