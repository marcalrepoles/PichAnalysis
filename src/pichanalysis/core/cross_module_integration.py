"""Rebuildable offline index of identities in recognized, frozen analysis runs.

The index contains navigation keys and provenance, never scientific tables. A
missing or damaged index can always be regenerated from run-owned artifacts.
"""
from __future__ import annotations

import hashlib
import json
import os
import sqlite3
from decimal import Decimal, InvalidOperation
from contextlib import closing
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

import pandas as pd

from .entity_identity import (CrossModuleMatch, EntityContext, IdentityType,
    InputLineage, MatchStatus, tokens)


SCHEMA_VERSION = 1
IDENTITY_COLUMNS = {
    IdentityType.UNIPROT: ("uniprot_accession", "uniprot_accessions", "original_uniprot", "canonical_uniprot", "UniProt"),
    IdentityType.GENE_SYMBOL: ("gene_symbol", "gene_symbols", "Gene_symbol"),
    IdentityType.NCBI_GENE: ("ncbi_gene_id", "ncbi_gene_ids", "NCBI_Gene_ID"),
    IdentityType.ENSEMBL_GENE: ("ensembl_gene_id",),
    IdentityType.ENSEMBL_PROTEIN: ("ensembl_protein_id",),
    IdentityType.REFSEQ: ("refseq_id",),
}
ORIGINAL_COLUMNS = ("original_identifier", "original_id", "display_identifier", "input_id")
FEATURE_COLUMNS = ("feature_id",)
ROW_COLUMNS = ("source_row",)
MODULE_ID_COLUMNS = ("GO_ID", "go_id", "kegg_gene_id", "pathway_id", "Reactome_ID",
    "reactome_entity_key", "string_protein_id", "complex_id", "entity_key",
    "interpro_id", "pfam_id", "MitoPathway")
SNAPSHOT_FIELDS = ("snapshot_id", "annotation_set_id", "database_snapshot_id")


@dataclass(frozen=True)
class CrossModuleAdapter:
    module_id: str
    list_runs: Callable
    load_run: Callable
    record_tables: tuple[str, ...]
    row_level: bool = False

    def extract_identity_records(self, outputs):
        tables = outputs.get("tables", {}) if isinstance(outputs, dict) else {}
        for key in self.record_tables:
            frame = tables.get(key) if isinstance(outputs, dict) else getattr(outputs, key, None)
            if isinstance(frame, pd.DataFrame):
                yield key, frame

    def resolve_matches(self, index, context, *, include_other_lineages=False):
        return index.matches(context, target_module=self.module_id,
            include_other_lineages=include_other_lineages)

    def load_detail(self, project, run_id, record_type, record_index):
        outputs = self.load_run(project, run_id)
        for name, frame in self.extract_identity_records(outputs):
            if name == record_type and 0 <= record_index < len(frame):
                return {key: str(value) for key, value in frame.iloc[record_index].items()}
        return {}

    def open_target(self, page, run_id, target_identity):
        history = getattr(page, "history", None) or getattr(page, "stat_history", None)
        if history is not None:
            index = history.findData(run_id)
            if index < 0:
                index = history.findText(run_id)
            if index >= 0:
                history.setCurrentIndex(index)
        return target_identity


class CrossModuleAdapterRegistry:
    def __init__(self):
        self._adapters: dict[str, CrossModuleAdapter] = {}

    def register(self, adapter: CrossModuleAdapter):
        if adapter.module_id in self._adapters:
            raise ValueError(f"Duplicate integration module: {adapter.module_id}")
        self._adapters[adapter.module_id] = adapter

    def get(self, module_id):
        return self._adapters[module_id]

    def values(self):
        return tuple(self._adapters.values())


def _historical_presence(project, run_id):
    root = Path(project.root) / "analyses/presence_absence/runs" / run_id
    required = {name: root / "tables" / f"{name}.csv" for name in
        ("classification", "quantitative_values", "presence_matrix")}
    if not (root / "metadata.json").is_file() or any(not path.is_file() for path in required.values()):
        raise FileNotFoundError("Incomplete Presence/Absence run")
    return {"metadata": json.loads((root / "metadata.json").read_text(encoding="utf-8")),
        "tables": {name: pd.read_csv(path, dtype=str, keep_default_na=False)
            for name, path in required.items()}, "run_root": root}


def _historical_go(project, run_id):
    root = Path(project.root) / "analyses/GO/runs" / run_id
    path = root / "annotation/go_annotations.csv"
    if not path.is_file() or not (root / "metadata.json").is_file():
        raise FileNotFoundError("Incomplete GO run")
    return {"metadata": json.loads((root / "metadata.json").read_text(encoding="utf-8")),
        "tables": {"annotations": pd.read_csv(path, dtype=str, keep_default_na=False)},
        "run_root": root}


def _mapping_runs(project):
    from .mapping_analysis import read_mapping_outputs
    raw = Path(project.root) / "mapping/raw"
    runs = {path.name for path in raw.iterdir() if path.is_dir() and
        (path / "metadata.json").is_file()} if raw.is_dir() else set()
    try:
        output = read_mapping_outputs(project)
        if output.metadata.get("run_id"):
            runs.add(str(output.metadata["run_id"]))
    except RuntimeError:
        pass
    return sorted(runs)


def _historical_mapping(project, run_id):
    from .mapping_analysis import read_mapping_outputs
    output = read_mapping_outputs(project)
    if str(output.metadata.get("run_id")) != str(run_id):
        raise FileNotFoundError("Historical mapping data are not available for safe cross-module resolution.")
    return output


def default_registry():
    from . import (complex_analysis, differential_analysis, go_analysis,
        interpro_pfam_analysis, kegg_analysis, mapping_analysis, mitocarta_analysis,
        mtdna_analysis, presence_analysis, proteomics_qc, reactome_analysis,
        string_analysis)
    registry = CrossModuleAdapterRegistry()
    for module_id, listing, loader, tables, row_level in (
        ("mapping", _mapping_runs, _historical_mapping, ("catalog",), True),
        ("presence_absence", presence_analysis.list_presence_runs, _historical_presence,
            ("classification",), True),
        ("go", go_analysis.list_go_runs, _historical_go, ("annotations",), False),
        ("kegg", kegg_analysis.list_kegg_runs, kegg_analysis.read_kegg_outputs,
            ("mapping",), False),
        ("reactome", reactome_analysis.list_reactome_runs, reactome_analysis.read_reactome_outputs,
            ("mapping", "membership"), False),
        ("mitocarta", mitocarta_analysis.list_mitocarta_runs, mitocarta_analysis.read_mitocarta_outputs,
            ("mapping", "membership"), False),
        ("domains", interpro_pfam_analysis.list_interpro_pfam_runs,
            interpro_pfam_analysis.read_interpro_pfam_outputs,
            ("mapping", "interpro_locations", "pfam_locations"), False),
        ("string", string_analysis.list_string_runs, string_analysis.read_string_outputs,
            ("mapping", "node_metrics"), False),
        ("complexes", complex_analysis.list_complex_runs, complex_analysis.read_complex_outputs,
            ("mapping", "protein_to_complexes"), False),
        ("mtdna_evidence", mtdna_analysis.list_mtdna_runs, mtdna_analysis.read_mtdna_outputs,
            ("mapping", "entity_to_evidence"), False),
        ("proteomics_qc", proteomics_qc.list_runs, proteomics_qc.load_run,
            ("feature_metadata", "feature_detection"), True),
        ("differential", differential_analysis.list_runs, differential_analysis.load_run,
            ("all_results", "qualitative_candidates"), True),
    ):
        registry.register(CrossModuleAdapter(module_id, listing, loader, tables, row_level))
    return registry


def _sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _first(row, names):
    for name in names:
        value = str(row.get(name, "")).strip()
        if value and value.casefold() not in {"na", "nan", "none", "null"}:
            return value
    return ""


def _project_id(project):
    return hashlib.sha256((str(Path(project.root).resolve()).casefold() + "|"
        + str(project.config.get("created_at", ""))).encode("utf-8")).hexdigest()


def _frozen_values_fingerprint(values, identifier_column, columns):
    """Fingerprint all original rows and the exact selected sample-column set."""
    columns = sorted(columns)
    digest = hashlib.sha256()
    for position, row in enumerate(values.to_dict("records"), 1):
        normalized = {}
        for name in columns:
            raw = str(row.get(name, "")).strip()
            if raw.casefold() in {"na", "nan", "null"}:
                raw = ""
            if raw:
                try:
                    raw = format(Decimal(raw).normalize(), "f")
                except InvalidOperation:
                    pass
            normalized[name] = raw
        payload = [position, str(row.get(identifier_column, "")), normalized]
        digest.update(json.dumps(payload, ensure_ascii=False, sort_keys=True).encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()
def _quantitative_lineage(outputs, module_id, project_id):
    if module_id not in {"proteomics_qc", "differential", "presence_absence"}:
        return InputLineage(project_id)
    root = Path(outputs["run_root"])
    meta = outputs["metadata"]
    if module_id == "proteomics_qc":
        matrix = root / "input/quantitative_matrix.csv"
        samples = pd.read_csv(root / "input/sample_metadata.csv", dtype=str, keep_default_na=False)
        features = pd.read_csv(root / "input/feature_metadata.csv", dtype=str, keep_default_na=False)
        values = pd.read_csv(matrix, dtype=str, keep_default_na=False)
        renamed = {row.sample_id: row.column_name for row in samples.itertuples(index=False)}
        values = values.rename(columns=renamed)
        values["original_identifier"] = features.original_identifier
        fingerprint = _frozen_values_fingerprint(values, "original_identifier", list(renamed.values()))
        return InputLineage(project_id, fingerprint, _sha(matrix), "Frozen QC original values")
    if module_id == "differential":
        parent_meta = json.loads((root / "input/preparation_parameters.json").read_text(encoding="utf-8"))
        raw_hash = parent_meta.get("input_sha256", {}).get("original_quantitative_matrix.csv")
        parent = Path(meta.get("preparation_run_path", ""))
        matrix = parent / "input/original_quantitative_matrix.csv"
        if matrix.is_file():
            samples = pd.read_csv(root / "input/sample_metadata.csv", dtype=str, keep_default_na=False)
            features = pd.read_csv(root / "input/feature_metadata.csv", dtype=str, keep_default_na=False)
            values = pd.read_csv(matrix, dtype=str, keep_default_na=False)
            renamed = {row.sample_id: row.column_name for row in samples.itertuples(index=False)}
            values = values.rename(columns=renamed)
            values["original_identifier"] = features.original_identifier
            fingerprint = _frozen_values_fingerprint(values, "original_identifier", list(renamed.values()))
            return InputLineage(project_id, fingerprint, raw_hash, "Frozen Differential parent original values")
        return InputLineage(project_id, raw_hash, None, "Frozen Differential parent original hash")
    if module_id == "presence_absence":
        values = outputs["tables"]["quantitative_values"]
        columns = [name for name in values if name not in
            {"source_row", "original_id", "identifier_type", "protein_group", "feature_id"}]
        fingerprint = _frozen_values_fingerprint(values, "original_id", columns)
        return InputLineage(project_id, fingerprint, None, "Frozen Presence/Absence original values")
    return InputLineage(project_id)

def _run_root(project, module_id, run_id, outputs):
    if isinstance(outputs, dict):
        return Path(outputs.get("run_root", project.root))
    return Path(getattr(outputs, "run_root", getattr(outputs, "root", project.root)))


def _outputs_metadata(outputs):
    return outputs["metadata"] if isinstance(outputs, dict) else outputs.metadata


def _frame_for(adapter, outputs, key):
    if isinstance(outputs, dict):
        if key == "feature_metadata":
            path = Path(outputs["run_root"]) / "input/feature_metadata.csv"
            return pd.read_csv(path, dtype=str, keep_default_na=False)
        return outputs.get("tables", {}).get(key)
    return getattr(outputs, key, None)


def _identities(row):
    found = []
    for identity_type, columns in IDENTITY_COLUMNS.items():
        for column in columns:
            if column in row:
                for token in tokens(row[column]):
                    if identity_type == IdentityType.NCBI_GENE and token.startswith("NCBI:"):
                        token = token[5:]
                    found.append((identity_type.value, token))
    for column in MODULE_ID_COLUMNS:
        if column in row:
            for token in tokens(row[column]):
                found.append((IdentityType.MODULE_SPECIFIC.value, f"{column}:{token}"))
    return tuple(dict.fromkeys(found))


def _artifact_hashes(root):
    candidates = ("metadata.json", "summary.csv", "input/parameters.json",
        "input/feature_metadata.csv", "results/differential_results_all.csv")
    return {relative: _sha(root / relative) for relative in candidates
        if (root / relative).is_file()}
def rebuild_cross_module_index(project, registry=None):
    registry = registry or default_registry()
    root = Path(project.root) / "analyses/integration"
    root.mkdir(parents=True, exist_ok=True)
    index_path = root / "cross_module_index.sqlite"
    temporary = root / "cross_module_index.sqlite.part"
    temporary.unlink(missing_ok=True)
    connection = sqlite3.connect(temporary)
    connection.executescript("""
        CREATE TABLE runs(module TEXT NOT NULL, run_id TEXT NOT NULL, status TEXT NOT NULL,
            run_path TEXT NOT NULL, created_at TEXT, snapshot TEXT, project_id TEXT NOT NULL,
            input_hash TEXT, matrix_hash TEXT, lineage_source TEXT,
            PRIMARY KEY(module, run_id));
        CREATE TABLE records(id INTEGER PRIMARY KEY, module TEXT NOT NULL, run_id TEXT NOT NULL,
            record_type TEXT NOT NULL, record_index INTEGER NOT NULL, feature_id TEXT,
            source_row INTEGER, original_identifier TEXT, target TEXT, artifact TEXT);
        CREATE TABLE identifiers(record_id INTEGER NOT NULL, type TEXT NOT NULL, value TEXT NOT NULL);
        CREATE INDEX by_identifier ON identifiers(type, value);
        CREATE INDEX by_feature ON records(module, run_id, feature_id, source_row);
    """)
    project_id = _project_id(project)
    manifest_runs = []
    for adapter in registry.values():
        for run_id in adapter.list_runs(project):
            try:
                outputs = adapter.load_run(project, run_id)
                meta = _outputs_metadata(outputs)
                run_root = _run_root(project, adapter.module_id, run_id, outputs)
                lineage = _quantitative_lineage(outputs, adapter.module_id, project_id)
                snapshot = _first(meta, SNAPSHOT_FIELDS)
                connection.execute("INSERT INTO runs VALUES (?,?,?,?,?,?,?,?,?,?)",
                    (adapter.module_id, run_id, "Ready", str(run_root),
                    str(meta.get("created_at", meta.get("started_at", ""))), snapshot,
                    project_id, lineage.frozen_input_hash, lineage.matrix_hash, lineage.source))
                manifest_runs.append({"module": adapter.module_id, "run_id": run_id,
                    "path": str(run_root), "status": "Ready", "snapshot": snapshot,
                    "artifact_hashes": _artifact_hashes(run_root)})
                for table_name in adapter.record_tables:
                    frame = _frame_for(adapter, outputs, table_name)
                    if not isinstance(frame, pd.DataFrame):
                        continue
                    for position, row in enumerate(frame.to_dict("records")):
                        feature = _first(row, FEATURE_COLUMNS)
                        raw_row = _first(row, ROW_COLUMNS)
                        source_row = int(raw_row) if raw_row.isdecimal() else None
                        original = _first(row, ORIGINAL_COLUMNS)
                        identities = _identities(row)
                        target = _first(row, ("feature_id", "entity_key", "uniprot_accession",
                            "ncbi_gene_id", "gene_symbol", "Reactome_ID", "pathway_id",
                            "complex_id", "string_protein_id", "GO_ID", "go_id", "kegg_gene_id",
                            "interpro_id", "pfam_id")) or original
                        if not (feature or source_row or identities or target):
                            continue
                        artifact = table_name
                        cursor = connection.execute("INSERT INTO records(module,run_id,record_type,record_index,feature_id,source_row,original_identifier,target,artifact) VALUES (?,?,?,?,?,?,?,?,?)",
                            (adapter.module_id, run_id, table_name, position, feature, source_row,
                             original, target, artifact))
                        connection.executemany("INSERT INTO identifiers VALUES (?,?,?)",
                            ((cursor.lastrowid, kind, value) for kind, value in identities))
            except (OSError, ValueError, RuntimeError, KeyError, sqlite3.Error) as error:
                run_root = Path(project.root) / "analyses" / adapter.module_id / "runs" / str(run_id)
                connection.execute("INSERT OR REPLACE INTO runs VALUES (?,?,?,?,?,?,?,?,?,?)",
                    (adapter.module_id, run_id, "Historical mapping unavailable" if adapter.module_id == "mapping" else "Incomplete", str(run_root), "", "",
                     project_id, None, None, str(error)))
                manifest_runs.append({"module": adapter.module_id, "run_id": run_id,
                    "path": str(run_root), "status": "Historical mapping unavailable" if adapter.module_id == "mapping" else "Incomplete", "reason": str(error)})
    connection.commit()
    if connection.execute("PRAGMA quick_check").fetchone()[0] != "ok":
        connection.close(); temporary.unlink(missing_ok=True)
        raise RuntimeError("Cross-module index integrity check failed")
    connection.close()
    os.replace(temporary, index_path)
    manifest = {"schema_version": SCHEMA_VERSION, "build_timestamp": datetime.now(timezone.utc).isoformat(),
        "project": str(project.root), "project_id": project_id, "runs": manifest_runs}
    manifest_path = root / "cross_module_index_manifest.json"
    staged = manifest_path.with_suffix(".json.part")
    staged.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    os.replace(staged, manifest_path)
    return CrossModuleIndex(index_path, registry)


class CrossModuleIndex:
    def __init__(self, path, registry=None):
        self.path = Path(path)
        self.registry = registry or default_registry()

    def healthy(self):
        try:
            with closing(sqlite3.connect(self.path)) as connection:
                return connection.execute("PRAGMA quick_check").fetchone()[0] == "ok" and \
                    connection.execute("SELECT COUNT(*) FROM runs").fetchone() is not None
        except (OSError, sqlite3.DatabaseError):
            return False

    def run_status(self, module_id, run_id):
        with closing(sqlite3.connect(self.path)) as connection:
            row = connection.execute("SELECT status,lineage_source FROM runs WHERE module=? AND run_id=?",
                (module_id, run_id)).fetchone()
        return tuple(row) if row else None
    def search(self, query):
        query = str(query).strip()
        if not query:
            return []
        with closing(sqlite3.connect(self.path)) as connection:
            connection.row_factory = sqlite3.Row
            return [dict(row) for row in connection.execute("""
                SELECT DISTINCT r.*, u.created_at, u.snapshot, u.status, u.input_hash
                FROM records r JOIN runs u USING(module,run_id)
                LEFT JOIN identifiers i ON i.record_id=r.id
                WHERE lower(r.feature_id)=lower(?) OR CAST(r.source_row AS TEXT)=?
                   OR lower(r.original_identifier)=lower(?) OR lower(i.value)=lower(?)
                   OR lower(r.target)=lower(?)
                ORDER BY r.module,r.run_id,r.id""", (query, query, query, query, query))]

    def context(self, record_id):
        with closing(sqlite3.connect(self.path)) as connection:
            connection.row_factory = sqlite3.Row
            row = connection.execute("SELECT r.*,u.project_id,u.input_hash,u.matrix_hash,u.lineage_source,u.snapshot FROM records r JOIN runs u USING(module,run_id) WHERE r.id=?", (record_id,)).fetchone()
            if row is None:
                raise KeyError(record_id)
            identifiers = connection.execute("SELECT type,value FROM identifiers WHERE record_id=?", (record_id,)).fetchall()
        grouped = {kind.value: tuple(item["value"] for item in identifiers if item["type"] == kind.value)
            for kind in IdentityType}
        return EntityContext(row["module"], row["run_id"], row["record_type"],
            row["feature_id"] or None, row["source_row"], row["original_identifier"] or None,
            uniprot_accessions=grouped[IdentityType.UNIPROT.value],
            gene_symbols=grouped[IdentityType.GENE_SYMBOL.value],
            ncbi_gene_ids=grouped[IdentityType.NCBI_GENE.value],
            ensembl_gene_ids=grouped[IdentityType.ENSEMBL_GENE.value],
            ensembl_protein_ids=grouped[IdentityType.ENSEMBL_PROTEIN.value],
            refseq_ids=grouped[IdentityType.REFSEQ.value],
            module_specific_ids=tuple(tuple(value.split(":", 1)) for value in
                grouped[IdentityType.MODULE_SPECIFIC.value]),
            protein_group_members=tokens(row["original_identifier"])
                if len(tokens(row["original_identifier"])) > 1 else (),
            input_lineage=InputLineage(row["project_id"], row["input_hash"],
                row["matrix_hash"], row["lineage_source"]),
            database_snapshot_provenance=row["snapshot"] or None)

    def matches(self, context, *, target_module=None, include_other_lineages=False):
        with closing(sqlite3.connect(self.path)) as connection:
            connection.row_factory = sqlite3.Row
            runs = [dict(row) for row in connection.execute("SELECT * FROM runs WHERE module<>?",
                (context.source_module,)) if target_module is None or row["module"] == target_module]
            results = []
            identity_pairs = [(kind.value, value) for kind, values in (
                (IdentityType.UNIPROT, context.uniprot_accessions),
                (IdentityType.NCBI_GENE, context.ncbi_gene_ids),
                (IdentityType.ENSEMBL_GENE, context.ensembl_gene_ids),
                (IdentityType.ENSEMBL_PROTEIN, context.ensembl_protein_ids),
                (IdentityType.REFSEQ, context.refseq_ids),
                (IdentityType.GENE_SYMBOL, context.gene_symbols)) for value in values]
            if any(kind != IdentityType.GENE_SYMBOL.value for kind, _ in identity_pairs):
                identity_pairs = [(kind, value) for kind, value in identity_pairs
                    if kind != IdentityType.GENE_SYMBOL.value]
            for run in runs:
                if run["status"] != "Ready":
                    results.append(CrossModuleMatch(run["module"], run["run_id"], "",
                        MatchStatus.INCOMPLETE_RUN, run["lineage_source"] or "Required artifact unavailable"))
                    continue
                lineage = InputLineage(run["project_id"], run["input_hash"],
                    run["matrix_hash"], run["lineage_source"])
                compatible = bool(context.input_lineage and context.input_lineage.compatible(lineage))
                if not compatible and not include_other_lineages:
                    continue
                hits = {}
                if compatible and context.source_row is not None:
                    for row in connection.execute("SELECT * FROM records WHERE module=? AND run_id=? AND source_row=? AND (feature_id=? OR feature_id IS NULL OR feature_id='' OR ?='')",
                        (run["module"], run["run_id"], context.source_row, context.feature_id or "", context.feature_id or "")):
                        hits[row["id"]] = (row, MatchStatus.EXACT, "Matched by exact source row", context.feature_id)
                for kind, value in identity_pairs:
                    for row in connection.execute("SELECT r.* FROM identifiers i JOIN records r ON r.id=i.record_id WHERE r.module=? AND r.run_id=? AND i.type=? AND i.value=?",
                        (run["module"], run["run_id"], kind, value)):
                        if row["id"] not in hits:
                            basis = {"uniprot": "Matched by exact UniProt accession",
                                "ncbi_gene": "Matched by NCBI Gene ID",
                                "gene_symbol": "Matched by persisted Gene Symbol"}.get(kind, f"Matched by exact {kind}")
                            hits[row["id"]] = (row, MatchStatus.MAPPED, basis, value)
                symbol_targets = {row["target"] for row, _, basis, _ in hits.values()
                    if "Gene Symbol" in basis}
                if len(symbol_targets) > 1:
                    hits = {key: (row, MatchStatus.AMBIGUOUS,
                        "Multiple target entities share a persisted Gene Symbol", value)
                        if "Gene Symbol" in basis else (row, status, basis, value)
                        for key, (row, status, basis, value) in hits.items()}
                if len(context.uniprot_accessions) > 1 or len(context.ncbi_gene_ids) > 1:
                    hits = {key: (row, MatchStatus.AMBIGUOUS, "Multiple source biological identities", value)
                        for key, (row, _, _, value) in hits.items()}
                for row, status, basis, value in hits.values():
                    if not compatible:
                        basis += " — biological identifier match in a different or unverified input lineage"
                    results.append(CrossModuleMatch(run["module"], run["run_id"], row["target"] or "",
                        status, basis, value, row["target"] or "", compatible,
                        run["snapshot"] or "", run["created_at"] or "", row["artifact"] or "",
                        {"record_id": str(row["id"]), "record_type": row["record_type"],
                         "source_row": str(row["source_row"] or ""), "feature_id": row["feature_id"] or ""}))
            return results


def open_cross_module_index(project, registry=None):
    registry = registry or default_registry()
    root = Path(project.root) / "analyses/integration"
    path = root / "cross_module_index.sqlite"
    manifest_path = root / "cross_module_index_manifest.json"
    index = CrossModuleIndex(path, registry)
    if path.is_file() and index.healthy() and manifest_path.is_file():
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            recorded = {(item["module"], item["run_id"]) for item in manifest["runs"]}
            current = {(adapter.module_id, run_id) for adapter in registry.values()
                for run_id in adapter.list_runs(project)}
            if manifest["schema_version"] == SCHEMA_VERSION and recorded == current:
                unchanged = all(item.get("artifact_hashes", {}) ==
                    _artifact_hashes(Path(item["path"])) for item in manifest["runs"]
                    if item.get("status") == "Ready")
                if unchanged:
                    return index
        except (OSError, ValueError, KeyError):
            pass
    return rebuild_cross_module_index(project, registry)
