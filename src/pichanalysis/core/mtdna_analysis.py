"""Offline gene-level preparation, provenance and history for mtDNA Evidence analysis.

R alone computes frequencies, contingency tests, FDR and visualizations.
"""
from __future__ import annotations
from .resources import r_script, r_scripts_dir

import json
import re
import shutil
import sqlite3
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from .r_runtime import RRuntime


class MtdnaAnalysisError(RuntimeError): pass
class MtdnaUnavailableError(MtdnaAnalysisError): pass
class MtdnaPresenceRequiredError(MtdnaAnalysisError): pass
class MtdnaNoEntitiesError(MtdnaAnalysisError): pass
class MtdnaMissingOutputError(MtdnaAnalysisError): pass
class MtdnaRExecutionError(MtdnaAnalysisError): pass


class MtdnaTargetOutsideBackgroundError(MtdnaAnalysisError):
    def __init__(self, details):
        self.details = details
        super().__init__(f"{len(details['entities_outside_background'])} target entities are outside the background.")


@dataclass(frozen=True)
class MtdnaParameters:
    target_selection: str = "All mapped entities"
    background_selection: str = "All mapped entities"
    manual_rows: tuple[int, ...] = ()
    background_manual_rows: tuple[int, ...] = ()
    minimum_overlap: int = 3
    fdr_cutoff: float = 0.05
    top_n: int = 20
    allow_target_outside_background: bool = False


TABLES = {
    "mapping": "mapping/experimental_entity_mapping.csv",
    "ambiguous": "mapping/ambiguous.csv", "unmapped": "mapping/unmapped.csv",
    "entity_evidence": "evidence/entity_evidence_summary.csv",
    "evidence_records": "evidence/target_evidence_records.csv",
    "category_frequency": "categories/category_frequency.csv",
    "source_frequency": "sources/source_frequency.csv",
    "source_agreement": "sources/source_agreement_target.csv",
    "source_count": "sources/source_count_distribution.csv",
    "enrichment_all": "enrichment/category_enrichment_all.csv",
    "enrichment_significant": "enrichment/category_enrichment_significant.csv",
    "enrichment_excluded": "enrichment/category_enrichment_excluded.csv",
    "comparison": "comparison/target_category_matrix.csv",
    "entity_to_evidence": "navigation/entity_to_evidence.csv",
    "category_to_entities": "navigation/category_to_entities.csv",
    "source_to_entities": "navigation/source_to_entities.csv",
    "summary": "summary.csv",
}
PLOTS = ("category_frequency", "source_frequency", "source_count_distribution",
         "enrichment_dot_plot", "evidence_records_per_entity", "target_category_comparison")


def _tokens(value):
    if value is None or pd.isna(value): return []
    return [part for part in re.split(r"[;,|\s]+", str(value).strip()) if part and part.upper() not in {"NA", "NAN"}]


def _text(value):
    return "" if value is None or pd.isna(value) else str(value).strip()


def _gene_id(value):
    value = _text(value)
    return value[:-2] if value.endswith(".0") and value[:-2].isdigit() else value


def _read_index(snapshot):
    with sqlite3.connect(Path(snapshot) / "mtdna_evidence.sqlite") as connection:
        return {name: pd.read_sql_query(f'SELECT * FROM "{name}"', connection).fillna("")
                for name in ("entities", "evidence_records", "category_membership", "categories", "source_agreement")}


def map_mtdna_experiment(project, database):
    """Resolve catalog rows without pretending an absent evidence hit is unmapped."""
    path = Path(project.root) / "mapping/tables/protein_catalog.csv"
    if not path.is_file(): raise MtdnaAnalysisError("Experimental protein catalog is unavailable.")
    catalog = pd.read_csv(path, dtype=str, keep_default_na=False)
    if "source_row" not in catalog: raise MtdnaAnalysisError("Protein catalog lacks source_row.")
    snapshot = database.active_snapshot()
    if snapshot is None: raise MtdnaUnavailableError("mtDNA Evidence database is not ready.")
    index = _read_index(snapshot)
    entities = index["entities"]
    by_key = {row.entity_key: row for row in entities.itertuples(index=False)}
    symbol_keys = {}
    uniprot_keys = {}
    for row in entities.itertuples(index=False):
        if row.gene_symbol: symbol_keys.setdefault(str(row.gene_symbol).upper(), set()).add(row.entity_key)
        for accession in _tokens(row.uniprot_accessions):
            uniprot_keys.setdefault(accession.upper(), set()).add(row.entity_key)
    records = []
    for source_row, frame in catalog.groupby("source_row", sort=False):
        original = sorted({_text(v) for v in frame.get("original_id", pd.Series(dtype=str)) if _text(v)})
        accessions = sorted({t for value in frame.get("uniprot_accession", pd.Series(dtype=str)) for t in _tokens(value)})
        ids = sorted({_gene_id(t) for value in frame.get("ncbi_gene_id", pd.Series(dtype=str)) for t in _tokens(value)})
        symbols = sorted({t.upper() for value in frame.get("gene_symbol", pd.Series(dtype=str)) for t in _tokens(value)})
        input_status = {_text(v).lower() for v in frame.get("mapping_status", pd.Series(dtype=str))}
        candidates = set()
        if len(ids) > 1:
            status, key = "ambiguous", ""
        elif len(ids) == 1:
            key = f"NCBI:{ids[0]}"
            # The experimental gene-level ID is authoritative. A contradictory
            # evidence-index UniProt hit never silently overrides it.
            for accession in accessions:
                candidates.update(uniprot_keys.get(accession.upper(), set()))
            for symbol in symbols:
                candidates.update(symbol_keys.get(symbol, set()))
            if any(candidate != key and candidate.startswith("NCBI:") for candidate in candidates):
                status, key = "ambiguous", ""
            else:
                status = "mapped_with_mtdna_evidence" if key in by_key else "mapped_no_mtdna_evidence"
        else:
            for symbol in symbols: candidates.update(symbol_keys.get(symbol, set()))
            for accession in accessions: candidates.update(uniprot_keys.get(accession.upper(), set()))
            if len(candidates) == 1:
                key = next(iter(candidates)); status = "mapped_with_mtdna_evidence"
            elif len(candidates) > 1:
                key = ""; status = "ambiguous"
            elif len(symbols) == 1 and "unmapped" not in input_status and "ambiguous" not in input_status:
                key = f"SYMBOL:{symbols[0]}"; status = "mapped_no_mtdna_evidence"
            else:
                key = ""; status = "ambiguous" if "ambiguous" in input_status else "unmapped"
        records.append({"source_row": str(source_row), "original_identifier": ";".join(original),
            "original_uniprot": ";".join(accessions), "gene_symbol": ";".join(symbols),
            "ncbi_gene_id": ";".join(ids), "entity_key": key, "mapping_status": status,
            "evidence_status": "has_mtdna_related_evidence" if key in by_key else
                               "no_mtdna_related_evidence" if key else "not_resolved",
            "source_rows": str(source_row)})
    return pd.DataFrame(records, columns=("source_row", "original_identifier", "original_uniprot",
        "gene_symbol", "ncbi_gene_id", "entity_key", "mapping_status", "evidence_status", "source_rows"))


def _selected_rows(project, selection, manual_rows):
    if selection in {"All mapped entities", "all_mapped"}: return None
    if selection in {"Manual selection", "manual"}: return {str(value) for value in manual_rows}
    path = Path(project.root) / "analyses/presence_absence/tables/classification.csv"
    if not path.is_file(): raise MtdnaPresenceRequiredError("Presence/Absence outputs are required for this target.")
    classes = pd.read_csv(path, dtype=str, keep_default_na=False)
    if not {"source_row", "classification"}.issubset(classes):
        raise MtdnaPresenceRequiredError("Presence/Absence classification is incomplete.")
    label = selection[6:] if selection.startswith("class:") else selection
    selected = ~classes.classification.isin({"Not reproducibly detected", "Sporadic"}) if label == "Reproducibly detected" else classes.classification.eq(label)
    return set(classes.loc[selected, "source_row"].astype(str))


def _entity_set(mapping, rows):
    frame = mapping[mapping.mapping_status.str.startswith("mapped_")]
    if rows is not None: frame = frame[frame.source_row.isin(rows)]
    return sorted(set(frame.entity_key) - {""})


def _write_csv(frame, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False)


def prepare_mtdna_run(project, manager, run_id, parameters=None):
    p = parameters or MtdnaParameters()
    db = manager.mtdna_evidence
    snapshot = db.active_snapshot()
    if snapshot is None: raise MtdnaUnavailableError("mtDNA Evidence database is not ready.")
    db.validate_snapshot(snapshot)
    if str(project.config.get("organism_tax_id") or "") != "9606":
        raise MtdnaAnalysisError("mtDNA Evidence analysis currently supports Homo sapiens only.")
    if p.minimum_overlap < 1 or not 0 <= p.fdr_cutoff <= 1 or p.top_n < 1:
        raise MtdnaAnalysisError("Minimum overlap, FDR cutoff or Top N is invalid.")
    root = Path(project.root) / "analyses/mtDNA"
    run = root / "runs" / run_id
    provenance = Path(project.root) / "scripts/runs" / f"{run_id}_mtdna"
    if run.exists() or provenance.exists(): raise MtdnaAnalysisError("This mtDNA run ID already exists.")
    mapping = map_mtdna_experiment(project, db)
    target_rows = _selected_rows(project, p.target_selection, p.manual_rows)
    background_rows = _selected_rows(project, p.background_selection, p.background_manual_rows)
    target = _entity_set(mapping, target_rows)
    background = _entity_set(mapping, background_rows)
    if not target or not background: raise MtdnaNoEntitiesError("No uniquely resolved gene-level entities are available.")
    outside = sorted(set(target) - set(background))
    if outside and not p.allow_target_outside_background:
        raise MtdnaTargetOutsideBackgroundError({"entities_outside_background": outside,
            "initial_target_size": len(target), "initial_background_size": len(background)})
    if outside: target = sorted(set(target) & set(background))
    if not target: raise MtdnaNoEntitiesError("The authorized target is empty after background adjustment.")
    index = _read_index(snapshot)
    manifest = db.manifest(snapshot)
    (run / "inputs").mkdir(parents=True)
    provenance.mkdir(parents=True)
    _write_csv(mapping, run / TABLES["mapping"])
    _write_csv(mapping[mapping.mapping_status.eq("ambiguous")], run / TABLES["ambiguous"])
    _write_csv(mapping[mapping.mapping_status.eq("unmapped")], run / TABLES["unmapped"])
    for name, keys in (("target", target), ("background", background)):
        _write_csv(pd.DataFrame({"entity_key": keys}), run / "inputs" / f"{name}.csv")
    for name, frame in index.items(): _write_csv(frame, run / "inputs" / f"{name}.csv")
    classification = Path(project.root) / "analyses/presence_absence/tables/classification.csv"
    if classification.is_file():
        classes = pd.read_csv(classification, dtype=str, keep_default_na=False)
        if {"source_row", "classification"}.issubset(classes):
            _write_csv(classes[["source_row", "classification"]], run / "inputs/classifications.csv")

    metadata = {"run_id": run_id, "created_at": datetime.now(timezone.utc).isoformat(),
        "snapshot_id": snapshot.name, "snapshot_path": str(snapshot),
        "snapshot_sqlite_sha256": manifest["sqlite_sha256"],
        "mitocarta_snapshot_id": manifest["mitocarta_snapshot_id"],
        "go_snapshot_id": manifest["go_snapshot_id"], "ncbi_accession": manifest["ncbi_accession"],
        "mitocarta_source_hashes": manifest["mitocarta_source_hashes"],
        "go_source_hashes": manifest["go_source_hashes"], "ncbi_raw_sha256": manifest["ncbi_raw_sha256"],
        "target_definition": p.target_selection, "background_definition": p.background_selection,
        "initial_target_size": len(target) + len(outside), "initial_background_size": len(background),
        "target_size": len(target), "background_size": len(background),
        "target_outside_background": outside, "target_adjustment_authorized": bool(outside and p.allow_target_outside_background),
        "minimum_overlap": p.minimum_overlap, "fdr_cutoff": p.fdr_cutoff, "top_n": p.top_n,
        "mapping_rule": "Unique experimental NCBI Gene ID first; compatible unique mtDNA entity or unambiguous symbol/UniProt bridge; unresolved gene with no evidence remains in experimental universe.",
        "enrichment_definition": "Gene-level category membership; one-sided Fisher greater, target versus experimental background minus target; BH over all eligible categories.",
        "network_access": False}
    (run / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    (provenance / "parameters.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    (provenance / "mtdna_evidence_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    scripts = r_scripts_dir()
    for source in (scripts / "10_mtdna_analysis.R", scripts / "lib/mtdna_analysis.R"):
        shutil.copy2(source, provenance / source.name)
    return run, provenance


def read_mtdna_outputs(project, run_id):
    run = Path(project.root) / "analyses/mtDNA/runs" / run_id
    mandatory = {name: run / relative for name, relative in TABLES.items()}
    missing = [str(path) for path in (*mandatory.values(), run / "metadata.json") if not path.is_file()]
    if missing: raise MtdnaMissingOutputError("Missing mtDNA run artifacts: " + ", ".join(missing))
    try:
        return {"run_root": run, "metadata": json.loads((run / "metadata.json").read_text(encoding="utf-8")),
            "tables": {name: pd.read_csv(path, dtype=str, keep_default_na=False) for name, path in mandatory.items()},
            "workbook": (run / "mtDNA_analysis.xlsx") if (run / "mtDNA_analysis.xlsx").is_file() else None,
            "plots": tuple(sorted((run / "plots").glob("*.png")))}
    except (OSError, ValueError, json.JSONDecodeError) as error:
        raise MtdnaMissingOutputError(f"Invalid historical mtDNA run: {error}") from error


def list_mtdna_runs(project):
    folder = Path(project.root) / "analyses/mtDNA/runs"
    return sorted((path.name for path in folder.iterdir() if path.is_dir()), reverse=True) if folder.is_dir() else []


def run_mtdna_analysis(project, manager, runtime: RRuntime, *, run_id, parameters=None, script=None, timeout=300):
    run, provenance = prepare_mtdna_run(project, manager, run_id, parameters)
    entry = script or r_script("10_mtdna_analysis.R")
    try: result = runtime.run(entry, "--run", str(run), "--provenance", str(provenance), timeout=timeout)
    except RuntimeError as error: raise MtdnaRExecutionError(str(error)) from error
    if result.returncode: raise MtdnaRExecutionError((result.stderr or result.stdout or "Rscript failed.").strip())
    outputs = read_mtdna_outputs(project, run_id)
    latest = Path(project.root) / "analyses/mtDNA"
    shutil.copy2(run / "summary.csv", latest / "summary.csv")
    shutil.copy2(run / "mtDNA_analysis.xlsx", latest / "mtDNA_analysis.xlsx")
    return outputs
