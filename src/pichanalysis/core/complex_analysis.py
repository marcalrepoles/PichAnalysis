"""Offline Complex Portal run preparation and persisted-result loading.

Scientific coverage, frequency and enrichment are calculated only in R.
"""
from __future__ import annotations

import json
import re
import shutil
import sqlite3
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from .complex_portal_database import TABLE_COLUMNS, UNIPROT, canonicalize_uniprot_for_complex_portal
from .r_runtime import RRuntime


class ComplexAnalysisError(RuntimeError): pass
class ComplexPortalUnavailableError(ComplexAnalysisError): pass
class ComplexUnsupportedOrganismError(ComplexAnalysisError): pass
class NoComplexProteinsError(ComplexAnalysisError): pass
class ComplexPresenceRequiredError(ComplexAnalysisError): pass
class InvalidComplexComponentError(ComplexAnalysisError): pass
class ComplexRExecutionError(ComplexAnalysisError): pass
class MissingComplexOutputError(ComplexAnalysisError): pass


class ComplexTargetOutsideBackgroundError(ComplexAnalysisError):
    def __init__(self, details):
        self.details = details
        super().__init__(f"Target outside background: {len(details['entities_outside_background'])} canonical protein(s) require authorization.")


@dataclass(frozen=True)
class ComplexParameters:
    target_selection: str = "All mapped proteins"
    background_selection: str = "All mapped proteins"
    manual_rows: tuple[int, ...] = ()
    background_manual_rows: tuple[int, ...] = ()
    minimum_overlap: int = 3
    fdr_cutoff: float = 0.05
    top_n: int = 20
    allow_target_outside_background: bool = False


TABLES = {
    "mapping": "mapping/protein_mapping.csv", "ambiguous": "mapping/ambiguous.csv", "unmapped": "mapping/unmapped.csv",
    "component_groups": "coverage/component_groups.csv", "alternative_component_groups": "coverage/alternative_component_groups.csv",
    "complex_coverage": "coverage/complex_coverage.csv", "coverage_class_counts": "coverage/coverage_class_counts.csv",
    "direct_participants": "coverage/direct_participants.csv", "complex_details": "coverage/complex_details.csv",
    "stoichiometry_summary": "coverage/stoichiometry_summary.csv", "nonprotein_summary": "coverage/nonprotein_summary.csv",
    "nested_complex_summary": "coverage/nested_complex_summary.csv", "complex_frequency": "frequency/complex_frequency.csv",
    "complex_enrichment_all": "enrichment/complex_enrichment_all.csv",
    "complex_enrichment_significant": "enrichment/complex_enrichment_significant.csv",
    "complex_enrichment_excluded": "enrichment/complex_enrichment_excluded.csv",
    "protein_to_complexes": "navigation/protein_to_complexes.csv", "complex_to_proteins": "navigation/complex_to_proteins.csv",
    "summary": "summary.csv",
}
PLOTS = ("coverage_classes", "top_complex_coverage", "top_complex_frequency", "enrichment_dot_plot",
         "complex_size_vs_coverage", "target_member_distribution")


def _tokens(value):
    if value is None or pd.isna(value):
        return []
    return [token for token in re.split(r"[;,|\s]+", str(value).strip()) if token and token.upper() != "NA"]


def map_complex_experiment(project, database):
    """Resolve each experimental entity once; retain every original accession."""
    path = Path(project.root) / "mapping/tables/protein_catalog.csv"
    if not path.is_file():
        raise ComplexAnalysisError("Experimental protein catalog is unavailable.")
    catalog = pd.read_csv(path, dtype=str).fillna("")
    if not {"source_row", "uniprot_accession"}.issubset(catalog):
        raise ComplexAnalysisError("Protein catalog lacks source_row or uniprot_accession.")
    records = []
    for source_row, frame in catalog.groupby("source_row", sort=False):
        accessions = sorted({token for value in frame.uniprot_accession for token in _tokens(value)})
        normalized = [(accession, canonicalize_uniprot_for_complex_portal(accession)) for accession in accessions]
        canonical = sorted({item["lookup_accession"] for _, item in normalized if UNIPROT.fullmatch(item["lookup_accession"])})
        status = "mapped_unique" if len(canonical) == 1 else "ambiguous" if len(canonical) > 1 else "unmapped"
        membership = ""
        if status == "mapped_unique":
            membership = "member_of_one_or_more_complexes" if database.get_complexes_for_uniprot(canonical[0]) else "no_curated_complex_membership"
        records.append({
            "source_row": str(source_row), "original_id": ";".join(sorted(set(frame.get("original_id", pd.Series(dtype=str)).astype(str)))),
            "original_uniprot": ";".join(accessions), "canonical_uniprot": canonical[0] if status == "mapped_unique" else "",
            "candidate_canonical_uniprot": ";".join(canonical),
            "gene_symbol": ";".join(sorted({value for value in frame.get("gene_symbol", pd.Series(dtype=str)).astype(str) if value})),
            "isoform_normalized": any(item["isoform_normalized"] for _, item in normalized),
            "mapping_status": status, "membership_status": membership, "source_rows": str(source_row),
        })
    return pd.DataFrame(records, columns=("source_row", "original_id", "original_uniprot", "canonical_uniprot",
        "candidate_canonical_uniprot", "gene_symbol", "isoform_normalized", "mapping_status", "membership_status", "source_rows"))


def _selected_rows(project, selection, manual_rows):
    if selection in {"All mapped proteins", "all_mapped"}:
        return None
    if selection in {"Manual selection", "manual"}:
        return {str(value) for value in manual_rows}
    path = Path(project.root) / "analyses/presence_absence/tables/classification.csv"
    if not path.is_file():
        raise ComplexPresenceRequiredError("Presence/Absence outputs are required for this target.")
    classes = pd.read_csv(path, dtype=str).fillna("")
    if not {"source_row", "classification"}.issubset(classes):
        raise ComplexPresenceRequiredError("Presence/Absence classification is incomplete.")
    label = selection[6:] if selection.startswith("class:") else selection
    selected = classes.classification != "Not reproducibly detected" if label == "Reproducibly detected" else classes.classification == label
    return set(classes.loc[selected, "source_row"].astype(str))


def _canonical_set(mapping, rows):
    frame = mapping[mapping.mapping_status == "mapped_unique"]
    if rows is not None:
        frame = frame[frame.source_row.isin(rows)]
    return sorted(set(frame.canonical_uniprot) - {""})


def _source_rows(mapping, canonical, selected_rows=None):
    frame = mapping[(mapping.mapping_status == "mapped_unique") & (mapping.canonical_uniprot.isin(canonical))]
    if selected_rows is not None:
        frame = frame[frame.source_row.isin(selected_rows)]
    return pd.DataFrame([{
        "canonical_uniprot": accession,
        "target_supporting_rows": ";".join(sorted(frame.loc[frame.canonical_uniprot == accession, "source_row"], key=lambda value: (len(value), value))),
    } for accession in canonical], columns=("canonical_uniprot", "target_supporting_rows"))


def prepare_complex_run(project, manager, run_id, parameters=None):
    p = parameters or ComplexParameters()
    database = manager.complex_portal
    snapshot = database.active_snapshot()
    if snapshot is None:
        raise ComplexPortalUnavailableError("Complex Portal database is unavailable.")
    try:
        database.validate_snapshot(snapshot)
    except (OSError, ValueError) as error:
        raise ComplexPortalUnavailableError(f"Complex Portal snapshot is invalid: {error}") from error
    if str(project.config.get("organism_tax_id") or "") != "9606":
        raise ComplexUnsupportedOrganismError("Complex Portal analysis currently supports Homo sapiens only.")
    if p.minimum_overlap < 1 or not 0 <= p.fdr_cutoff <= 1 or p.top_n < 1:
        raise ComplexAnalysisError("Minimum overlap, FDR cutoff or Top N is invalid.")
    root = Path(project.root) / "analyses/Complexes"
    run = root / "runs" / run_id
    provenance = Path(project.root) / "scripts/runs" / f"{run_id}_complexes"
    if run.exists() or provenance.exists():
        raise ComplexAnalysisError("This Complex Portal run ID already exists.")
    mapping = map_complex_experiment(project, database)
    target_rows = _selected_rows(project, p.target_selection, p.manual_rows)
    background_rows = _selected_rows(project, p.background_selection, p.background_manual_rows)
    target = _canonical_set(mapping, target_rows)
    background = _canonical_set(mapping, background_rows)
    if not target or not background:
        raise NoComplexProteinsError("No uniquely resolved canonical UniProt proteins are available in target or background.")
    outside = sorted(set(target) - set(background))
    if outside and not p.allow_target_outside_background:
        raise ComplexTargetOutsideBackgroundError({"entities_outside_background": outside,
            "initial_target_size": len(target), "initial_background_size": len(background)})
    if outside:
        target = sorted(set(target) & set(background))
    if not target:
        raise NoComplexProteinsError("The authorized target contains no canonical proteins after background adjustment.")
    manifest = database.manifest(snapshot)
    (run / "inputs").mkdir(parents=True)
    provenance.mkdir(parents=True)
    (run / "mapping").mkdir()
    mapping.to_csv(run / TABLES["mapping"], index=False)
    mapping[mapping.mapping_status == "ambiguous"].to_csv(run / TABLES["ambiguous"], index=False)
    mapping[mapping.mapping_status == "unmapped"].to_csv(run / TABLES["unmapped"], index=False)
    _source_rows(mapping, target, target_rows).to_csv(run / "inputs/target.csv", index=False)
    _source_rows(mapping, background, background_rows).rename(columns={"target_supporting_rows": "background_supporting_rows"}).to_csv(run / "inputs/background.csv", index=False)
    sqlite_path = snapshot / "complex_portal.sqlite"
    try:
        with sqlite3.connect(sqlite_path) as conn:
            for name in TABLE_COLUMNS:
                pd.read_sql_query(f'SELECT * FROM "{name}"', conn).to_csv(run / "inputs" / f"{name}.csv", index=False)
    except (sqlite3.Error, OSError) as error:
        raise ComplexPortalUnavailableError(f"Cannot read the local Complex Portal snapshot: {error}") from error
    metadata = {
        "run_id": run_id, "created_at": datetime.now(timezone.utc).isoformat(),
        "snapshot_id": snapshot.name, "snapshot_path": str(snapshot),
        "complex_portal_release": manifest.get("complex_portal_release", "unknown"),
        "raw_sha256": manifest["raw_file"]["sha256"], "sqlite_sha256": manifest["sqlite_sha256"],
        "database_scope": "manually_curated", "organism_tax_id": "9606",
        "target_definition": p.target_selection, "background_definition": p.background_selection,
        "initial_target_size": len(target) + len(outside), "initial_background_size": len(background),
        "target_size": len(target), "background_size": len(background),
        "target_outside_background": outside, "target_adjustment_authorized": bool(outside and p.allow_target_outside_background),
        "minimum_overlap": p.minimum_overlap, "fdr_cutoff": p.fdr_cutoff, "top_n": p.top_n,
        "canonicalization": "Valid UniProt isoform suffixes map to canonical accession only within this module; original input retained.",
        "component_group_rule": "One requirement per official expanded protein component, independent of stoichiometric copy count.",
        "alternative_component_rule": "A bracketed alternative set is one group; any detected option covers it.",
        "coverage_definition": "Detected target identity for every expected protein component group; does not demonstrate physical assembly, co-localization, activity or correct stoichiometry.",
        "enrichment_definition": "One-sided Fisher greater for unique canonical UniProt membership in target versus background minus target; BH across eligible curated complexes.",
        "network_access": False,
    }
    (run / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    (provenance / "parameters.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    (provenance / "complex_portal_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    scripts = Path(__file__).resolve().parents[3] / "r_scripts"
    for source in (scripts / "09_complex_analysis.R", scripts / "lib/complex_analysis.R"):
        shutil.copy2(source, provenance / source.name)
    return run, provenance


def read_complex_outputs(project, run_id):
    run = Path(project.root) / "analyses/Complexes/runs" / run_id
    mandatory = {name: run / relative for name, relative in TABLES.items()}
    mandatory["metadata"] = run / "metadata.json"
    # Tables and metadata define a valid historical run. Optional presentation
    # artifacts may be absent, but must never be borrowed from another run.
    missing = [str(path) for path in mandatory.values() if not path.is_file()]
    if missing:
        raise MissingComplexOutputError("Missing expected Complex Portal output: " + ", ".join(missing))
    try:
        return {"run_root": run, "metadata": json.loads(mandatory["metadata"].read_text(encoding="utf-8")),
            "tables": {name: pd.read_csv(path) for name, path in mandatory.items() if name in TABLES},
            "workbook": (run / "Complex_analysis.xlsx") if (run / "Complex_analysis.xlsx").is_file() else None,
            "plots": tuple(sorted((run / "plots").glob("*.*")))}
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as error:
        raise MissingComplexOutputError(f"Invalid Complex Portal run: {error}") from error


def list_complex_runs(project):
    folder = Path(project.root) / "analyses/Complexes/runs"
    return sorted(path.name for path in folder.iterdir() if path.is_dir()) if folder.is_dir() else []


def run_complex_analysis(project, manager, runtime: RRuntime, *, run_id, parameters=None, script=None, timeout=300):
    run, provenance = prepare_complex_run(project, manager, run_id, parameters)
    entry = script or Path(__file__).resolve().parents[3] / "r_scripts/09_complex_analysis.R"
    try:
        result = runtime.run(entry, "--run", str(run), "--provenance", str(provenance), timeout=timeout)
    except RuntimeError as error:
        raise ComplexRExecutionError(f"R execution failure: {error}") from error
    if result.returncode:
        message = (result.stderr or result.stdout or "Rscript exited with an error.").strip()
        if "Invalid Complex Portal component representation" in message:
            raise InvalidComplexComponentError(message)
        raise ComplexRExecutionError("R execution failure: " + message)
    outputs = read_complex_outputs(project, run_id)
    latest = Path(project.root) / "analyses/Complexes"
    latest.mkdir(parents=True, exist_ok=True)
    shutil.copy2(run / "summary.csv", latest / "summary.csv")
    shutil.copy2(run / "Complex_analysis.xlsx", latest / "Complex_analysis.xlsx")
    (latest / "mapping").mkdir(exist_ok=True)
    for name in ("protein_mapping.csv", "ambiguous.csv", "unmapped.csv"):
        shutil.copy2(run / "mapping" / name, latest / "mapping" / name)
    return outputs
