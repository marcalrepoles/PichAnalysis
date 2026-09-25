from __future__ import annotations
from .resources import r_script, r_scripts_dir

import json
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from .database_manager import DatabaseManager
from .databases.reactome import CORE_FILES
from .organism import get_organism
from .project import Project
from .r_runtime import RResult, RRuntime


class ReactomeAnalysisError(RuntimeError):
    pass


class ReactomeDatabaseUnavailableError(ReactomeAnalysisError):
    pass


class UnsupportedReactomeOrganismError(ReactomeAnalysisError):
    pass


class TargetOutsideBackgroundError(ReactomeAnalysisError):
    def __init__(self, details: dict[str, Any]):
        self.details = dict(details)
        outside = self.details.get("entities_outside_background", [])
        if isinstance(outside, str):
            outside = [outside]
        self.details["entities_outside_background"] = outside
        count = len(outside)
        super().__init__(f"Target outside background: {count} canonical Reactome entity/entities require a user decision.")


class ReactomeRExecutionError(ReactomeAnalysisError):
    pass


class MissingReactomeOutputError(ReactomeAnalysisError):
    pass


@dataclass(frozen=True)
class ReactomeReadiness:
    ready: bool
    reason: str


@dataclass(frozen=True)
class ReactomeParameters:
    target_selection: str = "All mapped entities"
    background_selection: str = "All mapped entities"
    manual_rows: tuple[int, ...] = ()
    minimum_overlap: int = 3
    fdr_cutoff: float = 0.05
    top_n: int = 20
    allow_target_outside_background: bool = False
    mapping_policy: str = "strict unique mapping; UniProt first; NCBI Gene fallback"


@dataclass(frozen=True)
class ReactomeOutputs:
    mapping: pd.DataFrame
    unmapped: pd.DataFrame
    ambiguous: pd.DataFrame
    membership: pd.DataFrame
    hierarchy: pd.DataFrame
    ancestry: pd.DataFrame
    frequency: pd.DataFrame
    enrichment: pd.DataFrame
    significant: pd.DataFrame
    excluded: pd.DataFrame
    summary: pd.DataFrame
    metadata: dict[str, Any]
    graphs: tuple[Path, ...]
    workbook: Path
    run_root: Path
    provenance_root: Path


def reactome_readiness(project: Project | None, manager: DatabaseManager) -> ReactomeReadiness:
    if project is None:
        return ReactomeReadiness(False, "Open a project.")
    organism = get_organism(project)
    if organism is None:
        return ReactomeReadiness(False, "Configure the project organism.")
    if organism.tax_id != "9606" or organism.name != "Homo sapiens":
        return ReactomeReadiness(False, "Unsupported organism: Reactome analysis currently supports only Homo sapiens.")
    if not manager.reactome.is_core_ready():
        return ReactomeReadiness(False, "Reactome database unavailable: install Reactome Core Data first.")
    if not (project.root / "mapping" / "tables" / "protein_catalog.csv").is_file():
        return ReactomeReadiness(False, "Run Identification and Annotation before Reactome analysis.")
    return ReactomeReadiness(True, "Ready for Reactome analysis.")


def _catalog(project: Project) -> pd.DataFrame:
    return pd.read_csv(project.root / "mapping" / "tables" / "protein_catalog.csv")


def _presence_classification(project: Project) -> pd.DataFrame:
    path = project.root / "analyses" / "presence_absence" / "tables" / "classification.csv"
    if not path.is_file():
        raise ReactomeAnalysisError("Presence/Absence outputs are required for the selected Reactome target.")
    return pd.read_csv(path)


def select_reactome_set(project: Project, selection: str, manual_rows: tuple[int, ...] | list[int] = ()) -> pd.DataFrame:
    catalog = _catalog(project)
    if selection in {"All mapped entities", "all_mapped"}:
        if "mapping_status" in catalog:
            return catalog[catalog["mapping_status"] != "unmapped"].copy()
        return catalog.copy()
    if selection in {"Manual selection", "manual"}:
        return catalog[catalog["source_row"].isin(set(manual_rows))].copy()
    classification = _presence_classification(project)
    if selection == "Reproducibly detected":
        rows = classification.loc[classification["classification"] != "Not reproducibly detected", "source_row"]
    elif selection in {"Shared", "Sporadic"}:
        rows = classification.loc[classification["classification"] == selection, "source_row"]
    else:
        label = selection[6:] if selection.startswith("class:") else selection
        rows = classification.loc[classification["classification"] == label, "source_row"]
    return catalog[catalog["source_row"].isin(rows)].copy()


def prepare_reactome_arguments(project: Project, manager: DatabaseManager, *, run_id: str,
    parameters: ReactomeParameters | None = None) -> list[str]:
    params = parameters or ReactomeParameters()
    state = reactome_readiness(project, manager)
    if not state.ready:
        if state.reason.startswith("Unsupported organism"):
            raise UnsupportedReactomeOrganismError(state.reason)
        if state.reason.startswith("Reactome database unavailable"):
            raise ReactomeDatabaseUnavailableError(state.reason)
        raise ReactomeAnalysisError(state.reason)
    if params.minimum_overlap < 1:
        raise ReactomeAnalysisError("Minimum overlap must be at least 1.")
    if not 0 <= params.fdr_cutoff <= 1:
        raise ReactomeAnalysisError("FDR threshold must be between 0 and 1.")
    if params.top_n < 1:
        raise ReactomeAnalysisError("Top N must be at least 1.")

    analysis_root = project.root / "analyses" / "Reactome"
    if (analysis_root / "runs" / run_id).exists() or (project.root / "scripts" / "runs" / f"{run_id}_reactome").exists():
        raise ReactomeAnalysisError("This Reactome run ID already exists.")
    raw = analysis_root / "raw" / run_id
    raw.mkdir(parents=True, exist_ok=False)
    target_file = raw / "target_catalog.csv"
    background_file = raw / "background_catalog.csv"
    select_reactome_set(project, params.target_selection, params.manual_rows).to_csv(target_file, index=False)
    select_reactome_set(project, params.background_selection, params.manual_rows).to_csv(background_file, index=False)

    snapshot = manager.reactome.active_snapshot()
    if snapshot is None:
        raise ReactomeDatabaseUnavailableError("Reactome database unavailable: no active Core Data snapshot.")
    manifest = manager.reactome.manifest(snapshot)
    paths = manager.reactome.core_paths()
    if set(paths) != set(CORE_FILES) or any(not paths[name].is_file() for name in CORE_FILES):
        raise ReactomeDatabaseUnavailableError("Reactome database unavailable: active snapshot Core Data files are incomplete.")
    organism = get_organism(project)
    assert organism is not None
    return [
        "--catalog", str(project.root / "mapping" / "tables" / "protein_catalog.csv"),
        "--target-file", str(target_file), "--background-file", str(background_file),
        "--output", str(analysis_root), "--project", str(project.root),
        "--pathways", str(paths["ReactomePathways.txt"]),
        "--relations", str(paths["ReactomePathwaysRelation.txt"]),
        "--uniprot-map", str(paths["UniProt2Reactome.txt"]),
        "--ncbi-map", str(paths["NCBI2Reactome.txt"]),
        "--diagrams", str(paths["humanPathwaysWithDiagrams.txt"]),
        "--summations", str(paths["pathway2summation.txt"]),
        "--run-id", run_id, "--snapshot-id", snapshot.name,
        "--snapshot-manifest", str(snapshot / "manifest.json"),
        "--snapshot-hashes", json.dumps(manager.reactome.source_hashes(), sort_keys=True),
        "--reactome-release", str(manifest.get("release_version") or "unknown"),
        "--organism", organism.name, "--tax-id", organism.tax_id,
        "--target-name", params.target_selection, "--background-name", params.background_selection,
        "--minimum-overlap", str(params.minimum_overlap), "--fdr-cutoff", str(params.fdr_cutoff),
        "--top-n", str(params.top_n),
        "--allow-target-adjustment", str(params.allow_target_outside_background).lower(),
        "--mapping-policy", params.mapping_policy,
    ]


def read_reactome_outputs(project: Project, run_id: str | None = None) -> ReactomeOutputs:
    base = project.root / "analyses" / "Reactome"
    root = base / "runs" / run_id if run_id else base
    metadata_path = root / "metadata.json" if run_id else base / "latest_metadata.json"
    expected = {
        "mapping": root / "mapping" / "reactome_entity_mapping.csv",
        "unmapped": root / "mapping" / "unmapped.csv",
        "ambiguous": root / "mapping" / "ambiguous.csv",
        "membership": root / "pathways" / "pathway_membership.csv",
        "hierarchy": root / "pathways" / "pathway_hierarchy.csv",
        "frequency": root / "frequency" / "reactome_frequency.csv",
        "enrichment": root / "enrichment" / "reactome_enrichment_all.csv",
        "significant": root / "enrichment" / "reactome_enrichment_significant.csv",
        "summary": root / "summary.csv", "metadata": metadata_path,
        "workbook": root / "Reactome_analysis.xlsx",
    }
    missing = [str(path) for path in expected.values() if not path.is_file()]
    graphs = tuple(sorted(path for path in (root / "graphs").glob("*.*") if path.suffix.lower() in {".png", ".pdf"}))
    if missing:
        detail = ", ".join(missing)
        raise MissingReactomeOutputError(f"Missing expected output: {detail}")
    try:
        metadata = json.loads(expected["metadata"].read_text(encoding="utf-8"))
        ancestry_path = root / "pathways" / "pathway_ancestry.csv"
        excluded_path = root / "enrichment" / "enrichment_excluded.csv"
        return ReactomeOutputs(
            pd.read_csv(expected["mapping"]), pd.read_csv(expected["unmapped"]),
            pd.read_csv(expected["ambiguous"]), pd.read_csv(expected["membership"]),
            pd.read_csv(expected["hierarchy"]),
            pd.read_csv(ancestry_path) if ancestry_path.is_file() else pd.DataFrame(),
            pd.read_csv(expected["frequency"]),
            pd.read_csv(expected["enrichment"]), pd.read_csv(expected["significant"]),
            pd.read_csv(excluded_path) if excluded_path.is_file() else pd.DataFrame(),
            pd.read_csv(expected["summary"]), metadata, graphs, expected["workbook"], root,
            project.root / "scripts" / "runs" / f"{metadata['run_id']}_reactome")
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as error:
        raise MissingReactomeOutputError(f"Missing expected output: invalid Reactome results ({error}).") from error


def run_reactome_analysis(project: Project, manager: DatabaseManager, runtime: RRuntime, *, run_id: str,
    parameters: ReactomeParameters | None = None, script: Path | None = None, timeout: int = 180) -> ReactomeOutputs:
    arguments = prepare_reactome_arguments(project, manager, run_id=run_id, parameters=parameters)
    entry = script or r_script("05_reactome_analysis.R")
    result: RResult = runtime.run(entry, *arguments, timeout=timeout)
    if result.returncode != 0:
        error_path = project.root / "analyses" / "Reactome" / "target_background_error.json"
        if error_path.is_file():
            try:
                raise TargetOutsideBackgroundError(json.loads(error_path.read_text(encoding="utf-8")))
            except json.JSONDecodeError:
                pass
        message = (result.stderr or result.stdout or "Rscript exited with an error.").strip()
        raise ReactomeRExecutionError(f"R execution failure: {message}")
    return read_reactome_outputs(project, run_id)


def list_reactome_runs(project: Project) -> list[str]:
    root = project.root / "analyses" / "Reactome" / "runs"
    return sorted(path.name for path in root.iterdir() if path.is_dir()) if root.is_dir() else []


def reactome_pathway_proteins(outputs: ReactomeOutputs, reactome_id: str) -> pd.DataFrame:
    members = outputs.membership[outputs.membership["Reactome_ID"].astype(str) == str(reactome_id)].copy()
    if members.empty:
        return members
    annotations = outputs.mapping[outputs.mapping["mapping_status"] == "mapped_unique"].copy()
    extra = [column for column in ("reactome_entity_key", "protein_name", "uniprot_function", "ncbi_summary") if column in annotations]
    if len(extra) > 1:
        annotations = annotations[extra].drop_duplicates("reactome_entity_key")
        members = members.merge(annotations, on="reactome_entity_key", how="left")
    return members.drop_duplicates(["reactome_entity_key", "Reactome_ID"]).reset_index(drop=True)


def reactome_protein_pathways(outputs: ReactomeOutputs, entity_key: str) -> pd.DataFrame:
    rows = outputs.membership[outputs.membership["reactome_entity_key"].astype(str) == str(entity_key)].copy()
    if rows.empty:
        return rows
    wanted = [column for column in ("Reactome_ID", "Pathway") if column in rows]
    result = rows[wanted].drop_duplicates()
    if not outputs.ancestry.empty and "pathway_id" in outputs.ancestry:
        result = result.merge(outputs.ancestry, left_on="Reactome_ID", right_on="pathway_id", how="left").drop(columns=["pathway_id"])
    return result.reset_index(drop=True)


def export_reactome_artifact(source: Path, destination: Path) -> Path:
    source, destination = Path(source), Path(destination)
    if not source.is_file():
        raise MissingReactomeOutputError(f"Missing expected output: {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    return destination
