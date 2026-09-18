from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from .database_manager import DatabaseManager
from .organism import get_organism
from .project import Project
from .r_runtime import RResult, RRuntime


class MitoCartaAnalysisError(RuntimeError): pass
class MitoCartaDatabaseUnavailableError(MitoCartaAnalysisError): pass
class UnsupportedMitoCartaOrganismError(MitoCartaAnalysisError): pass
class NoGeneResolvedEntitiesError(MitoCartaAnalysisError): pass
class PresenceAbsenceRequiredError(MitoCartaAnalysisError): pass
class NoMitoCartaGenesInTargetError(MitoCartaAnalysisError): pass
class NoMitoCartaGenesInBackgroundError(MitoCartaAnalysisError): pass
class MitoCartaRExecutionError(MitoCartaAnalysisError): pass
class MissingMitoCartaOutputError(MitoCartaAnalysisError): pass


class MitoCartaTargetOutsideBackgroundError(MitoCartaAnalysisError):
    def __init__(self, details: dict[str, Any]):
        self.details=dict(details);outside=self.details.get("entities_outside_background",[])
        if isinstance(outside,str):outside=[outside]
        self.details["entities_outside_background"]=outside
        super().__init__(f"Target outside background: {len(outside)} canonical gene entity/entities require a user decision.")


@dataclass(frozen=True)
class MitoCartaParameters:
    target_selection: str = "All mapped entities"
    background_selection: str = "All mapped entities"
    manual_rows: tuple[int,...] = ()
    minimum_overlap: int = 3
    fdr_cutoff: float = 0.05
    top_n: int = 20
    allow_target_outside_background: bool = False
    mapping_policy: str = "gene-level resolution: NCBI Gene ID, unique gene symbol, then unique UniProt"


@dataclass(frozen=True)
class MitoCartaOutputs:
    mapping: pd.DataFrame; ambiguous: pd.DataFrame; unmapped: pd.DataFrame; membership: pd.DataFrame
    coverage: pd.DataFrame; overall: pd.DataFrame; subcompartment_frequency: pd.DataFrame
    subcompartment_enrichment: pd.DataFrame; subcompartment_significant: pd.DataFrame
    pathway_frequency: pd.DataFrame; pathway_enrichment: pd.DataFrame; pathway_significant: pd.DataFrame
    pathway_excluded: pd.DataFrame; summary: pd.DataFrame; metadata: dict[str,Any]
    graphs: tuple[Path,...]; workbook: Path; run_root: Path; provenance_root: Path


def mitocarta_readiness(project: Project|None,manager: DatabaseManager)->tuple[bool,str]:
    if project is None:return False,"Open a project."
    organism=get_organism(project)
    if organism is None:return False,"Configure the project organism."
    if organism.name!="Homo sapiens" or organism.tax_id!="9606":return False,"Unsupported organism: MitoCarta analysis currently supports only Homo sapiens."
    if not manager.mitocarta.is_ready():return False,"MitoCarta database unavailable: install MitoCarta3.0 Core Data first."
    if not (project.root/"mapping"/"tables"/"protein_catalog.csv").is_file():return False,"Run Identification and Annotation before MitoCarta analysis."
    return True,"Ready for MitoCarta analysis."


def _catalog(project:Project)->pd.DataFrame:return pd.read_csv(project.root/"mapping"/"tables"/"protein_catalog.csv")


def select_mitocarta_set(project:Project,selection:str,manual_rows=())->pd.DataFrame:
    catalog=_catalog(project)
    if selection in {"All mapped entities","all_mapped"}:return catalog.copy()
    if selection in {"Manual selection","manual"}:return catalog[catalog.source_row.isin(set(manual_rows))].copy()
    path=project.root/"analyses"/"presence_absence"/"tables"/"classification.csv"
    if not path.is_file():raise PresenceAbsenceRequiredError("Presence/Absence outputs are required for the selected MitoCarta target.")
    classification=pd.read_csv(path);label=selection[6:] if selection.startswith("class:") else selection
    if selection=="Reproducibly detected":rows=classification.loc[classification.classification!="Not reproducibly detected","source_row"]
    else:rows=classification.loc[classification.classification==label,"source_row"]
    return catalog[catalog.source_row.isin(rows)].copy()


def prepare_mitocarta_arguments(project:Project,manager:DatabaseManager,*,run_id:str,parameters:MitoCartaParameters|None=None)->list[str]:
    params=parameters or MitoCartaParameters();ready,reason=mitocarta_readiness(project,manager)
    if not ready:
        if reason.startswith("Unsupported organism"):raise UnsupportedMitoCartaOrganismError(reason)
        if reason.startswith("MitoCarta database unavailable"):raise MitoCartaDatabaseUnavailableError(reason)
        raise MitoCartaAnalysisError(reason)
    if params.minimum_overlap<1:raise MitoCartaAnalysisError("Minimum overlap must be at least 1.")
    if not 0<=params.fdr_cutoff<=1:raise MitoCartaAnalysisError("FDR threshold must be between 0 and 1.")
    if params.top_n<1:raise MitoCartaAnalysisError("Top N must be at least 1.")
    root=project.root/"analyses"/"MitoCarta"
    if (root/"runs"/run_id).exists() or (project.root/"scripts"/"runs"/f"{run_id}_mitocarta").exists():raise MitoCartaAnalysisError("This MitoCarta run ID already exists.")
    raw=root/"raw"/run_id;raw.mkdir(parents=True,exist_ok=False);target=raw/"target_catalog.csv";background=raw/"background_catalog.csv"
    select_mitocarta_set(project,params.target_selection,params.manual_rows).to_csv(target,index=False);select_mitocarta_set(project,params.background_selection,params.manual_rows).to_csv(background,index=False)
    snapshot=manager.mitocarta.active_snapshot()
    if snapshot is None:raise MitoCartaDatabaseUnavailableError("MitoCarta database unavailable: no active snapshot.")
    paths={name:manager.mitocarta.table_path(name) for name in ("genes","subcompartments","pathways","hierarchy")}
    if any(path is None or not path.is_file() for path in paths.values()):raise MitoCartaDatabaseUnavailableError("MitoCarta database unavailable: active snapshot tables are incomplete.")
    organism=get_organism(project);manifest=manager.mitocarta.manifest(snapshot);assert organism
    return ["--catalog",str(project.root/"mapping"/"tables"/"protein_catalog.csv"),"--target-file",str(target),"--background-file",str(background),"--output",str(root),"--project",str(project.root),"--genes",str(paths["genes"]),"--subcompartments",str(paths["subcompartments"]),"--pathway-membership",str(paths["pathways"]),"--pathway-hierarchy",str(paths["hierarchy"]),"--run-id",run_id,"--snapshot-id",snapshot.name,"--snapshot-manifest",str(snapshot/"manifest.json"),"--snapshot-hashes",json.dumps(manager.mitocarta.source_hashes(),sort_keys=True),"--mitocarta-version",str(manifest.get("version","3.0")),"--organism",organism.name,"--tax-id",organism.tax_id,"--target-name",params.target_selection,"--background-name",params.background_selection,"--minimum-overlap",str(params.minimum_overlap),"--fdr-cutoff",str(params.fdr_cutoff),"--top-n",str(params.top_n),"--allow-target-adjustment",str(params.allow_target_outside_background).lower(),"--mapping-policy",params.mapping_policy]


def read_mitocarta_outputs(project:Project,run_id:str|None=None)->MitoCartaOutputs:
    base=project.root/"analyses"/"MitoCarta";root=base/"runs"/run_id if run_id else base;metadata_path=root/"metadata.json" if run_id else base/"latest_metadata.json"
    files={"mapping":root/"mapping/mitocarta_mapping.csv","ambiguous":root/"mapping/ambiguous.csv","unmapped":root/"mapping/unmapped.csv","membership":root/"mapping/mitocarta_membership.csv","coverage":root/"membership/mitocarta_summary.csv","overall":root/"enrichment/mitocarta_overall_enrichment.csv","sub_frequency":root/"subcompartments/subcompartment_frequency.csv","sub_all":root/"subcompartments/subcompartment_enrichment_all.csv","sub_sig":root/"subcompartments/subcompartment_enrichment_significant.csv","path_frequency":root/"mitopathways/mitopathway_frequency.csv","path_all":root/"mitopathways/mitopathway_enrichment_all.csv","path_sig":root/"mitopathways/mitopathway_enrichment_significant.csv","path_excluded":root/"mitopathways/mitopathway_enrichment_excluded.csv","summary":root/"summary.csv","workbook":root/"MitoCarta_analysis.xlsx","metadata":metadata_path}
    missing=[str(path) for path in files.values() if not path.is_file()]
    if missing:raise MissingMitoCartaOutputError("Missing expected output: "+", ".join(missing))
    try:
        frames=[pd.read_csv(files[name]) for name in ("mapping","ambiguous","unmapped","membership","coverage","overall","sub_frequency","sub_all","sub_sig","path_frequency","path_all","path_sig","path_excluded","summary")];metadata=json.loads(metadata_path.read_text(encoding="utf-8"));graphs=tuple(sorted((root/"graphs").glob("*.*")))
        return MitoCartaOutputs(*frames,metadata,graphs,files["workbook"],root,project.root/"scripts"/"runs"/f"{metadata['run_id']}_mitocarta")
    except (OSError,ValueError,KeyError,json.JSONDecodeError) as exc:raise MissingMitoCartaOutputError(f"Missing expected output: invalid MitoCarta results ({exc}).") from exc


def run_mitocarta_analysis(project:Project,manager:DatabaseManager,runtime:RRuntime,*,run_id:str,parameters:MitoCartaParameters|None=None,script:Path|None=None,timeout:int=180)->MitoCartaOutputs:
    arguments=prepare_mitocarta_arguments(project,manager,run_id=run_id,parameters=parameters);entry=script or Path(__file__).resolve().parents[3]/"r_scripts"/"06_mitocarta_analysis.R";result:RResult=runtime.run(entry,*arguments,timeout=timeout)
    if result.returncode!=0:
        error_path=project.root/"analyses"/"MitoCarta"/"target_background_error.json"
        if error_path.is_file():
            try:raise MitoCartaTargetOutsideBackgroundError(json.loads(error_path.read_text(encoding="utf-8")))
            except json.JSONDecodeError:pass
        message=(result.stderr or result.stdout or "Rscript exited with an error.").strip()
        if "No gene-resolved entities" in message:raise NoGeneResolvedEntitiesError(message)
        if "No MitoCarta genes in background" in message:raise NoMitoCartaGenesInBackgroundError(message)
        raise MitoCartaRExecutionError(f"R execution failure: {message}")
    return read_mitocarta_outputs(project,run_id)


def list_mitocarta_runs(project:Project)->list[str]:
    root=project.root/"analyses"/"MitoCarta"/"runs";return sorted(path.name for path in root.iterdir() if path.is_dir()) if root.is_dir() else []
