from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from .database_manager import DatabaseManager
from .interpro_database import InvalidAnnotationSet
from .project import Project
from .r_runtime import RResult, RRuntime


class InterProPfamAnalysisError(RuntimeError): pass
class AnnotationMetadataUnavailableError(InterProPfamAnalysisError): pass
class AnnotationSetUnavailableError(InterProPfamAnalysisError): pass
class AnnotationSetIncompleteError(InterProPfamAnalysisError): pass
class AnnotationSetCoverageError(InterProPfamAnalysisError): pass
class NoUniquelyResolvedProteinsError(InterProPfamAnalysisError): pass
class InterProPfamPresenceRequiredError(InterProPfamAnalysisError): pass
class InterProPfamRExecutionError(InterProPfamAnalysisError): pass
class MissingInterProPfamOutputError(InterProPfamAnalysisError): pass


class InterProPfamTargetOutsideBackgroundError(InterProPfamAnalysisError):
    def __init__(self,details):self.details=dict(details);super().__init__(f"Target outside background: {len(self.details.get('entities_outside_background',[]))} protein(s) require a user decision.")


@dataclass(frozen=True)
class InterProPfamParameters:
    annotation_set_id: str|None=None
    target_selection: str="All mapped proteins"
    background_selection: str="All mapped proteins"
    manual_rows: tuple[int,...]=()
    minimum_overlap: int=3
    fdr_cutoff: float=.05
    top_n: int=20
    allow_target_outside_background: bool=False


@dataclass(frozen=True)
class InterProPfamOutputs:
    mapping:pd.DataFrame;ambiguous:pd.DataFrame;unmapped:pd.DataFrame;coverage:pd.DataFrame
    interpro_frequency:pd.DataFrame;interpro_by_type:pd.DataFrame;pfam_frequency:pd.DataFrame
    interpro_enrichment:pd.DataFrame;interpro_significant:pd.DataFrame;interpro_excluded:pd.DataFrame
    pfam_enrichment:pd.DataFrame;pfam_significant:pd.DataFrame;pfam_excluded:pd.DataFrame
    repeated_features:pd.DataFrame;location_summary:pd.DataFrame;interpro_locations:pd.DataFrame;pfam_locations:pd.DataFrame
    pfam_architecture:pd.DataFrame;interpro_architecture:pd.DataFrame;pfam_architecture_frequency:pd.DataFrame;interpro_architecture_frequency:pd.DataFrame;integration:pd.DataFrame
    summary:pd.DataFrame;metadata:dict[str,Any];plots:tuple[Path,...];workbook:Path;run_root:Path;provenance_root:Path


def _catalog(project):return pd.read_csv(project.root/"mapping/tables/protein_catalog.csv")
def _tokens(value):
    if pd.isna(value):return []
    import re
    return [item for item in re.split(r"[;,|\s]+",str(value).strip()) if item and item.upper()!="NA"]
def _resolved_accessions(frame):
    result=[]
    if "source_row" not in frame or "uniprot_accession" not in frame:return result
    for _,group in frame.groupby("source_row",sort=False):
        accessions=sorted(set(item for value in group.uniprot_accession for item in _tokens(value)))
        if len(accessions)==1:result.append(accessions[0])
    return sorted(set(result))


def select_interpro_pfam_set(project,selection,manual_rows=()):
    catalog=_catalog(project)
    if selection in {"All mapped proteins","all_mapped"}:return catalog
    if selection in {"Manual selection","manual"}:return catalog[catalog.source_row.isin(set(manual_rows))]
    path=project.root/"analyses/presence_absence/tables/classification.csv"
    if not path.is_file():raise InterProPfamPresenceRequiredError("Presence/Absence outputs are required for the selected protein set.")
    classes=pd.read_csv(path);label=selection[6:] if selection.startswith("class:") else selection
    rows=classes.loc[classes.classification!="Not reproducibly detected","source_row"] if selection=="Reproducibly detected" else classes.loc[classes.classification==label,"source_row"]
    return catalog[catalog.source_row.isin(rows)]


def prepare_interpro_pfam_arguments(project:Project,manager:DatabaseManager,*,run_id:str,parameters:InterProPfamParameters|None=None):
    params=parameters or InterProPfamParameters();database=manager.interpro
    if not database.metadata_manifest():raise AnnotationMetadataUnavailableError("Annotation metadata unavailable.")
    root=database._set_root(project.root,params.annotation_set_id)
    if root is None or not root.is_dir():raise AnnotationSetUnavailableError("Annotation set unavailable.")
    try:annotation=database.load_annotation_set(project.root,params.annotation_set_id)
    except InvalidAnnotationSet as exc:raise AnnotationSetIncompleteError(f"Annotation set incomplete: {exc}") from exc
    if params.minimum_overlap<1:raise InterProPfamAnalysisError("Minimum overlap must be at least 1.")
    if not 0<=params.fdr_cutoff<=1:raise InterProPfamAnalysisError("FDR threshold must be between 0 and 1.")
    if params.top_n<1:raise InterProPfamAnalysisError("Top N must be at least 1.")
    output=project.root/"analyses/InterPro_Pfam";raw=output/"raw"/run_id
    if (output/"runs"/run_id).exists() or (project.root/"scripts/runs"/f"{run_id}_interpro_pfam").exists():raise InterProPfamAnalysisError("This InterPro/Pfam run ID already exists.")
    raw.mkdir(parents=True,exist_ok=False);target_catalog=select_interpro_pfam_set(project,params.target_selection,params.manual_rows);background_catalog=select_interpro_pfam_set(project,params.background_selection,params.manual_rows)
    background_accessions=_resolved_accessions(background_catalog)
    if not background_accessions:raise NoUniquelyResolvedProteinsError("No uniquely resolved proteins are available in the selected background.")
    covered=set(annotation.proteins.requested_accession.astype(str));missing=sorted(set(background_accessions)-covered)
    if missing:raise AnnotationSetCoverageError("Annotation set does not cover all required proteins: "+", ".join(missing))
    target_catalog.to_csv(raw/"target_catalog.csv",index=False);background_catalog.to_csv(raw/"background_catalog.csv",index=False)
    tables=root/"tables";manifest=annotation.manifest
    paths={"proteins":"proteins.csv","interpro-membership":"protein_interpro_membership.csv","interpro-locations":"interpro_locations.csv","pfam-membership":"protein_pfam_membership.csv","pfam-locations":"pfam_locations.csv","interpro-entries":"interpro_entries.csv","pfam-entries":"pfam_entries.csv"}
    arguments=["--catalog",str(project.root/"mapping/tables/protein_catalog.csv"),"--target-file",str(raw/"target_catalog.csv"),"--background-file",str(raw/"background_catalog.csv")]
    for flag,name in paths.items():arguments.extend([f"--{flag}",str(tables/name)])
    arguments.extend(["--annotation-manifest",str(root/"manifest.json"),"--annotation-hashes",json.dumps(manifest.get("table_hashes",{}),sort_keys=True),"--annotation-set-id",str(manifest["annotation_set_id"]),"--interpro-release",str(manifest.get("interpro_release","unknown")),"--pfam-release",str(manifest.get("pfam_release","unknown")),"--output",str(output),"--project",str(project.root),"--run-id",run_id,"--target-name",params.target_selection,"--background-name",params.background_selection,"--minimum-overlap",str(params.minimum_overlap),"--fdr-cutoff",str(params.fdr_cutoff),"--top-n",str(params.top_n),"--allow-target-adjustment",str(params.allow_target_outside_background).lower()])
    return arguments


def read_interpro_pfam_outputs(project:Project,run_id:str|None=None):
    base=project.root/"analyses/InterPro_Pfam";root=base/"runs"/run_id if run_id else base;metadata_path=root/"metadata.json" if run_id else base/"latest_metadata.json"
    files={"mapping":"mapping/protein_mapping.csv","ambiguous":"mapping/ambiguous.csv","unmapped":"mapping/unmapped.csv","coverage":"mapping/annotation_coverage.csv","interpro_frequency":"frequency/interpro_frequency.csv","interpro_by_type":"frequency/interpro_frequency_by_type.csv","pfam_frequency":"frequency/pfam_frequency.csv","interpro_enrichment":"enrichment/interpro_enrichment_all.csv","interpro_significant":"enrichment/interpro_enrichment_significant.csv","interpro_excluded":"enrichment/interpro_enrichment_excluded.csv","pfam_enrichment":"enrichment/pfam_enrichment_all.csv","pfam_significant":"enrichment/pfam_enrichment_significant.csv","pfam_excluded":"enrichment/pfam_enrichment_excluded.csv","repeated_features":"architecture/repeated_features.csv","location_summary":"architecture/location_summary.csv","interpro_locations":"architecture/interpro_locations.csv","pfam_locations":"architecture/pfam_locations.csv","pfam_architecture":"architecture/pfam_architecture.csv","interpro_architecture":"architecture/interpro_domain_architecture.csv","pfam_architecture_frequency":"architecture/pfam_architecture_frequency.csv","interpro_architecture_frequency":"architecture/interpro_architecture_frequency.csv","integration":"architecture/pfam_interpro_integration.csv","summary":"summary.csv"}
    mandatory={name:root/path for name,path in files.items()};mandatory.update(metadata=metadata_path,workbook=root/"InterPro_Pfam_analysis.xlsx");missing=[str(path) for path in mandatory.values() if not path.is_file()]
    if missing:raise MissingInterProPfamOutputError("Missing expected output: "+", ".join(missing))
    try:
        frames={name:pd.read_csv(root/path) for name,path in files.items()};metadata=json.loads(metadata_path.read_text(encoding="utf-8"));return InterProPfamOutputs(*[frames[name] for name in files],metadata,tuple(sorted((root/"plots").glob("*.*"))),mandatory["workbook"],root,project.root/"scripts/runs"/f"{metadata['run_id']}_interpro_pfam")
    except (OSError,ValueError,KeyError,json.JSONDecodeError) as exc:raise MissingInterProPfamOutputError(f"Invalid InterPro/Pfam outputs: {exc}") from exc


def run_interpro_pfam_analysis(project,manager,runtime:RRuntime,*,run_id,parameters=None,script=None,timeout=180):
    arguments=prepare_interpro_pfam_arguments(project,manager,run_id=run_id,parameters=parameters);entry=script or Path(__file__).resolve().parents[3]/"r_scripts/07_interpro_pfam_analysis.R";result:RResult=runtime.run(entry,*arguments,timeout=timeout)
    if result.returncode!=0:
        error=project.root/"analyses/InterPro_Pfam/target_background_error.json"
        if error.is_file():
            try:raise InterProPfamTargetOutsideBackgroundError(json.loads(error.read_text(encoding="utf-8")))
            except json.JSONDecodeError:pass
        message=(result.stderr or result.stdout or "Rscript exited with an error.").strip()
        if "No uniquely resolved proteins" in message:raise NoUniquelyResolvedProteinsError(message)
        raise InterProPfamRExecutionError(f"R execution failure: {message}")
    return read_interpro_pfam_outputs(project,run_id)


def list_interpro_pfam_runs(project):
    root=project.root/"analyses/InterPro_Pfam/runs";return sorted(path.name for path in root.iterdir() if path.is_dir()) if root.is_dir() else []


def interpro_pfam_experiment_accessions(project):return _resolved_accessions(_catalog(project))
def _protein_field(value):return set() if pd.isna(value) else {item.strip() for item in str(value).split(";") if item.strip()}
def interpro_for_protein(outputs,accession):return outputs.interpro_frequency[outputs.interpro_frequency.Proteins.map(lambda value:str(accession) in _protein_field(value))].copy()
def pfam_for_protein(outputs,accession):return outputs.pfam_frequency[outputs.pfam_frequency.Proteins.map(lambda value:str(accession) in _protein_field(value))].copy()
def proteins_for_interpro(outputs,feature):
    rows=outputs.interpro_frequency[outputs.interpro_frequency.InterPro_ID.astype(str)==str(feature)];values=set().union(*(_protein_field(value) for value in rows.Proteins)) if not rows.empty else set();return pd.DataFrame({"UniProt":sorted(values)})
def proteins_for_pfam(outputs,feature):
    rows=outputs.pfam_frequency[outputs.pfam_frequency.Pfam_ID.astype(str)==str(feature)];values=set().union(*(_protein_field(value) for value in rows.Proteins)) if not rows.empty else set();return pd.DataFrame({"UniProt":sorted(values)})
def interpro_locations_for_protein(outputs,accession):return outputs.interpro_locations[outputs.interpro_locations.uniprot_accession.astype(str)==str(accession)].copy()
def pfam_locations_for_protein(outputs,accession):return outputs.pfam_locations[outputs.pfam_locations.uniprot_accession.astype(str)==str(accession)].copy()
def architectures_for_protein(outputs,accession):
    pfam=outputs.pfam_architecture[outputs.pfam_architecture.UniProt.astype(str)==str(accession)].copy();interpro=outputs.interpro_architecture[outputs.interpro_architecture.UniProt.astype(str)==str(accession)].copy();return pfam,interpro
