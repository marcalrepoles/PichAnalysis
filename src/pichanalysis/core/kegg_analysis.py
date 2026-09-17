from __future__ import annotations
import hashlib,json,shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any
import pandas as pd
from .database_manager import DatabaseManager
from .go_analysis import available_sets,select_set
from .organism import get_organism
from .project import Project

@dataclass(frozen=True)
class KEGGReadiness: ready:bool; reason:str
@dataclass(frozen=True)
class KEGGOutputs:
    mapping:pd.DataFrame;frequency:pd.DataFrame;enrichment:pd.DataFrame;summary:pd.DataFrame
    metadata:dict[str,Any];tables:tuple[Path,...];graphs:tuple[Path,...];root:Path

def kegg_readiness(project:Project|None,manager:DatabaseManager)->KEGGReadiness:
    if project is None:return KEGGReadiness(False,"Open a project.")
    organism=get_organism(project)
    if organism is None:return KEGGReadiness(False,"Configure the project organism.")
    if organism.tax_id!="9606":return KEGGReadiness(False,"KEGG pathway analysis is currently available only for Homo sapiens.")
    compatible,reason=manager.pathway_analysis_compatibility()
    if not compatible:return KEGGReadiness(False,"This KEGG snapshot does not contain the identifier mapping tables required for pathway analysis." if reason=="Missing identifier mapping tables" else "KEGG database is not ready for pathway analysis.")
    if not (project.root/"mapping/tables/protein_catalog.csv").is_file():return KEGGReadiness(False,"Run Identification and Annotation first to resolve stable gene identifiers.")
    return KEGGReadiness(True,"Ready for KEGG pathway analysis.")

def prepare_kegg_arguments(project:Project,manager:DatabaseManager,*,target_selection:str,background_selection:str,
 manual_rows:list[int]|None,mapping_policy:str,fdr_cutoff:float,p_cutoff:float,min_count:int,top_n:int,run_id:str,
 allow_target_outside_background:bool=False)->list[str]:
    state=kegg_readiness(project,manager)
    if not state.ready:raise ValueError(state.reason)
    if mapping_policy not in {"unique-only","all-candidates"}:raise ValueError("Invalid KEGG mapping policy.")
    target=select_set(project,target_selection,manual_rows);background=select_set(project,background_selection,manual_rows)
    raw=project.root/"analyses/KEGG/raw"/run_id
    if raw.exists():raise ValueError("This KEGG run ID already exists.")
    raw.mkdir(parents=True);target_file=raw/"target_entities.csv";background_file=raw/"background_entities.csv"
    target.to_csv(target_file,index=False);background.to_csv(background_file,index=False)
    snapshot=manager.active_snapshot();assert snapshot
    manifest=manager.manifest(snapshot);manifest_path=snapshot/"manifest.json"
    organism=get_organism(project);assert organism
    return ["--catalog",str(project.root/"mapping/tables/protein_catalog.csv"),"--target-file",str(target_file),"--background-file",str(background_file),
      "--output",str(project.root/"analyses/KEGG"),"--project",str(project.root),"--pathways",str(snapshot/"tables/pathways.tsv"),
      "--genes",str(snapshot/"tables/genes.tsv"),"--ncbi-map",str(snapshot/"tables/ncbi_geneid_to_kegg.tsv"),"--uniprot-map",str(snapshot/"tables/uniprot_to_kegg.tsv"),
      "--gene-pathway",str(snapshot/"tables/gene_to_pathway.tsv"),"--pathway-gene",str(snapshot/"tables/pathway_to_gene.tsv"),"--kgml-dir",str(snapshot/"kgml"),
      "--mapping-policy",mapping_policy,"--fdr-cutoff",str(fdr_cutoff),"--p-cutoff",str(p_cutoff),"--min-count",str(min_count),"--top-n",str(top_n),
      "--run-id",run_id,"--snapshot-id",snapshot.name,"--snapshot-manifest",str(manifest_path),"--snapshot-sha256",hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
      "--database-schema-version",str(manifest.get("database_schema_version",1)),"--organism",organism.name,"--organism-code","hsa","--tax-id","9606",
      "--target-name",target_selection,"--background-name",background_selection,"--allow-target-adjustment",str(allow_target_outside_background).lower()]

def read_kegg_outputs(project:Project,run_id:str|None=None)->KEGGOutputs:
    base=project.root/"analyses/KEGG";root=base/"runs"/run_id if run_id else base
    try:return KEGGOutputs(pd.read_csv(root/"mapping/kegg_mapping.csv"),pd.read_csv(root/"frequency/kegg_pathway_frequency.csv"),pd.read_csv(root/"enrichment/kegg_pathway_all.csv"),pd.read_csv(root/"summary.csv"),json.loads((root/("metadata.json" if run_id else "latest_metadata.json")).read_text(encoding="utf-8")),tuple(sorted(root.glob("**/*.csv"))),tuple(sorted((root/"graphs").glob("*.png"))),root)
    except (OSError,ValueError,json.JSONDecodeError) as error:raise RuntimeError(f"Invalid KEGG analysis results: {error}") from error
def list_kegg_runs(project:Project)->list[str]:
    root=project.root/"analyses/KEGG/runs";return sorted(x.name for x in root.iterdir() if x.is_dir()) if root.is_dir() else []
def export_kegg(source:Path,destination:Path)->Path:destination.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,destination);return destination
