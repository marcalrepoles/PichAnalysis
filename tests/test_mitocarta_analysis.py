from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from pichanalysis.core.database_manager import DatabaseManager
from pichanalysis.core.database_registry import DatabaseState
from pichanalysis.core.mitocarta_analysis import (
    MitoCartaDatabaseUnavailableError, MitoCartaParameters, MitoCartaRExecutionError,
    MitoCartaTargetOutsideBackgroundError, PresenceAbsenceRequiredError,
    UnsupportedMitoCartaOrganismError, list_mitocarta_runs, prepare_mitocarta_arguments,
    read_mitocarta_outputs, run_mitocarta_analysis, select_mitocarta_set,
)
from pichanalysis.core.organism import set_organism
from pichanalysis.core.project import create_project
from pichanalysis.core.r_runtime import RResult, RRuntime


def project_fixture(root: Path, organism="Homo sapiens", tax_id="9606"):
    project=create_project(root,"MitoProject");set_organism(project,organism,tax_id)
    catalog=pd.DataFrame([
        {"source_row":1,"original_id":"A","ncbi_gene_id":"1","gene_symbol":"A","uniprot_accession":"U1","protein_name":"Protein A"},
        {"source_row":2,"original_id":"B","ncbi_gene_id":"","gene_symbol":"B","uniprot_accession":"U2","protein_name":"Protein B"},
        {"source_row":3,"original_id":"C","ncbi_gene_id":"","gene_symbol":"","uniprot_accession":"U3","protein_name":"Protein C"},
        {"source_row":4,"original_id":"D","ncbi_gene_id":"4","gene_symbol":"D","uniprot_accession":"U4","protein_name":"Protein D"},
        {"source_row":5,"original_id":"E","ncbi_gene_id":"5","gene_symbol":"E","uniprot_accession":"U5","protein_name":"Protein E"},
        {"source_row":6,"original_id":"U1;U2","ncbi_gene_id":"","gene_symbol":"","uniprot_accession":"U1;U2","protein_name":"Ambiguous group"},
        {"source_row":7,"original_id":"unknown","ncbi_gene_id":"999","gene_symbol":"Z","uniprot_accession":"UX","protein_name":"Unknown"},
        {"source_row":8,"original_id":"A duplicate","ncbi_gene_id":"1","gene_symbol":"A","uniprot_accession":"U1","protein_name":"Protein A duplicate"},
    ]);catalog.to_csv(project.root/"mapping/tables/protein_catalog.csv",index=False)
    tables=project.root/"analyses/presence_absence/tables";tables.mkdir(parents=True);pd.DataFrame({"source_row":range(1,9),"classification":["Shared","Shared","A-specific","Sporadic","A-specific","Shared","Not reproducibly detected","Shared"]}).to_csv(tables/"classification.csv",index=False)
    return project


def manager_fixture(root: Path):
    manager=DatabaseManager(root/"db");snap=manager.mitocarta.create_staging_snapshot();tables=snap/"tables"
    pd.DataFrame([
        {"gene_symbol":"A","gene_description":"A","ncbi_gene_id":"1","uniprot_accession":"U1","synonyms":"","maestro_score":10,"evidence":"literature","sub_compartment_raw":"matrix","mitopathways_raw":"Root > P1","tissues_raw":"all","is_mitocarta":True},
        {"gene_symbol":"B","gene_description":"B","ncbi_gene_id":"2","uniprot_accession":"U2","synonyms":"","maestro_score":9,"evidence":"MS/MS","sub_compartment_raw":"MIM|IMS","mitopathways_raw":"Root > P1; Root > P2","tissues_raw":"heart","is_mitocarta":True},
        {"gene_symbol":"C","gene_description":"C","ncbi_gene_id":"3","uniprot_accession":"U3","synonyms":"","maestro_score":1,"evidence":"","sub_compartment_raw":"","mitopathways_raw":"","tissues_raw":"","is_mitocarta":False},
        {"gene_symbol":"D","gene_description":"D","ncbi_gene_id":"4","uniprot_accession":"U4","synonyms":"","maestro_score":8,"evidence":"","sub_compartment_raw":"matrix","mitopathways_raw":"Root > P2","tissues_raw":"","is_mitocarta":True},
        {"gene_symbol":"E","gene_description":"E","ncbi_gene_id":"5","uniprot_accession":"U5","synonyms":"","maestro_score":0,"evidence":"","sub_compartment_raw":"","mitopathways_raw":"","tissues_raw":"","is_mitocarta":False},
    ]).to_csv(tables/"mitocarta_genes.csv",index=False)
    pd.DataFrame({"gene_symbol":["A","B","B","D"],"subcompartment":["matrix","MIM","IMS","matrix"]}).to_csv(tables/"subcompartment_membership.csv",index=False)
    pd.DataFrame({"mitopathway":["P1","P1","P2","P2"],"gene_symbol":["A","B","B","D"]}).to_csv(tables/"mitopathway_membership.csv",index=False)
    pd.DataFrame({"mitopathway":["P1","P2"],"pathway_full":["Root > P1","Root > P2"],"pathway_level_1":["Root","Root"],"pathway_level_2":["P1","P2"]}).to_csv(tables/"mitopathway_hierarchy.csv",index=False)
    manifest=manager.mitocarta.manifest(snap);manifest.update(status=DatabaseState.READY,version="3.0",files=[{"filename":"Human.MitoCarta3.0.xls","sha256":"a"},{"filename":"Human.MitoPathways3.0.gmx","sha256":"b"}]);manager.mitocarta._write_manifest(snap,manifest);manager.mitocarta._atomic_json(manager.mitocarta.active_pointer,{"snapshot_id":snap.name});return manager,snap


def test_parameters_active_snapshot_and_target_serialization(tmp_path):
    project=project_fixture(tmp_path/"projects");manager,snapshot=manager_fixture(tmp_path);params=MitoCartaParameters(target_selection="Manual selection",manual_rows=(1,2),minimum_overlap=1)
    args=prepare_mitocarta_arguments(project,manager,run_id="parameters",parameters=params)
    assert args[args.index("--snapshot-id")+1]==snapshot.name and args[args.index("--minimum-overlap")+1]=="1"
    assert len(pd.read_csv(args[args.index("--target-file")+1]))==2 and json.loads(args[args.index("--snapshot-hashes")+1])=={"Human.MitoCarta3.0.xls":"a","Human.MitoPathways3.0.gmx":"b"}


def test_missing_database_unsupported_organism_and_presence_requirement(tmp_path):
    project=project_fixture(tmp_path/"projects");manager=DatabaseManager(tmp_path/"empty")
    with pytest.raises(MitoCartaDatabaseUnavailableError):prepare_mitocarta_arguments(project,manager,run_id="missing")
    mouse=project_fixture(tmp_path/"mouse",organism="Mus musculus",tax_id="10090");ready,_=manager_fixture(tmp_path/"ready")
    with pytest.raises(UnsupportedMitoCartaOrganismError):prepare_mitocarta_arguments(mouse,ready,run_id="mouse")
    (project.root/"analyses/presence_absence/tables/classification.csv").unlink()
    with pytest.raises(PresenceAbsenceRequiredError):select_mitocarta_set(project,"Shared")


def test_structured_r_failure_and_target_background_error(tmp_path):
    project=project_fixture(tmp_path/"projects");manager,_=manager_fixture(tmp_path)
    class Failed:
        def run(self,*args,**kwargs):return RResult(tuple(),1,"","planned R failure")
    with pytest.raises(MitoCartaRExecutionError,match="planned R failure"):run_mitocarta_analysis(project,manager,Failed(),run_id="failed")
    error_path=project.root/"analyses/MitoCarta/target_background_error.json";error_path.parent.mkdir(parents=True,exist_ok=True);error_path.write_text(json.dumps({"entities_outside_background":["NCBI:1"],"initial_target_size":1,"initial_background_size":0}))
    with pytest.raises(MitoCartaTargetOutsideBackgroundError) as caught:run_mitocarta_analysis(project,manager,Failed(),run_id="outside")
    assert caught.value.details["initial_target_size"]==1


def test_python_to_r_offline_scientific_smoke_persistence_and_provenance(tmp_path):
    project=project_fixture(tmp_path/"projects");manager,snapshot=manager_fixture(tmp_path);runtime=RRuntime();assert runtime.available
    outputs=run_mitocarta_analysis(project,manager,runtime,run_id="smoke",parameters=MitoCartaParameters(target_selection="Shared",minimum_overlap=1,top_n=10),timeout=180)
    statuses=outputs.mapping.drop_duplicates("source_row").mapping_status.value_counts().to_dict();assert statuses=={"mapped_unique":6,"ambiguous":1,"unmapped":1}
    assert len(outputs.membership)==5 and len(outputs.overall)==1 and len(outputs.subcompartment_frequency)>=2 and len(outputs.pathway_frequency)==2
    assert outputs.workbook.is_file() and len(outputs.graphs)==10 and all(path.stat().st_size for path in outputs.graphs)
    assert outputs.metadata["snapshot_id"]==snapshot.name and outputs.metadata["network_access"] is False
    assert (outputs.provenance_root/"06_mitocarta_analysis.R").is_file() and (outputs.provenance_root/"lib/mitocarta_analysis.R").is_file() and (outputs.provenance_root/"database_hashes.json").is_file()
    assert list_mitocarta_runs(project)==["smoke"] and read_mitocarta_outputs(project,"smoke").metadata["run_id"]=="smoke"
    summary=dict(zip(outputs.summary.metric,outputs.summary.value.astype(str)));assert summary["MitoCarta version"]=="3.0" and summary["FDR method"]=="Benjamini-Hochberg"
    print("MITOCARTA_SCIENTIFIC_SMOKE "+json.dumps({"mapping_status":statuses,"membership_genes":len(outputs.membership),"overall_2x2":outputs.overall.iloc[0][["target_mitocarta","target_non_mitocarta","reference_mitocarta","reference_non_mitocarta"]].to_dict(),"odds_ratio":outputs.overall.iloc[0].odds_ratio,"ci":[outputs.overall.iloc[0].ci_low,outputs.overall.iloc[0].ci_high],"p_value":outputs.overall.iloc[0].p_value,"subcompartment_frequency_rows":len(outputs.subcompartment_frequency),"subcompartment_tests":len(outputs.subcompartment_enrichment),"mitopathway_frequency_rows":len(outputs.pathway_frequency),"mitopathway_tests":len(outputs.pathway_enrichment),"plots":len(outputs.graphs),"workbook":str(outputs.workbook),"snapshot_id":snapshot.name},default=str))


def test_output_validation(tmp_path):
    project=project_fixture(tmp_path/"projects")
    from pichanalysis.core.mitocarta_analysis import MissingMitoCartaOutputError
    with pytest.raises(MissingMitoCartaOutputError):read_mitocarta_outputs(project,"absent")
