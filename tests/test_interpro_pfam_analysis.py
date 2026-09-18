from __future__ import annotations

import json
from pathlib import Path
import pandas as pd
import pytest

from pichanalysis.core.database_manager import DatabaseManager
from pichanalysis.core.interpro_pfam_analysis import (
    AnnotationSetCoverageError, AnnotationSetUnavailableError, InterProPfamParameters,
    InterProPfamRExecutionError, InterProPfamTargetOutsideBackgroundError,
    MissingInterProPfamOutputError, list_interpro_pfam_runs, prepare_interpro_pfam_arguments,
    read_interpro_pfam_outputs, run_interpro_pfam_analysis,
)
from pichanalysis.core.project import create_project
from pichanalysis.core.r_runtime import RResult, RRuntime
from test_interpro_database import FakeProvider


def fixture(tmp_path):
    project=create_project(tmp_path/"projects","InterProProject");catalog=pd.DataFrame([
        {"source_row":1,"original_id":"P1","uniprot_accession":"P1"},{"source_row":2,"original_id":"P1 duplicate","uniprot_accession":"P1"},
        {"source_row":3,"original_id":"P2","uniprot_accession":"P2"},{"source_row":4,"original_id":"empty","uniprot_accession":"EMPTY"},
        {"source_row":5,"original_id":"group","uniprot_accession":"P1;P2"},{"source_row":6,"original_id":"unknown","uniprot_accession":""},
    ]);catalog.to_csv(project.root/"mapping/tables/protein_catalog.csv",index=False);presence=project.root/"analyses/presence_absence/tables";presence.mkdir(parents=True);pd.DataFrame({"source_row":range(1,7),"classification":["Shared","Shared","A-specific","Shared","Shared","Sporadic"]}).to_csv(presence/"classification.csv",index=False)
    manager=DatabaseManager(tmp_path/"db");manager.interpro.provider=FakeProvider();manager.interpro.download_metadata();manager.interpro.build_annotation_set(project.root,["P1","P2","EMPTY"],annotation_set_id="setA");return project,manager

def test_parameters_validation_coverage_and_structured_errors(tmp_path):
    project,manager=fixture(tmp_path);args=prepare_interpro_pfam_arguments(project,manager,run_id="params",parameters=InterProPfamParameters(annotation_set_id="setA",target_selection="Shared",minimum_overlap=1));assert args[args.index("--annotation-set-id")+1]=="setA" and args[args.index("--minimum-overlap")+1]=="1"
    with pytest.raises(AnnotationSetUnavailableError):prepare_interpro_pfam_arguments(project,manager,run_id="absent",parameters=InterProPfamParameters(annotation_set_id="missing"))
    catalog=pd.read_csv(project.root/"mapping/tables/protein_catalog.csv");catalog.loc[len(catalog)]={"source_row":7,"original_id":"outside","uniprot_accession":"P9"};catalog.to_csv(project.root/"mapping/tables/protein_catalog.csv",index=False)
    with pytest.raises(AnnotationSetCoverageError,match="does not cover"):prepare_interpro_pfam_arguments(project,manager,run_id="coverage",parameters=InterProPfamParameters(annotation_set_id="setA"))

def test_offline_scientific_smoke_runs_provenance_and_history(tmp_path):
    project,manager=fixture(tmp_path);manager.interpro.provider=object();runtime=RRuntime();assert runtime.available
    a=run_interpro_pfam_analysis(project,manager,runtime,run_id="A",parameters=InterProPfamParameters(annotation_set_id="setA",target_selection="Shared",minimum_overlap=1,top_n=10),timeout=180)
    assert len(a.mapping)>0 and len(a.coverage)==4 and len(a.interpro_frequency)>0 and len(a.pfam_frequency)>0;assert len(a.plots)==12 and a.workbook.is_file()
    repeated=a.repeated_features.query("Source_database == 'Pfam' and Feature_ID == 'PF00001'");assert repeated.iloc[0].Occurrence_count==2
    p1=a.pfam_architecture.query("UniProt == 'P1'").iloc[0];assert p1.Architecture.count("PF00001")==2 and bool(p1.Has_repeated_feature) and bool(p1.Has_overlap)
    assert len(a.pfam_locations.query("uniprot_accession == 'P1' and pfam_id == 'PF00001'"))==2
    assert len(a.pfam_frequency.query("Pfam_ID == 'PF00001'").iloc[0].Proteins.split(";"))>=1
    assert a.metadata["network_access"] is False and a.metadata["annotation_set_id"]=="setA";assert (a.provenance_root/"07_interpro_pfam_analysis.R").is_file() and (a.provenance_root/"annotation_set_manifest.json").is_file()
    b=run_interpro_pfam_analysis(project,manager,runtime,run_id="B",parameters=InterProPfamParameters(annotation_set_id="setA",target_selection="A-specific",minimum_overlap=1),timeout=180)
    again=read_interpro_pfam_outputs(project,"A");pd.testing.assert_frame_equal(a.summary,again.summary);assert b.metadata["run_id"]=="B" and list_interpro_pfam_runs(project)==["A","B"]

def test_target_safeguard_r_failure_and_missing_output(tmp_path):
    project,manager=fixture(tmp_path);runtime=RRuntime()
    with pytest.raises(InterProPfamTargetOutsideBackgroundError):run_interpro_pfam_analysis(project,manager,runtime,run_id="blocked",parameters=InterProPfamParameters(annotation_set_id="setA",target_selection="All mapped proteins",background_selection="A-specific",minimum_overlap=1))
    adjusted=run_interpro_pfam_analysis(project,manager,runtime,run_id="adjusted",parameters=InterProPfamParameters(annotation_set_id="setA",target_selection="All mapped proteins",background_selection="A-specific",minimum_overlap=1,allow_target_outside_background=True));assert adjusted.metadata["target_background_adjustment"]["user_decision"]=="Continue"
    class Failed:
        def run(self,*a,**k):return RResult(tuple(),1,"","planned R failure")
    with pytest.raises(InterProPfamRExecutionError,match="planned R failure"):run_interpro_pfam_analysis(project,manager,Failed(),run_id="failed",parameters=InterProPfamParameters(annotation_set_id="setA"))
    with pytest.raises(MissingInterProPfamOutputError):read_interpro_pfam_outputs(project,"missing")
