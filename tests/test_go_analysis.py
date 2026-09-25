import json

import pandas as pd
import pytest

from pichanalysis.core.go_analysis import (
    available_sets, export_go, go_readiness, list_go_runs, prepare_go_arguments,
    read_go_outputs, select_set,
)
from pichanalysis.core.organism import set_organism
from pichanalysis.core.project import create_project


def go_project(tmp_path):
    project=create_project(tmp_path,"GO")
    processed=project.root/"input"/"processed"/"data.csv";processed.write_text("ID\nTFAM\nUNKNOWN\n",encoding="utf-8")
    project.config["input"]["processed_file"]="input/processed/data.csv";project.save();set_organism(project,"Homo sapiens","9606")
    catalog=pd.DataFrame([
        {"source_row":1,"original_id":"TFAM","mapping_status":"mapped_unique","gene_symbol":"TFAM","ncbi_gene_id":"7019","uniprot_accession":"Q00059"},
        {"source_row":2,"original_id":"UNKNOWN","mapping_status":"unmapped","gene_symbol":"XUNKNOWN","ncbi_gene_id":None,"uniprot_accession":None},
    ])
    catalog.to_csv(project.root/"mapping"/"tables"/"protein_catalog.csv",index=False)
    return project


def base_parameters():
    return dict(target_selection="mapped",background_selection="mapped",manual_rows=None,
        ontologies=["BP","MF","CC"],evidence_filter="all",fdr_cutoff=.05,p_cutoff=1,
        min_count=1,top_n=20,simplify=False,simplify_cutoff=.7,run_id="run1")


def test_set_selection(tmp_path):
    project=go_project(tmp_path)
    assert "all_experiment" in available_sets(project)
    assert len(select_set(project,"all_experiment"))==2
    assert len(select_set(project,"mapped"))==1


def test_background_selection(tmp_path):
    project=go_project(tmp_path)
    assert select_set(project,"mapped").iloc[0]["gene_symbol"]=="TFAM"


def test_target_outside_background_is_reported(tmp_path):
    project=go_project(tmp_path);params=base_parameters();params["target_selection"]="all_experiment"
    with pytest.raises(ValueError,match="outside the selected background"):
        prepare_go_arguments(project,**params)


def test_parameter_construction(tmp_path):
    project=go_project(tmp_path);args=prepare_go_arguments(project,**base_parameters())
    assert json.loads(args[args.index("--ontologies")+1])==["BP","MF","CC"]
    assert args[args.index("--background-name")+1]=="mapped"


def test_read_outputs(tmp_path):
    project=go_project(tmp_path);root=project.root/"analyses"/"GO"
    for folder in ("annotation","frequency","enrichment","graphs"):(root/folder).mkdir(parents=True,exist_ok=True)
    pd.DataFrame({"GO_ID":["GO:1"]}).to_csv(root/"annotation"/"go_annotations.csv",index=False)
    pd.DataFrame({"input_id":["X"]}).to_csv(root/"annotation"/"unannotated_proteins.csv",index=False)
    pd.DataFrame({"terms_tested":[1],"significant_terms":[0]}).to_csv(root/"summary.csv",index=False)
    (root/"latest_metadata.json").write_text(json.dumps({"orgdb":"org.Hs.eg.db"}),encoding="utf-8")
    outputs=read_go_outputs(project);assert outputs.metadata["orgdb"]=="org.Hs.eg.db"


def test_history_and_export(tmp_path):
    project=go_project(tmp_path);runs=project.root/"analyses"/"GO"/"runs";(runs/"b").mkdir(parents=True);(runs/"a").mkdir()
    assert list_go_runs(project)==["a","b"]
    source=tmp_path/"x.csv";source.write_bytes(b"x\n1\n");destination=tmp_path/"out"/"x.csv";export_go(source,destination);assert destination.read_bytes()==source.read_bytes()


def test_blocked_without_organism(tmp_path):
    project=create_project(tmp_path,"NoOrg")
    assert not go_readiness(project).ready


def test_blocked_without_mappable_catalog(tmp_path):
    project=create_project(tmp_path,"NoCatalog");set_organism(project,"Homo sapiens","9606")
    assert not go_readiness(project).ready


def test_unsupported_organism(tmp_path):
    project=create_project(tmp_path,"Other");set_organism(project,"Danio rerio","7955")
    assert "not configured" in go_readiness(project).reason
