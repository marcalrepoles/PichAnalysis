import json

import pandas as pd

from pichanalysis.core.presence_analysis import (
    build_presence_arguments, export_presence, list_presence_runs, presence_readiness,
    quantification_types, quantitative_columns, read_presence_outputs, replicate_counts,
    suggested_minimum,
)
from pichanalysis.core.project import create_project


def project_with_quantification(tmp_path):
    project=create_project(tmp_path,"Presence")
    data=project.root/"input"/"processed"/"data.csv"
    data.write_text("ID,A1,A2,A3,B1,B2\nP1,1,1,0,1,0\n",encoding="utf-8")
    project.config["input"].update({"processed_file":"input/processed/data.csv","rows":1})
    project.config["columns"]={
        "ID":{"role":"identifier","identifier_type":"uniprot","primary_identifier":True},
        "A1":{"role":"quantification","condition":"A","replicate":"1","quantification_type":"lfq_intensity"},
        "A2":{"role":"quantification","condition":"A","replicate":"2","quantification_type":"lfq_intensity"},
        "A3":{"role":"quantification","condition":"A","replicate":"3","quantification_type":"lfq_intensity"},
        "B1":{"role":"quantification","condition":"B","replicate":"1","quantification_type":"raw_intensity"},
        "B2":{"role":"quantification","condition":"B","replicate":"2","quantification_type":"raw_intensity"},
    }
    project.save(); return project


def test_reads_quantitative_configuration(tmp_path):
    project=project_with_quantification(tmp_path)
    assert len(quantitative_columns(project))==5
    assert replicate_counts(project,"lfq_intensity")=={"A":3}


def test_quantification_type_selection_and_default(tmp_path):
    project=project_with_quantification(tmp_path)
    assert quantification_types(project)==["lfq_intensity","raw_intensity"]
    assert suggested_minimum(project,"lfq_intensity")==2


def test_validation_requires_conditions(tmp_path):
    project=project_with_quantification(tmp_path)
    project.config["columns"]["A1"]["condition"]=""
    assert not presence_readiness(project,"lfq_intensity").ready


def test_builds_safe_r_arguments(tmp_path):
    project=project_with_quantification(tmp_path)
    args=build_presence_arguments(project,quantification_type="lfq_intensity",selected_conditions=["A"],
        zero_is_missing=True,threshold=100,rule_mode="count",rule_value=2,predominant=True,run_id="run1")
    mapping=json.loads(args[args.index("--column-map")+1])
    assert [item["column"] for item in mapping]==["A1","A2","A3"]
    assert args[args.index("--threshold")+1]=="100.0"


def test_reads_outputs(tmp_path):
    project=project_with_quantification(tmp_path); root=project.root/"analyses"/"presence_absence"
    (root/"tables").mkdir(parents=True); (root/"graphs").mkdir()
    frame=pd.DataFrame({"source_row":[1],"original_id":["P1"]})
    for name in ("presence_matrix","condition_detection","classification","condition_summary"): frame.to_csv(root/"tables"/f"{name}.csv",index=False)
    (root/"latest_metadata.json").write_text(json.dumps({"total_entities":1}),encoding="utf-8")
    (root/"graphs"/"plot.png").write_bytes(b"png")
    outputs=read_presence_outputs(project)
    assert outputs.metadata["total_entities"]==1 and len(outputs.graphs)==1


def test_run_history(tmp_path):
    project=project_with_quantification(tmp_path); runs=project.root/"analyses"/"presence_absence"/"runs"
    (runs/"run2").mkdir(parents=True); (runs/"run1").mkdir()
    assert list_presence_runs(project)==["run1","run2"]


def test_export_preserves_source(tmp_path):
    source=tmp_path/"source.csv"; source.write_bytes(b"a\n1\n")
    destination=tmp_path/"export"/"copy.csv"; before=source.read_bytes()
    export_presence(source,destination)
    assert source.read_bytes()==before==destination.read_bytes()


def test_interface_disabled_without_valid_configuration(tmp_path, monkeypatch):
    monkeypatch.setenv("QT_QPA_PLATFORM","offscreen")
    from PySide6.QtWidgets import QApplication
    from pichanalysis.ui.presence_page import PresencePage
    app=QApplication.instance() or QApplication([])
    page=PresencePage(); page.set_project(create_project(tmp_path,"Empty"))
    assert not page.run_button.isEnabled()
