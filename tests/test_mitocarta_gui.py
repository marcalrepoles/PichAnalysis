from __future__ import annotations

import os
from dataclasses import replace
import pandas as pd
import pytest
from PySide6.QtWidgets import QApplication, QMessageBox

from pichanalysis.core.mitocarta_analysis import (
    MitoCartaParameters, MitoCartaTargetOutsideBackgroundError, list_mitocarta_runs,
    mitocarta_compartment_genes, mitocarta_gene_pathways, mitocarta_gene_subcompartments,
    mitocarta_pathway_genes, read_mitocarta_outputs, run_mitocarta_analysis,
)
from pichanalysis.core.r_runtime import RRuntime
from pichanalysis.ui.main_window import MainWindow
from pichanalysis.ui.mitocarta_page import MitoCartaAnalysisWorker, MitoCartaPage
from test_mitocarta_analysis import manager_fixture, project_fixture

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

@pytest.fixture(scope="module")
def app():
    instance=QApplication.instance() or QApplication([]);yield instance

@pytest.fixture()
def completed_runs(tmp_path):
    project=project_fixture(tmp_path/"projects");manager,_=manager_fixture(tmp_path);runtime=RRuntime();assert runtime.available
    first=run_mitocarta_analysis(project,manager,runtime,run_id="20260101_010101_A",parameters=MitoCartaParameters(target_selection="Shared",minimum_overlap=1),timeout=180)
    second=run_mitocarta_analysis(project,manager,runtime,run_id="20260101_010102_B",parameters=MitoCartaParameters(target_selection="A-specific",minimum_overlap=1),timeout=180)
    return project,manager,first,second

def test_loader_navigation_and_a_b_a(completed_runs):
    project,_manager,first,second=completed_runs
    assert list_mitocarta_runs(project)==["20260101_010101_A","20260101_010102_B"]
    a=read_mitocarta_outputs(project,"20260101_010101_A");b=read_mitocarta_outputs(project,"20260101_010102_B");again=read_mitocarta_outputs(project,"20260101_010101_A")
    for field in ("summary","membership","overall","subcompartment_frequency","subcompartment_enrichment","pathway_frequency","pathway_enrichment","mapping"):pd.testing.assert_frame_equal(getattr(a,field),getattr(again,field))
    assert a.metadata["snapshot_id"]==again.metadata["snapshot_id"]==first.metadata["snapshot_id"]
    assert a.metadata["run_id"]!=b.metadata["run_id"] and not a.summary.equals(b.summary)
    genes=mitocarta_compartment_genes(a,"matrix");assert not genes.empty
    key=str(genes.iloc[0].canonical_gene_key);assert not mitocarta_gene_subcompartments(a,key).empty;assert isinstance(mitocarta_gene_pathways(a,key),pd.DataFrame)
    assert not mitocarta_pathway_genes(a,str(a.pathway_frequency.iloc[0].MitoPathway)).empty
    assert mitocarta_gene_subcompartments(a,"not-a-gene").empty and mitocarta_compartment_genes(a,"not-a-compartment").empty
    assert mitocarta_gene_pathways(a,"not-a-gene").empty and mitocarta_pathway_genes(a,"not-a-pathway").empty

def test_page_readiness_results_filters_navigation_history_and_graphs(app,completed_runs):
    project,manager,first,second=completed_runs;page=MitoCartaPage(manager);page.set_project(project)
    assert "Ready" in page.database.text() and page.run_button.isEnabled();assert page.target.findText("Shared")>=0 and page.background.currentText()=="All gene-resolved entities in the experiment"
    assert page.fdr.value()==pytest.approx(.05) and page.minimum.value()==3 and page.top_n.value()==20
    page.show_outputs(first);assert page.summary.rowCount()>0 and page.membership.rowCount()>0 and page.overall.rowCount()==1
    assert page.sub_frequency.rowCount()>0 and page.path_frequency.rowCount()>0 and page.mapping.rowCount()>0;assert page.graph_choice.count()>0 and not page.graph_preview.pixmap().isNull()
    page.membership_filter.setCurrentIndex(page.membership_filter.findData("mito"));assert page.membership.rowCount()>0
    page.gene_choice.setCurrentIndex(0);assert page.gene_compartments.rowCount()+page.gene_pathways.rowCount()>0
    page.sub_choice.setCurrentIndex(0);assert page.sub_genes.rowCount()>0;page.path_choice.setCurrentIndex(0);assert page.path_genes.rowCount()>0
    ia=page.history.findData(first.metadata["run_id"]);ib=page.history.findData(second.metadata["run_id"])
    page.history.setCurrentIndex(ia);assert page.outputs.metadata["run_id"]==first.metadata["run_id"]
    page.history.setCurrentIndex(ib);assert page.outputs.metadata["run_id"]==second.metadata["run_id"]
    page.history.setCurrentIndex(ia);assert page.outputs.metadata["run_id"]==first.metadata["run_id"];page.close()

def test_page_not_installed_mapping_presence_and_incompatible_organism(app,tmp_path):
    from pichanalysis.core.database_manager import DatabaseManager
    from pichanalysis.core.organism import set_organism
    project=project_fixture(tmp_path/"projects");missing=MitoCartaPage(DatabaseManager(tmp_path/"empty"));missing.set_project(project)
    assert "not installed" in missing.readiness.text().lower() and not missing.run_button.isEnabled()
    manager,_=manager_fixture(tmp_path/"ready");(project.root/"mapping/tables/protein_catalog.csv").unlink();page=MitoCartaPage(manager);page.set_project(project)
    assert "Identification and Annotation" in page.readiness.text() and not page.run_button.isEnabled()
    project=project_fixture(tmp_path/"other");(project.root/"analyses/presence_absence/tables/classification.csv").unlink();page.set_project(project);assert not page.target.model().item(page.target.findText("Shared")).isEnabled()
    set_organism(project,"Mus musculus","10090");page.refresh_readiness();assert "Homo sapiens only" in page.readiness.text() and not page.run_button.isEnabled()

def test_worker_success_failure_and_structured_safeguard(app,completed_runs):
    project,manager,first,_=completed_runs;params=MitoCartaParameters();events=[]
    success=MitoCartaAnalysisWorker(project,manager,RRuntime(),"x",params,runner=lambda *a,**k:first);success.succeeded.connect(lambda value:events.append(("success",value)));success.run();assert events[-1][0]=="success"
    failure=MitoCartaAnalysisWorker(project,manager,RRuntime(),"x",params,runner=lambda *a,**k:(_ for _ in ()).throw(RuntimeError("planned")));failure.failed.connect(lambda value:events.append(("failed",value)));failure.run();assert events[-1]==("failed","planned")
    outside=MitoCartaAnalysisWorker(project,manager,RRuntime(),"x",params,runner=lambda *a,**k:(_ for _ in ()).throw(MitoCartaTargetOutsideBackgroundError({"entities_outside_background":["NCBI:1"],"initial_target_size":2,"initial_background_size":1})));outside.target_outside_background.connect(lambda value:events.append(("outside",value)));outside.run();assert events[-1][0]=="outside"

def test_python_r_worker_continue_and_gui_smoke(app,completed_runs):
    project,manager,_first,_second=completed_runs;runtime=RRuntime();events=[]
    original=MitoCartaParameters(target_selection="All mapped entities",background_selection="Shared",minimum_overlap=1)
    blocked=MitoCartaAnalysisWorker(project,manager,runtime,"worker_blocked",original);blocked.target_outside_background.connect(lambda details:events.append(("blocked",details)));blocked.run()
    assert events[-1][0]=="blocked" and events[-1][1]["initial_target_size"]>events[-1][1]["initial_background_size"]
    adjusted=replace(original,allow_target_outside_background=True);completed=MitoCartaAnalysisWorker(project,manager,runtime,"worker_adjusted",adjusted);completed.succeeded.connect(lambda outputs:events.append(("completed",outputs)));completed.run()
    assert events[-1][0]=="completed";page=MitoCartaPage(manager);page.set_project(project);page.show_outputs(events[-1][1])
    assert page.summary.rowCount()>0 and page.membership.rowCount()>0 and page.overall.rowCount()==1 and page.graph_choice.count()>0;page.close()

def test_missing_optional_graph_is_safe_and_mandatory_output_is_reported(completed_runs):
    from pichanalysis.core.mitocarta_analysis import MissingMitoCartaOutputError
    project,_manager,first,_second=completed_runs
    first.graphs[0].unlink();assert len(read_mitocarta_outputs(project,first.metadata["run_id"]).graphs)==len(first.graphs)-1
    (first.run_root/"summary.csv").unlink()
    with pytest.raises(MissingMitoCartaOutputError):read_mitocarta_outputs(project,first.metadata["run_id"])

def test_continue_and_cancel_decisions(app,monkeypatch):
    window=MainWindow();choices=iter([QMessageBox.ButtonRole.RejectRole,QMessageBox.ButtonRole.AcceptRole])
    class FakeBox:
        ButtonRole=QMessageBox.ButtonRole
        def __init__(self,*a):self.buttons=[];self.chosen=None
        def setWindowTitle(self,*a):pass
        def setText(self,text):assert "2" in text and "4" in text and "3" in text
        def addButton(self,label,role):self.buttons.append((label,role));return label
        def exec(self):
            role=next(choices);self.chosen=next(label for label,item_role in self.buttons if item_role==role)
        def clickedButton(self):return self.chosen
    monkeypatch.setattr("pichanalysis.ui.main_window.QMessageBox",FakeBox);details={"entities_outside_background":["A","B"],"initial_target_size":4,"initial_background_size":3}
    assert window._confirm_mitocarta_adjustment(details) is False;assert window._confirm_mitocarta_adjustment(details) is True;window.close()
