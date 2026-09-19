from __future__ import annotations

import os
from dataclasses import replace
import pandas as pd
import pytest
from PySide6.QtWidgets import QApplication, QMessageBox

from pichanalysis.core.interpro_pfam_analysis import InterProPfamParameters, read_interpro_pfam_outputs, run_interpro_pfam_analysis
from pichanalysis.core.r_runtime import RRuntime
from pichanalysis.ui.interpro_pfam_page import AnnotationSetBuildWorker, InterProPfamAnalysisWorker, InterProPfamPage
from pichanalysis.ui.main_window import MainWindow
from test_interpro_database import FakeProvider
from test_interpro_pfam_analysis import fixture

os.environ.setdefault("QT_QPA_PLATFORM","offscreen")

@pytest.fixture(scope="module")
def app():yield QApplication.instance() or QApplication([])

@pytest.fixture()
def gui_runs(tmp_path):
    project,manager=fixture(tmp_path);manager.interpro.build_annotation_set(project.root,["P1","P2","EMPTY"],annotation_set_id="setB");manager.interpro.provider=object();runtime=RRuntime()
    a=run_interpro_pfam_analysis(project,manager,runtime,run_id="A",parameters=InterProPfamParameters(annotation_set_id="setA",target_selection="Shared",minimum_overlap=1))
    b=run_interpro_pfam_analysis(project,manager,runtime,run_id="B",parameters=InterProPfamParameters(annotation_set_id="setB",target_selection="A-specific",minimum_overlap=1))
    return project,manager,a,b

def test_page_readiness_results_navigation_filters_graphs_and_history(app,gui_runs):
    project,manager,a,b=gui_runs;page=InterProPfamPage(manager);page.set_project(project)
    assert "Ready" in page.metadata.text() and page.annotation_set.count()==2 and page.run_button.isEnabled();assert page.background.currentText()=="All annotation-set-covered proteins in the experiment"
    assert page.fdr.value()==pytest.approx(.05) and page.minimum.value()==3 and page.top_n.value()==20
    page.show_outputs(a);assert page.summary.rowCount()>0 and page.coverage.rowCount()==3 and page.interpro_frequency.rowCount()>0 and page.pfam_frequency.rowCount()>0
    assert page.graph_choice.count()==6 and not page.graph_preview.pixmap().isNull();assert page.repeated.rowCount()>0
    page.interpro_type.setCurrentIndex(1);assert page.interpro_frequency.rowCount()>0
    page.protein.setCurrentIndex(page.protein.findText("P1"));assert page.protein_interpro.rowCount()>0 and page.protein_pfam.rowCount()>0 and page.protein_pfam_locations.rowCount()>=3
    assert page.pfam_architecture.rowCount()>0 and "PF00001" in page.pfam_architecture.item(0,1).text()
    ia=page.history.findData("A");ib=page.history.findData("B");page.history.setCurrentIndex(ia);assert page.outputs.metadata["annotation_set_id"]=="setA";page.history.setCurrentIndex(ib);assert page.outputs.metadata["annotation_set_id"]=="setB";page.history.setCurrentIndex(ia);assert page.outputs.metadata["annotation_set_id"]=="setA";page.close()

def test_readiness_missing_metadata_set_presence_and_coverage(app,tmp_path):
    project,manager=fixture(tmp_path);page=InterProPfamPage(manager);page.set_project(project);assert page.run_button.isEnabled()
    manager.interpro.active_metadata_pointer.unlink();page.refresh_readiness();assert "not installed" in page.readiness.text().lower() and not page.run_button.isEnabled()
    empty_manager=type(manager)(tmp_path/"empty");empty_manager.interpro.provider=FakeProvider();empty_manager.interpro.download_metadata();empty_project=type(project)(tmp_path/"empty_project",project.config);empty=InterProPfamPage(empty_manager);empty.set_project(empty_project);assert "Build or select" in empty.readiness.text()
    (project.root/"analyses/presence_absence/tables/classification.csv").unlink();empty._populate_targets();assert not empty.target.model().item(empty.target.findText("Shared")).isEnabled();page.close();empty.close()

def test_acquisition_and_analysis_workers(app,tmp_path):
    project,manager=fixture(tmp_path);events=[];build=AnnotationSetBuildWorker(manager.interpro,project.root,["P1" ]);build.succeeded.connect(lambda value:events.append("built"));build.run();assert events[-1]=="built"
    canceled=AnnotationSetBuildWorker(manager.interpro,project.root,["P1"]);canceled.cancel();canceled.canceled.connect(lambda value:events.append("cancelled"));canceled.run();assert events[-1]=="cancelled"
    failed=AnnotationSetBuildWorker(manager.interpro,project.root,["P1"]);failed.database=object();failed.failed.connect(lambda value:events.append("failed"));failed.run();assert events[-1]=="failed"
    output=read_interpro_pfam_outputs(project,"missing") if False else object();worker=InterProPfamAnalysisWorker(project,manager,RRuntime(),"x",InterProPfamParameters(),runner=lambda *a,**k:output);worker.succeeded.connect(lambda value:events.append("analysis"));worker.run();assert events[-1]=="analysis"
    error=InterProPfamAnalysisWorker(project,manager,RRuntime(),"x",InterProPfamParameters(),runner=lambda *a,**k:(_ for _ in ()).throw(RuntimeError("planned")));error.failed.connect(lambda value:events.append("error"));error.run();assert events[-1]=="error"

def test_continue_cancel_dialog(app,monkeypatch):
    window=MainWindow();choices=iter([QMessageBox.ButtonRole.RejectRole,QMessageBox.ButtonRole.AcceptRole])
    class Box:
        ButtonRole=QMessageBox.ButtonRole
        def __init__(self,*a):self.buttons=[]
        def setWindowTitle(self,*a):pass
        def setText(self,text):assert "2" in text and "4" in text and "3" in text
        def addButton(self,label,role):self.buttons.append((label,role));return label
        def exec(self):
            desired=next(choices);self.choice=next(label for label,role in self.buttons if role==desired)
        def clickedButton(self):return self.choice
    monkeypatch.setattr("pichanalysis.ui.main_window.QMessageBox",Box);details={"entities_outside_background":["P1","P2"],"initial_target_size":4,"initial_background_size":3};assert not window._confirm_interpro_pfam_adjustment(details);assert window._confirm_interpro_pfam_adjustment(details);window.close()

def test_missing_optional_plot_and_required_output(gui_runs):
    from pichanalysis.core.interpro_pfam_analysis import MissingInterProPfamOutputError
    project,_manager,a,_b=gui_runs;a.plots[0].unlink();assert len(read_interpro_pfam_outputs(project,"A").plots)==len(a.plots)-1;(a.run_root/"summary.csv").unlink()
    with pytest.raises(MissingInterProPfamOutputError):read_interpro_pfam_outputs(project,"A")
