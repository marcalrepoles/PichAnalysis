"""Offscreen offline application -> worker -> R -> results/history smoke."""
from __future__ import annotations

import os
import tempfile
from unittest.mock import patch
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtCore import QEventLoop, QTimer
from PySide6.QtWidgets import QApplication

from pichanalysis.core.mtdna_analysis import MtdnaParameters
from pichanalysis.core.r_runtime import RRuntime
from pichanalysis.ui.main_window import MainWindow
from pichanalysis.ui.mtdna_page import MtdnaAnalysisWorker
from smoke_mtdna_analysis import build_fixture
from test_mtdna_evidence_database import sources


def run_worker(project,manager,run_id,selection):
    worker=MtdnaAnalysisWorker(project,manager,RRuntime(),run_id,
        MtdnaParameters(target_selection=selection,minimum_overlap=2))
    loop=QEventLoop();outcome={}
    worker.succeeded.connect(lambda value:outcome.update(outputs=value))
    worker.failed.connect(lambda value:outcome.update(error=value))
    worker.target_outside.connect(lambda value:outcome.update(error=str(value)))
    worker.finished.connect(loop.quit)
    QTimer.singleShot(60000,loop.quit)
    worker.start();loop.exec();worker.wait(1000)
    assert not worker.isRunning(),"mtDNA worker did not finish"
    assert "error" not in outcome,outcome.get("error")
    assert "outputs" in outcome,"mtDNA worker produced no outputs"
    return outcome["outputs"]


def main():
    app=QApplication.instance() or QApplication([])
    root=Path(tempfile.mkdtemp(prefix="pichanalysis-mtdna-gui-"))
    project,manager=build_fixture(root)
    window=MainWindow()
    window.database_manager_page.manager=manager
    page=window.analyses_page.mtdna_page
    page.manager=manager
    window._activate_project(project,"Offline mtDNA GUI smoke")
    assert page.database.text().find("Status: Ready")>=0
    assert page.run_button.isEnabled()
    assert page.target.count()>2 and page.background.currentData()=="All mapped entities"
    assert page.minimum_overlap.value()==3 and page.fdr.value()==.05 and page.top_n.value()==20
    page.target.setCurrentIndex(page.target.findData("A-specific"))
    page.minimum_overlap.setValue(2)
    page.run_button.click()
    assert window.mtdna_worker is not None
    gui_loop=QEventLoop()
    window.mtdna_worker.finished.connect(gui_loop.quit)
    QTimer.singleShot(60000,gui_loop.quit)
    gui_loop.exec()
    assert window.mtdna_worker is None and page.outputs is not None
    assert page.outputs["tables"]["summary"].shape[0]>0
    outputs_a=run_worker(project,manager,"gui-A","A-specific")
    page.show_outputs(outputs_a)
    assert page.tables["summary"].rowCount()>0
    assert page.tables["entity_evidence"].rowCount()>0
    assert page.tables["category_frequency"].rowCount()>0
    assert page.tables["enrichment_all"].rowCount()>0
    assert page.tables["source_frequency"].rowCount()==3
    assert page.tables["evidence_records"].rowCount()>0
    assert page.tables["comparison"].rowCount()>0
    assert page.graph_choice.count()>=5
    page.tables["entity_evidence"].selectRow(0)
    assert page.entity_detail.rowCount()>0
    page.tables["category_frequency"].selectRow(0)
    assert page.category_entities.rowCount()>0
    page.tabs.setCurrentIndex(1)
    with patch("pichanalysis.ui.mtdna_page.QFileDialog.getSaveFileName",return_value=(str(root/"export_table.csv"),"")):
        page._export_current()
    with patch("pichanalysis.ui.mtdna_page.QFileDialog.getSaveFileName",return_value=(str(root/"export_workbook.xlsx"),"")):
        page._export_workbook()
    with patch("pichanalysis.ui.mtdna_page.QFileDialog.getSaveFileName",return_value=(str(root/"export_graph.png"),"")):
        page._export_graph()
    page.tabs.setCurrentIndex(page.tabs.indexOf(page.tables["entity_evidence"].parentWidget()))
    with patch("pichanalysis.ui.mtdna_page.QFileDialog.getSaveFileName",return_value=(str(root/"export_evidence.csv"),"")):
        page._export_evidence()
    assert all((root/name).is_file() for name in ("export_table.csv","export_workbook.xlsx",
        "export_graph.png","export_evidence.csv"))
    outputs_b=run_worker(project,manager,"gui-B","B-specific")
    page.show_outputs(outputs_b)
    page.history.setCurrentIndex(page.history.findData("gui-A"))
    assert page.outputs["metadata"]["run_id"]=="gui-A"
    assert all(page._frames[key].equals(outputs_a["tables"][key]) for key in outputs_a["tables"])
    snapshot_x=outputs_a["metadata"]["snapshot_id"]
    mito,go=sources(manager.root,"fixture-B")
    manager.mtdna_evidence.build(mitocarta=mito,go_snapshot=go,ncbi_xml=root/"ncbi.xml")
    page.set_project(project)
    outputs_y=run_worker(project,manager,"gui-Y","A-specific")
    assert outputs_y["metadata"]["snapshot_id"]!=snapshot_x
    page.show_outputs(outputs_y)
    page.history.setCurrentIndex(page.history.findData("gui-A"))
    assert page.outputs["metadata"]["snapshot_id"]==snapshot_x
    assert page.graph_choice.count()>=5
    window.close()
    assert not window.isVisible()
    print(f"mtDNA GUI smoke passed: A->B->A and snapshot {snapshot_x}->{outputs_y['metadata']['snapshot_id']}->{snapshot_x}")


if __name__=="__main__":main()
