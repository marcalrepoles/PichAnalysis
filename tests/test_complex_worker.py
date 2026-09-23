"""Complexes worker and GUI pipeline smoke without network access."""
import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtCore import QEventLoop, QTimer
from PySide6.QtWidgets import QApplication

from pichanalysis.core.complex_analysis import ComplexParameters, ComplexTargetOutsideBackgroundError
from pichanalysis.core.r_runtime import RRuntime
from pichanalysis.ui.complexes_page import ComplexAnalysisWorker, ComplexPage
from smoke_complex_analysis import build_fixture


def _finish(worker):
    loop = QEventLoop()
    worker.finished.connect(loop.quit)
    worker.start()
    QTimer.singleShot(30000, loop.quit)  # Failure watchdog, not a synchronization sleep.
    loop.exec()
    assert not worker.isRunning()
    QApplication.processEvents()


def test_worker_python_r_gui(tmp_path):
    app = QApplication.instance() or QApplication([])
    project, manager = build_fixture(tmp_path)
    page = ComplexPage(manager); page.set_project(project)
    runtime = RRuntime()
    if not runtime.available:
        import pytest
        pytest.skip("Rscript unavailable")
    page.target.setCurrentIndex(page.target.findData("Manual selection"))
    for row in range(6): page.manual.selectRow(row)
    parameters=ComplexParameters(**page.parameters())
    assert parameters.manual_rows==(1,2,3,4,5,6)
    worker = ComplexAnalysisWorker(project,manager,runtime,"worker_smoke",parameters)
    failures=[]
    worker.failed.connect(failures.append)
    worker.succeeded.connect(page.show_outputs)
    page.set_running(True)
    assert not page.run_button.isEnabled()
    _finish(worker)
    page.set_running(False)
    assert not failures and page.outputs is not None
    page.show_complex_details("CPX-1")
    assert page.coverage_progress.format()=="Protein-component coverage: 3 / 3 groups"
    assert page.graph_choice.count()==6
    assert page.tables["complex_enrichment_all"].rowCount()>0
    assert page.tables["protein_to_complexes"].rowCount()>0
    page.close(); app.processEvents()


def test_worker_failure_and_safeguard(tmp_path):
    app = QApplication.instance() or QApplication([])
    project, manager = build_fixture(tmp_path)
    def failure(*args,**kwargs): raise RuntimeError("controlled R failure")
    bad = ComplexAnalysisWorker(project,manager,RRuntime(),"bad",ComplexParameters(),runner=failure)
    errors=[]; bad.failed.connect(errors.append); _finish(bad)
    assert errors==["controlled R failure"]
    details={"entities_outside_background":["P12345"],"initial_target_size":2,"initial_background_size":1}
    def outside(*args,**kwargs): raise ComplexTargetOutsideBackgroundError(details)
    guarded=ComplexAnalysisWorker(project,manager,RRuntime(),"outside",ComplexParameters(),runner=outside)
    received=[]; guarded.target_outside.connect(received.append); _finish(guarded)
    assert received==[details]
    assert not (project.root/"analyses/Complexes/runs/outside").exists()
    app.processEvents()
