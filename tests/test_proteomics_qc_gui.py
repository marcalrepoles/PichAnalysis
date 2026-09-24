import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication, QMessageBox, QTableWidget

from pichanalysis.core.proteomics_qc import QCParameters, run_proteomics_qc
from pichanalysis.core.r_runtime import RRuntime
from pichanalysis.ui.proteomics_qc_page import ProteomicsQCPage
from tests.smoke_proteomics_qc import fixture_project
from tests.smoke_proteomics_qc_gui import smoke


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


def test_selection_defaults_and_minimum(app, tmp_path):
    project, columns, _ = fixture_project(tmp_path)
    page = ProteomicsQCPage()
    page.set_project(project)
    assert page.selected_columns() == tuple(columns)
    assert page.transformation.currentData() == "log2_positive"
    assert page.zero_missing.isChecked()
    for row in range(1, page.samples.rowCount()):
        page.samples.item(row, 0).setCheckState(Qt.CheckState.Unchecked)
    assert not page.run_button.isEnabled()
    assert "at least two" in page.readiness.text()
    page.close()


@pytest.mark.parametrize("family,transform", [("raw_intensity", "log2_positive"),
    ("spectral_count", "log2p1"), ("other_quantitative", "none")])
def test_family_defaults(app, tmp_path, family, transform):
    project, columns, _ = fixture_project(tmp_path)
    for col in columns:
        project.config["columns"][col]["quantification_type"] = family
    page = ProteomicsQCPage()
    page.set_project(project)
    assert page.transformation.currentData() == transform
    assert page.zero_missing.isChecked() == (family != "other_quantitative")
    page.close()


def test_mixed_family_and_missing_metadata(app, tmp_path):
    project, columns, _ = fixture_project(tmp_path)
    project.config["columns"]["B1"]["quantification_type"] = "spectral_count"
    page = ProteomicsQCPage()
    page.set_project(project)
    assert page.selected_columns() == tuple(col for col in columns if col != "B1")
    page.samples.item(3, 0).setCheckState(Qt.CheckState.Checked)
    assert not page.run_button.isEnabled()
    assert "incompatible" in page.readiness.text()
    page.close()
    project.config["columns"]["B1"]["quantification_type"] = "lfq_intensity"
    project.config["columns"]["A1"]["condition"] = ""
    page = ProteomicsQCPage()
    page.set_project(project)
    assert not page.run_button.isEnabled()
    assert "A1: missing condition" in page.readiness.text()
    page.close()


def test_na_rendering_and_no_mutation(app):
    page = ProteomicsQCPage()
    table = QTableWidget()
    import pandas as pd
    frame = pd.DataFrame({"pearson": ["", "0"], "cv_percent": ["", "0"]})
    page._fill(table, frame)
    assert table.item(0, 0).text() == "NA"
    assert table.item(0, 1).text() == "NA"
    assert table.item(1, 0).text() == "0"
    assert frame.iloc[0, 0] == ""
    page.close()


def test_r_failure_worker_and_repeated_run_guard(app, tmp_path, monkeypatch):
    project, columns, _ = fixture_project(tmp_path)
    page = ProteomicsQCPage()
    page.set_project(project)
    monkeypatch.setattr(QMessageBox, "warning", lambda *args: None)
    class Failure:
        code = "r_execution_failure"
        details = {}
        def __str__(self):
            return "R failed"
    page._running = True
    page._start_run()
    assert page._workers == []
    page._running = False
    page._run_failed(Failure())
    assert "R execution failure" in page.status.text()
    page.close()


def test_incomplete_run_and_optional_graph(app, tmp_path):
    project, columns, _ = fixture_project(tmp_path)
    outputs = run_proteomics_qc(project, RRuntime(), run_id="gui_unit",
        parameters=QCParameters(tuple(columns)))
    page = ProteomicsQCPage()
    page.set_project(project)
    png = outputs["run_root"] / "plots/pca.png"
    png.unlink()
    from pichanalysis.core.proteomics_qc import load_run
    page.show_outputs(load_run(project, "gui_unit"))
    page.graph_choice.setCurrentIndex(page.graph_choice.findText("PCA"))
    assert "Graph not available" in page.graph_preview.text()
    (outputs["run_root"] / "summary.csv").unlink()
    page._refresh_history()
    assert "Incomplete run" in page.history.itemText(0)
    page.history.setCurrentIndex(-1)
    page.history.setCurrentIndex(0)
    assert "Incomplete run" in page.status.text()
    page.close()


def test_python_r_gui_smoke(app):
    smoke()


def test_all_missing_sample_and_unavailable_pca(app, tmp_path):
    import pandas as pd
    project, columns, path = fixture_project(tmp_path)
    data = pd.read_csv(path)
    data["A3"] = float("nan")
    for row in range(len(data)):
        data.loc[row, columns[row % 6]] = float("nan")
    data.to_csv(path, index=False)
    outputs = run_proteomics_qc(project, RRuntime(), run_id="degenerate",
        parameters=QCParameters(tuple(columns)))
    page = ProteomicsQCPage()
    page.set_project(project)
    page.show_outputs(outputs)
    assert len(page._frames["sample_metadata"]) == 6
    assert page._frames["pca_summary"].iloc[0].status != "Ready"
    assert "PCA is not available" in page.pca_note.text()
    assert "no_detected_features" in page._frames["warnings"].warning_code.tolist()
    page.close()


def test_worker_failure_and_safe_close(app, tmp_path, monkeypatch):
    from PySide6.QtCore import QEventLoop, QTimer
    from pichanalysis.ui.main_window import MainWindow
    from pichanalysis.ui.proteomics_qc_page import ProteomicsQCWorker
    project, columns, _ = fixture_project(tmp_path)
    worker = ProteomicsQCWorker(project, RRuntime(), "fake", QCParameters(tuple(columns)),
        runner=lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError("synthetic R failure")))
    failure = []
    worker.failed.connect(failure.append)
    loop = QEventLoop()
    timer = QTimer()
    timer.setSingleShot(True)
    timer.timeout.connect(loop.quit)
    worker.finished.connect(loop.quit)
    timer.start(10000)
    worker.start()
    loop.exec()
    timer.stop()
    worker.wait()
    assert failure and "synthetic R failure" in str(failure[0])
    window = MainWindow()
    page = window.analyses_page.proteomics_qc_page
    window.analyses_page.set_project(project)
    window.show()
    page._running = True
    monkeypatch.setattr(QMessageBox, "information", lambda *args: None)
    window.close()
    assert window.isVisible()
    page._running = False
    window.close()
    assert not window.isVisible()

def test_worker_invalid_configuration(app, tmp_path):
    from PySide6.QtCore import QEventLoop, QTimer
    from pichanalysis.ui.proteomics_qc_page import ProteomicsQCWorker
    project, columns, _ = fixture_project(tmp_path)
    project.config["columns"]["B1"]["quantification_type"] = "spectral_count"
    worker = ProteomicsQCWorker(project, RRuntime(), "invalid", QCParameters(tuple(columns)))
    failures = []
    worker.failed.connect(failures.append)
    loop = QEventLoop()
    timer = QTimer()
    timer.setSingleShot(True)
    timer.timeout.connect(loop.quit)
    worker.finished.connect(loop.quit)
    timer.start(10000)
    worker.start()
    loop.exec()
    timer.stop()
    worker.wait()
    assert len(failures) == 1
    assert failures[0].code == "mixed_quantification_types"