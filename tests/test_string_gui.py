import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication

from pichanalysis.core.r_runtime import RRuntime
from pichanalysis.core.string_analysis import StringParameters, read_string_outputs, run_string_analysis
from pichanalysis.ui.main_window import MainWindow
from pichanalysis.ui.string_page import StringPage
from pichanalysis.ui.string_viewer import VIEWER_NODE_LIMIT
from test_string_analysis import IDS, scientific_fixture


def app():
    return QApplication.instance() or QApplication([])


def test_string_tab_and_readiness(tmp_path):
    app()
    window = MainWindow()
    assert window.analyses_page.string_page is not None
    page = window.analyses_page.string_page
    project, manager = scientific_fixture(tmp_path)
    page.manager = manager
    page.set_project(project)
    assert "Ready for STRING" in page.readiness.text()
    assert page.seed_count.text() == "Selected seeds: 3"
    assert page.parameters()["combined_score_threshold"] is None
    page.threshold_preset.setCurrentIndex(3)
    page.custom_threshold.setValue(555)
    assert page.parameters()["combined_score_threshold"] == 555
    page.network_type.setCurrentIndex(1)
    page.max_hop.setCurrentIndex(2)
    page.selection_mode.setCurrentIndex(1)
    assert page.parameters()["network_type"] == "physical"
    assert page.parameters()["max_hop"] == 2
    assert page.parameters()["degree1_selection_mode"] == "strict_common"
    window.close()


def test_python_r_gui_history_viewer_and_export(tmp_path):
    app()
    project, manager = scientific_fixture(tmp_path)
    runtime = RRuntime()
    if not runtime.available:
        import pytest
        pytest.skip("Rscript unavailable")
    a = run_string_analysis(project, manager, runtime, run_id="a", parameters=StringParameters(max_hop=2), timeout=120)
    b = run_string_analysis(project, manager, runtime, run_id="b", parameters=StringParameters(max_hop=1, network_type="physical"), timeout=120)
    page = StringPage(manager)
    page.set_project(project)
    page.show_outputs(a)
    assert len(page.viewer.node_items) == len(a["tables"]["expanded_nodes"])
    assert len(page.viewer.edge_items) == len(a["tables"]["expanded_edges"])
    assert page.viewer.search_node() is None
    page.viewer.search.setText("PA")
    assert page.viewer.search_node() == IDS["A"]
    assert "STRING ID" in page.node_details.text()
    page.viewer.edge_selected.emit(IDS["A"], IDS["X"])
    assert "confidence score" in page.edge_details.text()
    image = tmp_path / "network.png"
    page.viewer.export_image(image)
    assert image.is_file() and image.stat().st_size > 0
    page.show_outputs(b)
    assert "Physical" in page.viewer.title.text()
    page.show_outputs(read_string_outputs(project, "a"))
    assert "Functional" in page.viewer.title.text()
    first_snapshot = a["metadata"]["snapshot_id"]
    manager.string.install_from_directory(tmp_path / "raw")
    newer = run_string_analysis(project, manager, runtime, run_id="c", parameters=StringParameters(max_hop=0), timeout=120)
    page.show_outputs(newer)
    assert newer["metadata"]["snapshot_id"] != first_snapshot
    page.show_outputs(read_string_outputs(project, "a"))
    assert first_snapshot in page.history_status.text()
    assert len(page.viewer.edge_items) == len(a["tables"]["expanded_edges"])
    assert VIEWER_NODE_LIMIT >= len(page.viewer.node_items)
