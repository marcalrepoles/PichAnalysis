import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest
from PySide6.QtWidgets import QApplication

from pichanalysis.core.r_runtime import RRuntime
from pichanalysis.core.string_analysis import StringParameters
from pichanalysis.ui.string_page import StringAnalysisWorker, StringPage
from test_string_analysis import IDS, scientific_fixture


def test_thresholds_strict_and_worker_to_gui(tmp_path):
    app = QApplication.instance() or QApplication([])
    project, manager = scientific_fixture(tmp_path)
    runtime = RRuntime()
    if not runtime.available:
        pytest.skip("Rscript unavailable")
    page = StringPage(manager)
    page.set_project(project)
    assert "Ready for STRING" in page.readiness.text()
    for preset, score in (("Exploratory", 150), ("Functional", 400), ("Robust", 700)):
        worker = StringAnalysisWorker(project, manager, runtime, preset.lower(), StringParameters(threshold_preset=preset, max_hop=2))
        results = []
        errors = []
        worker.succeeded.connect(results.append)
        worker.failed.connect(errors.append)
        worker.run()
        assert not errors and len(results) == 1
        page.show_outputs(results[0])
        assert f"Threshold ≥ {score}" in page.viewer.title.text()
        assert page.outputs["metadata"]["threshold_preset"] == preset
        assert page.history.count() >= 1
        assert page.tables["summary"].rowCount() > 0
        assert page.tables["node_metrics"].rowCount() > 0
        assert page.graph_choice.count() > 0
    strict_worker = StringAnalysisWorker(project, manager, runtime, "strict", StringParameters(max_hop=2, degree1_selection_mode="strict_common"))
    strict_results = []
    strict_worker.succeeded.connect(strict_results.append)
    strict_worker.run()
    page.show_outputs(strict_results[0])
    assert page.tables["common_direct_neighbors"].rowCount() == 1
    assert page.outputs["tables"]["common_direct_neighbors"].string_protein_id.tolist() == [IDS["X"]]
    assert "Strict common" in page.viewer.title.text()
