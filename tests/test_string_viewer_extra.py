import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pandas as pd
import pytest
from PySide6.QtWidgets import QApplication

from pichanalysis.core.r_runtime import RRuntime
from pichanalysis.core.string_analysis import MissingStringOutputError, StringParameters, read_string_outputs, run_string_analysis
from pichanalysis.ui.string_page import StringPage
from pichanalysis.ui.string_viewer import VIEWER_NODE_LIMIT
from test_string_analysis import IDS, scientific_fixture


def test_viewer_filters_limit_exports_and_incomplete_run(tmp_path):
    app = QApplication.instance() or QApplication([])
    project, manager = scientific_fixture(tmp_path)
    runtime = RRuntime()
    if not runtime.available:
        pytest.skip("Rscript unavailable")
    outputs = run_string_analysis(project, manager, runtime, run_id="filters", parameters=StringParameters(max_hop=2, degree1_selection_mode="strict_common"), timeout=120)
    page = StringPage(manager)
    page.set_project(project)
    page.show_outputs(outputs)
    viewer = page.viewer
    viewer.highlight_common_nodes()
    assert viewer.highlighted == {IDS["X"]}
    assert viewer.node_items[IDS["X"]].pen().widthF() == 4
    viewer.layers.setCurrentIndex(1)
    assert len(viewer.node_items) == 3
    viewer.layers.setCurrentIndex(0)
    viewer.show_labels.setChecked(False)
    assert all(not item.label.isVisible() for item in viewer.node_items.values())
    node = viewer.node_items[IDS["A"]]
    original = node.pos()
    node.setPos(original.x() + 100, original.y())
    viewer.refresh_edges()
    viewer.reset_layout()
    assert node.pos() == original
    viewer.zoom_in.click()
    viewer.reset_view.click()
    svg = tmp_path / "view.svg"
    viewer.export_image(svg)
    assert svg.is_file() and svg.stat().st_size > 0
    large = {**outputs, "tables": dict(outputs["tables"])}
    template = outputs["tables"]["expanded_nodes"].iloc[0].to_dict()
    nodes = pd.DataFrame([{**template, "string_protein_id": f"synthetic-{i}", "hop_level": 0 if i == 0 else 1} for i in range(VIEWER_NODE_LIMIT + 1)])
    large["tables"]["expanded_nodes"] = nodes
    large["tables"]["node_metrics"] = pd.DataFrame({"string_protein_id": nodes.string_protein_id, "component_id": [1] * len(nodes)})
    large["tables"]["hubs"] = pd.DataFrame(columns=["string_protein_id"])
    large["tables"]["seeds"] = pd.DataFrame(columns=["string_protein_id"])
    viewer.load_outputs(large)
    assert not viewer.node_items
    assert str(VIEWER_NODE_LIMIT + 1) in viewer.message.text()
    viewer.layers.setCurrentIndex(1)
    assert len(viewer.node_items) == 1
    plot = next(outputs["run_root"].glob("plots/*.png"))
    plot.unlink()
    assert len(read_string_outputs(project, "filters")["plots"]) == 11
    edge_path = outputs["run_root"] / "networks/expanded_edges.csv"
    edge_path.unlink()
    with pytest.raises(MissingStringOutputError):
        read_string_outputs(project, "filters")
    page._load_history(page.history.findData("filters"))
    assert page.tables["summary"].rowCount() == 0
    assert not page.viewer.node_items
    assert "incomplete" in page.history_status.text()
