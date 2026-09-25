import os
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication

from pichanalysis.core.experiment_comparison import ComparisonSpec, TableSpec, compare
from pichanalysis.ui.experiment_comparison_page import ExperimentComparisonPage, ResultTable


os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


def test_sort_search_and_selected_row_detail():
    app = QApplication.instance() or QApplication([])
    view = ResultTable("Shared")
    view.set_frame(pd.DataFrame({"comparison_entity": ["P12345", "Q11111"],
        "A_source_rows": ["2", "1"]}))
    view.table.sortByColumn(0, Qt.DescendingOrder)
    view.table.selectRow(0)
    assert view.selected_record()["comparison_entity"] == "Q11111"
    assert "Q11111" in view.detail.toPlainText()
    view.search.setText("p12345")
    assert view.proxy.rowCount() == 1
    view.table.selectRow(0)
    assert view.selected_record()["comparison_entity"] == "P12345"
    view.close()


def test_frozen_table_picker_and_explicit_context_side(tmp_path, monkeypatch):
    app = QApplication.instance() or QApplication([])
    a, b = tmp_path / "a.csv", tmp_path / "b.csv"
    pd.DataFrame({"id": ["P12345"], "description": ["A"]}).to_csv(a, index=False)
    pd.DataFrame({"id": ["P12345"], "description": ["B"]}).to_csv(b, index=False)
    project = SimpleNamespace(root=tmp_path, config={"input": {"processed_file": None}})
    compare(project, ComparisonSpec(TableSpec(str(a), "id", "uniprot"),
        TableSpec(str(b), "id", "uniprot"), "uniprot"), run_id="frozen")
    page = ExperimentComparisonPage()
    page.set_project(project)
    from pichanalysis.ui import experiment_comparison_page as gui
    monkeypatch.setattr(gui.QInputDialog, "getItem", lambda *args: (
        "Comparison frozen / dataset A", True))
    page.source_a._frozen()
    assert page.source_a.path.is_file()
    assert page.source_a.source_kind == "comparison:frozen:a"
    page.loaded = {"config": {"a": {"source_kind": "differential:diff_a", "alias": "A"},
        "b": {"source_kind": "differential:diff_b", "alias": "B"}}}
    page.views["shared"].set_frame(pd.DataFrame({"comparison_entity": ["P12345"],
        "A_present": [True], "B_present": [True]}))
    page.tabs.setCurrentWidget(page.views["shared"])
    page.views["shared"].table.selectRow(0)
    monkeypatch.setattr(gui.QInputDialog, "getItem", lambda *args: ("Dataset B - B", True))
    emitted = []
    page.biological_context_requested.connect(lambda *args: emitted.append(args))
    page._open_context()
    assert emitted == [("differential", "diff_b", "P12345", None)]
    page.close()
