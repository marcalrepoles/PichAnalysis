import os
from types import SimpleNamespace

import pandas as pd
from PySide6.QtWidgets import QApplication

from pichanalysis.ui.experiment_comparison_page import ExperimentComparisonPage


os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


def test_gui_compare_history_and_export(tmp_path, monkeypatch):
    app = QApplication.instance() or QApplication([])
    a, b = tmp_path / "a.csv", tmp_path / "b.csv"
    pd.DataFrame({"accession": ["P12345", "P12346"], "description": ["A", "B"]}).to_csv(a, index=False)
    pd.DataFrame({"accession": ["P12345", "P12347"], "description": ["A", "C"]}).to_csv(b, index=False)
    project = SimpleNamespace(root=tmp_path, config={"input": {"processed_file": None}})
    page = ExperimentComparisonPage()
    page.set_project(project)
    page.source_a.set_source(a)
    page.source_b.set_source(b)
    page.source_a.kind.setCurrentIndex(page.source_a.kind.findData("uniprot"))
    page.source_b.kind.setCurrentIndex(page.source_b.kind.findData("uniprot"))
    page.source_a.confirmed.setChecked(True)
    page.source_b.confirmed.setChecked(True)
    page.comparison_type.setCurrentIndex(page.comparison_type.findData("uniprot"))
    page._compare()
    assert page.loaded is not None
    assert page.views["shared"].table.model().rowCount() == 1
    assert page.views["a_only"].table.model().rowCount() == 1
    assert page.views["b_only"].table.model().rowCount() == 1
    run_id = page.loaded["run_root"].name
    page.set_project(project)
    page.history.setCurrentIndex(page.history.findData(run_id))
    page._load_history()
    assert page.loaded["run_root"].name == run_id
    export = tmp_path / "export.xlsx"
    from pichanalysis.ui import experiment_comparison_page as gui
    monkeypatch.setattr(gui.QFileDialog, "getSaveFileName", lambda *args: (str(export), ""))
    page._export_workbook()
    assert "Summary" in pd.ExcelFile(export).sheet_names
    page.close()
