import os
from types import SimpleNamespace

import pytest
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (QApplication, QComboBox, QLineEdit,
    QListWidget, QListWidgetItem, QAbstractItemView)

from pichanalysis.ui.analyses_page import AnalysesPage

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


@pytest.mark.parametrize("destination,derived,rows", [
    ("go", "a_only", [1]), ("kegg", "b_only", [2]),
    ("reactome", "shared", [1, 3]), ("string", "opposite_direction", [3]),
])
def test_handoff_preloads_target_without_starting_analysis(tmp_path, monkeypatch,
        destination, derived, rows):
    QApplication.instance() or QApplication([])
    target = QComboBox()
    target.addItem("Manual selection", "manual" if destination in {"go", "kegg"}
        else "Manual selection")
    manual_list = QListWidget()
    manual_list.setSelectionMode(QAbstractItemView.SelectionMode.MultiSelection)
    for value in (1, 2, 3):
        item = QListWidgetItem(str(value))
        item.setData(Qt.ItemDataRole.UserRole, value)
        manual_list.addItem(item)
    page = SimpleNamespace(target=target, manual_rows=QLineEdit(),
        manual=QLineEdit() if destination == "kegg" else manual_list)
    tab_calls = []
    tabs = SimpleNamespace(setCurrentWidget=tab_calls.append)
    parent = SimpleNamespace(project=SimpleNamespace(root=tmp_path), module_tabs=tabs,
        go_page=page, kegg_page=page, reactome_page=page,
        mitocarta_page=page, mtdna_page=page, interpro_pfam_page=page,
        string_page=page, complex_page=page)
    from pichanalysis.ui import analyses_page
    messages = []
    monkeypatch.setattr(analyses_page.QMessageBox, "information",
        lambda *args: messages.append(args[-1]))
    monkeypatch.setattr(analyses_page.QMessageBox, "warning",
        lambda *args: pytest.fail(str(args[-1])))
    handoff = {"comparison_run_id": "comparison-a", "derived_set": derived,
        "destination": destination, "source_rows": rows}
    AnalysesPage.open_derived_target(parent, handoff)
    assert tab_calls == [page]
    assert target.currentData() in {"manual", "Manual selection"}
    if destination == "go":
        assert page.manual_rows.text() == "1"
    elif destination == "kegg":
        assert page.manual.text() == "2"
    else:
        assert sorted(item.data(Qt.ItemDataRole.UserRole)
            for item in page.manual.selectedItems()) == rows
    assert page.derived_target_provenance["derived_set"] == derived
    assert messages and "run the analysis manually" in messages[0]
    assert not (tmp_path / "analyses" / destination / "runs").exists()
