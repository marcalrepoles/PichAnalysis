"""Offline UI navigation contracts and legacy Mapping isolation."""
import json
from types import SimpleNamespace

import pandas as pd
import pytest
from PySide6.QtWidgets import QApplication, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget

from pichanalysis.core.cross_module_integration import (
    _historical_mapping, _mapping_runs, open_cross_module_index,
)
from pichanalysis.core.project import create_project
from pichanalysis.core.mapping_analysis import read_mapping_outputs
from pichanalysis.ui.analyses_page import AnalysesPage
from pichanalysis.ui.cross_module_actions import attach_explore_action
from pichanalysis.ui.cross_module_navigation import open_persisted_target


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


def test_legacy_mapping_never_borrows_latest(tmp_path):
    project = create_project(tmp_path, "integration")
    raw = project.root / "mapping/raw"
    for run_id in ("old", "new"):
        folder = raw / run_id
        folder.mkdir(parents=True)
        (folder / "metadata.json").write_text(json.dumps({"run_id": run_id}), encoding="utf-8")
    tables = project.root / "mapping/tables"
    tables.mkdir(parents=True, exist_ok=True)
    latest = pd.DataFrame({"source_row": [1], "original_id": ["LATEST"],
                           "uniprot_accession": ["P99999"]})
    for name in ("protein_catalog", "id_mapping", "unmapped", "ambiguous"):
        latest.to_csv(tables / f"{name}.csv", index=False)
    (project.root / "mapping/latest_metadata.json").write_text(
        json.dumps({"run_id": "new"}), encoding="utf-8")
    assert _mapping_runs(project) == ["new", "old"]
    with pytest.raises(FileNotFoundError, match="Historical mapping data"):
        _historical_mapping(project, "old")
    index = open_cross_module_index(project)
    status = index.run_status("mapping", "old")
    assert status[0] == "Historical mapping unavailable"
    assert "safe cross-module resolution" in status[1]
    assert not any(row["run_id"] == "old" for row in index.search("LATEST"))
    assert any(row["run_id"] == "new" for row in index.search("LATEST"))
    app = QApplication.instance() or QApplication([])
    page = AnalysesPage()
    page.set_project(project)
    page.show_outputs(read_mapping_outputs(project))
    page.show()
    page.preview.selectRow(0)
    page.explore_mapping.click()
    assert page.cross_module_explorer.context.source_module == "mapping"
    assert page.cross_module_explorer.context.source_run_id == "new"
    assert not page.cross_module_explorer.explore_feature("mapping", "old", "LATEST")
    assert "safe cross-module resolution" in page.cross_module_explorer.status.text()
    page.close()


@pytest.mark.parametrize("module_id", ["presence_absence", "go", "kegg", "reactome",
    "mitocarta", "domains", "string", "complexes", "mtdna_evidence"])
def test_shared_explore_action_uses_loaded_run_and_selection(app, module_id):
    page = QWidget()
    page.setLayout(QVBoxLayout())
    page.outputs = SimpleNamespace(metadata={"run_id": "historical_A"})
    table = QTableWidget(1, 1)
    table.setHorizontalHeaderLabels(["uniprot_accession"])
    table.setItem(0, 0, QTableWidgetItem("P12345-2"))
    page.layout().addWidget(table)
    seen = []
    button, _ = attach_explore_action(page, module_id, lambda *args: seen.append(args))
    page.show()
    table.selectRow(0)
    button.click()
    assert seen == [(module_id, "historical_A", "P12345-2")]
    page.close()


def test_open_target_loads_requested_run_not_latest_and_keeps_fallback(app):
    class Page(QWidget):
        def __init__(self):
            super().__init__()
            self.setLayout(QVBoxLayout())
            self.outputs = None

        def show_outputs(self, outputs):
            self.outputs = outputs

    page = Page()
    latest = SimpleNamespace(metadata={"run_id": "B"})
    old = SimpleNamespace(metadata={"run_id": "A"})
    page.outputs = latest
    adapter = SimpleNamespace(module_id="synthetic", list_runs=lambda project: ["A", "B"],
        load_run=lambda project, run_id: {"A": old, "B": latest}[run_id])
    focused = open_persisted_target(None, adapter, page, "A", "P12345-2")
    assert not focused
    assert page.outputs is old
    assert "Loaded run: A | Target: P12345-2" in page.cross_module_target_hint.text()
    page.close()
