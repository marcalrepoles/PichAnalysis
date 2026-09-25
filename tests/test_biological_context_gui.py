"""Offline GUI contracts for the single-run Biological Context Explorer."""
from __future__ import annotations

from types import SimpleNamespace
from pathlib import Path

import pandas as pd
import pytest
from PySide6.QtGui import QImage
from PySide6.QtWidgets import QApplication, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget

from pichanalysis.core.biological_context import ContextItem
from pichanalysis.core.entity_identity import EntityContext, InputLineage
from pichanalysis.ui import biological_context_page as gui
from pichanalysis.ui.cross_module_actions import attach_biological_context_action


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


def items(tmp_path):
    official = tmp_path / "official.png"
    image = QImage(12, 12, QImage.Format_RGB32)
    image.fill(0xFFFFFF)
    assert image.save(str(official), "PNG")
    return [
        ContextItem("kegg", "K1", "hsa00010", "KEGG pathway", "X", "hsa:1",
            pd.DataFrame({"kegg_gene_id": ["hsa:1"], "gene_symbol": ["GENE"]}), pd.DataFrame(), official),
        ContextItem("reactome", "R1", "R-HSA-1", "Reactome pathway", "Y", "U:P1",
            pd.DataFrame({"reactome_entity_key": ["U:P1"], "Gene_symbol": ["GENE"]}), pd.DataFrame(), official),
        ContextItem("string", "S1", "9606.P1", "STRING network", "Z", "9606.P1",
            pd.DataFrame({"string_protein_id": ["9606.P1"], "gene_symbol": ["GENE"]}),
            pd.DataFrame({"protein_a": ["9606.P1"], "protein_b": ["9606.P1"]}), network_type="physical"),
    ]


def test_shared_action_requires_selected_persisted_row(app):
    page = QWidget()
    page.setLayout(QVBoxLayout())
    page.outputs = SimpleNamespace(metadata={"run_id": "historical-A"})
    table = QTableWidget(1, 2)
    table.setHorizontalHeaderLabels(["uniprot_accession", "source_row"])
    table.setItem(0, 0, QTableWidgetItem("P12345-2"))
    table.setItem(0, 1, QTableWidgetItem("9"))
    page.layout().addWidget(table)
    seen = []
    button = attach_biological_context_action(page, "reactome", lambda *args: seen.append(args))
    page.show()
    assert not button.isEnabled()
    table.selectRow(0)
    assert button.isEnabled()
    button.click()
    assert seen == [("reactome", "historical-A", "P12345-2", 9)]
    page.close()


def test_explorer_tables_official_view_generation_export_history(app, tmp_path, monkeypatch):
    source = EntityContext("differential", "D1", "all_results", "F1", 1, "P1",
        uniprot_accessions=("P1",), gene_symbols=("GENE",),
        input_lineage=InputLineage("project", "same-input"))
    index = object()
    prepared = items(tmp_path)
    monkeypatch.setattr(gui, "resolve_entity", lambda *args: (index, source))
    monkeypatch.setattr(gui, "find_contexts", lambda *args, **kwargs: prepared)
    monkeypatch.setattr(gui, "differential_choices", lambda *args: ["D1", "D2"])
    monkeypatch.setattr(gui, "other_differential_runs", lambda *args: [])
    image = QImage(100, 80, QImage.Format_RGB32)
    image.fill(0xFFFFFF)
    monkeypatch.setattr(gui, "render_context", lambda *args: image)
    opened = []
    monkeypatch.setattr(gui.QDesktopServices, "openUrl", lambda url: opened.append(url.toLocalFile()) or True)
    page = gui.BiologicalContextPage()
    page.set_project(SimpleNamespace(root=tmp_path))
    page.show()
    page.open_entity("differential", "D1", "F1", 1)
    assert page.differential.currentData() == "D1"
    assert [page.tables[m].rowCount() for m in ("kegg", "reactome", "string")] == [1, 1, 1]
    page.tables["kegg"].selectRow(0)
    assert page.official_button.isEnabled()
    page.open_official()
    assert [Path(value) for value in opened] == [prepared[0].official_image]
    page.build_view()
    assert page.current_root and (page.current_root / "figures/context.png").is_file()
    assert page.current_figure is not None
    export = tmp_path / "export.png"
    monkeypatch.setattr(gui.QFileDialog, "getSaveFileName",
        lambda *args: (str(export), "PNG (*.png)"))
    page.export_png()
    assert export.is_file()
    pdf = tmp_path / "export.pdf"
    monkeypatch.setattr(gui.QFileDialog, "getSaveFileName",
        lambda *args: (str(pdf), "PDF (*.pdf)"))
    page.export_pdf()
    assert pdf.read_bytes().startswith(b"%PDF")
    page.tables["reactome"].selectRow(0)
    page.build_view()
    assert page.current_root and page.current_root != prepared[0].official_image
    page.tables["string"].selectRow(0)
    assert not page.official_button.isEnabled()
    page.build_view()
    assert page.current_root
    history_root = page.current_root
    page.set_project(SimpleNamespace(root=tmp_path))
    position = page.history.findData(str(history_root))
    assert position >= 0
    page.history.setCurrentIndex(position)
    assert "snapshot Z" in page.notice.text()
    assert page.current_figure is not None
    assert not page.export_table_button.isEnabled()
    assert not page.official_button.isEnabled()
    page.close()


def test_incomplete_source_clears_previous_figure(app, tmp_path, monkeypatch):
    page = gui.BiologicalContextPage()
    page.set_project(SimpleNamespace(root=tmp_path))
    page.current_figure = QImage(10, 10, QImage.Format_RGB32)
    monkeypatch.setattr(gui, "resolve_entity", lambda *args: (_ for _ in ()).throw(
        RuntimeError("Source run is incomplete")))
    with pytest.raises(RuntimeError, match="incomplete"):
        page.open_entity("differential", "broken", "F1")
    assert page.current_figure is None and page.figure.pixmap().isNull()
    page.close()

def test_mapping_source_action_carries_frozen_row(app, tmp_path):
    import json
    from pichanalysis.core.project import create_project
    from pichanalysis.core.mapping_analysis import read_mapping_outputs
    from pichanalysis.ui.analyses_page import AnalysesPage
    project = create_project(tmp_path, "mapping-context")
    tables = project.root / "mapping/tables"
    tables.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame({"source_row": [9], "original_id": ["P12345-2"],
                          "uniprot_accession": ["P12345-2"], "gene_symbol": ["GENE"]})
    for name in ("protein_catalog", "id_mapping", "unmapped", "ambiguous"):
        frame.to_csv(tables / f"{name}.csv", index=False)
    (project.root / "mapping/latest_metadata.json").write_text(
        json.dumps({"run_id": "mapping-A"}), encoding="utf-8")
    page = AnalysesPage()
    page.set_project(project)
    page.show_outputs(read_mapping_outputs(project))
    assert not page.mapping_biological_context.isEnabled()
    page.show()
    page.preview.selectRow(0)
    assert page.mapping_biological_context.isEnabled()
    seen = []
    page.open_biological_context = lambda *args: seen.append(args)
    page._mapping_biological_context()
    assert seen == [("mapping", "mapping-A", "P12345-2", 9)]
    page.close()


def test_differential_source_action_carries_feature_and_row(app):
    from pichanalysis.ui.differential_analysis_page import DifferentialAnalysisPage
    from PySide6.QtCore import Qt
    page = DifferentialAnalysisPage()
    page.statistics = {"metadata": {"run_id": "D1"}}
    page.results_table.setRowCount(1)
    page.results_table.setColumnCount(1)
    cell = QTableWidgetItem("F1")
    cell.setData(Qt.UserRole, "F1")
    cell.setData(Qt.UserRole + 1, "9")
    page.results_table.setItem(0, 0, cell)
    assert not page.biological_context_button.isEnabled()
    page.show()
    page.results_table.selectRow(0)
    assert page.biological_context_button.isEnabled()
    seen = []
    page.biological_context_requested.connect(lambda *args: seen.append(args))
    page.biological_context_button.click()
    assert seen == [("differential", "D1", "F1", 9)]
    page.close()
