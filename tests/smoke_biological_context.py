"""Offline single-run Biological Context GUI smoke using prepared frozen fixtures."""
from __future__ import annotations

import tempfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pandas as pd
from PySide6.QtGui import QImage
from PySide6.QtWidgets import QApplication

from pichanalysis.core.biological_context import ContextItem
from pichanalysis.core.entity_identity import EntityContext, InputLineage
from pichanalysis.core.project import create_project
from pichanalysis.ui import biological_context_page as gui


def main() -> int:
    app = QApplication.instance() or QApplication([])
    with tempfile.TemporaryDirectory(prefix="pich_context_smoke_") as temporary:
        project = create_project(Path(temporary), "Biological context fixture")
        snapshot = project.root / "fixture-snapshot"
        snapshot.mkdir()
        official = snapshot / "official.png"
        image = QImage(20, 20, QImage.Format_RGB32)
        image.fill(0xFFFFFF)
        assert image.save(str(official), "PNG")
        source = EntityContext("differential", "D1", "all_results", "F1", 1, "P1",
            uniprot_accessions=("P1",), gene_symbols=("GENE",),
            input_lineage=InputLineage("fixture", "frozen-input"))
        items = [
            ContextItem("kegg", "K1", "hsa00010", "Fixture KEGG pathway", "X", "hsa:1",
                pd.DataFrame({"kegg_gene_id": ["hsa:1"], "gene_symbol": ["GENE"]}),
                pd.DataFrame(), official),
            ContextItem("reactome", "R1", "R-HSA-1", "Fixture Reactome pathway", "Y", "U:P1",
                pd.DataFrame({"reactome_entity_key": ["U:P1"], "Gene_symbol": ["GENE"]}),
                pd.DataFrame(), official),
            ContextItem("string", "S1", "9606.P1", "Fixture STRING network", "Z", "9606.P1",
                pd.DataFrame({"string_protein_id": ["9606.P1", "9606.P2"],
                    "gene_symbol": ["GENE", "OTHER"]}),
                pd.DataFrame({"protein_a": ["9606.P1"], "protein_b": ["9606.P2"]})),
        ]
        with patch.object(gui, "resolve_entity", return_value=(SimpleNamespace(), source)), \
             patch.object(gui, "find_contexts", return_value=items), \
             patch.object(gui, "differential_choices", return_value=["D1"]), \
             patch.object(gui, "other_differential_runs", return_value=[]):
            page = gui.BiologicalContextPage()
            page.set_project(project)
            page.show()
            page.open_entity("differential", "D1", "F1", 1)
            assert page.differential.currentData() == "D1"
            # This fixture has no real statistics; neutral rendering is explicit.
            page.differential.setCurrentIndex(0)
            assert [page.tables[module].rowCount() for module in ("kegg", "reactome", "string")] == [1, 1, 1]
            for module in ("kegg", "reactome", "string"):
                page.tables[module].selectRow(0)
                if module != "string":
                    assert page.official_button.isEnabled()
                    assert page.current_item.official_image.read_bytes() == official.read_bytes()
                page.build_view()
                assert page.current_root and (page.current_root / "figures/context.png").is_file()
            destination = project.root / "export.png"
            with patch.object(gui.QFileDialog, "getSaveFileName", return_value=(str(destination), "PNG (*.png)")):
                page.export_png()
            assert destination.is_file()
            history = page.current_root
            page.set_project(project)
            page.history.setCurrentIndex(page.history.findData(str(history)))
            assert page.current_figure is not None
            assert "snapshot Z" in page.notice.text()
            page.close()
        app.processEvents()
    print("Biological Context offline GUI smoke: passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())