"""Offline Database Manager Gene Ontology card checks."""
from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication

from pichanalysis.core.database_manager import DatabaseManager
from pichanalysis.ui.gene_ontology_card import GeneOntologyCard


def test_gene_ontology_card_absent(tmp_path):
    app = QApplication.instance() or QApplication([])
    manager = DatabaseManager(tmp_path / "databases")
    card = GeneOntologyCard(manager)
    assert "Not installed" in card.state_label.text()
    assert card.primary_button.text() == "Download"
    assert card.primary_button.isEnabled()
    assert not card.cancel_button.isEnabled()
    card.close()
