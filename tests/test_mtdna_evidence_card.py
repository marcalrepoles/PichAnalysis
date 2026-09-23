"""Database Manager mtDNA evidence dependency/readiness UI tests."""
import os

os.environ.setdefault("QT_QPA_PLATFORM","offscreen")

from PySide6.QtWidgets import QApplication

from pichanalysis.core.database_manager import DatabaseManager
from pichanalysis.ui.database_manager_page import DatabaseManagerPage


def test_card_reports_missing_dependencies_and_does_not_build(tmp_path):
    app=QApplication.instance() or QApplication([])
    page=DatabaseManagerPage(DatabaseManager(tmp_path))
    card=page.mtdna_evidence_card
    assert "MitoCarta database is required." in card.details.text()
    assert "Gene Ontology database is required." in card.details.text()
    assert not card.primary_button.isEnabled()
    card.start()
    assert not card.is_running() and not page.manager.mtdna_evidence.is_ready()
    page.close();app.processEvents()
