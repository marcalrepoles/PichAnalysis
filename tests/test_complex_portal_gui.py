import os
import shutil

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication, QMessageBox

from pichanalysis.core.database_manager import DatabaseManager
from pichanalysis.core.databases.complex_portal import ComplexPortalCancelled
from pichanalysis.ui.complex_portal_database_card import ComplexPortalDownloadWorker
from pichanalysis.ui.database_manager_page import DatabaseManagerPage
from test_complex_portal_database import fixture


class FixtureProvider:
    def __init__(self, source, outcome="ok"):
        self.source = source
        self.outcome = outcome

    def download(self, target, progress=None, cancel_requested=None):
        if self.outcome == "cancel" or (cancel_requested and cancel_requested()):
            raise ComplexPortalCancelled("Complex Portal download was canceled.")
        if self.outcome == "error":
            raise RuntimeError("Fixture download failed")
        shutil.copy2(self.source, target)
        if progress:
            progress(target.stat().st_size, target.stat().st_size)
        import hashlib
        digest = hashlib.sha256(target.read_bytes()).hexdigest()
        return {"filename": target.name, "source_url": "fixture://curated-human",
                "resolved_url": "fixture://curated-human", "last_modified": None, "etag": None,
                "content_length": target.stat().st_size, "size": target.stat().st_size, "sha256": digest}


def test_card_worker_update_cancel_and_progress(tmp_path, monkeypatch):
    app = QApplication.instance() or QApplication([])
    manager = DatabaseManager(tmp_path / "db")
    source = fixture(tmp_path / "9606.tsv")
    page = DatabaseManagerPage(manager)
    card = page.complex_portal_card
    assert "Not installed" in card.state_label.text()
    assert "Manually curated" in card.details.text()
    assert card.primary_button.text() == "Download"
    worker = ComplexPortalDownloadWorker(manager, FixtureProvider(source))
    results, errors, progress = [], [], []
    worker.succeeded.connect(results.append)
    worker.failed.connect(errors.append)
    worker.progress.connect(progress.append)
    worker.run()
    assert results and not errors
    assert [event["stage"] for event in progress] == ["download", "download", "validate", "parse", "parse", "membership", "index", "validate"]
    card.refresh()
    assert "Ready" in card.state_label.text()
    assert card.primary_button.text() == "Update"
    first = manager.complex_portal.active_snapshot()
    worker = ComplexPortalDownloadWorker(manager, FixtureProvider(source))
    worker.run()
    assert manager.complex_portal.active_snapshot() != first
    active = manager.complex_portal.active_snapshot()
    worker = ComplexPortalDownloadWorker(manager, FixtureProvider(source, "cancel"))
    canceled = []
    worker.canceled.connect(canceled.append)
    worker.run()
    assert canceled and manager.complex_portal.active_snapshot() == active
    worker = ComplexPortalDownloadWorker(manager, FixtureProvider(source, "error"))
    errors = []
    worker.failed.connect(errors.append)
    worker.run()
    assert errors and manager.complex_portal.active_snapshot() == active
    monkeypatch.setattr(QMessageBox, "warning", lambda *args: None)
    card.failed("Fixture download failed")
    assert "previous snapshot remains active" in card.progress_text.text()
    card.canceled("Complex Portal download was canceled.")
    assert "canceled" in card.progress_text.text()
    card.show_progress({"stage": "download", "message": "Downloading Complex Portal data...", "bytes_downloaded": 10, "bytes_total": 100})
    assert card.progress.value() == 10
    card.show_progress({"stage": "index", "message": "Indexing database..."})
    assert "Indexing" in card.progress_text.text()
    assert not page.is_running()
