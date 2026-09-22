import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtCore import QEventLoop, QTimer
from PySide6.QtWidgets import QApplication

from pichanalysis.core.database_manager import DatabaseManager
from pichanalysis.ui.database_manager_page import DatabaseManagerPage
from test_complex_portal_database import fixture
from test_complex_portal_gui import FixtureProvider


def test_card_runs_download_on_qthread_and_finishes(tmp_path):
    app = QApplication.instance() or QApplication([])
    source = fixture(tmp_path / "9606.tsv")
    manager = DatabaseManager(tmp_path / "db")
    manager.complex_portal.provider = FixtureProvider(source)
    page = DatabaseManagerPage(manager)
    card = page.complex_portal_card
    loop = QEventLoop()
    timed_out = []
    def timeout():
        timed_out.append(True)
        loop.quit()
    timer = QTimer()
    timer.setSingleShot(True)
    timer.timeout.connect(timeout)
    card.start()
    worker = card.worker
    assert worker is not None
    assert page.is_running()
    worker.finished.connect(loop.quit)
    timer.start(10000)
    loop.exec()
    timer.stop()
    app.processEvents()
    assert not timed_out
    assert manager.complex_portal.validate_snapshot()
    assert not page.is_running()
    assert "Ready" in card.state_label.text()
