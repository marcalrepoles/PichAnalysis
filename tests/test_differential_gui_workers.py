"""Differential GUI worker lifecycle checks without arbitrary sleeps."""
import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest
from PySide6.QtCore import QEventLoop, QTimer
from PySide6.QtGui import QCloseEvent
from PySide6.QtWidgets import QApplication, QMessageBox

from pichanalysis.ui.differential_analysis_page import DifferentialAnalysisPage, DifferentialWorker
from pichanalysis.ui.main_window import MainWindow
from tests.smoke_differential_preparation import fixture_project


def _wait(worker):
    loop = QEventLoop()
    timer = QTimer()
    timer.setSingleShot(True)
    timer.timeout.connect(loop.quit)
    worker.finished.connect(loop.quit)
    timer.start(10000)
    loop.exec()
    assert timer.isActive()
    QApplication.processEvents()


@pytest.fixture
def app():
    return QApplication.instance() or QApplication([])


def test_worker_success_and_failure(app):
    success, failure = [], []
    good = DifferentialWorker(lambda *args, **kwargs: "ready", None, None, "good", None)
    good.succeeded.connect(success.append)
    good.start(); _wait(good)
    assert success == ["ready"]

    def broken(*args, **kwargs):
        raise RuntimeError("simulated R failure")

    bad = DifferentialWorker(broken, None, None, "bad", None)
    bad.failed.connect(failure.append)
    bad.start(); _wait(bad)
    assert len(failure) == 1 and "simulated R failure" in str(failure[0])


def test_repeated_execution_is_blocked(app, tmp_path):
    project, _, _ = fixture_project(tmp_path)
    page = DifferentialAnalysisPage()
    page.set_project(project)
    assert page.prepare_button.isEnabled()
    page._running = True
    page._update_state()
    page._start_preparation()
    assert page.worker is None
    assert not page.prepare_button.isEnabled()
    page._running = False
    page.close()


def test_main_window_refuses_close_during_differential_worker(app, monkeypatch):
    window = MainWindow()
    window.analyses_page.differential_analysis_page._running = True
    monkeypatch.setattr(QMessageBox, "information", lambda *args, **kwargs: None)
    event = QCloseEvent()
    window.closeEvent(event)
    assert not event.isAccepted()
    window.analyses_page.differential_analysis_page._running = False
    window.close()
