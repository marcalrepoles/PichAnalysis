import os
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest
from PySide6.QtWidgets import QApplication, QGroupBox, QMessageBox

from pichanalysis.core.database_manager import DatabaseManager
from pichanalysis.core.database_registry import DatabaseState
from pichanalysis.core.databases.reactome import CORE_FILES, ReactomeDownloadCancelled, ReactomeProvider
from pichanalysis.ui.database_manager_page import DatabaseManagerPage, ReactomeCoreDownloadWorker


CORE_DATA = {
    "ReactomePathways.txt": b"R-HSA-1\tPathway one\tHomo sapiens\n",
    "ReactomePathwaysRelation.txt": b"R-HSA-1\tR-HSA-2\n",
    "UniProt2Reactome.txt": b"P04637\tR-HSA-1\tTP53\tEvidence\tHomo sapiens\n",
    "NCBI2Reactome.txt": b"7157\tR-HSA-1\tTP53\tEvidence\tHomo sapiens\n",
    "humanPathwaysWithDiagrams.txt": b"R-HSA-1\n",
    "pathway2summation.txt": b"R-HSA-1\tA valid summary\n",
}


class FakeReactomeProvider:
    def __init__(self, fail_at: str | None = None, cancel_on: str | None = None):
        self.fail_at = fail_at
        self.cancel_on = cancel_on
        self.cancel_callback = None
        self.calls = []

    def download(self, name, destination, progress=None, cancel_requested=None):
        self.calls.append(name)
        if name == self.fail_at:
            raise RuntimeError("planned transport failure")
        data = CORE_DATA[name]
        halfway = max(len(data) // 2, 1)
        if progress:
            progress(name, halfway, len(data))
        if name == self.cancel_on:
            assert self.cancel_callback is not None
            self.cancel_callback()
        if cancel_requested and cancel_requested():
            raise ReactomeDownloadCancelled("Reactome Core Data download was canceled.")
        Path(destination).write_bytes(data)
        if progress:
            progress(name, len(data), len(data))
        return {"name": name, "file_size": len(data)}


def write_core(directory: Path, label: str = "Pathway one") -> Path:
    directory.mkdir()
    for name, data in CORE_DATA.items():
        if name == "ReactomePathways.txt":
            data = f"R-HSA-1\t{label}\tHomo sapiens\n".encode()
        (directory / name).write_bytes(data)
    return directory


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


def run_worker(manager, provider):
    worker = ReactomeCoreDownloadWorker(manager, provider)
    events = {"progress": [], "success": [], "failure": [], "canceled": []}
    worker.progress.connect(events["progress"].append)
    worker.succeeded.connect(events["success"].append)
    worker.failed.connect(events["failure"].append)
    worker.canceled.connect(events["canceled"].append)
    worker.run()
    return worker, events


def test_reactome_worker_success_activates_ready_snapshot_and_emits_progress(tmp_path):
    manager = DatabaseManager(tmp_path / "db")
    _, events = run_worker(manager, FakeReactomeProvider())
    active = manager.reactome.active_snapshot()
    assert active is not None
    assert manager.reactome.manifest(active)["status"] == DatabaseState.READY
    assert len(events["success"]) == 1 and not events["failure"] and not events["canceled"]
    assert {item["filename"] for item in events["progress"]} == set(CORE_FILES)
    last = events["progress"][-1]
    assert last["file_index"] == 6 and last["file_total"] == 6
    assert last["bytes_downloaded"] == last["bytes_total"] > 0


def test_reactome_worker_failure_does_not_activate_snapshot(tmp_path):
    manager = DatabaseManager(tmp_path / "db")
    _, events = run_worker(manager, FakeReactomeProvider(fail_at=CORE_FILES[2]))
    assert manager.reactome.active_snapshot() is None
    assert manager.reactome.state() == DatabaseState.ERROR
    assert events["failure"] == ["planned transport failure"] and not events["success"]


def test_reactome_worker_cancel_is_cooperative_and_does_not_activate(tmp_path):
    manager = DatabaseManager(tmp_path / "db")
    provider = FakeReactomeProvider(cancel_on=CORE_FILES[0])
    worker = ReactomeCoreDownloadWorker(manager, provider)
    provider.cancel_callback = worker.cancel
    canceled = []
    worker.canceled.connect(canceled.append)
    worker.run()
    assert canceled == ["Reactome Core Data download was canceled."]
    assert manager.reactome.active_snapshot() is None
    assert manager.reactome.state() == DatabaseState.INCOMPLETE
    assert provider.calls == [CORE_FILES[0]]


def test_provider_cancel_during_stream_removes_part_file(tmp_path, monkeypatch):
    class Response:
        headers = {"Content-Length": "4"}

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def read(self, _size):
            return b"data"

    monkeypatch.setattr("urllib.request.urlopen", lambda *args, **kwargs: Response())
    checks = iter((False, True))
    destination = tmp_path / CORE_FILES[0]
    with pytest.raises(ReactomeDownloadCancelled):
        ReactomeProvider(retries=0).download(destination.name, destination, cancel_requested=lambda: next(checks))
    assert not destination.exists()
    assert not destination.with_name(destination.name + ".part").exists()


def test_valid_update_activates_b_and_preserves_a(tmp_path):
    manager = DatabaseManager(tmp_path / "db")
    first = manager.reactome.install_from_directory(write_core(tmp_path / "a", "A"))
    _, events = run_worker(manager, FakeReactomeProvider())
    second = manager.reactome.active_snapshot()
    assert events["success"] and second != first
    assert first.is_dir() and second.is_dir()


def test_failed_update_preserves_active_a(tmp_path):
    manager = DatabaseManager(tmp_path / "db")
    first = manager.reactome.install_from_directory(write_core(tmp_path / "a", "A"))
    _, events = run_worker(manager, FakeReactomeProvider(fail_at=CORE_FILES[1]))
    assert events["failure"]
    assert manager.reactome.active_snapshot() == first
    assert manager.reactome.is_core_ready()


def test_canceled_update_preserves_active_a(tmp_path):
    manager = DatabaseManager(tmp_path / "db")
    first = manager.reactome.install_from_directory(write_core(tmp_path / "a", "A"))
    provider = FakeReactomeProvider(cancel_on=CORE_FILES[0])
    worker = ReactomeCoreDownloadWorker(manager, provider)
    provider.cancel_callback = worker.cancel
    worker.run()
    assert manager.reactome.active_snapshot() == first
    assert manager.reactome.is_core_ready()


def test_reactome_card_not_installed_downloading_ready_and_failed_update(app, tmp_path, monkeypatch):
    manager = DatabaseManager(tmp_path / "db")
    page = DatabaseManagerPage(manager)
    page.show()
    app.processEvents()
    assert page.findChild(type(page.reactome_state_label), None) is not None
    assert page.reactome_state_label.text() == "Core Data: Not installed"
    assert page.reactome_primary_button.text() == "Download"
    assert page.reactome_primary_button.isEnabled()
    assert not page.reactome_cancel_button.isEnabled()

    class RunningWorker:
        @staticmethod
        def isRunning():
            return True

    page.reactome_worker = RunningWorker()
    page._refresh_reactome()
    page._show_reactome_progress({
        "file_index": 3, "file_total": 6, "filename": "UniProt2Reactome.txt",
        "bytes_downloaded": 45, "bytes_total": 90,
    })
    assert page.reactome_state_label.text() == "Core Data: Downloading"
    assert not page.reactome_primary_button.isEnabled() and page.reactome_cancel_button.isEnabled()
    assert "File 3 of 6" in page.reactome_progress_text.text()
    assert "UniProt2Reactome.txt" in page.reactome_progress_text.text()

    page.reactome_worker = None
    active = manager.reactome.install_from_directory(write_core(tmp_path / "ready", "Ready"))
    manifest = manager.reactome.manifest(active)
    manifest["release_version"] = "97"
    manager.reactome._write_manifest(active, manifest)
    page.refresh()
    assert page.reactome_state_label.text() == "Core Data: Ready"
    assert page.reactome_primary_button.text() == "Update"
    assert "Release: 97" in page.reactome_details.text()
    assert f"Snapshot: {active.name}" in page.reactome_details.text()
    assert "Downloaded: —" not in page.reactome_details.text()

    monkeypatch.setattr(QMessageBox, "warning", lambda *args: QMessageBox.StandardButton.Ok)
    page._reactome_failed("planned update failure")
    assert page.reactome_state_label.text() == "Core Data: Ready"
    assert "planned update failure" in page.reactome_error_label.text()
    page.close()


def test_reactome_card_uses_indeterminate_progress_when_size_is_unknown(app, tmp_path):
    page = DatabaseManagerPage(DatabaseManager(tmp_path / "db"))
    page._show_reactome_progress({
        "file_index": 1, "file_total": 6, "filename": CORE_FILES[0],
        "bytes_downloaded": 0, "bytes_total": 0,
    })
    assert page.reactome_progress.minimum() == 0 and page.reactome_progress.maximum() == 0
    assert "Size unknown" in page.reactome_progress_text.text()
    page.close()


def test_main_window_reactome_gui_smoke_offscreen(app, tmp_path, monkeypatch):
    import pichanalysis.ui.database_manager_page as database_page_module
    from pichanalysis.ui.main_window import MainWindow

    manager = DatabaseManager(tmp_path / "db")
    monkeypatch.setattr(database_page_module, "DatabaseManager", lambda: manager)
    window = MainWindow()
    window.show()
    window.navigation.setCurrentRow(3)
    app.processEvents()
    page = window.database_manager_page
    group_titles = {group.title() for group in page.findChildren(QGroupBox)}
    assert "Reactome — Homo sapiens" in group_titles
    assert page.reactome_state_label.text() == "Core Data: Not installed"

    class RunningWorker:
        @staticmethod
        def isRunning():
            return True

    page.reactome_worker = RunningWorker()
    page._refresh_reactome()
    page._show_reactome_progress({
        "file_index": 2, "file_total": 6, "filename": CORE_FILES[1],
        "bytes_downloaded": 10, "bytes_total": 20,
    })
    assert page.reactome_cancel_button.isEnabled()
    assert not page.reactome_primary_button.isEnabled()

    page.reactome_worker = None
    manager.reactome.install_from_directory(write_core(tmp_path / "gui-ready", "GUI"))
    page.refresh()
    assert page.reactome_state_label.text() == "Core Data: Ready"
    assert page.reactome_primary_button.text() == "Update"
    assert window.close()
    app.processEvents()
