"""Database Manager card for curated human Complex Portal snapshots."""
from __future__ import annotations

import threading

from PySide6.QtCore import QThread, QUrl, Signal
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import QGroupBox, QHBoxLayout, QLabel, QMessageBox, QProgressBar, QPushButton, QVBoxLayout

from ..core.database_registry import DatabaseState
from ..core.databases.complex_portal import ComplexPortalCancelled


class ComplexPortalDownloadWorker(QThread):
    progress = Signal(dict)
    succeeded = Signal(str)
    failed = Signal(str)
    canceled = Signal(str)

    def __init__(self, manager, provider=None):
        super().__init__()
        self.manager = manager
        self._cancel = threading.Event()
        if provider is not None:
            self.manager.complex_portal.provider = provider

    def run(self):
        try:
            snapshot = self.manager.complex_portal.download_and_build(self.progress.emit, self._cancel.is_set)
            self.succeeded.emit(str(snapshot))
        except ComplexPortalCancelled as error:
            self.canceled.emit(str(error))
        except Exception as error:
            self.failed.emit(str(error))

    def cancel(self):
        self._cancel.set()


class ComplexPortalDatabaseCard(QGroupBox):
    def __init__(self, manager):
        super().__init__("Complex Portal — Homo sapiens")
        self.manager = manager
        self.worker = None
        self._workers = []
        self.last_error = ""
        self.state_label = QLabel()
        self.details = QLabel()
        self.details.setWordWrap(True)
        self.progress = QProgressBar()
        self.progress.setRange(0, 1)
        self.progress_text = QLabel("No download is running.")
        self.progress_text.setWordWrap(True)
        self.error_label = QLabel()
        self.error_label.setWordWrap(True)
        self.primary_button = QPushButton("Download")
        self.cancel_button = QPushButton("Cancel")
        self.folder_button = QPushButton("Open database folder")
        layout = QVBoxLayout(self)
        for widget in (self.state_label, self.details, self.progress, self.progress_text, self.error_label):
            layout.addWidget(widget)
        actions = QHBoxLayout()
        for button in (self.primary_button, self.cancel_button, self.folder_button):
            actions.addWidget(button)
        layout.addLayout(actions)
        self.primary_button.clicked.connect(self.start)
        self.cancel_button.clicked.connect(self.cancel)
        self.folder_button.clicked.connect(self.open_folder)
        self.refresh()

    def refresh(self):
        database = self.manager.complex_portal
        active = database.active_snapshot()
        manifest = database.manifest(active) if active else {}
        running = self.is_running()
        state = DatabaseState.DOWNLOADING if running else database.state()
        self.state_label.setText(f"Status: {state}")
        self.details.setText(
            "Scope: Manually curated complexes (predicted hu.MAP/MuSIC excluded)\n"
            f"Release: {manifest.get('complex_portal_release', 'Unknown')}\n"
            f"Snapshot: {active.name if active else '—'}\n"
            f"Downloaded: {manifest.get('completed_at', '—')}\n"
            f"Curated complexes: {manifest.get('complex_count', '—')}\n"
            f"Unique proteins: {manifest.get('unique_protein_count', '—')}\n"
            f"Protein memberships: {manifest.get('protein_membership_count', '—')}\n"
            f"Non-protein participants: {manifest.get('nonprotein_participant_count', '—')}\n"
            f"Location: {database.root}"
        )
        self.primary_button.setText("Update" if active else "Download")
        self.primary_button.setEnabled(not running)
        self.cancel_button.setEnabled(running)
        self.error_label.setText(f"Last attempt failed: {self.last_error}" if self.last_error else "")

    def start(self):
        if self.is_running():
            return
        self.last_error = ""
        worker = ComplexPortalDownloadWorker(self.manager)
        self.worker = worker
        self._workers.append(worker)
        worker.progress.connect(self.show_progress)
        worker.succeeded.connect(self.finished_ok)
        worker.failed.connect(self.failed)
        worker.canceled.connect(self.canceled)
        worker.finished.connect(lambda w=worker: self._worker_finished(w))
        worker.finished.connect(worker.deleteLater)
        worker.start()
        self.refresh()

    def _worker_finished(self, worker):
        if worker in self._workers:
            self._workers.remove(worker)
        self.refresh()

    def show_progress(self, payload):
        if payload.get("stage") == "download":
            done = int(payload.get("bytes_downloaded") or 0)
            total = int(payload.get("bytes_total") or 0)
            if total:
                self.progress.setRange(0, total)
                self.progress.setValue(done)
            else:
                self.progress.setRange(0, 0)
        else:
            self.progress.setRange(0, 0)
        self.progress_text.setText(str(payload.get("message", "Processing Complex Portal data...")))

    def finished_ok(self, _snapshot):
        self.worker = None
        self.progress.setRange(0, 1)
        self.progress.setValue(1)
        self.progress_text.setText("Curated human Complex Portal database is ready.")
        self.refresh()

    def failed(self, message):
        self.worker = None
        self.last_error = message
        self.progress_text.setText("Complex Portal update failed; the previous snapshot remains active.")
        self.refresh()
        QMessageBox.warning(self, "Complex Portal download failed", message)

    def canceled(self, message):
        self.worker = None
        self.progress.setRange(0, 1)
        self.progress.setValue(0)
        self.progress_text.setText(message)
        self.refresh()

    def cancel(self):
        if self.worker:
            self.progress_text.setText("Cancellation requested. The current safe operation will stop first.")
            self.worker.cancel()

    def is_running(self):
        return bool(self._workers)

    def open_folder(self):
        self.manager.complex_portal.root.mkdir(parents=True, exist_ok=True)
        if not QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.manager.complex_portal.root))):
            QMessageBox.warning(self, "Database folder", "Could not open the database folder.")
