"""Database Manager controls for official human Gene Ontology snapshots."""
from __future__ import annotations

import threading

from PySide6.QtCore import QThread, QUrl, Signal
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import QGroupBox, QHBoxLayout, QLabel, QMessageBox, QPushButton, QVBoxLayout

from ..core.gene_ontology_database import GeneOntologyCancelled


class GeneOntologyWorker(QThread):
    progress = Signal(dict)
    succeeded = Signal(str)
    failed = Signal(str)
    canceled = Signal(str)

    def __init__(self, database):
        super().__init__()
        self.database = database
        self._cancel = threading.Event()

    def run(self):
        try:
            snapshot = self.database.download_and_build(progress=self.progress.emit,
                cancel_requested=self._cancel.is_set)
            self.succeeded.emit(str(snapshot))
        except GeneOntologyCancelled as error:
            self.canceled.emit(str(error))
        except Exception as error:
            self.failed.emit(str(error))

    def cancel(self):
        self._cancel.set()


class GeneOntologyCard(QGroupBox):
    snapshot_changed = Signal()
    def __init__(self, manager):
        super().__init__("Gene Ontology — Homo sapiens")
        self.manager = manager
        self.worker = None
        self._workers = []
        self.last_error = ""
        self.state_label = QLabel()
        self.details = QLabel()
        self.details.setWordWrap(True)
        self.progress_label = QLabel("No download is running.")
        self.progress_label.setWordWrap(True)
        self.error_label = QLabel()
        self.error_label.setWordWrap(True)
        self.primary_button = QPushButton("Download")
        self.cancel_button = QPushButton("Cancel")
        self.folder_button = QPushButton("Open database folder")
        layout = QVBoxLayout(self)
        for widget in (self.state_label, self.details, self.progress_label, self.error_label):
            layout.addWidget(widget)
        buttons = QHBoxLayout()
        for widget in (self.primary_button, self.cancel_button, self.folder_button):
            buttons.addWidget(widget)
        layout.addLayout(buttons)
        self.primary_button.clicked.connect(self.start)
        self.cancel_button.clicked.connect(self.cancel)
        self.folder_button.clicked.connect(self.open_folder)
        self.refresh()

    def is_running(self):
        return bool(self._workers)

    def refresh(self):
        db = self.manager.gene_ontology
        active = db.active_snapshot()
        manifest = db.manifest(active) if active else {}
        counts = manifest.get("counts", {})
        self.state_label.setText(f"Status: {'Downloading' if self.is_running() else 'Ready' if active else 'Not installed'}")
        self.details.setText(f"Snapshot: {active.name if active else '—'}\n"
            f"Ontology version: {manifest.get('ontology_version', '—')}\n"
            f"Human annotations: {counts.get('annotations', '—')}\n"
            f"Downloaded: {manifest.get('completed_at', '—')}\n"
            f"Annotation version: {manifest.get('annotation_version', '—')}\n"
            f"Terms: {counts.get('terms', '—')}\n"
            f"Relations: {counts.get('relations', '—')}\nLocation: {db.root}")
        self.primary_button.setText("Update" if active else "Download")
        self.primary_button.setEnabled(not self.is_running())
        self.cancel_button.setEnabled(self.is_running())
        self.error_label.setText(f"Last attempt failed: {self.last_error}" if self.last_error else "")

    def start(self):
        if self.is_running():
            return
        self.last_error = ""
        worker = GeneOntologyWorker(self.manager.gene_ontology)
        self.worker = worker
        self._workers.append(worker)
        worker.progress.connect(lambda payload: self.progress_label.setText(payload.get("message", "Working...")))
        worker.succeeded.connect(self._success)
        worker.failed.connect(self._failure)
        worker.canceled.connect(self._canceled)
        worker.finished.connect(lambda w=worker: self._finished(w))
        worker.finished.connect(worker.deleteLater)
        worker.start()
        self.refresh()

    def _finished(self, worker):
        if worker in self._workers:
            self._workers.remove(worker)
        self.refresh()

    def _success(self, _snapshot):
        self.worker = None
        self.progress_label.setText("Gene Ontology snapshot is ready.")
        self.refresh()
        self.snapshot_changed.emit()

    def _failure(self, message):
        self.worker = None
        self.last_error = message
        self.progress_label.setText("Update failed; the previous snapshot remains active.")
        self.refresh()
        QMessageBox.warning(self, "Gene Ontology update failed", message)

    def _canceled(self, message):
        self.worker = None
        self.progress_label.setText(message)
        self.refresh()

    def cancel(self):
        if self.worker:
            self.worker.cancel()

    def open_folder(self):
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.manager.gene_ontology.root)))
