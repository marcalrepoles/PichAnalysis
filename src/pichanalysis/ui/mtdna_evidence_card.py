"""Database Manager card for derived human mtDNA evidence snapshots."""
from __future__ import annotations

import threading

from PySide6.QtCore import QThread, QUrl, Signal
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import QGroupBox, QHBoxLayout, QLabel, QMessageBox, QProgressBar, QPushButton, QVBoxLayout

from ..core.mtdna_evidence_database import MtdnaBuildCancelled


class MtdnaEvidenceWorker(QThread):
    progress=Signal(dict)
    succeeded=Signal(str)
    failed=Signal(str)
    canceled=Signal(str)

    def __init__(self,manager,ncbi_xml=None):
        super().__init__();self.manager=manager;self.ncbi_xml=ncbi_xml;self._cancel=threading.Event()

    def run(self):
        try:
            snapshot=self.manager.mtdna_evidence.build(mitocarta=self.manager.mitocarta,
                ncbi_xml=self.ncbi_xml,progress=self.progress.emit,cancel_requested=self._cancel.is_set)
            self.succeeded.emit(str(snapshot))
        except MtdnaBuildCancelled as error:self.canceled.emit(str(error))
        except Exception as error:self.failed.emit(str(error))

    def cancel(self):self._cancel.set()


class MtdnaEvidenceCard(QGroupBox):
    def __init__(self,manager):
        super().__init__("mtDNA Evidence — Homo sapiens")
        self.manager=manager;self.worker=None;self._workers=[];self.last_error=""
        self.state_label=QLabel();self.details=QLabel();self.details.setWordWrap(True)
        self.progress=QProgressBar();self.progress.setRange(0,1)
        self.progress_text=QLabel("No build is running.");self.progress_text.setWordWrap(True)
        self.error_label=QLabel();self.error_label.setWordWrap(True)
        self.primary_button=QPushButton("Build");self.cancel_button=QPushButton("Cancel")
        self.folder_button=QPushButton("Open database folder")
        layout=QVBoxLayout(self)
        for widget in (self.state_label,self.details,self.progress,self.progress_text,self.error_label):layout.addWidget(widget)
        actions=QHBoxLayout()
        for widget in (self.primary_button,self.cancel_button,self.folder_button):actions.addWidget(widget)
        layout.addLayout(actions)
        self.primary_button.clicked.connect(self.start);self.cancel_button.clicked.connect(self.cancel)
        self.folder_button.clicked.connect(self.open_folder)
        self.refresh()

    def is_running(self):return bool(self._workers)

    def refresh(self):
        db=self.manager.mtdna_evidence;active=db.active_snapshot();manifest=db.manifest(active) if active else {}
        dependencies=db.dependencies(self.manager.mitocarta)
        running=self.is_running()
        self.state_label.setText(f"Status: {'Building' if running else 'Ready' if active else 'Not installed'}")
        missing=[]
        if not dependencies["mitocarta_ready"]:missing.append("MitoCarta database is required.")
        if not dependencies["go_ready"]:missing.append("Gene Ontology database is required.")
        counts=manifest.get("counts",{})
        self.details.setText("\n".join([
            f"Snapshot: {active.name if active else '—'}",
            f"MitoCarta source: {'Ready' if dependencies['mitocarta_ready'] else 'Not installed'}",
            f"MitoCarta snapshot: {dependencies['mitocarta_snapshot'].name if dependencies['mitocarta_snapshot'] else '—'}",
            f"Gene Ontology source: {'Ready' if dependencies['go_ready'] else 'Not installed'}",
            f"GO snapshot: {dependencies['go_snapshot'].name if dependencies['go_snapshot'] else '—'}",
            f"NCBI mtDNA reference: {manifest.get('ncbi_accession','NC_012920.1')}",
            f"Evidence entities: {counts.get('entities','—')}",f"Evidence records: {counts.get('evidence_records','—')}",
            f"MitoCarta-supported genes: {counts.get('mitocarta_supported_genes','—')}",
            f"GO-supported genes: {counts.get('go_supported_genes','—')}",
            f"mtDNA-encoded proteins: {counts.get('mtdna_encoded_genes','—')}",
            f"Categories: {counts.get('categories','—')}",*missing]))
        self.primary_button.setText("Update" if active else "Build")
        self.primary_button.setEnabled(not running and not missing)
        self.cancel_button.setEnabled(running)
        self.error_label.setText(f"Last attempt failed: {self.last_error}" if self.last_error else "")

    def start(self):
        if self.is_running():return
        deps=self.manager.mtdna_evidence.dependencies(self.manager.mitocarta)
        if not deps["mitocarta_ready"] or not deps["go_ready"]:self.refresh();return
        self.last_error="";worker=MtdnaEvidenceWorker(self.manager)
        self.worker=worker;self._workers.append(worker)
        worker.progress.connect(lambda payload:self.progress_text.setText(payload.get("message","Building evidence index...")))
        worker.succeeded.connect(self.finished_ok);worker.failed.connect(self.failed);worker.canceled.connect(self.canceled)
        worker.finished.connect(lambda w=worker:self._worker_finished(w));worker.finished.connect(worker.deleteLater)
        worker.start();self.progress.setRange(0,0);self.refresh()

    def _worker_finished(self,worker):
        if worker in self._workers:self._workers.remove(worker)
        self.refresh()

    def finished_ok(self,_snapshot):
        self.worker=None;self.progress.setRange(0,1);self.progress.setValue(1)
        self.progress_text.setText("mtDNA evidence index is ready.");self.refresh()

    def failed(self,message):
        self.worker=None;self.last_error=message;self.progress.setRange(0,1);self.progress.setValue(0)
        self.progress_text.setText("mtDNA evidence build failed; the previous snapshot remains active.")
        self.refresh();QMessageBox.warning(self,"mtDNA evidence build failed",message)

    def canceled(self,message):
        self.worker=None;self.progress.setRange(0,1);self.progress.setValue(0)
        self.progress_text.setText(message);self.refresh()

    def cancel(self):
        if self.worker:self.worker.cancel()

    def open_folder(self):
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.manager.mtdna_evidence.root)))
