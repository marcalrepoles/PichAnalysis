from __future__ import annotations
import threading
from PySide6.QtCore import QThread, QUrl, Signal
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import QGroupBox, QHBoxLayout, QLabel, QMessageBox, QProgressBar, QPushButton, QVBoxLayout
from ..core.database_registry import DatabaseState
from ..core.databases.string import StringDownloadCancelled, StringProvider

class StringDownloadWorker(QThread):
    progress=Signal(dict);succeeded=Signal(str);failed=Signal(str);canceled=Signal(str)
    def __init__(self,manager,provider: StringProvider|None=None):
        super().__init__();self.manager=manager;self._cancel=threading.Event()
        if provider is not None:self.manager.string.provider=provider
    def run(self):
        try:self.succeeded.emit(str(self.manager.string.download_and_build(self.progress.emit,self._cancel.is_set)))
        except StringDownloadCancelled as error:self.canceled.emit(str(error))
        except Exception as error:self.failed.emit(str(error))
    def cancel(self):self._cancel.set()

class StringDatabaseCard(QGroupBox):
    def __init__(self,manager):
        super().__init__("STRING — Homo sapiens");self.manager=manager;self.worker=None;self.last_error="";layout=QVBoxLayout(self)
        self.state_label=QLabel();self.details=QLabel();self.details.setWordWrap(True);self.progress=QProgressBar();self.progress.setRange(0,1);self.progress_text=QLabel("No download is running.");self.progress_text.setWordWrap(True);self.error_label=QLabel();self.error_label.setWordWrap(True)
        self.primary_button=QPushButton("Download");self.cancel_button=QPushButton("Cancel");self.folder_button=QPushButton("Open database folder");actions=QHBoxLayout()
        for button in (self.primary_button,self.cancel_button,self.folder_button):actions.addWidget(button)
        for widget in (self.state_label,self.details,self.progress,self.progress_text,self.error_label):layout.addWidget(widget)
        layout.addLayout(actions);self.primary_button.clicked.connect(self.start);self.cancel_button.clicked.connect(self.cancel);self.folder_button.clicked.connect(self.open_folder);self.refresh()
    def refresh(self):
        db=self.manager.string;active=db.active_snapshot();m=db.manifest(active) if active else {};running=self.is_running();state=DatabaseState.DOWNLOADING if running else db.state();self.state_label.setText(f"Status: {state}")
        self.details.setText(f"STRING version: {m.get('string_version','Unknown')}\nSnapshot: {active.name if active else '—'}\nDownloaded: {m.get('completed_at','—')}\nFunctional network: {m.get('functional_edge_count','—')} edges\nPhysical network: {m.get('physical_edge_count','—')} edges\nProteins: {m.get('protein_count','—')}\nAliases: {m.get('alias_count','—')}\nLicense: CC BY 4.0\nLocation: {db.root}")
        self.primary_button.setText("Update" if active else "Download");self.primary_button.setEnabled(not running);self.cancel_button.setEnabled(running);self.error_label.setText(f"Last attempt failed: {self.last_error}" if self.last_error else "")
    def start(self):
        if self.is_running():return
        self.last_error="";self.worker=StringDownloadWorker(self.manager);self.worker.progress.connect(self.show_progress);self.worker.succeeded.connect(self.finished_ok);self.worker.failed.connect(self.failed);self.worker.canceled.connect(self.canceled);self.worker.finished.connect(self.worker.deleteLater);self.worker.start();self.refresh()
    def show_progress(self,payload):
        if payload.get("stage")=="Downloading":
            done,total=int(payload.get("bytes_downloaded",0)),int(payload.get("bytes_total",0))
            if total:self.progress.setRange(0,total);self.progress.setValue(done);size=f"{done:,} / {total:,} bytes"
            else:self.progress.setRange(0,0);size=f"{done:,} bytes" if done else "Size unknown"
            self.progress_text.setText(f"Downloading STRING — Homo sapiens\nFile {payload.get('file_index')} of {payload.get('file_total')}: {payload.get('filename')}\n{size}")
        else:self.progress.setRange(0,0);self.progress_text.setText(str(payload.get("message","Indexing STRING database...")))
    def finished_ok(self,_):self.worker=None;self.progress.setRange(0,1);self.progress.setValue(1);self.progress_text.setText("STRING database is ready.");self.refresh()
    def failed(self,message):self.worker=None;self.last_error=message;self.progress_text.setText("STRING update failed; the previous snapshot remains active.");self.refresh();QMessageBox.warning(self,"STRING database download failed",message)
    def canceled(self,message):self.worker=None;self.progress.setRange(0,1);self.progress.setValue(0);self.progress_text.setText(message);self.refresh()
    def cancel(self):
        if self.worker:self.progress_text.setText("Cancellation requested. The current safe operation will stop first.");self.worker.cancel()
    def is_running(self):return bool(self.worker and self.worker.isRunning())
    def open_folder(self):
        self.manager.string.root.mkdir(parents=True,exist_ok=True)
        if not QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.manager.string.root))):QMessageBox.warning(self,"Database folder","Could not open the database folder.")
