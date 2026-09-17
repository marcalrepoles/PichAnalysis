from __future__ import annotations

from PySide6.QtCore import QThread, QUrl, Signal
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import (QCheckBox, QDialog, QDialogButtonBox, QFormLayout, QGroupBox,
    QHBoxLayout, QLabel, QMessageBox, QProgressBar, QPushButton, QVBoxLayout, QWidget)

from ..core.database_manager import DatabaseManager
from ..core.database_registry import DatabaseState
from ..core.databases.kegg import KEGG_USAGE_URL

ACADEMIC_NOTICE = ("The KEGG REST API is provided for academic use by academic users. "
    "Non-academic use requires the appropriate KEGG license. Please review the official KEGG API terms before downloading.")

class DatabaseDownloadWorker(QThread):
    progress = Signal(dict); succeeded = Signal(str); failed = Signal(str)
    def __init__(self, manager: DatabaseManager, components: list[str], resume: bool = False):
        super().__init__(); self.manager=manager; self.components=components; self.resume=resume
    def run(self) -> None:
        try: self.succeeded.emit(str(self.manager.download(self.components, self.progress.emit, resume=self.resume)))
        except Exception as error: self.failed.emit(str(error))
    def cancel(self) -> None: self.manager.cancel()

class ComponentDialog(QDialog):
    def __init__(self, manager: DatabaseManager, parent=None):
        super().__init__(parent); self.setWindowTitle("Configure KEGG download")
        self.core=QCheckBox("Core tables (required)"); self.core.setChecked(True); self.core.setEnabled(False)
        self.entries=QCheckBox("Pathway entries"); self.kgml=QCheckBox("KGML pathway files"); self.images=QCheckBox("Pathway PNG images")
        _,_,free=manager.disk_info(); destination=QLabel(str(manager.database_root)); destination.setWordWrap(True)
        form=QFormLayout(); form.addRow("Destination",destination); form.addRow("Free disk space",QLabel(f"{free/(1024**3):.1f} GiB"))
        layout=QVBoxLayout(self); layout.addLayout(form); group=QGroupBox("Components"); gl=QVBoxLayout(group)
        for box in (self.core,self.entries,self.kgml,self.images): gl.addWidget(box)
        layout.addWidget(group); buttons=QDialogButtonBox(QDialogButtonBox.StandardButton.Ok|QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.accept); buttons.rejected.connect(self.reject); layout.addWidget(buttons)
    def components(self)->list[str]:
        return [name for name,box in (("entries",self.entries),("kgml",self.kgml),("images",self.images)) if box.isChecked()]

class DatabaseManagerPage(QWidget):
    def __init__(self, manager: DatabaseManager|None=None)->None:
        super().__init__(); self.manager=manager or DatabaseManager(); self.worker:DatabaseDownloadWorker|None=None
        title=QLabel("KEGG Homo sapiens"); title.setStyleSheet("font-size: 20px; font-weight: 600")
        self.state_label=QLabel(); self.details=QLabel(); self.details.setWordWrap(True)
        self.progress=QProgressBar(); self.progress.setRange(0,1); self.progress_text=QLabel("No download is running.")
        self.download_button=QPushButton("Download"); self.update_button=QPushButton("Update"); self.resume_button=QPushButton("Resume")
        self.cancel_button=QPushButton("Cancel"); self.folder_button=QPushButton("Open database folder")
        actions=QHBoxLayout()
        for button in (self.download_button,self.update_button,self.resume_button,self.cancel_button,self.folder_button): actions.addWidget(button)
        layout=QVBoxLayout(self)
        for widget in (title,self.state_label,self.details,self.progress,self.progress_text): layout.addWidget(widget)
        layout.addLayout(actions); layout.addStretch()
        self.download_button.clicked.connect(self._configure); self.update_button.clicked.connect(self._configure)
        self.resume_button.clicked.connect(lambda:self._start([],True)); self.cancel_button.clicked.connect(self.cancel); self.folder_button.clicked.connect(self._open_folder)
        self.refresh()
    def refresh(self)->None:
        state=self.manager.state(); active=self.manager.active_snapshot(); manifest=self.manager.manifest(active) if active else {}
        self.state_label.setText(f"Status: {state}")
        self.details.setText((f"Active snapshot: {active.name}\nLocation: {active}\nPathways: {manifest.get('counts',{}).get('pathways','—')}" if active
            else f"Location: {self.manager.database_root}\nNo complete snapshot is active."))
        running=bool(self.worker and self.worker.isRunning()); self.download_button.setEnabled(not running and state==DatabaseState.NOT_INSTALLED)
        self.update_button.setEnabled(not running and active is not None); self.resume_button.setEnabled(not running and self.manager.latest_incomplete_snapshot() is not None); self.cancel_button.setEnabled(running)
    def _configure(self)->None:
        if self.manager.active_snapshot() is None:
            notice=QMessageBox(self); notice.setWindowTitle("KEGG academic-use notice"); notice.setText(ACADEMIC_NOTICE); notice.setInformativeText(KEGG_USAGE_URL)
            notice.setStandardButtons(QMessageBox.StandardButton.Ok|QMessageBox.StandardButton.Cancel)
            if notice.exec()!=QMessageBox.StandardButton.Ok:return
        dialog=ComponentDialog(self.manager,self)
        if dialog.exec()==QDialog.DialogCode.Accepted:self._start(dialog.components(),False)
    def _start(self,components:list[str],resume:bool)->None:
        if self.worker and self.worker.isRunning():return
        self.worker=DatabaseDownloadWorker(self.manager,components,resume); self.worker.progress.connect(self._show_progress)
        self.worker.succeeded.connect(self._finished); self.worker.failed.connect(self._failed); self.worker.finished.connect(self.worker.deleteLater); self.worker.start(); self.refresh()
    def _show_progress(self,payload:dict)->None:
        total=max(int(payload.get("total",1)),1); self.progress.setRange(0,total); self.progress.setValue(int(payload.get("done",0)))
        self.progress_text.setText(str(payload.get("message",""))); self.state_label.setText(f"Status: {payload.get('status',DatabaseState.DOWNLOADING)}")
    def _finished(self,_snapshot:str)->None:self.worker=None;self.refresh()
    def _failed(self,message:str)->None:self.worker=None;self.refresh();QMessageBox.warning(self,"Database download failed",message)
    def cancel(self)->None:
        if self.worker:self.progress_text.setText("Cancellation requested. The current safe operation will finish first.");self.worker.cancel()
    def is_running(self)->bool:return bool(self.worker and self.worker.isRunning())
    def _open_folder(self)->None:
        self.manager.database_root.mkdir(parents=True,exist_ok=True);QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.manager.database_root)))
