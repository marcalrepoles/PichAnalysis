from __future__ import annotations

import threading
from datetime import datetime

from PySide6.QtCore import QThread, QUrl, Signal
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import (
    QCheckBox, QDialog, QDialogButtonBox, QFormLayout, QGroupBox, QHBoxLayout,
    QLabel, QMessageBox, QProgressBar, QPushButton, QVBoxLayout, QWidget,
)

from ..core.database_manager import DatabaseManager
from ..core.database_registry import DatabaseState
from ..core.databases.kegg import KEGG_USAGE_URL
from ..core.databases.reactome import CORE_FILES, DIAGRAM_ARCHIVE, ReactomeDownloadCancelled, ReactomeProvider
from ..core.databases.mitocarta import FILES as MITOCARTA_FILES, MitoCartaDownloadCancelled, MitoCartaProvider
from ..core.interpro_database import InterProCancelled, InterProProvider
from .string_database_card import StringDatabaseCard, StringDownloadWorker

ACADEMIC_NOTICE = (
    "The KEGG REST API is provided for academic use by academic users. "
    "Non-academic use requires the appropriate KEGG license. Please review the official KEGG API terms before downloading."
)


class DatabaseDownloadWorker(QThread):
    progress = Signal(dict)
    succeeded = Signal(str)
    failed = Signal(str)

    def __init__(self, manager: DatabaseManager, components: list[str], resume: bool = False):
        super().__init__()
        self.manager = manager
        self.components = components
        self.resume = resume

    def run(self) -> None:
        try:
            self.succeeded.emit(str(self.manager.download(self.components, self.progress.emit, resume=self.resume)))
        except Exception as error:
            self.failed.emit(str(error))

    def cancel(self) -> None:
        self.manager.cancel()


class ReactomeCoreDownloadWorker(QThread):
    progress = Signal(dict)
    succeeded = Signal(str)
    failed = Signal(str)
    canceled = Signal(str)

    def __init__(self, manager: DatabaseManager, provider: ReactomeProvider | None = None):
        super().__init__()
        self.manager = manager
        self.provider = provider or ReactomeProvider()
        self._cancel = threading.Event()

    def run(self) -> None:
        snapshot = None
        try:
            snapshot = self.manager.reactome.create_staging_snapshot()
            self.manager.reactome.mark_downloading(snapshot)
            raw = snapshot / "raw"
            for index, name in enumerate(CORE_FILES, 1):
                if self._cancel.is_set():
                    raise ReactomeDownloadCancelled("Reactome Core Data download was canceled.")
                self.progress.emit(self._payload(index, name, 0, 0))
                self.provider.download(
                    name,
                    raw / name,
                    progress=lambda filename, done, total, i=index: self.progress.emit(
                        self._payload(i, filename, done, total)
                    ),
                    cancel_requested=self._cancel.is_set,
                )
            completed = self.manager.reactome.finalize_snapshot(snapshot)
            self.succeeded.emit(str(completed))
        except ReactomeDownloadCancelled as error:
            if snapshot is not None:
                try:
                    self.manager.reactome.mark_incomplete(snapshot)
                except Exception:
                    self.manager.logger.exception("Could not record the incomplete Reactome snapshot")
            self.canceled.emit(str(error))
        except Exception as error:
            if snapshot is not None and self.manager.reactome.manifest(snapshot).get("status") != DatabaseState.ERROR:
                try:
                    self.manager.reactome.mark_error(snapshot, str(error))
                except Exception:
                    self.manager.logger.exception("Could not record the failed Reactome snapshot")
            self.manager.logger.exception("Reactome Core Data download failed")
            self.failed.emit(str(error))

    @staticmethod
    def _payload(index: int, filename: str, done: int, total: int) -> dict:
        return {
            "status": DatabaseState.DOWNLOADING,
            "file_index": index,
            "file_total": len(CORE_FILES),
            "filename": filename,
            "bytes_downloaded": done,
            "bytes_total": total,
        }

    def cancel(self) -> None:
        self._cancel.set()


class ReactomeDiagramDownloadWorker(QThread):
    progress = Signal(dict)
    succeeded = Signal(str)
    failed = Signal(str)
    canceled = Signal(str)

    def __init__(self, manager: DatabaseManager, provider: ReactomeProvider | None = None):
        super().__init__(); self.manager = manager; self.provider = provider or ReactomeProvider(); self._cancel = threading.Event()

    def run(self) -> None:
        snapshot = None
        try:
            snapshot = self.manager.reactome.create_diagram_staging()
            self.manager.reactome.mark_downloading(snapshot, "diagrams")
            archive = snapshot / "archives" / DIAGRAM_ARCHIVE
            info = self.provider.download(DIAGRAM_ARCHIVE, archive,
                progress=lambda name, done, total: self.progress.emit({"stage":"download","filename":name,"bytes_downloaded":done,"bytes_total":total}),
                cancel_requested=self._cancel.is_set)
            completed = self.manager.reactome.install_diagrams(snapshot, archive, info,
                progress=lambda done,total,name:self.progress.emit({"stage":"index" if name.startswith("__INDEX") else "extract","filename":name,"items_done":done,"items_total":total}),
                cancel_requested=self._cancel.is_set)
            self.succeeded.emit(str(completed))
        except ReactomeDownloadCancelled as error:
            if snapshot is not None: self.manager.reactome.mark_incomplete(snapshot, "diagrams")
            self.canceled.emit(str(error))
        except Exception as error:
            if snapshot is not None and self.manager.reactome.manifest(snapshot).get("status") != DatabaseState.ERROR:
                self.manager.reactome.mark_component_error(snapshot, "diagrams", str(error))
            self.manager.logger.exception("Reactome diagram download failed"); self.failed.emit(str(error))

    def cancel(self) -> None: self._cancel.set()


class MitoCartaDownloadWorker(QThread):
    progress = Signal(dict); succeeded = Signal(str); failed = Signal(str); canceled = Signal(str)
    def __init__(self, manager: DatabaseManager, provider: MitoCartaProvider | None = None):
        super().__init__(); self.manager=manager; self.provider=provider or MitoCartaProvider(); self._cancel=threading.Event()
    def run(self) -> None:
        snapshot=None
        try:
            snapshot=self.manager.mitocarta.create_staging_snapshot();self.manager.mitocarta.mark_downloading(snapshot);information=[]
            for index,name in enumerate(MITOCARTA_FILES,1):
                if self._cancel.is_set():raise MitoCartaDownloadCancelled("MitoCarta3.0 download was canceled.")
                self.progress.emit(self._payload(index,name,0,0,"download"))
                information.append(self.provider.download(name,snapshot/"raw"/name,
                    progress=lambda filename,done,total,i=index:self.progress.emit(self._payload(i,filename,done,total,"download")),cancel_requested=self._cancel.is_set))
            if self._cancel.is_set():raise MitoCartaDownloadCancelled("MitoCarta3.0 download was canceled.")
            self.progress.emit(self._payload(2,"Validating and normalizing MitoCarta3.0",0,0,"processing"))
            completed=self.manager.mitocarta.finalize_snapshot(snapshot,information);self.succeeded.emit(str(completed))
        except MitoCartaDownloadCancelled as error:
            if snapshot is not None:self.manager.mitocarta.mark_incomplete(snapshot)
            self.canceled.emit(str(error))
        except Exception as error:
            if snapshot is not None and self.manager.mitocarta.manifest(snapshot).get("status")!=DatabaseState.ERROR:self.manager.mitocarta.mark_error(snapshot,str(error))
            self.manager.logger.exception("MitoCarta3.0 download failed");self.failed.emit(str(error))
    @staticmethod
    def _payload(index,name,done,total,stage):return {"stage":stage,"file_index":index,"file_total":len(MITOCARTA_FILES),"filename":name,"bytes_downloaded":done,"bytes_total":total}
    def cancel(self):self._cancel.set()


class InterProMetadataWorker(QThread):
    progress=Signal(dict);succeeded=Signal(str);failed=Signal(str);canceled=Signal(str)
    def __init__(self,manager:DatabaseManager,provider:InterProProvider|None=None):
        super().__init__();self.manager=manager;self._cancel=threading.Event()
        if provider is not None:self.manager.interpro.provider=provider
    def run(self):
        try:
            self.progress.emit({"message":"Downloading official InterPro/Pfam metadata"});path=self.manager.interpro.download_metadata(self._cancel.is_set);self.succeeded.emit(str(path))
        except InterProCancelled as error:self.canceled.emit(str(error))
        except Exception as error:self.failed.emit(str(error))
    def cancel(self):self._cancel.set()


class InterProAnnotationWorker(QThread):
    progress=Signal(dict);succeeded=Signal(str);failed=Signal(str);canceled=Signal(str)
    def __init__(self,database,project_path,accessions,annotation_set_id=None):
        super().__init__();self.database=database;self.project_path=project_path;self.accessions=tuple(accessions);self.annotation_set_id=annotation_set_id;self._cancel=threading.Event()
    def run(self):
        try:self.succeeded.emit(str(self.database.build_annotation_set(self.project_path,self.accessions,annotation_set_id=self.annotation_set_id,progress=self.progress.emit,cancel_requested=self._cancel.is_set)))
        except InterProCancelled as error:self.canceled.emit(str(error))
        except Exception as error:self.failed.emit(str(error))
    def cancel(self):self._cancel.set()


class ComponentDialog(QDialog):
    def __init__(self, manager: DatabaseManager, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Configure KEGG download")
        self.core = QCheckBox("Core tables (required)")
        self.core.setChecked(True)
        self.core.setEnabled(False)
        self.entries = QCheckBox("Pathway entries")
        self.kgml = QCheckBox("KGML pathway files")
        self.images = QCheckBox("Pathway PNG images")
        _, _, free = manager.disk_info()
        destination = QLabel(str(manager.database_root))
        destination.setWordWrap(True)
        form = QFormLayout()
        form.addRow("Destination", destination)
        form.addRow("Free disk space", QLabel(f"{free / (1024 ** 3):.1f} GiB"))
        layout = QVBoxLayout(self)
        layout.addLayout(form)
        group = QGroupBox("Components")
        group_layout = QVBoxLayout(group)
        for box in (self.core, self.entries, self.kgml, self.images):
            group_layout.addWidget(box)
        layout.addWidget(group)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def components(self) -> list[str]:
        return [
            name for name, box in (("entries", self.entries), ("kgml", self.kgml), ("images", self.images))
            if box.isChecked()
        ]


class DatabaseManagerPage(QWidget):
    def __init__(self, manager: DatabaseManager | None = None) -> None:
        super().__init__()
        self.manager = manager or DatabaseManager()
        self.worker: DatabaseDownloadWorker | None = None
        self.reactome_worker: ReactomeCoreDownloadWorker | None = None
        self.reactome_diagram_worker: ReactomeDiagramDownloadWorker | None = None
        self._reactome_attempt_error = ""
        self.mitocarta_worker: MitoCartaDownloadWorker | None = None
        self._mitocarta_attempt_error = ""
        self.interpro_worker: InterProMetadataWorker | None = None
        self._interpro_attempt_error = ""


        layout = QVBoxLayout(self)
        layout.addWidget(self._build_kegg_card())
        layout.addWidget(self._build_reactome_card())
        layout.addWidget(self._build_mitocarta_card())
        layout.addWidget(self._build_interpro_card())
        self.string_card = StringDatabaseCard(self.manager)
        layout.addWidget(self.string_card)
        layout.addStretch()
        self.refresh()

    def _build_kegg_card(self) -> QGroupBox:
        card = QGroupBox("KEGG Homo sapiens")
        layout = QVBoxLayout(card)
        self.state_label = QLabel()
        self.details = QLabel()
        self.details.setWordWrap(True)
        self.progress = QProgressBar()
        self.progress.setRange(0, 1)
        self.progress_text = QLabel("No download is running.")
        self.download_button = QPushButton("Download")
        self.update_button = QPushButton("Update")
        self.resume_button = QPushButton("Resume")
        self.cancel_button = QPushButton("Cancel")
        self.folder_button = QPushButton("Open database folder")
        actions = QHBoxLayout()
        for button in (self.download_button, self.update_button, self.resume_button, self.cancel_button, self.folder_button):
            actions.addWidget(button)
        for widget in (self.state_label, self.details, self.progress, self.progress_text):
            layout.addWidget(widget)
        layout.addLayout(actions)
        self.download_button.clicked.connect(self._configure)
        self.update_button.clicked.connect(self._configure)
        self.resume_button.clicked.connect(lambda: self._start([], True))
        self.cancel_button.clicked.connect(self.cancel)
        self.folder_button.clicked.connect(self._open_folder)
        return card

    def _build_reactome_card(self) -> QGroupBox:
        card = QGroupBox("Reactome — Homo sapiens")
        layout = QVBoxLayout(card)
        self.reactome_state_label = QLabel()
        self.reactome_diagram_state_label = QLabel()
        self.reactome_details = QLabel()
        self.reactome_details.setWordWrap(True)
        self.reactome_progress = QProgressBar()
        self.reactome_progress.setRange(0, 1)
        self.reactome_progress_text = QLabel("No download is running.")
        self.reactome_error_label = QLabel()
        self.reactome_error_label.setWordWrap(True)
        self.reactome_primary_button = QPushButton("Download")
        self.reactome_diagram_button = QPushButton("Download diagrams")
        self.reactome_cancel_button = QPushButton("Cancel")
        self.reactome_folder_button = QPushButton("Open database folder")
        actions = QHBoxLayout()
        for button in (self.reactome_primary_button, self.reactome_diagram_button, self.reactome_cancel_button, self.reactome_folder_button):
            actions.addWidget(button)
        for widget in (
            self.reactome_state_label, self.reactome_diagram_state_label, self.reactome_details, self.reactome_progress,
            self.reactome_progress_text, self.reactome_error_label,
        ):
            layout.addWidget(widget)
        layout.addLayout(actions)
        self.reactome_primary_button.clicked.connect(self._start_reactome)
        self.reactome_diagram_button.clicked.connect(self._start_reactome_diagrams)
        self.reactome_cancel_button.clicked.connect(self.cancel_reactome)
        self.reactome_folder_button.clicked.connect(self._open_reactome_folder)
        return card

    def _build_mitocarta_card(self) -> QGroupBox:
        card=QGroupBox("MitoCarta3.0 — Homo sapiens");layout=QVBoxLayout(card)
        self.mitocarta_state_label=QLabel();self.mitocarta_details=QLabel();self.mitocarta_details.setWordWrap(True)
        self.mitocarta_progress=QProgressBar();self.mitocarta_progress.setRange(0,1);self.mitocarta_progress_text=QLabel("No download is running.")
        self.mitocarta_error_label=QLabel();self.mitocarta_error_label.setWordWrap(True)
        self.mitocarta_primary_button=QPushButton("Download");self.mitocarta_cancel_button=QPushButton("Cancel");self.mitocarta_folder_button=QPushButton("Open database folder")
        actions=QHBoxLayout()
        for button in (self.mitocarta_primary_button,self.mitocarta_cancel_button,self.mitocarta_folder_button):actions.addWidget(button)
        for widget in (self.mitocarta_state_label,self.mitocarta_details,self.mitocarta_progress,self.mitocarta_progress_text,self.mitocarta_error_label):layout.addWidget(widget)
        layout.addLayout(actions);self.mitocarta_primary_button.clicked.connect(self._start_mitocarta);self.mitocarta_cancel_button.clicked.connect(self.cancel_mitocarta);self.mitocarta_folder_button.clicked.connect(self._open_mitocarta_folder);return card

    def _build_interpro_card(self)->QGroupBox:
        card=QGroupBox("InterPro / Pfam");layout=QVBoxLayout(card);self.interpro_state_label=QLabel();self.interpro_details=QLabel();self.interpro_details.setWordWrap(True);self.interpro_progress=QProgressBar();self.interpro_progress.setRange(0,1);self.interpro_progress_text=QLabel("No download is running.");self.interpro_error_label=QLabel();self.interpro_error_label.setWordWrap(True)
        self.interpro_primary_button=QPushButton("Download metadata");self.interpro_cancel_button=QPushButton("Cancel");self.interpro_folder_button=QPushButton("Open database folder");actions=QHBoxLayout()
        for button in (self.interpro_primary_button,self.interpro_cancel_button,self.interpro_folder_button):actions.addWidget(button)
        for widget in (self.interpro_state_label,self.interpro_details,self.interpro_progress,self.interpro_progress_text,self.interpro_error_label):layout.addWidget(widget)
        layout.addLayout(actions);self.interpro_primary_button.clicked.connect(self._start_interpro);self.interpro_cancel_button.clicked.connect(self.cancel_interpro);self.interpro_folder_button.clicked.connect(self._open_interpro_folder);return card

    def refresh(self) -> None:
        state = self.manager.state()
        active = self.manager.active_snapshot()
        manifest = self.manager.manifest(active) if active else {}
        self.state_label.setText(f"Status: {state}")
        self.details.setText(
            f"Active snapshot: {active.name}\nLocation: {active}\nPathways: {manifest.get('counts', {}).get('pathways', '—')}"
            if active else f"Location: {self.manager.database_root}\nNo complete snapshot is active."
        )
        running = bool(self.worker and self.worker.isRunning())
        self.download_button.setEnabled(not running and state == DatabaseState.NOT_INSTALLED)
        self.update_button.setEnabled(not running and active is not None)
        self.resume_button.setEnabled(not running and self.manager.latest_incomplete_snapshot() is not None)
        self.cancel_button.setEnabled(running)
        self._refresh_reactome()
        self._refresh_mitocarta()
        self._refresh_interpro()
        self.string_card.refresh()

    def _refresh_interpro(self):
        database=self.manager.interpro;manifest=database.metadata_manifest();running=bool(self.interpro_worker and self.interpro_worker.isRunning());state=DatabaseState.DOWNLOADING if running else database.state()
        self.interpro_state_label.setText(f"Status: {state}");self.interpro_details.setText(f"InterPro release: {manifest.get('interpro_release','Unknown')}\nPfam release: {manifest.get('pfam_release','Unknown')}\nMetadata: {'Ready' if manifest else 'Not installed'}\nCached proteins: {database.cached_protein_count()}\nDownloaded: {self._format_timestamp(manifest.get('completed_at'))}\nLocation: {database.root}")
        self.interpro_primary_button.setText("Update metadata" if manifest else "Download metadata");self.interpro_primary_button.setEnabled(not running);self.interpro_cancel_button.setEnabled(running);self.interpro_error_label.setText(f"Last attempt failed: {self._interpro_attempt_error}" if self._interpro_attempt_error else "")

    def _refresh_mitocarta(self) -> None:
        database=self.manager.mitocarta;active=database.active_snapshot();manifest=database.manifest(active) if active else {};running=bool(self.mitocarta_worker and self.mitocarta_worker.isRunning());state=DatabaseState.DOWNLOADING if running else database.state()
        self.mitocarta_state_label.setText(f"Status: {state}");self.mitocarta_details.setText(f"Version: {manifest.get('version','3.0')}\nSnapshot: {active.name if active else '—'}\nDownloaded: {self._format_timestamp(manifest.get('download_completed_at')) if active else '—'}\nCore Data: {'Ready' if active else 'Not installed'}\nLocation: {database.root}")
        self.mitocarta_primary_button.setText("Update" if active else "Download");self.mitocarta_primary_button.setEnabled(not running);self.mitocarta_cancel_button.setEnabled(running);self.mitocarta_error_label.setText(f"Last attempt failed: {self._mitocarta_attempt_error}" if self._mitocarta_attempt_error else "")

    def _refresh_reactome(self) -> None:
        database = self.manager.reactome
        active = database.active_snapshot()
        manifest = database.manifest(active) if active else {}
        core_running = bool(self.reactome_worker and self.reactome_worker.isRunning())
        diagram_running = bool(self.reactome_diagram_worker and self.reactome_diagram_worker.isRunning())
        running = core_running or diagram_running
        state = DatabaseState.DOWNLOADING if running else database.state()
        diagram_info = database.diagram_manifest()
        diagram_state = "Downloading" if diagram_running else database.diagram_state()
        self.reactome_state_label.setText(f"Core Data: {state}")
        self.reactome_diagram_state_label.setText(f"Pathway Diagrams: {diagram_state}")
        release = manifest.get("release_version")
        release_text = str(release) if release and str(release).lower() != "unknown" else "Unknown"
        snapshot_text = active.name if active else "—"
        downloaded = self._format_timestamp(manifest.get("download_completed_at")) if active else "—"
        self.reactome_details.setText(
            f"Release: {release_text}\nSnapshot: {snapshot_text}\nDownloaded: {downloaded}\nLocation: {database.root}"
        )
        self.reactome_primary_button.setText("Update" if active else "Download")
        self.reactome_primary_button.setEnabled(not running)
        self.reactome_diagram_button.setText("Update diagrams" if database.is_diagrams_ready() else "Download diagrams")
        self.reactome_diagram_button.setEnabled(not running and active is not None)
        self.reactome_cancel_button.setEnabled(running)
        self.reactome_error_label.setText(
            f"Last attempt failed: {self._reactome_attempt_error}" if self._reactome_attempt_error else ""
        )

    @staticmethod
    def _format_timestamp(value) -> str:
        if not value:
            return "—"
        try:
            return datetime.fromisoformat(str(value).replace("Z", "+00:00")).astimezone().strftime("%Y-%m-%d %H:%M:%S %Z")
        except (TypeError, ValueError):
            return str(value)

    def _configure(self) -> None:
        if self.manager.active_snapshot() is None:
            notice = QMessageBox(self)
            notice.setWindowTitle("KEGG academic-use notice")
            notice.setText(ACADEMIC_NOTICE)
            notice.setInformativeText(KEGG_USAGE_URL)
            notice.setStandardButtons(QMessageBox.StandardButton.Ok | QMessageBox.StandardButton.Cancel)
            if notice.exec() != QMessageBox.StandardButton.Ok:
                return
        dialog = ComponentDialog(self.manager, self)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            self._start(dialog.components(), False)

    def _start(self, components: list[str], resume: bool) -> None:
        if self.worker and self.worker.isRunning():
            return
        self.worker = DatabaseDownloadWorker(self.manager, components, resume)
        self.worker.progress.connect(self._show_progress)
        self.worker.succeeded.connect(self._finished)
        self.worker.failed.connect(self._failed)
        self.worker.finished.connect(self.worker.deleteLater)
        self.worker.start()
        self.refresh()

    def _show_progress(self, payload: dict) -> None:
        total = max(int(payload.get("total", 1)), 1)
        self.progress.setRange(0, total)
        self.progress.setValue(int(payload.get("done", 0)))
        self.progress_text.setText(str(payload.get("message", "")))
        self.state_label.setText(f"Status: {payload.get('status', DatabaseState.DOWNLOADING)}")

    def _finished(self, _snapshot: str) -> None:
        self.worker = None
        self.refresh()

    def _failed(self, message: str) -> None:
        self.worker = None
        self.refresh()
        QMessageBox.warning(self, "Database download failed", message)

    def cancel(self) -> None:
        if self.worker:
            self.progress_text.setText("Cancellation requested. The current safe operation will finish first.")
            self.worker.cancel()

    def _start_reactome(self) -> None:
        if self.reactome_worker and self.reactome_worker.isRunning():
            return
        self._reactome_attempt_error = ""
        self.reactome_worker = ReactomeCoreDownloadWorker(self.manager)
        self.reactome_worker.progress.connect(self._show_reactome_progress)
        self.reactome_worker.succeeded.connect(self._reactome_finished)
        self.reactome_worker.failed.connect(self._reactome_failed)
        self.reactome_worker.canceled.connect(self._reactome_canceled)
        self.reactome_worker.finished.connect(self.reactome_worker.deleteLater)
        self.reactome_worker.start()
        self._refresh_reactome()

    def _start_mitocarta(self) -> None:
        if self.mitocarta_worker and self.mitocarta_worker.isRunning():return
        self._mitocarta_attempt_error="";self.mitocarta_worker=MitoCartaDownloadWorker(self.manager);self.mitocarta_worker.progress.connect(self._show_mitocarta_progress);self.mitocarta_worker.succeeded.connect(self._mitocarta_finished);self.mitocarta_worker.failed.connect(self._mitocarta_failed);self.mitocarta_worker.canceled.connect(self._mitocarta_canceled);self.mitocarta_worker.finished.connect(self.mitocarta_worker.deleteLater);self.mitocarta_worker.start();self._refresh_mitocarta()

    def _show_mitocarta_progress(self,payload:dict)->None:
        if payload.get("stage")=="processing":self.mitocarta_progress.setRange(0,0);self.mitocarta_progress_text.setText(str(payload.get("filename")));return
        done,total=int(payload.get("bytes_downloaded",0)),int(payload.get("bytes_total",0))
        if total:self.mitocarta_progress.setRange(0,total);self.mitocarta_progress.setValue(done);size=f"{self._format_bytes(done)} / {self._format_bytes(total)}"
        else:self.mitocarta_progress.setRange(0,0);size=self._format_bytes(done) if done else "Size unknown"
        self.mitocarta_progress_text.setText(f"Downloading MitoCarta3.0\n{payload.get('filename')}\n{size}\nFile {payload.get('file_index')} of {payload.get('file_total')}")

    def _mitocarta_finished(self,_snapshot:str)->None:self.mitocarta_worker=None;self.mitocarta_progress.setRange(0,1);self.mitocarta_progress.setValue(1);self.mitocarta_progress_text.setText("MitoCarta3.0 Core Data is ready.");self.refresh()
    def _mitocarta_failed(self,message:str)->None:self.mitocarta_worker=None;self._mitocarta_attempt_error=message;self.mitocarta_progress_text.setText("MitoCarta3.0 update failed; the previous snapshot remains active.");self.refresh();QMessageBox.warning(self,"MitoCarta3.0 download failed",message)
    def _mitocarta_canceled(self,message:str)->None:self.mitocarta_worker=None;self.mitocarta_progress.setRange(0,1);self.mitocarta_progress.setValue(0);self.mitocarta_progress_text.setText(message);self.refresh()
    def cancel_mitocarta(self)->None:
        if self.mitocarta_worker:self.mitocarta_progress_text.setText("Cancellation requested. The current download will stop safely.");self.mitocarta_worker.cancel()

    def _start_interpro(self):
        if self.interpro_worker and self.interpro_worker.isRunning():return
        self._interpro_attempt_error="";self.interpro_worker=InterProMetadataWorker(self.manager);self.interpro_worker.progress.connect(lambda payload:self.interpro_progress_text.setText(str(payload.get("message",""))));self.interpro_worker.succeeded.connect(self._interpro_finished);self.interpro_worker.failed.connect(self._interpro_failed);self.interpro_worker.canceled.connect(self._interpro_canceled);self.interpro_worker.finished.connect(self.interpro_worker.deleteLater);self.interpro_worker.start();self._refresh_interpro()
    def _interpro_finished(self,_path):self.interpro_worker=None;self.interpro_progress.setValue(1);self.interpro_progress_text.setText("InterPro/Pfam metadata is ready.");self.refresh()
    def _interpro_failed(self,message):self.interpro_worker=None;self._interpro_attempt_error=message;self.interpro_progress_text.setText("InterPro/Pfam metadata update failed; the previous snapshot remains active.");self.refresh();QMessageBox.warning(self,"InterPro/Pfam metadata download failed",message)
    def _interpro_canceled(self,message):self.interpro_worker=None;self.interpro_progress.setValue(0);self.interpro_progress_text.setText(message);self.refresh()
    def cancel_interpro(self):
        if self.interpro_worker:self.interpro_progress_text.setText("Cancellation requested. The current request will stop safely.");self.interpro_worker.cancel()

    def _start_reactome_diagrams(self) -> None:
        if self.reactome_diagram_worker and self.reactome_diagram_worker.isRunning(): return
        self._reactome_attempt_error = ""; self.reactome_diagram_worker = ReactomeDiagramDownloadWorker(self.manager)
        self.reactome_diagram_worker.progress.connect(self._show_reactome_diagram_progress)
        self.reactome_diagram_worker.succeeded.connect(self._reactome_diagram_finished)
        self.reactome_diagram_worker.failed.connect(self._reactome_diagram_failed)
        self.reactome_diagram_worker.canceled.connect(self._reactome_diagram_canceled)
        self.reactome_diagram_worker.finished.connect(self.reactome_diagram_worker.deleteLater)
        self.reactome_diagram_worker.start(); self._refresh_reactome()

    def _show_reactome_diagram_progress(self, payload: dict) -> None:
        if payload.get("stage") == "download":
            done,total=int(payload.get("bytes_downloaded",0)),int(payload.get("bytes_total",0))
            if total:self.reactome_progress.setRange(0,total);self.reactome_progress.setValue(done)
            else:self.reactome_progress.setRange(0,0)
            self.reactome_progress_text.setText(f"Downloading Reactome Pathway Diagrams\n{payload.get('filename')}\n{self._format_bytes(done)}" + (f" / {self._format_bytes(total)}" if total else ""))
        elif payload.get("stage") == "extract":
            done,total=int(payload.get("items_done",0)),max(int(payload.get("items_total",1)),1);self.reactome_progress.setRange(0,total);self.reactome_progress.setValue(done);self.reactome_progress_text.setText(f"Validating and indexing diagrams\n{done} / {total}: {payload.get('filename','')}")
        else:
            self.reactome_progress.setRange(0,0);self.reactome_progress_text.setText("Indexing Reactome Pathway Diagrams")

    def _reactome_diagram_finished(self, _snapshot: str) -> None:
        self.reactome_diagram_worker=None;self.reactome_progress.setRange(0,1);self.reactome_progress.setValue(1);self.reactome_progress_text.setText("Reactome diagrams are ready.");self.refresh()

    def _reactome_diagram_failed(self, message: str) -> None:
        self.reactome_diagram_worker=None;self._reactome_attempt_error=message;self.reactome_progress_text.setText("Reactome diagram update failed; the previous snapshot remains active.");self.refresh();QMessageBox.warning(self,"Reactome diagram download failed",message)

    def _reactome_diagram_canceled(self, message: str) -> None:
        self.reactome_diagram_worker=None;self.reactome_progress.setRange(0,1);self.reactome_progress.setValue(0);self.reactome_progress_text.setText(message);self.refresh()

    def _show_reactome_progress(self, payload: dict) -> None:
        downloaded = int(payload.get("bytes_downloaded", 0))
        total = int(payload.get("bytes_total", 0))
        if total > 0:
            self.reactome_progress.setRange(0, total)
            self.reactome_progress.setValue(downloaded)
            byte_text = f"{self._format_bytes(downloaded)} / {self._format_bytes(total)}"
        else:
            self.reactome_progress.setRange(0, 0)
            byte_text = self._format_bytes(downloaded) if downloaded else "Size unknown"
        index = int(payload.get("file_index", 0))
        count = int(payload.get("file_total", len(CORE_FILES)))
        filename = str(payload.get("filename", ""))
        self.reactome_state_label.setText("Core Data: Downloading")
        self.reactome_progress_text.setText(
            f"Downloading Reactome Core Data\nFile {index} of {count}: {filename}\n{byte_text}"
        )

    @staticmethod
    def _format_bytes(value: int) -> str:
        size = float(value)
        for unit in ("B", "KB", "MB", "GB"):
            if size < 1024 or unit == "GB":
                return f"{size:.1f} {unit}"
            size /= 1024
        return f"{size:.1f} GB"

    def _reactome_finished(self, _snapshot: str) -> None:
        self.reactome_worker = None
        self.reactome_progress.setRange(0, 1)
        self.reactome_progress.setValue(1)
        self.reactome_progress_text.setText("Reactome Core Data is ready.")
        self.refresh()

    def _reactome_failed(self, message: str) -> None:
        self.reactome_worker = None
        self._reactome_attempt_error = message
        self.reactome_progress_text.setText("Reactome Core Data download failed.")
        self.refresh()
        QMessageBox.warning(self, "Reactome download failed", message)

    def _reactome_canceled(self, message: str) -> None:
        self.reactome_worker = None
        self.reactome_progress.setRange(0, 1)
        self.reactome_progress.setValue(0)
        self.reactome_progress_text.setText(message)
        self.refresh()

    def cancel_reactome(self) -> None:
        if self.reactome_worker:
            self.reactome_progress_text.setText("Cancellation requested. The current download will stop safely.")
            self.reactome_worker.cancel()
        if self.reactome_diagram_worker:
            self.reactome_progress_text.setText("Cancellation requested. The diagram update will stop safely.")
            self.reactome_diagram_worker.cancel()

    def is_running(self) -> bool:
        return bool(
            (self.worker and self.worker.isRunning())
            or (self.reactome_worker and self.reactome_worker.isRunning())
            or (self.reactome_diagram_worker and self.reactome_diagram_worker.isRunning())
            or (self.mitocarta_worker and self.mitocarta_worker.isRunning())
            or (self.interpro_worker and self.interpro_worker.isRunning())
            or self.string_card.is_running()
        )

    def _open_folder(self) -> None:
        self.manager.database_root.mkdir(parents=True, exist_ok=True)
        if not QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.manager.database_root))):
            QMessageBox.warning(self, "Database folder", "Could not open the database folder.")

    def _open_reactome_folder(self) -> None:
        self.manager.reactome.root.mkdir(parents=True, exist_ok=True)
        if not QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.manager.reactome.root))):
            QMessageBox.warning(self, "Database folder", "Could not open the database folder.")

    def _open_mitocarta_folder(self) -> None:
        self.manager.mitocarta.root.mkdir(parents=True,exist_ok=True)
        if not QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.manager.mitocarta.root))):QMessageBox.warning(self,"Database folder","Could not open the database folder.")

    def _open_interpro_folder(self)->None:
        self.manager.interpro.root.mkdir(parents=True,exist_ok=True)
        if not QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.manager.interpro.root))):QMessageBox.warning(self,"Database folder","Could not open the database folder.")
