"""Explicit, offline builder and history viewer for Consolidated Reports."""
from __future__ import annotations

from PySide6.QtCore import QThread, Qt, QUrl, Signal
from PySide6.QtGui import QCloseEvent, QDesktopServices
from PySide6.QtWidgets import (QCheckBox, QComboBox, QFileDialog, QFormLayout,
    QHBoxLayout, QLabel, QLineEdit, QListWidget, QMessageBox, QPlainTextEdit,
    QPushButton, QTableWidget, QTableWidgetItem, QTreeWidget, QTreeWidgetItem,
    QVBoxLayout, QWidget)

from ..core.consolidated_report import (ReportConfig, ReportError, ReportSection,
    discover_runs, export_report, generate_report, list_reports, load_report)


class ReportWorker(QThread):
    succeeded = Signal(str)
    failed = Signal(str)

    def __init__(self, project, config, sources=None):
        super().__init__()
        self.project, self.config, self.sources = project, config, sources

    def run(self):
        try:
            path = generate_report(self.project, self.config, self.sources)
            self.succeeded.emit(path.name)
        except Exception as error:
            self.failed.emit(str(error))


class ConsolidatedReportPage(QWidget):
    def __init__(self):
        super().__init__()
        self.project = None
        self.sources = None
        self.runs = {}
        self.worker = None
        self._notes = {}
        self._current_note_key = None
        layout = QVBoxLayout(self)
        title = QLabel("Consolidated Report")
        title.setStyleSheet("font-size: 20px; font-weight: 600;")
        layout.addWidget(title)
        help_text = QLabel("Consolidated Reports organize existing PichAnalysis results. Generating a report does not rerun or reinterpret the selected analyses.")
        help_text.setWordWrap(True)
        layout.addWidget(help_text)
        self.title = QLineEdit("PichAnalysis Consolidated Report")
        self.subtitle = QLineEdit()
        self.author = QLineEdit()
        self.description = QPlainTextEdit()
        self.description.setMaximumHeight(60)
        form = QFormLayout()
        form.addRow("Report title", self.title)
        form.addRow("Subtitle", self.subtitle)
        form.addRow("Author", self.author)
        form.addRow("Description", self.description)
        layout.addLayout(form)
        source_layout = QHBoxLayout()
        self.sources_tree = QTreeWidget()
        self.sources_tree.setHeaderLabels(["Module / run / artifact", "Date or type", "Snapshot", "Status"])
        source_layout.addWidget(self.sources_tree, 3)
        side = QVBoxLayout()
        side.addWidget(QLabel("Section order - select runs explicitly"))
        self.order = QListWidget()
        side.addWidget(self.order, 1)
        row = QHBoxLayout()
        self.move_up = QPushButton("Move up")
        self.move_down = QPushButton("Move down")
        row.addWidget(self.move_up); row.addWidget(self.move_down)
        side.addLayout(row)
        side.addWidget(QLabel("User notes for selected run"))
        self.notes = QPlainTextEdit()
        self.notes.setPlaceholderText("Optional user-provided notes; not generated interpretation.")
        side.addWidget(self.notes, 1)
        self.run_metadata = QPlainTextEdit()
        self.run_metadata.setReadOnly(True)
        side.addWidget(self.run_metadata, 1)
        source_layout.addLayout(side, 2)
        layout.addLayout(source_layout, 1)
        preview_row = QHBoxLayout()
        preview_row.addWidget(QLabel("Table preview rows"))
        self.preview_rows = QComboBox()
        for number in (10, 20, 50, 100):
            self.preview_rows.addItem(str(number), number)
        self.preview_rows.setCurrentIndex(1)
        preview_row.addWidget(self.preview_rows)
        self.provenance = QCheckBox("Include provenance appendix")
        self.provenance.setChecked(True)
        self.attachments = QCheckBox("Include complete selected source files as attachments")
        self.attachments.setChecked(True)
        preview_row.addWidget(self.provenance); preview_row.addWidget(self.attachments)
        layout.addLayout(preview_row)
        caveat = QLabel("Large tables are shown as previews. The complete selected source table can be included as an attachment.")
        caveat.setWordWrap(True)
        layout.addWidget(caveat)
        layout.addWidget(QLabel("User discussion (optional, clearly identified as user-provided text)"))
        self.discussion = QPlainTextEdit()
        self.discussion.setMaximumHeight(70)
        layout.addWidget(self.discussion)
        self.preview = QLabel("Select runs and artifacts to preview report contents.")
        self.preview.setWordWrap(True)
        layout.addWidget(self.preview)
        buttons = QHBoxLayout()
        self.refresh = QPushButton("Refresh source runs")
        self.generate = QPushButton("Generate report")
        self.generate.setEnabled(False)
        buttons.addWidget(self.refresh); buttons.addWidget(self.generate)
        layout.addLayout(buttons)
        layout.addWidget(QLabel("Historical reports"))
        self.history = QTableWidget(0, 6)
        self.history.setHorizontalHeaderLabels(["Date/time", "Title", "Modules", "Selected runs", "Report ID", "Status"])
        layout.addWidget(self.history, 1)
        history_buttons = QHBoxLayout()
        self.open_report = QPushButton("Open report")
        self.open_folder = QPushButton("Open report folder")
        self.export_html = QPushButton("Export HTML...")
        self.export_pdf = QPushButton("Export PDF...")
        self.export_bundle = QPushButton("Export bundle...")
        for button in (self.open_report, self.open_folder, self.export_html,
                       self.export_pdf, self.export_bundle):
            history_buttons.addWidget(button)
        layout.addLayout(history_buttons)
        self.status = QLabel("Open a project to create a report.")
        self.status.setWordWrap(True)
        layout.addWidget(self.status)
        self.sources_tree.itemChanged.connect(self._selection_changed)
        self.order.currentRowChanged.connect(self._order_selected)
        self.notes.textChanged.connect(self._save_note)
        self.move_up.clicked.connect(lambda: self._move(-1))
        self.move_down.clicked.connect(lambda: self._move(1))
        self.refresh.clicked.connect(self.refresh_sources)
        self.generate.clicked.connect(self._generate)
        self.open_report.clicked.connect(self._open_report)
        self.open_folder.clicked.connect(self._open_folder)
        self.export_html.clicked.connect(lambda: self._export("html"))
        self.export_pdf.clicked.connect(lambda: self._export("pdf"))
        self.export_bundle.clicked.connect(lambda: self._export("bundle"))

    def set_project(self, project):
        self.project = project
        self._notes.clear()
        self._current_note_key = None
        self.order.clear()
        self.sources_tree.clear()
        self.refresh_sources()
        self.refresh_history()

    def is_running(self):
        return bool(self.worker and self.worker.isRunning())

    def refresh_sources(self):
        self.sources_tree.blockSignals(True)
        self.sources_tree.clear()
        self.order.clear()
        self.runs.clear()
        if not self.project:
            self.sources_tree.blockSignals(False)
            self.status.setText("Open a project to create a report.")
            self._update_preview()
            return
        try:
            found = discover_runs(self.project, self.sources)
        except Exception as error:
            self.sources_tree.blockSignals(False)
            self.status.setText(f"Could not discover report sources: {error}")
            return
        groups = {}
        for run in found:
            key = (run.module, run.run_id)
            self.runs[key] = run
            group = groups.get(run.module)
            if group is None:
                group = QTreeWidgetItem([run.label, "", "", ""])
                groups[run.module] = group
                self.sources_tree.addTopLevelItem(group)
            item = QTreeWidgetItem([run.run_id, run.run_date, run.snapshot, run.status])
            item.setData(0, Qt.ItemDataRole.UserRole, key)
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(0, Qt.CheckState.Unchecked)
            if run.status != "Ready":
                item.setDisabled(True)
                item.setToolTip(0, run.reason)
            group.addChild(item)
            for artifact in run.artifacts:
                child = QTreeWidgetItem([artifact.label, artifact.kind, artifact.artifact_id, "Available"])
                child.setData(0, Qt.ItemDataRole.UserRole, artifact.artifact_id)
                child.setFlags(child.flags() | Qt.ItemFlag.ItemIsUserCheckable)
                child.setCheckState(0, Qt.CheckState.Unchecked)
                item.addChild(child)
        self.sources_tree.blockSignals(False)
        self.status.setText(f"{len(found)} persisted runs discovered. Select runs and artifacts explicitly.")
        self._update_preview()

    def _run_item(self, key):
        for i in range(self.sources_tree.topLevelItemCount()):
            group = self.sources_tree.topLevelItem(i)
            for j in range(group.childCount()):
                child = group.child(j)
                if child.data(0, Qt.ItemDataRole.UserRole) == key:
                    return child
        return None

    def _selection_changed(self, item, _column):
        key = item.data(0, Qt.ItemDataRole.UserRole)
        if not isinstance(key, tuple):
            self._update_preview()
            return
        if item.checkState(0) == Qt.CheckState.Checked:
            if not any(self.order.item(i).data(Qt.ItemDataRole.UserRole) == key
                       for i in range(self.order.count())):
                from PySide6.QtWidgets import QListWidgetItem
                row = QListWidgetItem(f"{self.runs[key].label} - {key[1]}")
                row.setData(Qt.ItemDataRole.UserRole, key)
                self.order.addItem(row)
        else:
            for i in range(self.order.count()):
                if self.order.item(i).data(Qt.ItemDataRole.UserRole) == key:
                    self.order.takeItem(i)
                    break
        self._update_preview()

    def _order_selected(self, row):
        self._save_note()
        self._current_note_key = self.order.item(row).data(Qt.ItemDataRole.UserRole) if row >= 0 else None
        self.notes.blockSignals(True)
        self.notes.setPlainText(self._notes.get(self._current_note_key, ""))
        self.notes.blockSignals(False)
        run = self.runs.get(self._current_note_key)
        self.run_metadata.setPlainText("" if run is None else
            f"Module: {run.label}\nRun ID: {run.run_id}\nDate: {run.run_date}\n"
            f"Input lineage: {run.input_lineage}\nSnapshot/version: {run.snapshot}\n"
            f"Status: {run.status}\nSource: {run.source_path}")

    def _save_note(self):
        if self._current_note_key is not None:
            self._notes[self._current_note_key] = self.notes.toPlainText()

    def _move(self, offset):
        index = self.order.currentRow()
        target = index + offset
        if index < 0 or target < 0 or target >= self.order.count():
            return
        item = self.order.takeItem(index)
        self.order.insertItem(target, item)
        self.order.setCurrentRow(target)
        self._update_preview()

    def build_config(self):
        self._save_note()
        sections = []
        for i in range(self.order.count()):
            key = self.order.item(i).data(Qt.ItemDataRole.UserRole)
            run_item = self._run_item(key)
            artifacts = tuple(run_item.child(j).data(0, Qt.ItemDataRole.UserRole)
                for j in range(run_item.childCount())
                if run_item.child(j).checkState(0) == Qt.CheckState.Checked)
            sections.append(ReportSection(key[0], key[1], artifacts,
                int(self.preview_rows.currentData()), self._notes.get(key, "")))
        return ReportConfig(self.title.text(), self.subtitle.text(), self.author.text(),
            self.description.toPlainText(), tuple(sections), self.discussion.toPlainText(),
            self.provenance.isChecked(), self.attachments.isChecked())

    def _update_preview(self):
        config = self.build_config()
        selected = [self.runs[(section.module, section.run_id)] for section in config.sections]
        figures = sum(artifact.kind == "figure" for section, run in zip(config.sections, selected)
            for artifact in run.artifacts if artifact.artifact_id in section.artifact_ids)
        tables = sum(artifact.kind == "table" for section, run in zip(config.sections, selected)
            for artifact in run.artifacts if artifact.artifact_id in section.artifact_ids)
        lineages = len({run.input_lineage for run in selected})
        snapshots = len({run.snapshot for run in selected})
        self.preview.setText(f"Preview: {len(selected)} sections / runs; {figures} figures; {tables} tables; "
            f"{figures + tables if config.include_attachments else 0} selected-file attachments; "
            f"{lineages} input lineages; {snapshots} database snapshots.")
        self.generate.setEnabled(bool(self.project and selected and all(section.artifact_ids
            for section in config.sections) and not (self.worker and self.worker.isRunning())))

    def _generate(self):
        if not self.project or (self.worker and self.worker.isRunning()):
            return
        config = self.build_config()
        if not config.sections or any(not section.artifact_ids for section in config.sections):
            self.status.setText("Select at least one artifact for each selected run.")
            return
        self.worker = ReportWorker(self.project, config, self.sources)
        self.worker.succeeded.connect(self._generation_done)
        self.worker.failed.connect(self._generation_failed)
        self.worker.finished.connect(self._update_preview)
        self.status.setText("Generating consolidated report...")
        self.generate.setEnabled(False)
        self.worker.start()

    def _generation_done(self, report_id):
        self.status.setText(f"Report Ready: {report_id}")
        self.refresh_history()

    def _generation_failed(self, message):
        self.status.setText(f"Report Incomplete: {message}")
        self.refresh_history()

    def refresh_history(self):
        rows = list_reports(self.project) if self.project else []
        self.history.setRowCount(len(rows))
        for row, report in enumerate(rows):
            values = (report["date"], report["title"], report["modules"],
                report["runs"], report["report_id"], report["status"])
            for column, value in enumerate(values):
                self.history.setItem(row, column, QTableWidgetItem(str(value)))

    def _selected_report_id(self):
        row = self.history.currentRow()
        return self.history.item(row, 4).text() if row >= 0 and self.history.item(row, 4) else None

    def _open_report(self):
        report_id = self._selected_report_id()
        if not report_id or not self.project:
            return
        try:
            report = load_report(self.project, report_id)
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(report["html"])))
        except (OSError, ValueError, ReportError) as error:
            self.status.setText(f"Could not open historical report: {error}")

    def _open_folder(self):
        report_id = self._selected_report_id()
        if not report_id or not self.project:
            return
        try:
            report = load_report(self.project, report_id)
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(report["path"])))
        except (OSError, ValueError, ReportError) as error:
            self.status.setText(f"Could not open report folder: {error}")

    def _export(self, kind):
        report_id = self._selected_report_id()
        if not report_id or not self.project:
            return
        if kind == "bundle":
            parent = QFileDialog.getExistingDirectory(self, "Choose export parent folder")
            if not parent:
                return
            from pathlib import Path
            destination = Path(parent) / report_id
        else:
            name, _ = QFileDialog.getSaveFileName(self, f"Export report {kind.upper()}",
                f"{report_id}.{kind}", f"{kind.upper()} files (*.{kind})")
            if not name:
                return
            from pathlib import Path
            destination = Path(name)
        try:
            export_report(self.project, report_id, destination, kind)
            self.status.setText(f"Exported {kind.upper()}: {destination}")
        except (OSError, ValueError, ReportError) as error:
            self.status.setText(f"Export failed: {error}")

    def closeEvent(self, event: QCloseEvent):
        if self.worker and self.worker.isRunning():
            self.status.setText("Report generation is still running. Close after it completes.")
            event.ignore()
            return
        super().closeEvent(event)
