"""On-demand biological exploration of one frozen source result."""
from __future__ import annotations

import shutil
from pathlib import Path

from PySide6.QtCore import Qt, QUrl
from PySide6.QtGui import QDesktopServices, QPixmap, QPainter, QPdfWriter
from PySide6.QtWidgets import (QCheckBox, QComboBox, QFileDialog, QHBoxLayout, QLabel, QMessageBox,
    QPushButton, QScrollArea, QTableWidget, QTableWidgetItem, QTabWidget,
    QVBoxLayout, QWidget)

from ..core.biological_context import (BiologicalContextError, differential_choices, frozen_direction,
    find_contexts, list_saved_contexts, load_saved_context, other_differential_runs, resolve_entity,
    save_context_figure)
from ..core.biological_context_render import render_context
from ..core import differential_analysis


class BiologicalContextPage(QWidget):
    def __init__(self):
        super().__init__()
        self.project = None
        self.context = None
        self.index = None
        self.items = []
        self.current_item = None
        self.current_figure = None
        self.current_root = None
        self._unverified_differential = set()
        layout = QVBoxLayout(self)
        self.header = QLabel("Select one persisted result to explore its biological context.")
        self.header.setWordWrap(True)
        layout.addWidget(self.header)
        direction = QHBoxLayout()
        direction.addWidget(QLabel("Differential context:"))
        self.differential = QComboBox()
        self.differential.addItem("None", "")
        self.differential.currentIndexChanged.connect(self._direction_changed)
        direction.addWidget(self.differential)
        self.direction_note = QLabel("No compatible Differential Analysis run is selected. Directional status is not available.")
        self.direction_note.setWordWrap(True)
        direction.addWidget(self.direction_note, 1)
        layout.addLayout(direction)
        self.cross_lineage = QCheckBox("Include biological identifier links from other or unverified input lineages")
        self.cross_lineage.setToolTip("These links do not equate experimental features across input lineages.")
        self.cross_lineage.toggled.connect(self._refresh_contexts)
        layout.addWidget(self.cross_lineage)
        self.tabs = QTabWidget()
        self.tables = {}
        for module, title in (("kegg", "KEGG"), ("reactome", "Reactome"), ("string", "STRING")):
            page = QWidget()
            box = QVBoxLayout(page)
            table = QTableWidget(0, 7)
            table.setHorizontalHeaderLabels(["ID", "Title / type", "Run", "Detected members",
                "Members with direction", "Snapshot", "Official image"])
            table.itemSelectionChanged.connect(lambda m=module: self._selected(m))
            box.addWidget(table)
            self.tables[module] = table
            self.tabs.addTab(page, title)
        self.overview = QLabel("No context is loaded.")
        self.overview.setWordWrap(True)
        overview_page = QWidget()
        overview_layout = QVBoxLayout(overview_page)
        overview_layout.addWidget(self.overview)
        overview_layout.addStretch()
        self.tabs.insertTab(0, overview_page, "Overview")
        viewer = QWidget()
        viewer_layout = QVBoxLayout(viewer)
        actions = QHBoxLayout()
        self.official_button = QPushButton("Open official local image / diagram")
        self.official_button.clicked.connect(self.open_official)
        self.build_button = QPushButton("Build PichAnalysis view")
        self.build_button.clicked.connect(self.build_view)
        actions.addWidget(self.official_button)
        actions.addWidget(self.build_button)
        viewer_layout.addLayout(actions)
        self.notice = QLabel("Select a pathway or network first. Official images are never modified.")
        self.notice.setWordWrap(True)
        viewer_layout.addWidget(self.notice)
        self.figure = QLabel()
        self.figure.setAlignment(Qt.AlignCenter)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(self.figure)
        self.detail_tabs = QTabWidget()
        self.detail_tabs.addTab(scroll, "Figure")
        self.member_table = QTableWidget()
        self.edge_table = QTableWidget()
        self.detail_tabs.addTab(self.member_table, "Detected members")
        self.detail_tabs.addTab(self.edge_table, "Persisted network edges")
        viewer_layout.addWidget(self.detail_tabs, 1)
        self.tabs.addTab(viewer, "Visualization")
        exports = QWidget()
        export_layout = QVBoxLayout(exports)
        self.export_table_button = QPushButton("Export selected table as CSV...")
        self.export_table_button.clicked.connect(self.export_table)
        self.export_workbook_button = QPushButton("Export selected table as XLSX...")
        self.export_workbook_button.clicked.connect(self.export_workbook)
        self.export_png_button = QPushButton("Export generated PNG...")
        self.export_png_button.clicked.connect(self.export_png)
        self.export_pdf_button = QPushButton("Export generated PDF...")
        self.export_pdf_button.clicked.connect(self.export_pdf)
        self.export_official_button = QPushButton("Copy official image...")
        self.export_official_button.clicked.connect(self.export_official)
        self.history = QComboBox()
        self.history.addItem("Saved visualizations", None)
        self.history.currentIndexChanged.connect(self._load_history)
        for widget in (self.export_table_button, self.export_workbook_button, self.export_png_button,
                       self.export_pdf_button, self.export_official_button,
                       QLabel("Historical generated figures"), self.history):
            export_layout.addWidget(widget)
        export_layout.addStretch()
        self.tabs.addTab(exports, "Exports / History")
        layout.addWidget(self.tabs, 1)
        self._update_actions()

    def set_project(self, project):
        self.project = project
        self._clear()
        self.history.blockSignals(True)
        self.history.clear()
        self.history.addItem("Saved visualizations", None)
        if project:
            for root in list_saved_contexts(project):
                self.history.addItem(root.name, str(root))
        self.history.blockSignals(False)

    def _clear(self):
        self.context = self.index = self.current_item = self.current_figure = self.current_root = None
        self.items = []
        self._unverified_differential.clear()
        for table in self.tables.values():
            table.setRowCount(0)
        self.figure.clear()
        self._fill_detail_table(self.member_table, None)
        self._fill_detail_table(self.edge_table, None)
        self.overview.setText("No context is loaded.")
        self.notice.setText("Select a pathway or network first.")
        self.cross_lineage.blockSignals(True)
        self.cross_lineage.setChecked(False)
        self.cross_lineage.blockSignals(False)
        self.differential.blockSignals(True)
        self.differential.clear()
        self.differential.addItem("None", "")
        self.differential.blockSignals(False)
        self._update_actions()

    def open_entity(self, module, run_id, identifier, source_row=None):
        self._clear()
        if not self.project:
            raise BiologicalContextError("Open a project first.")
        try:
            self.index, self.context = resolve_entity(self.project, module, run_id, identifier, source_row)
            self.items = find_contexts(self.project, self.index, self.context)
            choices = differential_choices(self.index, self.context)
            self.differential.blockSignals(True)
            for choice in choices:
                self.differential.addItem(choice, choice)
            self._unverified_differential = set(other_differential_runs(self.index, choices))
            for choice in sorted(self._unverified_differential):
                self.differential.addItem(f"{choice} — different/unverified input lineage", choice)
            if module == "differential" and run_id in choices:
                self.differential.setCurrentIndex(self.differential.findData(run_id))
            self.differential.blockSignals(False)
            self._direction_changed()
            self._fill_tables()
            context = self.context
            self.overview.setText("\n".join((
                f"Original identifier: {context.original_identifier or 'Not recorded'}",
                f"Feature ID: {context.feature_id or 'Not recorded'}",
                f"Source row: {context.source_row if context.source_row is not None else 'Not recorded'}",
                f"UniProt: {', '.join(context.uniprot_accessions) or 'Not recorded'}",
                f"Gene Symbol: {', '.join(context.gene_symbols) or 'Not recorded'}",
                f"NCBI Gene: {', '.join(context.ncbi_gene_ids) or 'Not recorded'}",
                f"Source: {context.source_module} / {context.source_run_id}",
                f"Input lineage: {context.input_lineage.frozen_input_hash if context.input_lineage else 'Not recorded'}",
                f"Source snapshot: {context.database_snapshot_provenance or 'Not recorded'}",
                "Biological identifier links do not equate experimental features across input lineages.",
                "Up/down is not pathway activity.",
            )))
            self.tabs.setCurrentIndex(0)
        except Exception:
            self._clear()
            raise

    def _refresh_contexts(self):
        if not self.project or not self.context or not self.index:
            return
        self.current_item = self.current_figure = self.current_root = None
        self.figure.clear()
        try:
            self.items = find_contexts(self.project, self.index, self.context,
                include_other_lineages=self.cross_lineage.isChecked())
            self._fill_tables()
            self.notice.setText("Cross-lineage links are biological identifier matches, not the same experimental feature."
                if self.cross_lineage.isChecked() else "Select a pathway or network first.")
        except Exception as error:
            self.items = []
            self._fill_tables()
            self.notice.setText(str(error))
        self._update_actions()
    def _direction_changed(self):
        run = self.differential.currentData() or ""
        if self.context:
            self.current_item = self.current_figure = self.current_root = None
            self.figure.clear()
            self._fill_tables()
            self._update_actions()
        self.direction_note.setText(
            (f"Direction from frozen Differential run {run}: Condition A - Condition B. "
             "Different/unverified lineage: biological identifier matching only; not the same experimental feature."
             if run in self._unverified_differential else
             f"Direction from frozen Differential run {run}: Condition A - Condition B.")
            if run else "No compatible Differential Analysis run is selected. Directional status is not available.")

    def _fill_tables(self):
        for module, table in self.tables.items():
            items = [item for item in self.items if item.module == module]
            table.setRowCount(len(items))
            for row, item in enumerate(items):
                values = (item.item_id, item.title if module != "string" else item.network_type,
                    item.run_id, str(len(item.members)), "On selection" if self.differential.currentData() else "0",
                    item.snapshot or "Not recorded", "Available" if item.official_image else "Unavailable")
                for column, value in enumerate(values):
                    cell = QTableWidgetItem(value)
                    cell.setData(Qt.UserRole, next(index for index, candidate in enumerate(self.items) if candidate is item))
                    table.setItem(row, column, cell)

    @staticmethod
    def _fill_detail_table(table, frame):
        table.setRowCount(0)
        if frame is None:
            table.setColumnCount(0)
            return
        table.setColumnCount(len(frame.columns))
        table.setHorizontalHeaderLabels([str(name) for name in frame.columns])
        table.setRowCount(min(len(frame), 500))
        for row, values in enumerate(frame.head(500).itertuples(index=False, name=None)):
            for column, value in enumerate(values):
                table.setItem(row, column, QTableWidgetItem(str(value)))
    def _selected(self, module):
        table = self.tables[module]
        row = table.currentRow()
        if row < 0 or not table.item(row, 0):
            return
        self.current_item = self.items[table.item(row, 0).data(Qt.UserRole)]
        self._fill_detail_table(self.member_table, self.current_item.members)
        self._fill_detail_table(self.edge_table, self.current_item.edges)
        self.current_figure = self.current_root = None
        self.figure.clear()
        self.notice.setText("Official pathway image is not available in this local snapshot."
            if module != "string" and not self.current_item.official_image
            else "Choose an action to view the official resource or build a separate PichAnalysis view.")
        run = self.differential.currentData() or None
        if run and self.context:
            try:
                results = differential_analysis.load_run(self.project, run)["tables"]["all_results"]
                count = sum(frozen_direction(self.project, run, self.context, member, results) != "none"
                    for member in self.current_item.members.fillna("").to_dict("records"))
                table.item(row, 4).setText(str(count))
            except Exception as error:
                table.item(row, 4).setText("Unavailable")
                self.notice.setText(f"Differential context unavailable: {error}")
        self._update_actions()

    def _update_actions(self):
        item = self.current_item
        self.official_button.setEnabled(bool(item and item.official_image))
        self.build_button.setEnabled(bool(item and self.context))
        self.export_table_button.setEnabled(bool(item))
        self.export_workbook_button.setEnabled(bool(item))
        self.export_official_button.setEnabled(bool(item and item.official_image))
        self.export_png_button.setEnabled(self.current_figure is not None)
        self.export_pdf_button.setEnabled(self.current_figure is not None)

    def open_official(self):
        if self.current_item and self.current_item.official_image:
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.current_item.official_image)))

    def build_view(self):
        if not self.current_item or not self.context:
            return
        try:
            run = self.differential.currentData() or None
            self.current_figure = render_context(self.project, self.context, self.current_item, run)
            self.current_root = save_context_figure(self.project, self.context, self.current_item,
                run, self.current_figure)
            self.figure.setPixmap(QPixmap.fromImage(self.current_figure))
            self.notice.setText(f"Separate PichAnalysis view saved in {self.current_root}. "
                "Colors show detected member status, not pathway activity.")
            self.history.addItem(self.current_root.name, str(self.current_root))
            self.tabs.setCurrentIndex(4)
        except Exception as error:
            self.current_figure = self.current_root = None
            self.figure.clear()
            self.notice.setText(str(error))
        self._update_actions()

    def _save_path(self, title, default, pattern):
        value, _ = QFileDialog.getSaveFileName(self, title, default, pattern)
        return Path(value) if value else None

    def export_table(self):
        if self.current_item:
            path = self._save_path("Export context table", "biological_context.csv", "CSV (*.csv)")
            if path:
                self.current_item.members.to_csv(path, index=False)

    def export_workbook(self):
        if self.current_item:
            path = self._save_path("Export context table", "biological_context.xlsx", "Excel (*.xlsx)")
            if path:
                with __import__("pandas").ExcelWriter(path) as writer:
                    self.current_item.members.to_excel(writer, sheet_name="Members", index=False)
                    if not self.current_item.edges.empty:
                        self.current_item.edges.to_excel(writer, sheet_name="Edges", index=False)
    def export_png(self):
        if self.current_figure is not None:
            path = self._save_path("Export PichAnalysis figure", "biological_context.png", "PNG (*.png)")
            if path and not self.current_figure.save(str(path), "PNG"):
                QMessageBox.warning(self, "Export failed", "Could not save the PNG figure.")

    def export_pdf(self):
        if self.current_figure is not None:
            path = self._save_path("Export PichAnalysis figure", "biological_context.pdf", "PDF (*.pdf)")
            if path:
                writer = QPdfWriter(str(path))
                painter = QPainter(writer)
                painter.drawImage(writer.pageLayout().paintRectPixels(writer.resolution()), self.current_figure)
                painter.end()

    def export_official(self):
        if self.current_item and self.current_item.official_image:
            path = self._save_path("Copy official image", self.current_item.official_image.name, "PNG (*.png)")
            if path:
                shutil.copy2(self.current_item.official_image, path)

    def _load_history(self):
        root = self.history.currentData()
        if not root:
            return
        try:
            root = Path(root)
            metadata, figure = load_saved_context(root)
            self.current_item = None
            self.context = None
            self.current_figure = QPixmap(str(figure)).toImage()
            import pandas as pd
            self._fill_detail_table(self.member_table, pd.read_csv(root / "tables/members.csv", dtype=str))
            edges = root / "tables/edges.csv"
            self._fill_detail_table(self.edge_table, pd.read_csv(edges, dtype=str) if edges.is_file() and edges.stat().st_size else None)
            self.current_root = root
            self.figure.setPixmap(QPixmap.fromImage(self.current_figure))
            self.notice.setText(f"Historical figure: {metadata['source_module']} / {metadata['source_run']} "
                f"→ {metadata['item_module']} / {metadata['item_run']} / {metadata['item_id']} "
                f"(snapshot {metadata.get('snapshot') or 'not recorded'}).")
            self.tabs.setCurrentIndex(4)
            self._update_actions()
        except (OSError, ValueError, KeyError, BiologicalContextError) as error:
            self.figure.clear()
            self.current_figure = None
            self.notice.setText(str(error))
