from __future__ import annotations

from typing import Any

from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QFormLayout, QHBoxLayout, QLabel, QLineEdit,
    QPlainTextEdit, QProgressBar, QPushButton, QTableWidget, QTableWidgetItem,
    QTabWidget, QVBoxLayout, QWidget,
)

from ..core.column_mapping import experimental_design
from ..core.mapping_analysis import MappingOutputs, mapping_readiness
from ..core.organism import COMMON_ORGANISMS, get_organism
from ..core.project import Project
from .presence_page import PresencePage
from .go_page import GOPage


class AnalysesPage(QWidget):
    organism_requested = Signal(str, str)
    run_requested = Signal(bool)
    export_table_requested = Signal()
    export_workbook_requested = Signal()
    open_results_requested = Signal()

    def __init__(self) -> None:
        super().__init__()
        self.project: Project | None = None
        self.organism = QComboBox()
        for item in COMMON_ORGANISMS:
            self.organism.addItem(f"{item.name} — {item.tax_id}", (item.name, item.tax_id))
        self.organism.addItem("Outro...", None)
        self.custom_name = QLineEdit()
        self.custom_name.setPlaceholderText("Nome científico")
        self.custom_tax_id = QLineEdit()
        self.custom_tax_id.setPlaceholderText("NCBI Taxonomy ID")
        self.save_organism = QPushButton("Save organism")
        self.identifier = QLabel("—")
        self.identifier_type = QLabel("—")
        self.record_count = QLabel("0")
        form = QFormLayout()
        form.addRow("Organism", self.organism)
        form.addRow("Nome personalizado", self.custom_name)
        form.addRow("Taxonomy ID", self.custom_tax_id)
        form.addRow("Identificador principal", self.identifier)
        form.addRow("Identifier type", self.identifier_type)
        form.addRow("Número de registros", self.record_count)
        self.readiness = QLabel("Abra um projeto.")
        self.readiness.setWordWrap(True)
        self.refresh = QCheckBox("Atualizar anotações online")
        self.refresh.setToolTip("Desmarcado: usar dados já disponíveis no cache quando possível.")
        self.run_button = QPushButton("Mapear e anotar proteínas")
        self.run_button.setEnabled(False)
        self.progress = QProgressBar()
        self.progress.setRange(0, 0)
        self.progress.hide()
        self.summary = QLabel("No results available.")
        self.summary.setWordWrap(True)
        self.preview = QTableWidget()
        self.detail = QPlainTextEdit()
        self.detail.setReadOnly(True)
        result_tabs = QTabWidget()
        result_tabs.addTab(self.preview, "Catálogo")
        result_tabs.addTab(self.detail, "Detalhes da linha")
        self.export_table = QPushButton("Export table...")
        self.export_workbook = QPushButton("Export workbook...")
        self.open_results = QPushButton("Open results folder")
        actions = QHBoxLayout()
        for button in (self.export_table, self.export_workbook, self.open_results):
            button.setEnabled(False)
            actions.addWidget(button)
        mapping_page = QWidget()
        layout = QVBoxLayout(mapping_page)
        layout.addWidget(QLabel("Identificação e anotação"))
        layout.addLayout(form)
        layout.addWidget(self.save_organism)
        layout.addWidget(self.readiness)
        layout.addWidget(self.refresh)
        layout.addWidget(self.run_button)
        layout.addWidget(self.progress)
        layout.addWidget(self.summary)
        layout.addWidget(result_tabs, 1)
        layout.addLayout(actions)
        self.organism.currentIndexChanged.connect(self._organism_mode)
        self.save_organism.clicked.connect(self._save_organism)
        self.run_button.clicked.connect(lambda: self.run_requested.emit(self.refresh.isChecked()))
        self.export_table.clicked.connect(self.export_table_requested)
        self.export_workbook.clicked.connect(self.export_workbook_requested)
        self.open_results.clicked.connect(self.open_results_requested)
        self.preview.currentCellChanged.connect(self._show_detail)
        self._organism_mode()
        self.presence_page = PresencePage()
        self.go_page = GOPage()
        module_tabs = QTabWidget()
        module_tabs.addTab(mapping_page, "Identificação e anotação")
        module_tabs.addTab(self.presence_page, "Presence / absence")
        module_tabs.addTab(self.go_page, "Gene Ontology")
        outer_layout = QVBoxLayout(self)
        outer_layout.addWidget(module_tabs)

    def _organism_mode(self) -> None:
        custom = self.organism.currentData() is None
        self.custom_name.setEnabled(custom)
        self.custom_tax_id.setEnabled(custom)

    def _save_organism(self) -> None:
        selected = self.organism.currentData()
        if selected is None:
            self.organism_requested.emit(self.custom_name.text(), self.custom_tax_id.text())
        else:
            self.organism_requested.emit(*selected)

    def set_project(self, project: Project | None) -> None:
        self.project = project
        self.presence_page.set_project(project)
        self.go_page.set_project(project)
        if project is None:
            self.run_button.setEnabled(False)
            return
        current = get_organism(project)
        if current:
            found = self.organism.findData((current.name, current.tax_id))
            if found >= 0:
                self.organism.setCurrentIndex(found)
            else:
                self.organism.setCurrentIndex(self.organism.count()-1)
                self.custom_name.setText(current.name)
                self.custom_tax_id.setText(current.tax_id)
        design = experimental_design(project)
        self.identifier.setText(str(design["primary_identifier_column"] or "Not configured"))
        self.identifier_type.setText(str(design["primary_identifier_type"] or "Not configured"))
        self.record_count.setText(str(project.config.get("input", {}).get("rows") or 0))
        state = mapping_readiness(project)
        self.readiness.setText(state.reason)
        self.run_button.setEnabled(state.ready)

    def set_running(self, running: bool) -> None:
        self.progress.setVisible(running)
        self.run_button.setEnabled(not running and mapping_readiness(self.project).ready)
        self.readiness.setText("Em execução..." if running else mapping_readiness(self.project).reason)

    def show_outputs(self, outputs: MappingOutputs) -> None:
        meta = outputs.metadata
        self.summary.setText(
            f"IDs de entrada: {meta.get('input_count', 0)} | IDs únicos: {meta.get('unique_id_count', 0)} | "
            f"Mapeados unicamente: {meta.get('mapped_unique_count', 0)} | Ambíguos: {meta.get('ambiguous_count', 0)} | "
            f"Unmapped: {meta.get('unmapped_count', 0)} | Organism mismatches: {meta.get('organism_mismatch_count', 0)}"
        )
        frame = outputs.catalog.head(200)
        self.preview.setRowCount(len(frame)); self.preview.setColumnCount(len(frame.columns))
        self.preview.setHorizontalHeaderLabels([str(value) for value in frame.columns])
        for row, values in enumerate(frame.itertuples(index=False, name=None)):
            for column, value in enumerate(values):
                self.preview.setItem(row, column, QTableWidgetItem("" if value is None else str(value)))
        for button in (self.export_table, self.export_workbook, self.open_results):
            button.setEnabled(True)

    def _show_detail(self, row: int, _column: int, _old_row: int, _old_column: int) -> None:
        if row < 0:
            return
        wanted = ("uniprot_accession", "gene_symbol", "protein_name", "uniprot_function",
            "uniprot_subcellular_location", "ncbi_summary", "ncbi_summary_source",
            "ncbi_summary_date", "mapping_status")
        headers = {self.preview.horizontalHeaderItem(i).text(): i for i in range(self.preview.columnCount())}
        lines = []
        for name in wanted:
            if name in headers:
                item = self.preview.item(row, headers[name])
                lines.append(f"{name}: {item.text() if item else ''}")
        self.detail.setPlainText("\n\n".join(lines))
