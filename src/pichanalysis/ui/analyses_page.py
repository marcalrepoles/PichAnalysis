from __future__ import annotations

from typing import Any

from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QFormLayout, QHBoxLayout, QLabel, QLineEdit,
    QPlainTextEdit, QProgressBar, QPushButton, QTableWidget, QTableWidgetItem, QMessageBox,
    QTabWidget, QVBoxLayout, QWidget,
)

from ..core.column_mapping import experimental_design
from ..core.mapping_analysis import MappingOutputs, mapping_readiness
from ..core.organism import COMMON_ORGANISMS, get_organism
from ..core.project import Project
from .presence_page import PresencePage
from .go_page import GOPage
from .kegg_page import KEGGPage
from .reactome_page import ReactomePage
from .mitocarta_page import MitoCartaPage
from .interpro_pfam_page import InterProPfamPage
from .string_page import StringPage
from .complexes_page import ComplexPage
from .mtdna_page import MtdnaPage
from .proteomics_qc_page import ProteomicsQCPage
from .differential_analysis_page import DifferentialAnalysisPage
from .cross_module_explorer import CrossModuleExplorer
from .biological_context_page import BiologicalContextPage
from .experiment_comparison_page import ExperimentComparisonPage
from .cross_module_actions import attach_explore_action, attach_biological_context_action
from .cross_module_navigation import open_persisted_target
from ..core.cross_module_integration import default_registry
from ..core.database_manager import DatabaseManager


class AnalysesPage(QWidget):
    organism_requested = Signal(str, str)
    run_requested = Signal(bool)
    export_table_requested = Signal()
    export_workbook_requested = Signal()
    open_results_requested = Signal()

    def __init__(self, database_manager: DatabaseManager | None = None) -> None:
        super().__init__()
        self.project: Project | None = None
        self.organism = QComboBox()
        for item in COMMON_ORGANISMS:
            self.organism.addItem(f"{item.name} — {item.tax_id}", (item.name, item.tax_id))
        self.organism.addItem("Other...", None)
        self.custom_name = QLineEdit()
        self.custom_name.setPlaceholderText("Scientific name")
        self.custom_tax_id = QLineEdit()
        self.custom_tax_id.setPlaceholderText("NCBI Taxonomy ID")
        self.save_organism = QPushButton("Save organism")
        self.identifier = QLabel("—")
        self.identifier_type = QLabel("—")
        self.record_count = QLabel("0")
        form = QFormLayout()
        form.addRow("Organism", self.organism)
        form.addRow("Custom name", self.custom_name)
        form.addRow("Taxonomy ID", self.custom_tax_id)
        form.addRow("Primary identifier", self.identifier)
        form.addRow("Identifier type", self.identifier_type)
        form.addRow("Record count", self.record_count)
        self.readiness = QLabel("Open a project.")
        self.readiness.setWordWrap(True)
        self.refresh = QCheckBox("Refresh online annotations")
        self.refresh.setToolTip("Unchecked: use data already available in the cache when possible.")
        self.run_button = QPushButton("Map and annotate proteins")
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
        result_tabs.addTab(self.preview, "Catalog")
        result_tabs.addTab(self.detail, "Row details")
        self.export_table = QPushButton("Export table...")
        self.export_workbook = QPushButton("Export workbook...")
        self.open_results = QPushButton("Open results folder")
        self.explore_mapping = QPushButton("Explore across analyses...")
        self.explore_mapping.clicked.connect(self._explore_mapping)
        self.mapping_biological_context = QPushButton("Biological context...")
        self.mapping_biological_context.clicked.connect(self._mapping_biological_context)
        self.preview.itemSelectionChanged.connect(lambda: self.mapping_biological_context.setEnabled(bool(self.project and self.preview.selectedItems())))
        actions = QHBoxLayout()
        for button in (self.export_table, self.export_workbook, self.open_results, self.explore_mapping, self.mapping_biological_context):
            button.setEnabled(False)
            actions.addWidget(button)
        mapping_page = QWidget()
        layout = QVBoxLayout(mapping_page)
        layout.addWidget(QLabel("Identification and Annotation"))
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
        self.kegg_page = KEGGPage(database_manager or DatabaseManager())
        self.reactome_page = ReactomePage(database_manager or self.kegg_page.manager)
        self.mitocarta_page = MitoCartaPage(database_manager or self.kegg_page.manager)
        self.interpro_pfam_page = InterProPfamPage(database_manager or self.kegg_page.manager)
        self.string_page = StringPage(database_manager or self.kegg_page.manager)
        self.complex_page = ComplexPage(database_manager or self.kegg_page.manager)
        self.mtdna_page = MtdnaPage(database_manager or self.kegg_page.manager)
        self.proteomics_qc_page = ProteomicsQCPage()
        self.differential_analysis_page = DifferentialAnalysisPage()
        self.cross_module_explorer = CrossModuleExplorer()
        self.biological_context_page = BiologicalContextPage()
        self.experiment_comparison_page = ExperimentComparisonPage()
        self.cross_module_explorer.biological_context_requested.connect(self._open_entity_context)
        self.differential_analysis_page.biological_context_requested.connect(self.open_biological_context)
        self.experiment_comparison_page.biological_context_requested.connect(self.open_biological_context)
        self.experiment_comparison_page.derived_target_requested.connect(self.open_derived_target)
        self.cross_module_explorer.open_target_requested.connect(self.open_analysis_target)
        self.differential_analysis_page.explore_requested.connect(self.explore_context)
        self.proteomics_qc_page.explore_requested.connect(self.explore_context)
        for module_id, page in (
            ("presence_absence", self.presence_page), ("go", self.go_page),
            ("kegg", self.kegg_page), ("reactome", self.reactome_page),
            ("mitocarta", self.mitocarta_page), ("domains", self.interpro_pfam_page),
            ("string", self.string_page), ("complexes", self.complex_page),
            ("mtdna_evidence", self.mtdna_page)):
            attach_explore_action(page, module_id, self.explore_context)
            if module_id in {"kegg", "reactome", "string"}:
                attach_biological_context_action(page, module_id, self.open_biological_context)
        module_tabs = QTabWidget()
        self.module_tabs = module_tabs
        module_tabs.addTab(mapping_page, "Identification and Annotation")
        module_tabs.addTab(self.presence_page, "Presence / absence")
        module_tabs.addTab(self.go_page, "Gene Ontology")
        module_tabs.addTab(self.kegg_page, "KEGG Pathways")
        module_tabs.addTab(self.reactome_page, "Reactome")
        module_tabs.addTab(self.mitocarta_page, "MitoCarta")
        module_tabs.addTab(self.interpro_pfam_page, "InterPro / Pfam")
        module_tabs.addTab(self.string_page, "STRING")
        module_tabs.addTab(self.complex_page, "Complexes")
        module_tabs.addTab(self.mtdna_page, "mtDNA Evidence")
        module_tabs.addTab(self.proteomics_qc_page, "Proteomics QC")
        module_tabs.addTab(self.differential_analysis_page, "Differential Analysis")
        module_tabs.addTab(self.cross_module_explorer, "Cross-module Explorer")
        module_tabs.addTab(self.biological_context_page, "Biological Context")
        module_tabs.addTab(self.experiment_comparison_page, "Experiment Comparison")
        outer_layout = QVBoxLayout(self)
        outer_layout.addWidget(module_tabs)

    def explore_context(self, module_id: str, run_id: str, feature_id: str) -> None:
        self.module_tabs.setCurrentWidget(self.cross_module_explorer)
        self.cross_module_explorer.explore_feature(module_id, run_id, feature_id)

    def open_biological_context(self, module_id: str, run_id: str,
                                identifier: str, source_row=None) -> None:
        if not self.project:
            return
        self.module_tabs.setCurrentWidget(self.biological_context_page)
        try:
            self.biological_context_page.open_entity(module_id, run_id, identifier, source_row)
        except Exception as error:
            self.biological_context_page.notice.setText(str(error))
            QMessageBox.warning(self, "Biological context unavailable", str(error))

    def _open_entity_context(self, context) -> None:
        identifier = context.feature_id or context.original_identifier or (
            context.uniprot_accessions[0] if len(context.uniprot_accessions) == 1 else "")
        self.open_biological_context(context.source_module, context.source_run_id,
            identifier, context.source_row)

    def _mapping_biological_context(self) -> None:
        if not self.project or self.preview.currentRow() < 0:
            return
        from ..core.mapping_analysis import read_mapping_outputs
        try:
            outputs = read_mapping_outputs(self.project)
            headers = {self.preview.horizontalHeaderItem(column).text(): column
                for column in range(self.preview.columnCount())}
            row = self.preview.currentRow()
            source = self.preview.item(row, headers["source_row"]).text() if "source_row" in headers else ""
            identifier = next((self.preview.item(row, headers[name]).text() for name in
                ("original_id", "uniprot_accession", "gene_symbol")
                if name in headers and self.preview.item(row, headers[name]) and
                self.preview.item(row, headers[name]).text()), "")
            self.open_biological_context("mapping", str(outputs.metadata.get("run_id", "")),
                identifier or source, int(source) if source.isdigit() else None)
        except (OSError, RuntimeError, KeyError) as error:
            QMessageBox.warning(self, "Biological context unavailable", str(error))
    def _explore_mapping(self) -> None:
        row = self.preview.currentRow()
        if row < 0 or not self.project:
            return
        headers = {self.preview.horizontalHeaderItem(index).text(): index
            for index in range(self.preview.columnCount())}
        column = headers.get("original_id", headers.get("uniprot_accession"))
        if column is None or not self.preview.item(row, column):
            return
        value = self.preview.item(row, column).text()
        run_id = str(self.project.config.get("mapping_run_id", ""))
        from ..core.mapping_analysis import read_mapping_outputs
        try:
            run_id = str(read_mapping_outputs(self.project).metadata.get("run_id", run_id))
        except RuntimeError:
            pass
        self.module_tabs.setCurrentWidget(self.cross_module_explorer)
        self.cross_module_explorer.explore_feature("mapping", run_id, value)

    def open_derived_target(self, handoff: dict) -> None:
        """Preload a manual target, preserving every destination safeguard."""
        from PySide6.QtCore import Qt, QItemSelectionModel
        from ..core.derived_set_handoff import record
        pages = {"go": self.go_page, "kegg": self.kegg_page,
            "reactome": self.reactome_page, "mitocarta": self.mitocarta_page,
            "mtdna": self.mtdna_page, "domains": self.interpro_pfam_page,
            "string": self.string_page, "complexes": self.complex_page}
        page = pages.get(handoff.get("destination"))
        if page is None or self.project is None:
            return
        rows = set(handoff["source_rows"])
        if handoff["destination"] in {"go", "kegg"}:
            manual_index = page.target.findData("manual")
            if manual_index < 0:
                QMessageBox.warning(self, "Target unavailable", "The destination is not ready for manual target selection.")
                return
            page.target.setCurrentIndex(manual_index)
            control = page.manual_rows if handoff["destination"] == "go" else page.manual
            control.setText(", ".join(str(value) for value in sorted(rows)))
            found = rows
        else:
            manual_index = page.target.findData("Manual selection")
            if manual_index < 0:
                QMessageBox.warning(self, "Target unavailable", "The destination is not ready for manual target selection.")
                return
            page.target.setCurrentIndex(manual_index)
            control = page.manual
            if handoff["destination"] in {"mtdna", "complexes"}:
                frame = page.mapping_frame
                found = set()
                control.clearSelection()
                if not frame.empty:
                    for index, source in enumerate(frame.source_row):
                        if int(source) in rows:
                            control.selectionModel().select(control.model().index(index, 0),
                                QItemSelectionModel.SelectionFlag.Select | QItemSelectionModel.SelectionFlag.Rows)
                            found.add(int(source))
            else:
                found = set()
                control.clearSelection()
                for index in range(control.count()):
                    item = control.item(index)
                    value = int(item.data(Qt.ItemDataRole.UserRole))
                    if value in rows:
                        item.setSelected(True)
                        found.add(value)
        if found != rows:
            QMessageBox.warning(self, "Target unavailable",
                "The destination does not expose every selected source row. No analysis was started.")
            return
        page.derived_target_provenance = handoff.copy()
        record(self.project, handoff)
        self.module_tabs.setCurrentWidget(page)
        QMessageBox.information(self, "Review target",
            "The derived target is preloaded. Review its background, parameters and snapshot, then run the analysis manually.")
    def open_analysis_target(self, module_id: str, run_id: str,
                             target_identity: str, record_type: str = "", source_row: str = "") -> None:
        if not self.project:
            self.cross_module_explorer.status.setText("Open a project before navigating to a result.")
            return
        if module_id == "mapping":
            from ..core.mapping_analysis import read_mapping_outputs
            try:
                outputs = read_mapping_outputs(self.project)
                if str(outputs.metadata.get("run_id")) != str(run_id):
                    raise FileNotFoundError("Historical mapping data are not available for safe cross-module resolution.")
                self.show_outputs(outputs)
                self.module_tabs.setCurrentIndex(0)
                self.mapping_target_hint.setText(f"Loaded run: {run_id} | Target: {target_identity}")
            except (OSError, RuntimeError) as error:
                self.cross_module_explorer.status.setText(str(error))
            return
        pages = {"presence_absence": self.presence_page, "go": self.go_page,
            "kegg": self.kegg_page, "reactome": self.reactome_page,
            "mitocarta": self.mitocarta_page, "domains": self.interpro_pfam_page,
            "string": self.string_page, "complexes": self.complex_page,
            "mtdna_evidence": self.mtdna_page, "proteomics_qc": self.proteomics_qc_page,
            "differential": self.differential_analysis_page}
        page = pages.get(module_id)
        if page is None:
            self.cross_module_explorer.status.setText(f"No analysis page for {module_id}.")
            return
        previous_page = self.module_tabs.currentWidget()
        self.module_tabs.setCurrentWidget(page)
        try:
            adapter = default_registry().get(module_id)
            open_persisted_target(self.project, adapter, page, run_id, target_identity, source_row=source_row)
        except (OSError, KeyError, RuntimeError, ValueError) as error:
            self.module_tabs.setCurrentWidget(previous_page)
            self.cross_module_explorer.status.setText(f"Could not open {module_id} run {run_id}: {error}")
            return
        self.module_tabs.setCurrentWidget(page)
        self.cross_module_explorer.status.setText(
            f"Opened {module_id} run {run_id} for target {target_identity}.")
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
        self.kegg_page.set_project(project)
        self.reactome_page.set_project(project)
        self.mitocarta_page.set_project(project)
        self.interpro_pfam_page.set_project(project)
        self.string_page.set_project(project)
        self.complex_page.set_project(project)
        self.mtdna_page.set_project(project)
        self.proteomics_qc_page.set_project(project)
        self.differential_analysis_page.set_project(project)
        self.cross_module_explorer.set_project(project)
        self.biological_context_page.set_project(project)
        self.experiment_comparison_page.set_project(project)
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
        self.readiness.setText("Running..." if running else mapping_readiness(self.project).reason)

    def show_outputs(self, outputs: MappingOutputs) -> None:
        meta = outputs.metadata
        self.summary.setText(
            f"Input IDs: {meta.get('input_count', 0)} | Unique IDs: {meta.get('unique_id_count', 0)} | "
            f"Uniquely mapped: {meta.get('mapped_unique_count', 0)} | Ambiguous: {meta.get('ambiguous_count', 0)} | "
            f"Unmapped: {meta.get('unmapped_count', 0)} | Organism mismatches: {meta.get('organism_mismatch_count', 0)}"
        )
        frame = outputs.catalog.head(200)
        self.preview.setRowCount(len(frame)); self.preview.setColumnCount(len(frame.columns))
        self.preview.setHorizontalHeaderLabels([str(value) for value in frame.columns])
        for row, values in enumerate(frame.itertuples(index=False, name=None)):
            for column, value in enumerate(values):
                self.preview.setItem(row, column, QTableWidgetItem("" if value is None else str(value)))
        for button in (self.export_table, self.export_workbook, self.open_results, self.explore_mapping, self.mapping_biological_context):
            button.setEnabled(True)
        self.mapping_biological_context.setEnabled(bool(self.preview.selectedItems()))

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
