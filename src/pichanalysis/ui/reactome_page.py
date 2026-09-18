from __future__ import annotations

from pathlib import Path
import shutil
from typing import Callable

import pandas as pd
from PySide6.QtCore import QThread, Qt, QUrl, Signal
from PySide6.QtGui import QDesktopServices, QPixmap
from PySide6.QtWidgets import (
    QAbstractItemView, QComboBox, QDoubleSpinBox, QFileDialog, QFormLayout, QGroupBox,
    QHBoxLayout, QLabel, QListWidget, QListWidgetItem, QMessageBox, QProgressBar,
    QPushButton, QScrollArea, QSpinBox, QTableWidget, QTableWidgetItem, QTabWidget, QVBoxLayout, QWidget,
)

from ..core.database_manager import DatabaseManager
from ..core.organism import get_organism
from ..core.project import Project
from ..core.r_runtime import RRuntime
from ..core.reactome_analysis import (
    MissingReactomeOutputError, ReactomeOutputs, ReactomeParameters,
    TargetOutsideBackgroundError, export_reactome_artifact, list_reactome_runs,
    reactome_pathway_proteins, reactome_protein_pathways, reactome_readiness,
    read_reactome_outputs, run_reactome_analysis,
)


class ReactomeAnalysisWorker(QThread):
    succeeded = Signal(object)
    target_outside_background = Signal(dict)
    failed = Signal(str)

    def __init__(self, project: Project, manager: DatabaseManager, runtime: RRuntime,
        run_id: str, parameters: ReactomeParameters,
        runner: Callable = run_reactome_analysis) -> None:
        super().__init__()
        self.project, self.manager, self.runtime = project, manager, runtime
        self.run_id, self.parameters, self.runner = run_id, parameters, runner

    def run(self) -> None:
        try:
            self.succeeded.emit(self.runner(self.project, self.manager, self.runtime,
                run_id=self.run_id, parameters=self.parameters))
        except TargetOutsideBackgroundError as error:
            self.target_outside_background.emit(error.details)
        except Exception as error:
            self.failed.emit(str(error))


class ReactomePage(QWidget):
    run_requested = Signal(dict)
    open_database_requested = Signal()

    def __init__(self, manager: DatabaseManager) -> None:
        super().__init__()
        self.manager = manager
        self.project: Project | None = None
        self.outputs: ReactomeOutputs | None = None
        self._current_frames: dict[QTableWidget, pd.DataFrame] = {}
        self._diagram_pixmap = QPixmap()
        self._diagram_scale = 1.0

        self.database = QLabel()
        self.database.setWordWrap(True)
        self.open_database = QPushButton("Open Database Manager")
        self.target = QComboBox()
        self.background = QComboBox()
        self.background_help = QLabel("Background represents the proteins that could have been observed in this experiment.")
        self.background_help.setWordWrap(True)
        self.manual = QListWidget()
        self.manual.setSelectionMode(QAbstractItemView.SelectionMode.MultiSelection)
        self.manual.setMinimumHeight(100)
        self.manual_help = QLabel("Select experimental entities by their stable source row. Used when Manual selection is chosen.")
        self.manual_help.setWordWrap(True)
        self.presence_help = QLabel()
        self.presence_help.setWordWrap(True)

        self.advanced = QGroupBox("Advanced options")
        self.advanced.setCheckable(True)
        self.advanced.setChecked(False)
        self.fdr = QDoubleSpinBox()
        self.fdr.setDecimals(4)
        self.fdr.setRange(0.0001, 1.0)
        self.fdr.setValue(ReactomeParameters().fdr_cutoff)
        self.minimum = QSpinBox()
        self.minimum.setRange(1, 100000)
        self.minimum.setValue(ReactomeParameters().minimum_overlap)
        self.top_n = QSpinBox()
        self.top_n.setRange(1, 1000)
        self.top_n.setValue(ReactomeParameters().top_n)
        advanced_form = QFormLayout(self.advanced)
        advanced_form.addRow("FDR threshold", self.fdr)
        advanced_form.addRow("Minimum overlap", self.minimum)
        advanced_form.addRow("Top N pathways in plots", self.top_n)
        self.advanced.toggled.connect(self._toggle_advanced)
        self._toggle_advanced(False)

        self.method_help = QLabel(
            "Frequency: Shows how many proteins from the selected set are associated with each Reactome pathway.\n"
            "Enrichment: Tests whether a pathway is represented more often than expected relative to the selected experimental background.\n"
            "FDR: Corrects for testing many pathways simultaneously."
        )
        self.method_help.setWordWrap(True)
        self.readiness = QLabel()
        self.readiness.setWordWrap(True)
        self.run_button = QPushButton("5. Run Reactome analysis")
        self.progress = QProgressBar()
        self.progress.setRange(0, 0)
        self.progress.hide()
        self.running_text = QLabel()

        config = QWidget()
        form = QFormLayout()
        form.addRow("2. Target set", self.target)
        form.addRow("Manual selection", self.manual)
        form.addRow("3. Background", self.background)
        config_layout = QVBoxLayout(config)
        config_layout.addWidget(QLabel("1. Database"))
        config_layout.addWidget(self.database)
        config_layout.addWidget(self.open_database)
        config_layout.addLayout(form)
        config_layout.addWidget(self.manual_help)
        config_layout.addWidget(self.presence_help)
        config_layout.addWidget(self.background_help)
        config_layout.addWidget(QLabel("4. Analysis settings"))
        config_layout.addWidget(self.advanced)
        config_layout.addWidget(self.method_help)
        config_layout.addWidget(self.readiness)
        config_layout.addWidget(self.run_button)
        config_layout.addWidget(self.progress)
        config_layout.addWidget(self.running_text)

        self.summary = QTableWidget()
        self.frequency = QTableWidget()
        self.enrichment = QTableWidget()
        self.enrichment_filter = QComboBox()
        self.enrichment_filter.addItem("All tested", "all")
        self.enrichment_filter.addItem("Significant only", "significant")
        enrichment_page = QWidget()
        enrichment_layout = QVBoxLayout(enrichment_page)
        enrichment_layout.addWidget(self.enrichment_filter)
        self.no_significant = QLabel("No pathways met the selected FDR threshold.")
        self.no_significant.hide()
        enrichment_layout.addWidget(self.no_significant)
        enrichment_layout.addWidget(self.enrichment)

        self.pathway = QComboBox()
        self.pathway_members = QTableWidget()
        pathway_page = QWidget()
        pathway_layout = QVBoxLayout(pathway_page)
        pathway_layout.addWidget(QLabel("Select pathway by Reactome ID"))
        pathway_layout.addWidget(self.pathway)
        pathway_layout.addWidget(self.pathway_members)
        self.diagram_status = QLabel("No pathway selected.")
        self.diagram_status.setWordWrap(True)
        self.diagram_image = QLabel("No local diagram is available.")
        self.diagram_image.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.diagram_scroll = QScrollArea(); self.diagram_scroll.setWidget(self.diagram_image); self.diagram_scroll.setWidgetResizable(True); self.diagram_scroll.setMinimumHeight(280)
        self.diagram_open = QPushButton("Open pathway diagram"); self.diagram_fit = QPushButton("Fit to Window"); self.diagram_actual = QPushButton("Actual Size"); self.diagram_zoom_out = QPushButton("Zoom Out"); self.diagram_zoom_in = QPushButton("Zoom In"); self.diagram_export = QPushButton("Export Reactome Diagram...")
        diagram_actions = QHBoxLayout()
        for button in (self.diagram_open,self.diagram_fit,self.diagram_actual,self.diagram_zoom_out,self.diagram_zoom_in,self.diagram_export): diagram_actions.addWidget(button)
        pathway_layout.addWidget(self.diagram_status); pathway_layout.addLayout(diagram_actions); pathway_layout.addWidget(self.diagram_scroll)
        attribution = QLabel("Reactome pathway diagrams — © Reactome, licensed under CC BY 4.0."); attribution.setWordWrap(True); pathway_layout.addWidget(attribution)

        self.protein = QComboBox()
        self.protein_pathways = QTableWidget()
        protein_page = QWidget()
        protein_layout = QVBoxLayout(protein_page)
        protein_layout.addWidget(QLabel("Select experimental entity"))
        protein_layout.addWidget(self.protein)
        protein_layout.addWidget(self.protein_pathways)
        self.open_protein_diagram = QPushButton("Open pathway diagram"); protein_layout.addWidget(self.open_protein_diagram)

        self.mapping = QTableWidget()
        self.unmapped = QTableWidget()
        self.ambiguous = QTableWidget()
        self.hierarchy = QTableWidget()

        self.graph_choice = QComboBox()
        self.graph_preview = QLabel("No plots are available for this run.")
        self.graph_preview.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.graph_preview.setMinimumHeight(300)
        graph_page = QWidget()
        graph_layout = QVBoxLayout(graph_page)
        graph_layout.addWidget(self.graph_choice)
        graph_layout.addWidget(self.graph_preview, 1)

        self.history = QComboBox()
        self.history_status = QLabel()
        history_page = QWidget()
        history_layout = QVBoxLayout(history_page)
        history_layout.addWidget(QLabel("Previous Reactome runs"))
        history_layout.addWidget(self.history)
        history_layout.addWidget(self.history_status)
        history_layout.addStretch()

        self.tabs = QTabWidget()
        self.tabs.addTab(config, "Configuration")
        self.tabs.addTab(self.summary, "Summary")
        self.tabs.addTab(self.frequency, "Frequency")
        self.tabs.addTab(enrichment_page, "Enrichment")
        self.tabs.addTab(pathway_page, "Pathways")
        self.tabs.addTab(protein_page, "Proteins")
        self.tabs.addTab(self.mapping, "Mapping")
        self.tabs.addTab(self.unmapped, "Unmapped")
        self.tabs.addTab(self.ambiguous, "Ambiguous")
        self.tabs.addTab(self.hierarchy, "Hierarchy")
        self.tabs.addTab(graph_page, "Graphs")
        self.tabs.addTab(history_page, "History")

        self.export_table = QPushButton("Export current table...")
        self.export_workbook = QPushButton("Export workbook...")
        self.export_graph = QPushButton("Export graph...")
        self.open_folder = QPushButton("Open Reactome results folder")
        actions = QHBoxLayout()
        for button in (self.export_table, self.export_workbook, self.export_graph, self.open_folder):
            button.setEnabled(False)
            actions.addWidget(button)
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("Reactome — Homo sapiens"))
        layout.addWidget(self.tabs, 1)
        layout.addLayout(actions)

        self.open_database.clicked.connect(self.open_database_requested)
        self.run_button.clicked.connect(self._emit_run)
        self.target.currentIndexChanged.connect(self._selection_changed)
        self.background.currentIndexChanged.connect(self._selection_changed)
        self.enrichment_filter.currentIndexChanged.connect(self._show_enrichment)
        self.pathway.currentIndexChanged.connect(self._show_pathway_members)
        self.protein.currentIndexChanged.connect(self._show_protein_pathways)
        self.graph_choice.currentIndexChanged.connect(self._show_graph)
        self.history.currentIndexChanged.connect(self._load_history)
        self.export_table.clicked.connect(self._export_current_table)
        self.export_workbook.clicked.connect(self._export_workbook)
        self.export_graph.clicked.connect(self._export_graph)
        self.open_folder.clicked.connect(self._open_results_folder)
        self.diagram_open.clicked.connect(lambda:self._show_diagram(self.pathway.currentData()))
        self.diagram_fit.clicked.connect(self._fit_diagram); self.diagram_actual.clicked.connect(lambda:self._set_diagram_scale(1.0)); self.diagram_zoom_out.clicked.connect(lambda:self._set_diagram_scale(self._diagram_scale/1.25)); self.diagram_zoom_in.clicked.connect(lambda:self._set_diagram_scale(self._diagram_scale*1.25)); self.diagram_export.clicked.connect(self._export_diagram)
        self.protein_pathways.cellDoubleClicked.connect(self._open_protein_pathway)
        self.open_protein_diagram.clicked.connect(lambda:self._open_protein_pathway(self.protein_pathways.currentRow(),0))

    def set_project(self, project: Project | None) -> None:
        self.project = project
        self.outputs = None
        self._populate_sets()
        self._populate_manual()
        self.refresh_readiness()
        self._refresh_history()
        if project and (project.root / "analyses" / "Reactome" / "latest_metadata.json").is_file():
            try:
                self.show_outputs(read_reactome_outputs(project))
            except MissingReactomeOutputError as error:
                self.history_status.setText(f"Latest Reactome run is incomplete: {error}")

    def _toggle_advanced(self, expanded: bool) -> None:
        for child in self.advanced.findChildren(QWidget):
            child.setVisible(expanded)

    def refresh_readiness(self) -> None:
        snapshot = self.manager.reactome.active_snapshot()
        manifest = self.manager.reactome.manifest(snapshot) if snapshot else {}
        release = manifest.get("release_version")
        release_text = str(release) if release and str(release).lower() != "unknown" else "Unknown"
        status = "Ready" if self.manager.reactome.is_core_ready() else "Not installed"
        self.database.setText(
            f"Core Data status: {status}\nDiagrams: {'Ready' if self.manager.reactome.is_diagrams_ready() else 'Not installed'}\nRelease: {release_text}\nSnapshot: {snapshot.name if snapshot else '—'}"
        )
        state = reactome_readiness(self.project, self.manager)
        message = state.reason
        if not self.manager.reactome.is_core_ready():
            message = "Reactome Core Data is not installed."
        elif self.project and get_organism(self.project) and get_organism(self.project).tax_id != "9606":
            message = "Reactome analysis currently supports Homo sapiens only."
        self.readiness.setText(message)
        self.run_button.setEnabled(state.ready)

    def _populate_sets(self) -> None:
        self.target.clear()
        self.background.clear()
        base = [("All mapped entities", "All mapped entities")]
        for label, value in base:
            self.target.addItem(label, value)
        self.background.addItem("All Reactome-mapped entities in the experiment", "All mapped entities")
        presence_path = self.project.root / "analyses" / "presence_absence" / "tables" / "classification.csv" if self.project else Path()
        classifications: list[str] = []
        if presence_path.is_file():
            try:
                frame = pd.read_csv(presence_path)
                classifications = sorted(str(value) for value in frame.get("classification", pd.Series(dtype=str)).dropna().unique())
            except (OSError, ValueError):
                classifications = []
        dependent = ["Reproducibly detected", "Shared"]
        dependent.extend(value for value in classifications if value.endswith("-specific"))
        dependent.append("Sporadic")
        for label in dict.fromkeys(dependent):
            self.target.addItem(label, label)
            self.background.addItem(label, label)
            if not classifications:
                self.target.model().item(self.target.count()-1).setEnabled(False)
                self.background.model().item(self.background.count()-1).setEnabled(False)
        self.target.addItem("Manual selection", "Manual selection")
        self.background.addItem("Manual selection", "Manual selection")
        self.presence_help.setText(
            "Presence/Absence-derived sets use the existing classification results."
            if classifications else "Run Presence / absence to enable reproducibility, shared, condition-specific, and sporadic sets."
        )

    def _populate_manual(self) -> None:
        self.manual.clear()
        if not self.project:
            return
        path = self.project.root / "mapping" / "tables" / "protein_catalog.csv"
        if not path.is_file():
            return
        try:
            frame = pd.read_csv(path).drop_duplicates("source_row")
        except (OSError, ValueError, KeyError):
            return
        for row in frame.itertuples(index=False):
            source = int(getattr(row, "source_row"))
            original = str(getattr(row, "original_id", ""))
            gene = str(getattr(row, "gene_symbol", ""))
            item = QListWidgetItem(f"{source}: {gene} — {original}")
            item.setData(Qt.ItemDataRole.UserRole, source)
            self.manual.addItem(item)
        self._selection_changed()

    def _selection_changed(self) -> None:
        manual = self.target.currentData() == "Manual selection" or self.background.currentData() == "Manual selection"
        self.manual.setEnabled(manual)
        self.manual_help.setEnabled(manual)

    def parameters(self) -> dict:
        rows = tuple(int(item.data(Qt.ItemDataRole.UserRole)) for item in self.manual.selectedItems())
        return {
            "target_selection": self.target.currentData(),
            "background_selection": self.background.currentData(),
            "manual_rows": rows,
            "minimum_overlap": self.minimum.value(),
            "fdr_cutoff": self.fdr.value(),
            "top_n": self.top_n.value(),
        }

    def _emit_run(self) -> None:
        if (self.target.currentData() == "Manual selection" or self.background.currentData() == "Manual selection") and not self.manual.selectedItems():
            self.readiness.setText("Select at least one experimental entity for Manual selection.")
            return
        self.run_requested.emit(self.parameters())

    def set_running(self, running: bool) -> None:
        self.progress.setVisible(running)
        self.running_text.setText("Running Reactome analysis..." if running else "")
        self.run_button.setEnabled(not running and reactome_readiness(self.project, self.manager).ready)

    def _fill(self, table: QTableWidget, frame: pd.DataFrame) -> None:
        self._current_frames[table] = frame.copy()
        table.clear()
        table.setRowCount(len(frame))
        table.setColumnCount(len(frame.columns))
        table.setHorizontalHeaderLabels([str(column) for column in frame.columns])
        for row, values in enumerate(frame.itertuples(index=False, name=None)):
            for column, value in enumerate(values):
                table.setItem(row, column, QTableWidgetItem("" if pd.isna(value) else str(value)))

    def _clear_results(self) -> None:
        self.outputs = None
        self._current_frames.clear()
        for table in (self.summary, self.frequency, self.enrichment, self.pathway_members,
            self.protein_pathways, self.mapping, self.unmapped, self.ambiguous, self.hierarchy):
            table.clear(); table.setRowCount(0); table.setColumnCount(0)
        self.pathway.clear(); self.protein.clear(); self.graph_choice.clear()
        self.graph_preview.setText("No plots are available for this run.")

    def show_outputs(self, outputs: ReactomeOutputs) -> None:
        self.outputs = outputs
        self._fill(self.summary, outputs.summary)
        self._fill(self.frequency, outputs.frequency)
        self._fill(self.mapping, outputs.mapping)
        self._fill(self.unmapped, outputs.unmapped)
        self._fill(self.ambiguous, outputs.ambiguous)
        hierarchy = outputs.hierarchy
        if not outputs.ancestry.empty:
            hierarchy = hierarchy.merge(outputs.ancestry, left_on="child_pathway_id", right_on="pathway_id", how="left").drop(columns=["pathway_id"])
        self._fill(self.hierarchy, hierarchy)
        self._show_enrichment()
        self.no_significant.setVisible(outputs.significant.empty)

        self.pathway.blockSignals(True)
        self.pathway.clear()
        pathway_rows = outputs.membership[["Reactome_ID", "Pathway"]].drop_duplicates().sort_values("Reactome_ID")
        for row in pathway_rows.itertuples(index=False):
            self.pathway.addItem(f"{row.Reactome_ID} — {row.Pathway}", str(row.Reactome_ID))
        self.pathway.blockSignals(False)

        self.protein.blockSignals(True)
        self.protein.clear()
        entity_rows = outputs.membership.drop_duplicates("reactome_entity_key")
        for row in entity_rows.itertuples(index=False):
            gene = getattr(row, "Gene_symbol", "")
            self.protein.addItem(f"{gene} — {row.reactome_entity_key}", str(row.reactome_entity_key))
        self.protein.blockSignals(False)

        self.graph_choice.blockSignals(True)
        self.graph_choice.clear()
        pngs = {path.stem: path for path in outputs.graphs if path.suffix.lower() == ".png"}
        pdfs = {path.stem: path for path in outputs.graphs if path.suffix.lower() == ".pdf"}
        labels = {
            "reactome_frequency": "Frequency",
            "reactome_enrichment_dot": "Enrichment dot plot",
            "reactome_enrichment_fdr": "-log10(FDR) bar plot",
        }
        for stem, path in sorted(pngs.items()):
            self.graph_choice.addItem(labels.get(stem, stem), {"png": path, "pdf": pdfs.get(stem)})
        self.graph_choice.blockSignals(False)
        self._show_pathway_members()
        self._show_protein_pathways()
        self._show_graph()
        self.history_status.setText(f"Loaded Reactome run: {outputs.metadata.get('run_id', 'unknown')}")
        for button in (self.export_table, self.export_workbook, self.open_folder):
            button.setEnabled(True)
        self.export_graph.setEnabled(self.graph_choice.count() > 0)
        self._refresh_history(outputs.metadata.get("run_id"))

    def _show_enrichment(self) -> None:
        if not self.outputs:
            self._fill(self.enrichment, pd.DataFrame())
            return
        frame = self.outputs.significant if self.enrichment_filter.currentData() == "significant" else self.outputs.enrichment
        self._fill(self.enrichment, frame)

    def _show_pathway_members(self) -> None:
        if not self.outputs or not self.pathway.currentData():
            self._fill(self.pathway_members, pd.DataFrame())
            self._show_diagram(None); return
        frame=reactome_pathway_proteins(self.outputs, self.pathway.currentData()).copy();frame["Diagram_available"]=self.manager.reactome.has_diagram(self.pathway.currentData());self._fill(self.pathway_members, frame)
        self._show_diagram(self.pathway.currentData())

    def _show_protein_pathways(self) -> None:
        if not self.outputs or not self.protein.currentData():
            self._fill(self.protein_pathways, pd.DataFrame())
            return
        frame=reactome_protein_pathways(self.outputs, self.protein.currentData()).copy()
        if "Reactome_ID" in frame: frame["Diagram_available"]=frame["Reactome_ID"].map(self.manager.reactome.has_diagram)
        self._fill(self.protein_pathways, frame)

    def _show_diagram(self, reactome_id) -> None:
        path=self.manager.reactome.diagram_path(str(reactome_id)) if reactome_id else None
        self._diagram_pixmap=QPixmap(str(path)) if path else QPixmap(); self.diagram_export.setEnabled(bool(path));self.diagram_open.setEnabled(bool(path))
        if self._diagram_pixmap.isNull(): self.diagram_image.clear();self.diagram_image.setText("No local diagram is available for this pathway.");self.diagram_status.setText(f"{reactome_id or 'Pathway'}: diagram not available locally.");return
        self.diagram_status.setText(f"{reactome_id}: local diagram available.");self._fit_diagram()

    def _set_diagram_scale(self, scale: float) -> None:
        if self._diagram_pixmap.isNull(): return
        self._diagram_scale=max(0.1,min(8.0,scale));size=self._diagram_pixmap.size()*self._diagram_scale;self.diagram_image.setPixmap(self._diagram_pixmap.scaled(size,Qt.AspectRatioMode.KeepAspectRatio,Qt.TransformationMode.SmoothTransformation));self.diagram_image.resize(self.diagram_image.pixmap().size())

    def _fit_diagram(self) -> None:
        if self._diagram_pixmap.isNull(): return
        viewport=self.diagram_scroll.viewport().size();self._set_diagram_scale(min(viewport.width()/self._diagram_pixmap.width(),viewport.height()/self._diagram_pixmap.height()))

    def _export_diagram(self) -> None:
        source=self.manager.reactome.diagram_path(str(self.pathway.currentData()))
        if not source:return
        filename,_=QFileDialog.getSaveFileName(self,"Export original Reactome diagram",source.name,"PNG (*.png)")
        if filename:
            destination=Path(filename);destination=destination if destination.suffix else destination.with_suffix(".png")
            try:shutil.copy2(source,destination)
            except OSError as error:QMessageBox.warning(self,"Export failed",str(error))

    def _open_protein_pathway(self,row:int,_column:int) -> None:
        frame=self._current_frames.get(self.protein_pathways)
        if frame is None or row>=len(frame) or "Reactome_ID" not in frame:return
        identifier=str(frame.iloc[row]["Reactome_ID"]);index=self.pathway.findData(identifier)
        if index>=0:self.pathway.setCurrentIndex(index);self.tabs.setCurrentIndex(4)

    def _show_graph(self) -> None:
        data = self.graph_choice.currentData()
        if not data or not data.get("png") or not Path(data["png"]).is_file():
            self.graph_preview.clear()
            self.graph_preview.setText("No plots are available for this run.")
            return
        pixmap = QPixmap(str(data["png"]))
        self.graph_preview.setPixmap(pixmap.scaled(900, 600, Qt.AspectRatioMode.KeepAspectRatio, Qt.TransformationMode.SmoothTransformation))

    def _refresh_history(self, selected: str | None = None) -> None:
        self.history.blockSignals(True)
        self.history.clear()
        if self.project:
            for run_id in list_reactome_runs(self.project):
                metadata_path = self.project.root / "analyses" / "Reactome" / "runs" / run_id / "metadata.json"
                label = run_id
                try:
                    import json
                    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
                    label = (f"{run_id} | {metadata.get('target_definition', 'Unknown target')} | "
                        f"{metadata.get('background_definition', 'Unknown background')} | "
                        f"Release {metadata.get('reactome_release', 'Unknown')} | Snapshot {metadata.get('snapshot_id', 'Unknown')}")
                except (OSError, ValueError):
                    label = f"{run_id} | Incomplete metadata"
                self.history.addItem(label, run_id)
        if selected:
            index = self.history.findData(selected)
            if index >= 0:
                self.history.setCurrentIndex(index)
        else:
            self.history.setCurrentIndex(-1)
        self.history.blockSignals(False)

    def _load_history(self, index: int) -> None:
        if not self.project or index < 0:
            return
        run_id = self.history.itemData(index)
        if not run_id:
            return
        try:
            self.show_outputs(read_reactome_outputs(self.project, run_id))
        except MissingReactomeOutputError as error:
            self._clear_results()
            self.history_status.setText(f"Historical run is incomplete: {error}")

    def _export_current_table(self) -> None:
        widget = self.tabs.currentWidget()
        table = widget if isinstance(widget, QTableWidget) else {
            3: self.enrichment, 4: self.pathway_members, 5: self.protein_pathways,
        }.get(self.tabs.currentIndex())
        frame = self._current_frames.get(table) if table else None
        if frame is None:
            QMessageBox.information(self, "Export table", "The current panel does not contain an exportable table.")
            return
        filename, _ = QFileDialog.getSaveFileName(self, "Export current table", "reactome_table.csv", "CSV (*.csv)")
        if filename:
            try:
                frame.to_csv(filename, index=False)
            except OSError as error:
                QMessageBox.warning(self, "Export failed", str(error))

    def _export_workbook(self) -> None:
        if not self.outputs:
            return
        filename, _ = QFileDialog.getSaveFileName(self, "Export workbook", "Reactome_analysis.xlsx", "Excel workbook (*.xlsx)")
        if filename:
            try:
                export_reactome_artifact(self.outputs.workbook, Path(filename))
            except OSError as error:
                QMessageBox.warning(self, "Export failed", str(error))

    def _export_graph(self) -> None:
        data = self.graph_choice.currentData()
        if not data:
            return
        filename, selected_filter = QFileDialog.getSaveFileName(self, "Export graph", "reactome_plot.png", "PNG (*.png);;PDF (*.pdf)")
        if not filename:
            return
        suffix = ".pdf" if "PDF" in selected_filter or Path(filename).suffix.lower() == ".pdf" else ".png"
        source = data.get(suffix[1:])
        if not source:
            QMessageBox.warning(self, "Export failed", f"The {suffix.upper()} graph is not available for this run.")
            return
        destination = Path(filename)
        if not destination.suffix:
            destination = destination.with_suffix(suffix)
        try:
            export_reactome_artifact(source, destination)
        except OSError as error:
            QMessageBox.warning(self, "Export failed", str(error))

    def _open_results_folder(self) -> None:
        if self.outputs and not QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.outputs.run_root))):
            QMessageBox.warning(self, "Reactome results", "Could not open the Reactome results folder.")
