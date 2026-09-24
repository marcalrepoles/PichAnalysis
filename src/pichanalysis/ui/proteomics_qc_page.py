"""Offline Proteomics QC explorer. Scientific values come only from persisted R runs."""
from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd
from PySide6.QtCore import QThread, Qt, QUrl, Signal
from PySide6.QtGui import QDesktopServices, QPixmap
from PySide6.QtWidgets import (QAbstractItemView, QCheckBox, QComboBox, QFileDialog,
    QFormLayout, QGroupBox, QHBoxLayout, QLabel, QLineEdit, QMessageBox,
    QProgressBar, QPushButton, QScrollArea, QTableWidget, QTableWidgetItem,
    QTabWidget, QVBoxLayout, QWidget)

from ..core.column_mapping import experimental_design
from ..core.mapping_analysis import new_run_id
from ..core.proteomics_qc import (DEFAULT_TRANSFORMATION, DEFAULT_ZERO_MISSING,
    QCParameters, ProteomicsQCError, load_run, list_runs, run_proteomics_qc)
from ..core.r_runtime import RRuntime


FAMILIES = {"lfq_intensity": "LFQ intensity", "raw_intensity": "Raw intensity",
    "spectral_count": "Spectral count", "other_quantitative": "Other quantitative"}
TRANSFORMS = (("Log2 positive values", "log2_positive"),
    ("Log2(x + 1)", "log2p1"), ("No transformation", "none"))
GRAPH_LABELS = {"detected_features": "Detected features per sample",
    "missing_fraction": "Missing fraction per sample", "sample_distributions": "Sample distributions",
    "sample_boxplots": "Sample boxplots", "pearson_heatmap": "Pearson correlation",
    "spearman_heatmap": "Spearman correlation", "shared_feature_heatmap": "Shared feature counts",
    "jaccard_heatmap": "Detection Jaccard", "pca": "PCA", "condition_cv": "CV distribution",
    "feature_detection_frequency": "Feature detection frequency"}
GROUPS = {
    "Summary": (("summary", "Summary"),),
    "Samples": (("sample_metrics", "Sample metrics"), ("sample_metadata", "Sample metadata"),
        ("sample_diagnostics", "Sample diagnostics")),
    "Detection": (("sample_detection", "Sample detection"), ("feature_detection", "Feature detection"),
        ("feature_condition_detection", "Condition detection"), ("replicate_consistency", "Replicate consistency")),
    "Missingness": (("sample_missingness", "Sample missingness"), ("feature_missingness", "Feature missingness"),
        ("condition_missingness", "Condition missingness"), ("sample_distribution", "Distributions")),
    "Distributions": (("sample_distribution", "Sample distributions"),),
    "Correlations": (("pairwise_correlations", "Pairwise correlations"), ("pearson_matrix", "Pearson matrix"),
        ("spearman_matrix", "Spearman matrix"), ("shared_feature_matrix", "Shared feature counts"),
        ("within_condition_correlations", "Within-condition summary")),
    "Detection overlap": (("detection_overlap", "Detection overlap / Jaccard"), ("jaccard_matrix", "Jaccard matrix")),
    "Variability": (("feature_condition_cv", "Feature-condition CV"), ("condition_cv_summary", "Condition CV summary")),
    "PCA": (("pca_summary", "PCA status"), ("pca_scores", "Sample scores"),
        ("pca_variance", "Variance explained"), ("pca_loadings", "Feature loadings"),
        ("sample_distance_matrix", "Sample distances")),
    "Diagnostics": (("sample_diagnostics", "Sample diagnostics"), ("warnings", "Structural warnings")),
}
ERROR_TITLES = {"no_quantitative_columns": "No quantitative columns selected",
    "fewer_than_two_samples": "Fewer than two samples selected",
    "mixed_quantification_types": "Mixed incompatible quantification types",
    "missing_condition_metadata": "Missing condition metadata",
    "missing_replicate_metadata": "Missing replicate metadata",
    "invalid_numeric_value": "Invalid numeric value",
    "negative_value_incompatible": "Negative value incompatible with selected scale",
    "unsupported_transformation": "Unsupported transformation",
    "no_usable_quantitative_data": "No usable quantitative data",
    "r_execution_failure": "R execution failure",
    "missing_required_output": "Missing required output"}


class ProteomicsQCWorker(QThread):
    succeeded = Signal(object)
    failed = Signal(object)

    def __init__(self, project, runtime, run_id, parameters, runner=run_proteomics_qc):
        super().__init__()
        self.project, self.runtime, self.run_id, self.parameters = project, runtime, run_id, parameters
        self.runner = runner

    def run(self):
        try:
            self.succeeded.emit(self.runner(self.project, self.runtime, run_id=self.run_id,
                parameters=self.parameters))
        except Exception as error:
            self.failed.emit(error)


class ProteomicsQCPage(QWidget):
    explore_requested = Signal(str, str, str)

    def __init__(self, runtime=None):
        super().__init__()
        self.runtime = runtime or RRuntime()
        self.project = None
        self.outputs = None
        self._frames = {}
        self._visible_frames = {}
        self._workers = []
        self._running = False
        self._selected_graph = None
        self._fit = True
        self.samples = QTableWidget(0, 6)
        self.samples.setHorizontalHeaderLabels(["Use", "Sample", "Column", "Condition", "Replicate", "Quantification type"])
        self.samples.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.samples.itemChanged.connect(self._selection_changed)
        self.family = QLabel("No quantitative family selected")
        self.readiness = QLabel("Open a project to configure Proteomics QC.")
        self.readiness.setWordWrap(True)
        self.transformation = QComboBox()
        for label, value in TRANSFORMS:
            self.transformation.addItem(label, value)
        self.transformation.setToolTip("Transformation is used for quantitative QC metrics and plots. The original imported values are never overwritten.")
        self.zero_missing = QCheckBox("Treat zero as missing")
        self.zero_missing.setToolTip("Zero handling defines detection/missingness and is separate from quantitative transformation.")
        self.run_button = QPushButton("Run QC")
        self.run_button.clicked.connect(self._start_run)
        self.progress = QProgressBar()
        self.progress.setRange(0, 0)
        self.progress.hide()
        self.status = QLabel()
        config = QGroupBox("1. Quantitative samples")
        config_layout = QVBoxLayout(config)
        for widget in (self.samples, self.family, self.readiness):
            config_layout.addWidget(widget)
        handling = QGroupBox("2. Data handling")
        form = QFormLayout(handling)
        form.addRow("Transformation", self.transformation)
        form.addRow(self.zero_missing)
        form.addRow(QLabel("Proteomics QC does not impute missing values."))
        form.addRow(QLabel("QC metrics are descriptive. Samples are not automatically excluded."))
        run_group = QGroupBox("3. Run quality control")
        run_layout = QVBoxLayout(run_group)
        for widget in (self.run_button, self.progress, self.status):
            run_layout.addWidget(widget)
        self.history = QComboBox()
        self.history.currentIndexChanged.connect(self._load_history)
        self.history_metadata = QLabel("No run loaded.")
        self.history_metadata.setWordWrap(True)
        history_row = QHBoxLayout()
        history_row.addWidget(QLabel("Historical run"))
        history_row.addWidget(self.history, 1)
        result_box = QGroupBox("4. Results")
        result_layout = QVBoxLayout(result_box)
        result_layout.addLayout(history_row)
        result_layout.addWidget(self.history_metadata)
        self.tabs = QTabWidget()
        self.tables = {}
        self.table_tabs = {}
        self.searches = {}
        self.conditions = {}
        self.sample_detail = QLabel("Select a sample to inspect persisted QC details.")
        self.sample_detail.setWordWrap(True)
        for group, entries in GROUPS.items():
            nested = QTabWidget()
            for key, label in entries:
                page = QWidget()
                layout = QVBoxLayout(page)
                search = QLineEdit()
                search.setPlaceholderText("Search sample, condition, feature ID or identifier")
                condition = QComboBox()
                condition.addItem("All conditions", "")
                row = QHBoxLayout()
                row.addWidget(search, 1)
                row.addWidget(condition)
                layout.addLayout(row)
                table = QTableWidget()
                table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
                layout.addWidget(table, 1)
                self.tables[(group, key)] = table
                if group == "Samples" and key == "sample_metrics":
                    table.itemSelectionChanged.connect(self._show_sample_detail)
                self.searches[(group, key)] = search
                self.conditions[(group, key)] = condition
                search.textChanged.connect(lambda _text, g=group, k=key: self._filter(g, k))
                condition.currentIndexChanged.connect(lambda _index, g=group, k=key: self._filter(g, k))
                nested.addTab(page, label)
            self.table_tabs[group] = nested
            panel = QWidget()
            panel_layout = QVBoxLayout(panel)
            if group == "Samples":
                panel_layout.addWidget(self.sample_detail)
            if group == "Detection overlap":
                panel_layout.addWidget(QLabel("Jaccard similarity compares the sets of detected features and is distinct from quantitative correlation."))
            if group == "Variability":
                panel_layout.addWidget(QLabel("CV is calculated from original linear-scale values using observed replicates only. Missing values are not imputed."))
            if group == "PCA":
                self.pca_note = QLabel("PCA uses only features with quantitative values available across all selected samples. Missing values are not imputed.")
                self.pca_note.setWordWrap(True)
                panel_layout.addWidget(self.pca_note)
            if group == "Diagnostics":
                panel_layout.addWidget(QLabel("Diagnostics summarize sample-level QC metrics and do not automatically identify samples for exclusion."))
            panel_layout.addWidget(nested, 1)
            self.tabs.addTab(panel, group)
        self.graph_choice = QComboBox()
        self.graph_choice.currentIndexChanged.connect(self._show_graph)
        self.graph_preview = QLabel("Graph not available for this run.")
        self.graph_preview.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.graph_preview.setMinimumHeight(260)
        self.graph_scroll = QScrollArea()
        self.graph_scroll.setWidgetResizable(True)
        self.graph_scroll.setWidget(self.graph_preview)
        self.previous = QPushButton("Previous")
        self.next = QPushButton("Next")
        self.fit = QPushButton("Fit")
        self.actual = QPushButton("Actual size")
        self.previous.clicked.connect(lambda: self.graph_choice.setCurrentIndex(max(0, self.graph_choice.currentIndex()-1)))
        self.next.clicked.connect(lambda: self.graph_choice.setCurrentIndex(min(self.graph_choice.count()-1, self.graph_choice.currentIndex()+1)))
        self.fit.clicked.connect(lambda: self._set_fit(True))
        self.actual.clicked.connect(lambda: self._set_fit(False))
        graph_page = QWidget()
        graph_layout = QVBoxLayout(graph_page)
        graph_layout.addWidget(self.graph_choice)
        graph_layout.addWidget(self.graph_scroll, 1)
        navigation = QHBoxLayout()
        for widget in (self.previous, self.next, self.fit, self.actual):
            navigation.addWidget(widget)
        graph_layout.addLayout(navigation)
        self.tabs.addTab(graph_page, "Graphs")
        result_layout.addWidget(self.tabs, 1)
        self.export_table = QPushButton("Export current table...")
        self.export_workbook = QPushButton("Export workbook...")
        self.export_graph = QPushButton("Export graph...")
        self.open_results = QPushButton("Open Proteomics QC results folder")
        self.export_table.clicked.connect(self._export_current)
        self.export_workbook.clicked.connect(self._export_workbook)
        self.export_graph.clicked.connect(self._export_graph)
        self.open_results.clicked.connect(self._open_results)
        self.explore_button = QPushButton("Explore across analyses...")
        self.explore_button.clicked.connect(self._explore_selected)
        actions = QHBoxLayout()
        for widget in (self.export_table, self.export_workbook, self.export_graph, self.open_results, self.explore_button):
            actions.addWidget(widget)
        main = QVBoxLayout(self)
        for widget in (config, handling, run_group, result_box):
            main.addWidget(widget)
        main.addLayout(actions)
        self._update_actions()

    def _explore_selected(self):
        if not self.outputs:
            return
        table = self.tables[("Detection", "feature_detection")]
        row = table.currentRow()
        headers = {table.horizontalHeaderItem(column).text(): column
            for column in range(table.columnCount())}
        column = headers.get("Feature Id")
        if row < 0 or column is None or not table.item(row, column):
            return
        self.explore_requested.emit("proteomics_qc", self.outputs["metadata"]["run_id"],
            table.item(row, column).text())
    def set_project(self, project):
        self.project = project
        self.outputs = None
        self._frames.clear()
        self.samples.blockSignals(True)
        self.samples.setRowCount(0)
        if project:
            configured = experimental_design(project)["quantification_columns"]
            families = [entry["quantification_type"] for entry in configured]
            predominant = max(dict.fromkeys(families), key=families.count) if families else None
            self.samples.setRowCount(len(configured))
            for row, entry in enumerate(configured):
                family = str(entry["quantification_type"])
                use = QTableWidgetItem()
                use.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsUserCheckable)
                use.setCheckState(Qt.CheckState.Checked if family == predominant else Qt.CheckState.Unchecked)
                self.samples.setItem(row, 0, use)
                for col, value in enumerate((f"S{row+1:03d}", entry["column"], entry["condition"],
                    entry["replicate"], FAMILIES.get(family, family)), 1):
                    item = QTableWidgetItem(str(value or "Missing"))
                    item.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
                    self.samples.setItem(row, col, item)
                self.samples.item(row, 0).setData(Qt.ItemDataRole.UserRole, family)
        self.samples.blockSignals(False)
        self._selection_changed()
        self._refresh_history()
        self._update_actions()

    def selected_columns(self):
        return tuple(self.samples.item(row, 2).text() for row in range(self.samples.rowCount())
            if self.samples.item(row, 0).checkState() == Qt.CheckState.Checked)

    def _selected_family(self):
        return {self.samples.item(row, 0).data(Qt.ItemDataRole.UserRole)
            for row in range(self.samples.rowCount())
            if self.samples.item(row, 0).checkState() == Qt.CheckState.Checked}

    def _selection_changed(self):
        families = self._selected_family()
        family = next(iter(families)) if len(families) == 1 else None
        self.family.setText("Quantification family: " + (FAMILIES.get(family, family) if family else "not available"))
        if family:
            self.transformation.setCurrentIndex(self.transformation.findData(DEFAULT_TRANSFORMATION.get(family, "none")))
            self.zero_missing.setChecked(DEFAULT_ZERO_MISSING.get(family, False))
        if len(self.selected_columns()) < 2:
            message = "Select at least two quantitative samples."
        elif len(families) != 1 or family not in FAMILIES:
            message = "Selected samples use incompatible quantification types."
        else:
            incomplete = []
            for row in range(self.samples.rowCount()):
                if self.samples.item(row, 0).checkState() != Qt.CheckState.Checked:
                    continue
                for col, name in ((3, "condition"), (4, "replicate")):
                    if self.samples.item(row, col).text() == "Missing":
                        incomplete.append(f"{self.samples.item(row, 2).text()}: missing {name}")
            message = "; ".join(incomplete) if incomplete else "Ready for descriptive, offline quality control."
        self.readiness.setText(message)
        self.run_button.setEnabled(bool(self.project) and not self._running and message.startswith("Ready"))

    def parameters(self):
        return QCParameters(self.selected_columns(), self.transformation.currentData(), self.zero_missing.isChecked())

    def is_running(self):
        return self._running or any(worker.isRunning() for worker in self._workers)

    def _start_run(self):
        if not self.project or self.is_running() or not self.run_button.isEnabled():
            return
        worker = ProteomicsQCWorker(self.project, self.runtime, new_run_id(), self.parameters())
        self._workers.append(worker)
        self._running = True
        self.progress.show()
        self.status.setText("Running proteomics quality control...")
        self._selection_changed()
        worker.succeeded.connect(self._run_succeeded)
        worker.failed.connect(self._run_failed)
        worker.finished.connect(lambda w=worker: self._worker_finished(w))
        worker.finished.connect(worker.deleteLater)
        worker.start()

    def _worker_finished(self, worker):
        if worker in self._workers:
            self._workers.remove(worker)
        self._running = False
        self.progress.hide()
        self._selection_changed()

    def _run_succeeded(self, outputs):
        self.show_outputs(outputs)
        self.status.setText("Proteomics QC run completed.")

    def _run_failed(self, error):
        code = getattr(error, "code", "proteomics_qc_error")
        detail = getattr(error, "details", {}) or {}
        message = str(error)
        if code == "invalid_numeric_value":
            examples = detail.get("values", [])
            if examples:
                first = examples[0]
                message += f"\nSample: {first.get('sample')}\nSource row: {first.get('source_row')}\nInvalid value count: {detail.get('invalid_value_count')}"
        self.status.setText(f"{ERROR_TITLES.get(code, 'Proteomics QC error')}: {message}")
        QMessageBox.warning(self, ERROR_TITLES.get(code, "Proteomics QC error"), message)

    @staticmethod
    def _fill(table, frame):
        table.clear()
        table.setRowCount(len(frame))
        table.setColumnCount(len(frame.columns))
        table.setHorizontalHeaderLabels([str(column).replace("_", " ").title() for column in frame.columns])
        for row, values in enumerate(frame.itertuples(index=False, name=None)):
            for col, value in enumerate(values):
                table.setItem(row, col, QTableWidgetItem("NA" if value is None or str(value) == "" else str(value)))

    def _filter(self, group, key):
        if key not in self._frames:
            return
        frame = self._frames[key]
        term = self.searches[(group, key)].text().strip()
        condition = self.conditions[(group, key)].currentData()
        if condition and "condition" in frame:
            frame = frame[frame["condition"].astype(str) == condition]
        if term:
            mask = frame.astype(str).apply(lambda column: column.str.contains(term, case=False, regex=False)).any(axis=1)
            frame = frame[mask]
        self._visible_frames[(group, key)] = frame
        self._fill(self.tables[(group, key)], frame)

    def show_outputs(self, outputs):
        self.outputs = outputs
        self._frames = dict(outputs["tables"])
        run = outputs["run_root"]
        metadata = pd.read_csv(run / "input/sample_metadata.csv", dtype=str, keep_default_na=False)
        self._frames["sample_metadata"] = metadata
        diagnostics = self._frames["sample_diagnostics"]
        self._frames["sample_metrics"] = metadata.merge(diagnostics, left_on="sample_id", right_on="sample", how="left")
        for (group, key), combo in self.conditions.items():
            combo.blockSignals(True)
            combo.clear()
            combo.addItem("All conditions", "")
            for value in metadata.condition.drop_duplicates():
                combo.addItem(value, value)
            combo.blockSignals(False)
            self.searches[(group, key)].clear()
            self._filter(group, key)
        pca = self._frames["pca_summary"].iloc[0]
        status = str(pca.get("status", ""))
        self.pca_note.setText("PCA uses only features with quantitative values available across all selected samples. Missing values are not imputed. "
            + f"Centered: {'Yes' if str(pca.get('center')) == 'TRUE' else 'No'}. Unit-variance scaled: {'Yes' if str(pca.get('scale')) == 'TRUE' else 'No'}. "
            + ("" if status == "Ready" else f"PCA is not available for this run. {status}"))
        self.graph_choice.blockSignals(True)
        self.graph_choice.clear()
        plots = {path.stem: path for path in outputs["plots"]}
        for name, label in GRAPH_LABELS.items():
            self.graph_choice.addItem(label, plots.get(name))
        self.graph_choice.blockSignals(False)
        self._show_graph()
        meta = outputs["metadata"]
        self.history_metadata.setText(f"Run ID: {meta.get('run_id')} | Date/time: {meta.get('created_at')} | "
            f"Quantification type: {FAMILIES.get(meta.get('quantification_type'), meta.get('quantification_type'))} | "
            f"Transformation: {meta.get('transformation')} | Zero treated as missing: {meta.get('zero_is_missing')} | "
            f"Samples: {meta.get('sample_count')} | Conditions: {', '.join(meta.get('conditions', []))}")
        self._refresh_history(meta.get("run_id"))
        self._update_actions()

    def _show_sample_detail(self):
        table = self.tables[("Samples", "sample_metrics")]
        rows = table.selectionModel().selectedRows()
        frame = self._visible_frames.get(("Samples", "sample_metrics"))
        if not rows or frame is None or rows[0].row() >= len(frame):
            self.sample_detail.setText("Select a sample to inspect persisted QC details.")
            return
        record = frame.iloc[rows[0].row()]
        sample_id = str(record.get("sample_id", ""))
        detail = dict(record)
        for name in ("sample_detection", "sample_missingness"):
            source = self._frames.get(name)
            if source is not None and "sample" in source:
                found = source[source["sample"] == sample_id]
                if not found.empty:
                    detail.update(found.iloc[0].to_dict())
        fields = ("sample_id", "condition", "replicate", "detected_features", "missing_features",
            "missing_fraction", "zero_count", "na_count", "median_transformed_value",
            "IQR_transformed_value", "median_within_condition_pearson",
            "median_within_condition_spearman", "median_within_condition_jaccard", "PCA_PC1", "PCA_PC2")
        self.sample_detail.setText(" | ".join(f"{field.replace('_', ' ').title()}: {detail.get(field) or 'NA'}" for field in fields))
    def _refresh_history(self, selected=None):
        self.history.blockSignals(True)
        self.history.clear()
        if self.project:
            for run_id in list_runs(self.project):
                try:
                    meta = load_run(self.project, run_id)["metadata"]
                    label = (f"{meta.get('created_at')} | {FAMILIES.get(meta.get('quantification_type'), meta.get('quantification_type'))} | "
                        f"{meta.get('transformation')} | zero missing: {meta.get('zero_is_missing')} | "
                        f"{meta.get('sample_count')} samples | {', '.join(meta.get('conditions', []))} | {run_id}")
                except Exception:
                    label = f"{run_id} — Incomplete run"
                self.history.addItem(label, run_id)
        if selected:
            index = self.history.findData(selected)
            if index >= 0:
                self.history.setCurrentIndex(index)
        self.history.blockSignals(False)

    def _load_history(self):
        if not self.project or self.is_running():
            return
        run_id = self.history.currentData()
        if not run_id:
            return
        try:
            self.show_outputs(load_run(self.project, run_id))
            self.status.setText(f"Loaded historical run {run_id} offline.")
        except Exception as error:
            self.status.setText(f"Incomplete run: {error}")

    def _set_fit(self, fit):
        self._fit = fit
        self._show_graph()

    def _show_graph(self):
        source = self.graph_choice.currentData()
        if source and Path(source).is_file():
            self._selected_graph = Path(source)
            pixmap = QPixmap(str(source))
            if self._fit:
                pixmap = pixmap.scaled(1000, 600, Qt.AspectRatioMode.KeepAspectRatio,
                    Qt.TransformationMode.SmoothTransformation)
            self.graph_preview.setPixmap(pixmap)
        else:
            self._selected_graph = None
            self.graph_preview.clear()
            self.graph_preview.setText("Graph not available for this run.")
        self._update_actions()

    def _current_frame(self):
        group = self.tabs.tabText(self.tabs.currentIndex())
        nested = self.table_tabs.get(group)
        if not nested:
            return None
        key = GROUPS[group][nested.currentIndex()][0]
        return self._visible_frames.get((group, key))

    def _export_current(self):
        frame = self._current_frame()
        if frame is None:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export current table", "Proteomics_QC_table.csv", "CSV (*.csv)")
        if path:
            frame.to_csv(path, index=False)

    def _export_workbook(self):
        source = self.outputs["workbook"] if self.outputs else None
        if not source:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export workbook", "Proteomics_QC.xlsx", "Excel (*.xlsx)")
        if path:
            shutil.copy2(source, path)

    def _export_graph(self):
        source = self._selected_graph
        if not source:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export graph", source.name,
            "PNG (*.png);;PDF (*.pdf)")
        if path:
            chosen = Path(path)
            original = source.with_suffix(chosen.suffix.lower())
            if original.suffix.lower() not in {".png", ".pdf"} or not original.is_file():
                QMessageBox.warning(self, "Export graph", "The requested graph format is not available for this run.")
                return
            shutil.copy2(original, chosen)

    def _open_results(self):
        if self.project:
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.project.root / "analyses/Proteomics_QC")))

    def _update_actions(self):
        has_run = bool(self.outputs)
        self.export_table.setEnabled(has_run)
        self.export_workbook.setEnabled(has_run and bool(self.outputs["workbook"]))
        self.export_graph.setEnabled(has_run and bool(self._selected_graph))
        self.open_results.setEnabled(bool(self.project))
