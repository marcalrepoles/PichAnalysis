"""Two-stage, offline Differential Analysis interface over frozen R artifacts."""
from __future__ import annotations

import json
import shutil
import subprocess
from ..core.external_process import run_external
from pathlib import Path

import pandas as pd
from PySide6.QtCore import QThread, Qt, QUrl, Signal
from PySide6.QtGui import QDesktopServices, QPixmap
from PySide6.QtWidgets import (QAbstractItemView, QCheckBox, QComboBox, QDoubleSpinBox,
    QFileDialog, QFormLayout, QGroupBox, QHBoxLayout, QLabel, QLineEdit,
    QMessageBox, QPlainTextEdit, QProgressBar, QPushButton, QScrollArea,
    QSpinBox, QTableWidget, QTableWidgetItem, QTabWidget, QVBoxLayout, QWidget)

from ..core.column_mapping import experimental_design
from ..core.differential_analysis import (DifferentialParameters, list_preparation_runs,
    list_runs as statistics_runs, load_run as load_statistics, run_differential_analysis)
from ..core.differential_preparation import (PreparationParameters, list_runs as preparation_runs,
    load_run as load_preparation, run_differential_preparation)
from ..core.mapping_analysis import new_run_id
from ..core.r_runtime import RRuntime


FAMILIES = {"lfq_intensity": "LFQ intensity", "raw_intensity": "Raw intensity",
    "other_quantitative": "Other continuous quantitative", "spectral_count": "Spectral count"}
PREP_TABLES = {"Eligibility": "feature_eligibility", "Missingness patterns": "feature_missingness_patterns",
    "Qualitative candidates": "qualitative_detection_candidates", "Excluded features": "excluded_features",
    "Normalization": "sample_normalization", "Imputation": "imputed_cells"}
RESULT_FIELDS = [("Feature ID", "feature_id"), ("Identifier", "display_identifier"),
    ("Gene", "gene_symbol"), ("UniProt", "uniprot_accession"),
    ("Condition A", "condition_A"), ("Condition B", "condition_B"),
    ("n A", "n_A_model"), ("n B", "n_B_model"),
    ("Mean A", "mean_A_prepared"), ("Mean B", "mean_B_prepared"),
    ("Effect", "effect"), ("Log2FC", "log2FC"), ("Fold change", "fold_change"),
    ("Percent change", "percent_change"), ("Moderated t", "moderated_t"),
    ("P-value", "P.Value"), ("FDR", "adj.P.Val"), ("B", "B"),
    ("95% CI low", "CI_95_low"), ("95% CI high", "CI_95_high"),
    ("Residual df", "df_residual"), ("Imputed cells", "imputed_cell_count"),
    ("Result class", "result_class"), ("Exclusion reason", "test_exclusion_reason")]


class DifferentialWorker(QThread):
    succeeded = Signal(object)
    failed = Signal(object)

    def __init__(self, runner, project, runtime, run_id, parameters):
        super().__init__()
        self.runner, self.project, self.runtime = runner, project, runtime
        self.run_id, self.parameters = run_id, parameters

    def run(self):
        try:
            self.succeeded.emit(self.runner(self.project, self.runtime, run_id=self.run_id,
                parameters=self.parameters))
        except Exception as error:
            self.failed.emit(error)


class DifferentialAnalysisPage(QWidget):
    explore_requested = Signal(str, str, str)

    def __init__(self, runtime=None):
        super().__init__()
        self.runtime = runtime or RRuntime()
        self.project = None
        self.preparation = None
        self.statistics = None
        self.worker = None
        self._running = False
        self._fit = True
        self._frames = {}
        self._current_table_key = None
        self.stage_tabs = QTabWidget()
        self._build_preparation()
        self._build_statistics()
        self._build_results()
        layout = QVBoxLayout(self)
        layout.addWidget(self.stage_tabs)
        self._update_state()

    @staticmethod
    def _note(text):
        label = QLabel(text)
        label.setWordWrap(True)
        return label

    @staticmethod
    def _fill(table, frame, fields=None):
        fields = fields or [(column.replace("_", " ").title(), column) for column in frame.columns]
        table.setSortingEnabled(False)
        table.setRowCount(len(frame))
        table.setColumnCount(len(fields))
        table.setHorizontalHeaderLabels([label for label, _ in fields])
        for row, values in enumerate(frame.to_dict("records")):
            for column, (_, key) in enumerate(fields):
                value = str(values.get(key, ""))
                item = QTableWidgetItem(value if value and value.lower() not in ("nan", "na") else "NA")
                item.setData(Qt.ItemDataRole.UserRole, str(values.get("feature_id", "")))
                table.setItem(row, column, item)
        table.setSortingEnabled(True)

    def _build_preparation(self):
        page = QWidget()
        layout = QVBoxLayout(page)
        layout.addWidget(self._note("Differential Preparation — Prepare quantitative data, define missing-value handling, and create a frozen input matrix for differential statistics."))
        conditions = QHBoxLayout()
        self.condition_a, self.condition_b = QComboBox(), QComboBox()
        for name, widget in (("Condition A", self.condition_a), ("Condition B", self.condition_b)):
            conditions.addWidget(QLabel(name)); conditions.addWidget(widget)
            widget.currentIndexChanged.connect(self._conditions_changed)
        layout.addLayout(conditions)
        self.direction = self._note("Effect direction: Condition A − Condition B")
        layout.addWidget(self.direction)
        self.samples = QTableWidget(0, 6)
        self.samples.setHorizontalHeaderLabels(["Use", "Sample", "Column", "Condition", "Replicate", "Quantification type"])
        self.samples.itemChanged.connect(self._update_state)
        layout.addWidget(self.samples)
        self.family = self._note("Quantification family: not selected")
        self.preparation_readiness = self._note("Open a project to prepare differential data.")
        layout.addWidget(self.family); layout.addWidget(self.preparation_readiness)
        layout.addWidget(self._note("Sample inclusion is explicit. Proteomics QC does not automatically exclude samples."))
        handling = QGroupBox("Data handling and feature eligibility")
        form = QFormLayout(handling)
        self.transformation = QComboBox()
        self.transformation.addItem("Log2 positive values", "log2_positive")
        self.transformation.addItem("No transformation", "none")
        self.transformation.activated.connect(lambda *_: self.transformation.setProperty("user_changed", True))
        self.zero_missing = QCheckBox("Treat zero as missing")
        self.zero_missing.clicked.connect(lambda *_: self.zero_missing.setProperty("user_changed", True))
        self.normalization = QComboBox()
        self.normalization.addItem("No additional normalization", "none")
        self.normalization.addItem("Median centering", "median_center")
        self.minimum_observed = QSpinBox(); self.minimum_observed.setRange(1, 1000); self.minimum_observed.setValue(2)
        self.other_confirmed = QCheckBox("I confirm that Other quantitative values are continuous")
        self.other_confirmed.toggled.connect(self._update_state)
        form.addRow("Transformation", self.transformation)
        form.addRow(self.zero_missing)
        form.addRow(self._note("Zero handling defines detection and missingness before feature eligibility."))
        form.addRow("Normalization", self.normalization)
        form.addRow(self._note("Median centering is applied after transformation and before imputation using observed values only."))
        form.addRow("Minimum observed replicates per condition", self.minimum_observed)
        form.addRow(self._note("Feature eligibility is determined before imputation."))
        form.addRow(self.other_confirmed)
        layout.addWidget(handling)
        imputation = QGroupBox("Optional imputation")
        form = QFormLayout(imputation)
        self.imputation = QComboBox()
        for label, value in (("No imputation", "none"), ("MinProb", "MinProb"),
            ("QRILC", "QRILC"), ("KNN", "KNN")):
            self.imputation.addItem(label, value)
        self.imputation.currentIndexChanged.connect(self._imputation_changed)
        self.seed = QSpinBox(); self.seed.setRange(0, 2147483647); self.seed.setValue(12345)
        self.q = QDoubleSpinBox(); self.q.setRange(0.000001, .999999); self.q.setDecimals(6); self.q.setValue(.01)
        self.sigma = QDoubleSpinBox(); self.sigma.setRange(.000001, 100000); self.sigma.setDecimals(6); self.sigma.setValue(1)
        self.margin = QComboBox(); self.margin.addItem("Rows", 1); self.margin.addItem("Columns", 2)
        self.margin.setCurrentIndex(1)
        self.k = QSpinBox(); self.k.setRange(1, 100000); self.k.setValue(10)
        form.addRow("Imputation method", self.imputation)
        for key, label, widget in (("seed", "Random seed", self.seed), ("q", "q", self.q),
            ("sigma", "sigma", self.sigma), ("margin", "MARGIN", self.margin), ("k", "k", self.k)):
            form.addRow(label, widget)
        form.addRow(self._note("The software does not automatically infer whether missing values are MAR, MNAR, below the detection limit, or biological absence."))
        form.addRow(self._note("Imputation is applied only to partial missing values in features already eligible for continuous analysis. A condition with no observed values is never imputed."))
        layout.addWidget(imputation)
        self.prepare_button = QPushButton("Prepare differential data")
        self.prepare_button.clicked.connect(self._start_preparation)
        self.progress = QProgressBar(); self.progress.setRange(0, 0); self.progress.hide()
        self.status = self._note("")
        layout.addWidget(self.prepare_button); layout.addWidget(self.progress); layout.addWidget(self.status)
        self.prep_history = QComboBox(); self.prep_history.currentIndexChanged.connect(self._load_preparation_history)
        layout.addWidget(QLabel("Preparation history")); layout.addWidget(self.prep_history)
        self.prep_summary = QPlainTextEdit(); self.prep_summary.setReadOnly(True)
        layout.addWidget(self.prep_summary)
        self.prep_audit = QTabWidget(); self.prep_tables = {}
        for label, key in PREP_TABLES.items():
            table = QTableWidget(); self.prep_tables[key] = table
            self.prep_audit.addTab(table, label)
        layout.addWidget(self.prep_audit)
        self.stage_tabs.addTab(page, "Prepare data")
        self._imputation_changed()

    def _build_statistics(self):
        page = QWidget(); layout = QVBoxLayout(page)
        layout.addWidget(QLabel("Differential Statistics"))
        self.parent_choice = QComboBox(); self.parent_choice.currentIndexChanged.connect(self._parent_changed)
        layout.addWidget(QLabel("Ready Preparation run")); layout.addWidget(self.parent_choice)
        self.parent_summary = self._note("Select a Ready Preparation run.")
        layout.addWidget(self.parent_summary)
        self.scale = QComboBox(); self.scale.addItem("Choose prepared scale", None)
        self.scale.addItem("Log2", "log2"); self.scale.addItem("Continuous", "continuous")
        self.scale.currentIndexChanged.connect(self._scale_changed)
        form = QFormLayout(); form.addRow("Prepared scale", self.scale)
        form.addRow(self._note("Choose Log2 only if the prepared quantitative values are already on a log2 scale."))
        self.trend = QCheckBox("Use intensity trend in empirical Bayes moderation")
        self.robust = QCheckBox("Use robust empirical Bayes moderation")
        form.addRow(self.trend); form.addRow(self.robust)
        self.fdr = QDoubleSpinBox(); self.fdr.setRange(.000001, 1); self.fdr.setDecimals(6); self.fdr.setValue(.05)
        self.effect = QDoubleSpinBox(); self.effect.setRange(0, 1000000); self.effect.setDecimals(4)
        self.effect_label = QLabel("Minimum absolute effect")
        form.addRow("FDR threshold", self.fdr); form.addRow(self.effect_label, self.effect)
        form.addRow(self._note("The effect threshold changes result classification only. It does not change p-values or BH-adjusted FDR values."))
        self.top_n = QSpinBox(); self.top_n.setRange(1, 100000); self.top_n.setValue(20)
        self.volcano_labels = QSpinBox(); self.volcano_labels.setRange(0, 100000); self.volcano_labels.setValue(20)
        form.addRow("Top N in plots", self.top_n); form.addRow("Volcano labels", self.volcano_labels)
        layout.addLayout(form)
        self.statistics_button = QPushButton("Run differential statistics")
        self.statistics_button.clicked.connect(self._start_statistics)
        layout.addWidget(self.statistics_button); layout.addStretch()
        self.stage_tabs.addTab(page, "Run differential statistics")

    def _build_results(self):
        page = QWidget(); layout = QVBoxLayout(page)
        self.stat_history = QComboBox(); self.stat_history.currentIndexChanged.connect(self._load_statistics_history)
        layout.addWidget(QLabel("Statistics history")); layout.addWidget(self.stat_history)
        self.result_tabs = QTabWidget(); layout.addWidget(self.result_tabs, 1)
        self.result_tables = {}
        self.summary_table = QTableWidget()
        self.result_tabs.addTab(self.summary_table, "Summary")
        results_page = QWidget(); results_layout = QVBoxLayout(results_page)
        filters = QHBoxLayout()
        self.result_filter = QComboBox()
        for value in ("All", "Tested", "FDR significant", "Combined significant", "Higher in Condition A",
            "Higher in Condition B", "Not significant", "Not tested"):
            self.result_filter.addItem(value)
        self.effect_filter = QComboBox()
        for value in ("Any effect", "Passes effect threshold", "Does not pass effect threshold"):
            self.effect_filter.addItem(value)
        self.search = QLineEdit(); self.search.setPlaceholderText("Search feature ID, identifier, gene or UniProt")
        self.sort = QComboBox()
        for label, key in (("Original order", None), ("FDR", "adj.P.Val"), ("P-value", "P.Value"),
            ("Effect", "effect"), ("Absolute effect", "abs_effect"), ("B", "B"),
            ("Residual df", "df_residual")):
            self.sort.addItem(label, key)
        for widget in (self.result_filter, self.effect_filter, self.search, self.sort):
            filters.addWidget(widget)
        results_layout.addLayout(filters)
        self.results_table = QTableWidget(); self.results_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.results_table.itemSelectionChanged.connect(self._show_feature)
        results_layout.addWidget(self.results_table, 1)
        self.feature_detail = QPlainTextEdit(); self.feature_detail.setReadOnly(True)
        results_layout.addWidget(self.feature_detail)
        self.result_tabs.addTab(results_page, "Results")
        for signal in (self.result_filter.currentIndexChanged, self.effect_filter.currentIndexChanged,
            self.search.textChanged, self.sort.currentIndexChanged):
            signal.connect(self._filter_results)
        significant = QWidget(); significant_layout = QVBoxLayout(significant)
        significant_layout.addWidget(QLabel("Combined significant"))
        self.combined_table = QTableWidget(); significant_layout.addWidget(self.combined_table)
        significant_layout.addWidget(QLabel("FDR significant"))
        self.fdr_table = QTableWidget(); significant_layout.addWidget(self.fdr_table)
        self.result_tabs.addTab(significant, "Significant")
        qualitative = QWidget(); qualitative_layout = QVBoxLayout(qualitative)
        qualitative_layout.addWidget(self._note("These features were not included in the continuous differential model and do not have limma p-values or FDR values."))
        self.qualitative_table = QTableWidget(); qualitative_layout.addWidget(self.qualitative_table)
        self.result_tabs.addTab(qualitative, "Qualitative candidates")
        audit = QWidget(); audit_layout = QVBoxLayout(audit)
        self.model_table = QTableWidget(); audit_layout.addWidget(self.model_table)
        self.ebayes_table = QTableWidget(); audit_layout.addWidget(self.ebayes_table)
        self.result_tabs.addTab(audit, "Model audit")
        parent = QWidget(); parent_layout = QVBoxLayout(parent)
        self.parent_provenance = self._note(""); parent_layout.addWidget(self.parent_provenance)
        self.parent_audit_tabs = QTabWidget(); self.parent_audit_tables = {}
        for label, key in (("Feature eligibility", "feature_eligibility"),
            ("Missingness patterns", "feature_missingness_patterns"),
            ("Normalization", "sample_normalization"),
            ("Imputation summary", "feature_imputation_summary"), ("Imputed cells", "imputed_cells")):
            table = QTableWidget(); self.parent_audit_tables[key] = table
            self.parent_audit_tabs.addTab(table, label)
        parent_layout.addWidget(self.parent_audit_tabs)
        self.result_tabs.addTab(parent, "Preparation audit")
        graphs = QWidget(); graphs_layout = QVBoxLayout(graphs)
        self.graph_choice = QComboBox(); self.graph_choice.currentIndexChanged.connect(self._show_graph)
        self.graph_label = QLabel("Graph not available for this run.")
        self.graph_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.graph_scroll = QScrollArea(); self.graph_scroll.setWidgetResizable(True); self.graph_scroll.setWidget(self.graph_label)
        graphs_layout.addWidget(self.graph_choice); graphs_layout.addWidget(self.graph_scroll)
        controls = QHBoxLayout()
        for label, callback in (("Previous", lambda: self.graph_choice.setCurrentIndex(max(0, self.graph_choice.currentIndex()-1))),
            ("Next", lambda: self.graph_choice.setCurrentIndex(min(self.graph_choice.count()-1, self.graph_choice.currentIndex()+1))),
            ("Fit", lambda: self._set_fit(True)), ("Actual size", lambda: self._set_fit(False))):
            button = QPushButton(label); button.clicked.connect(callback); controls.addWidget(button)
        graphs_layout.addLayout(controls); self.result_tabs.addTab(graphs, "Graphs")
        self.history_table = QTableWidget(); self.result_tabs.addTab(self.history_table, "History")
        actions = QHBoxLayout()
        for label, method in (("Export current table...", self._export_table),
            ("Export workbook...", self._export_workbook), ("Export graph...", self._export_graph),
            ("Export selected feature details...", self._export_feature),
            ("Open Differential Analysis results folder", self._open_results),
            ("Open Preparation run folder", self._open_preparation)):
            button = QPushButton(label); button.clicked.connect(method); actions.addWidget(button)
        self.explore_button = QPushButton("Explore across analyses...")
        self.explore_button.clicked.connect(self._explore_selected)
        actions.addWidget(self.explore_button)
        layout.addLayout(actions)
        self.stage_tabs.addTab(page, "Results & history")

    def set_project(self, project):
        if self.is_running():
            raise RuntimeError("Wait for the current analysis before changing projects.")
        self.project, self.preparation, self.statistics = project, None, None
        self.transformation.setProperty("user_changed", False)
        self.zero_missing.setProperty("user_changed", False)
        self._check_imputation_packages()
        design = experimental_design(project) if project else {"quantification_columns": []}
        entries = design["quantification_columns"]
        names = list(dict.fromkeys(str(entry["condition"]) for entry in entries if entry["condition"]))
        for choice in (self.condition_a, self.condition_b):
            choice.blockSignals(True); choice.clear(); choice.addItems(names); choice.blockSignals(False)
        if len(names) > 1:
            self.condition_b.setCurrentIndex(1)
        self._conditions_changed()
        self._refresh_preparation_history()
        self._refresh_statistics_history()
        self._update_state()

    def _conditions_changed(self):
        a, b = self.condition_a.currentText(), self.condition_b.currentText()
        self.direction.setText(f"Effect direction: {a or 'Condition A'} − {b or 'Condition B'}")
        entries = experimental_design(self.project)["quantification_columns"] if self.project else []
        relevant = [entry for entry in entries if entry["condition"] in (a, b)]
        self.samples.blockSignals(True); self.samples.setRowCount(len(relevant))
        for row, entry in enumerate(relevant):
            use = QTableWidgetItem(); use.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsUserCheckable)
            use.setCheckState(Qt.CheckState.Checked); use.setData(Qt.ItemDataRole.UserRole, entry["quantification_type"])
            self.samples.setItem(row, 0, use)
            for column, value in enumerate((f"S{row+1:03d}", entry["column"], entry["condition"],
                entry["replicate"], FAMILIES.get(entry["quantification_type"], entry["quantification_type"])), 1):
                item = QTableWidgetItem(str(value or "Missing")); item.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
                self.samples.setItem(row, column, item)
        self.samples.blockSignals(False)
        self._update_state()

    def _check_imputation_packages(self):
        """Disable only methods whose installed R dependencies are unavailable."""
        available = {"none": True, "MinProb": False, "QRILC": False, "KNN": False}
        if self.runtime.executable:
            expression = 'cat(as.integer(vapply(c("MsCoreUtils","imputeLCMD","impute"), requireNamespace, logical(1), quietly=TRUE)), sep=" ")'
            try:
                result = run_external((self.runtime.executable, "-e", expression),
                    capture_output=True, text=True, timeout=15, shell=False)
                flags = [part == "1" for part in result.stdout.strip().split()]
                if result.returncode == 0 and len(flags) == 3:
                    available.update({"MinProb": flags[0] and flags[1],
                        "QRILC": flags[0] and flags[1], "KNN": flags[0] and flags[2]})
            except (OSError, subprocess.TimeoutExpired):
                pass
        for index in range(self.imputation.count()):
            self.imputation.model().item(index).setEnabled(available[self.imputation.itemData(index)])
        if not available.get(self.imputation.currentData(), False):
            self.imputation.setCurrentIndex(0)
    def selected_columns(self):
        return tuple(self.samples.item(row, 2).text() for row in range(self.samples.rowCount())
            if self.samples.item(row, 0).checkState() == Qt.CheckState.Checked)

    def _selected_family(self):
        return {self.samples.item(row, 0).data(Qt.ItemDataRole.UserRole)
            for row in range(self.samples.rowCount())
            if self.samples.item(row, 0).checkState() == Qt.CheckState.Checked}

    def _update_state(self, *_):
        if not hasattr(self, "prepare_button"):
            return
        a, b = self.condition_a.currentText(), self.condition_b.currentText()
        counts = {name: sum(self.samples.item(row, 0).checkState() == Qt.CheckState.Checked and
            self.samples.item(row, 3).text() == name for row in range(self.samples.rowCount())) for name in (a, b)}
        families = self._selected_family()
        family = next(iter(families)) if len(families) == 1 else None
        self.family.setText("Quantification family: " + FAMILIES.get(family, "not available"))
        if not self.project:
            message = "Open a project to prepare differential data."
        elif not a or not b or a == b:
            message = "Condition A and Condition B must differ."
        elif min(counts.values()) < 2:
            message = "Select at least two samples in each condition."
        elif family == "spectral_count":
            message = "Spectral-count differential analysis requires a count-based model and is not supported by the current intensity-based pipeline."
        elif family not in ("lfq_intensity", "raw_intensity", "other_quantitative"):
            message = "Selected samples use incompatible quantification types."
        elif family == "other_quantitative" and not self.other_confirmed.isChecked():
            message = "Confirm that Other quantitative values are continuous."
        else:
            message = "Ready to prepare quantitative data offline."
        self.preparation_readiness.setText(message)
        self.prepare_button.setEnabled(message.startswith("Ready") and not self.is_running())
        if family in ("lfq_intensity", "raw_intensity", "other_quantitative"):
            if not self.transformation.property("user_changed"):
                default = "none" if family == "other_quantitative" else "log2_positive"
                self.transformation.setCurrentIndex(self.transformation.findData(default))
            if not self.zero_missing.property("user_changed"):
                self.zero_missing.setChecked(family != "other_quantitative")
        self.statistics_button.setEnabled(bool(self.project and self.parent_choice.currentData()
            and self.scale.currentData() and not self.is_running()))

    def _imputation_changed(self, *_):
        method = self.imputation.currentData()
        for widget, visible in ((self.seed, method in ("MinProb", "QRILC")),
            (self.q, method == "MinProb"), (self.sigma, method in ("MinProb", "QRILC")),
            (self.margin, method in ("MinProb", "QRILC")), (self.k, method == "KNN")):
            widget.setVisible(visible)

    def preparation_parameters(self):
        return PreparationParameters(self.condition_a.currentText(), self.condition_b.currentText(),
            self.selected_columns(), self.other_confirmed.isChecked(), self.transformation.currentData(),
            self.zero_missing.isChecked(), self.normalization.currentData(), self.minimum_observed.value(),
            self.imputation.currentData(), self.seed.value(), self.q.value(), self.sigma.value(),
            self.margin.currentData(), self.k.value())

    def statistics_parameters(self):
        return DifferentialParameters(self.parent_choice.currentData(), self.scale.currentData(),
            self.trend.isChecked(), self.robust.isChecked(), self.fdr.value(), self.effect.value(),
            self.top_n.value(), self.volcano_labels.value())

    def is_running(self):
        return self._running or bool(self.worker and self.worker.isRunning())

    def _start_preparation(self):
        if not self.project or self.is_running() or not self.prepare_button.isEnabled():
            return
        self._start_worker(run_differential_preparation, self.preparation_parameters(), "preparation")

    def _start_statistics(self):
        if not self.project or self.is_running() or not self.statistics_button.isEnabled():
            return
        self._start_worker(run_differential_analysis, self.statistics_parameters(), "statistics")

    def _start_worker(self, runner, parameters, kind):
        self._running = True; self._update_state(); self.progress.show()
        self.status.setText("Preparing differential-analysis data..." if kind == "preparation" else "Running differential statistics...")
        worker = DifferentialWorker(runner, self.project, self.runtime, new_run_id(), parameters)
        self.worker = worker
        worker.succeeded.connect(lambda output: self._succeeded(kind, output))
        worker.failed.connect(self._failed)
        worker.finished.connect(lambda: self._worker_finished(worker))
        worker.start()

    def _succeeded(self, kind, output):
        self._running = False; self.progress.hide(); self.status.setText("Run completed.")
        if kind == "preparation":
            self._refresh_preparation_history(output["metadata"]["run_id"])
            self._show_preparation(output)
            self.stage_tabs.setCurrentIndex(1)
        else:
            self._refresh_statistics_history(output["metadata"]["run_id"])
            self._show_statistics(output)
            self.stage_tabs.setCurrentIndex(2)
        self._update_state()

    def _failed(self, error):
        self._running = False; self.progress.hide()
        details = getattr(error, "details", None)
        message = str(error) + ("\n" + json.dumps(details, indent=2) if details else "")
        self.status.setText(message)
        self._refresh_preparation_history(); self._refresh_statistics_history()
        self._update_state()
        QMessageBox.warning(self, "Differential Analysis failed", message)

    def _worker_finished(self, worker):
        if self.worker is worker:
            self.worker = None
        worker.deleteLater()
        self._update_state()

    def _refresh_preparation_history(self, selected=None):
        self.prep_history.blockSignals(True); self.prep_history.clear()
        self.parent_choice.blockSignals(True); self.parent_choice.clear()
        if self.project:
            ready = set(list_preparation_runs(self.project))
            for run_id in preparation_runs(self.project):
                root = Path(self.project.root) / "analyses/Differential/Preparation/runs" / run_id
                try:
                    meta = json.loads((root / "input/parameters.json").read_text(encoding="utf-8"))
                    summary = f"{meta.get('created_at','')} | {meta.get('condition_A','')} vs {meta.get('condition_B','')} | {meta.get('transformation','')} | {meta.get('normalization','')} | {meta.get('imputation_method','')} | {run_id}"
                except (OSError, ValueError):
                    summary = run_id
                label = summary if run_id in ready else f"Incomplete preparation run | {summary}"
                self.prep_history.addItem(label, run_id)
                if run_id in ready:
                    self.parent_choice.addItem(summary, run_id)
        self.prep_history.blockSignals(False); self.parent_choice.blockSignals(False)
        if selected:
            index = self.prep_history.findData(selected)
            if index >= 0: self.prep_history.setCurrentIndex(index)
            index = self.parent_choice.findData(selected)
            if index >= 0: self.parent_choice.setCurrentIndex(index)
        self._load_preparation_history(); self._parent_changed()

    def _load_preparation_history(self, *_):
        run_id = self.prep_history.currentData()
        if not self.project or not run_id:
            self.preparation = None; self.prep_summary.clear()
            for table in self.prep_tables.values(): table.clear(); table.setRowCount(0)
            return
        try:
            self._show_preparation(load_preparation(self.project, run_id))
        except Exception as error:
            self.preparation = None; self.prep_summary.setPlainText(f"Incomplete preparation run: {error}")
            for table in self.prep_tables.values(): table.clear(); table.setRowCount(0)

    def _show_preparation(self, output):
        self.preparation = output
        meta = output["metadata"]
        self.prep_summary.setPlainText("\n".join(f"{key.replace('_', ' ').title()}: {value}" for key, value in meta.items()
            if key in ("run_id", "condition_A", "condition_B", "sample_count_A", "sample_count_B",
                "quantification_type", "transformation", "normalization", "zero_is_missing",
                "minimum_observed_per_condition", "imputation_method", "imputation_seed")))
        summary = output["tables"]["summary"]
        if {"metric", "value"}.issubset(summary.columns):
            self.prep_summary.appendPlainText("\n" + "\n".join(f"{row.metric}: {row.value}" for row in summary.itertuples(index=False)))
        for key, table in self.prep_tables.items():
            self._fill(table, output["tables"][key])
        if meta.get("imputation_method") == "none":
            self.prep_audit.setTabText(self.prep_audit.indexOf(self.prep_tables["imputed_cells"]), "Imputation — No imputation was applied")
        else:
            self.prep_audit.setTabText(self.prep_audit.indexOf(self.prep_tables["imputed_cells"]), "Imputation")

    def _parent_changed(self, *_):
        run_id = self.parent_choice.currentData()
        if not self.project or not run_id:
            self.parent_summary.setText("No Ready Preparation run selected.")
            self._update_state(); return
        try:
            output = load_preparation(self.project, run_id)
            meta = output["metadata"]
            self.parent_summary.setText(f"{run_id} | {meta['condition_A']} vs {meta['condition_B']} | "
                f"{meta['transformation']} | {meta['normalization']} | {meta['imputation_method']} | "
                f"{len(output['tables']['prepared_matrix'])} eligible features")
            auto = meta["transformation"] == "log2_positive"
            self.scale.setCurrentIndex(self.scale.findData("log2") if auto else 0)
            self.scale.setEnabled(not auto)
        except Exception as error:
            self.parent_summary.setText(f"Incomplete preparation run: {error}")
            self.scale.setCurrentIndex(0)
        self._update_state()

    def _scale_changed(self, *_):
        self.effect_label.setText("Minimum absolute log2 fold change" if self.scale.currentData() == "log2"
            else "Minimum absolute difference" if self.scale.currentData() == "continuous" else "Minimum absolute effect")
        self._update_state()

    def _refresh_statistics_history(self, selected=None):
        self.stat_history.blockSignals(True); self.stat_history.clear()
        records = []
        if self.project:
            for run_id in statistics_runs(self.project):
                root = Path(self.project.root) / "analyses/Differential/Statistics/runs" / run_id
                try:
                    meta = json.loads((root / "input/parameters.json").read_text(encoding="utf-8"))
                    label = f"{meta.get('created_at','')} | {meta.get('condition_A','')} vs {meta.get('condition_B','')} | {meta.get('parent_preparation_run_id','')} | {meta.get('prepared_scale','')} | {run_id}"
                    load_statistics(self.project, run_id)
                except Exception:
                    label = f"Incomplete run | {run_id}"
                self.stat_history.addItem(label, run_id)
                records.append({"Run ID": run_id, "Status": "Incomplete run" if label.startswith("Incomplete") else "Ready", "Description": label})
        self.stat_history.blockSignals(False)
        self._fill(self.history_table, pd.DataFrame(records, columns=["Run ID", "Status", "Description"]))
        if selected:
            index = self.stat_history.findData(selected)
            if index >= 0: self.stat_history.setCurrentIndex(index)
        self._load_statistics_history()

    def _load_statistics_history(self, *_):
        run_id = self.stat_history.currentData()
        if not self.project or not run_id:
            self._clear_statistics(); return
        try:
            self._show_statistics(load_statistics(self.project, run_id))
        except Exception as error:
            self._clear_statistics()
            self.status.setText(f"Incomplete run: {error}")

    def _clear_statistics(self):
        self.statistics = None
        self._frames = {}
        for table in (self.summary_table, self.results_table, self.combined_table,
            self.fdr_table, self.qualitative_table, self.model_table, self.ebayes_table,
            *self.parent_audit_tables.values()):
            table.clear(); table.setRowCount(0)
        self.feature_detail.clear()
        self.parent_provenance.clear()
        self.graph_choice.clear()
        self.graph_label.setPixmap(QPixmap())
        self.graph_label.setText("Graph not available for this run.")
    def _show_statistics(self, output):
        self.statistics = output
        meta, tables = output["metadata"], output["tables"]
        self._frames = tables
        self._fill(self.summary_table, tables["summary"])
        self._fill(self.combined_table, tables["significant_results"], RESULT_FIELDS)
        self._fill(self.fdr_table, tables["fdr_significant_results"], RESULT_FIELDS)
        self._fill(self.qualitative_table, tables["qualitative_candidates"])
        self._fill(self.model_table, tables["model_audit"])
        self._fill(self.ebayes_table, tables["ebayes_audit"])
        self._filter_results()
        a, b = meta.get("condition_A", "Condition A"), meta.get("condition_B", "Condition B")
        self.result_tabs.setTabText(1, f"Results — positive: higher in {a}; negative: higher in {b}")
        self.parent_provenance.setText(" | ".join(f"{key.replace('_',' ').title()}: {meta.get(key,'NA')}" for key in
            ("parent_preparation_run_id", "preparation_transformation", "preparation_normalization",
                "preparation_imputation", "preparation_imputation_seed", "preparation_minimum_observed")))
        parent_id = meta.get("parent_preparation_run_id")
        try:
            parent = load_preparation(self.project, parent_id)
            for key, table in self.parent_audit_tables.items():
                self._fill(table, parent["tables"][key])
        except Exception:
            self.parent_provenance.setText(self.parent_provenance.text() + " | Parent Preparation run unavailable; frozen Statistics results remain available.")
            copied = {"feature_eligibility": "feature_eligibility.csv",
                "feature_imputation_summary": "feature_imputation_summary.csv"}
            for key, table in self.parent_audit_tables.items():
                path = output["run_root"] / "input" / copied.get(key, "")
                if key in copied and path.is_file():
                    self._fill(table, pd.read_csv(path, dtype=str, keep_default_na=False))
                else:
                    table.clear(); table.setRowCount(0)
        self.graph_choice.blockSignals(True); self.graph_choice.clear()
        available = {graph.stem: graph for graph in output["plots"]}
        for stem in ("volcano", "ma", "mean_a_vs_b", "p_value_histogram", "fdr_histogram",
            "effect_distribution", "residual_df_distribution", "top_effects", "imputation_diagnostic"):
            if stem == "imputation_diagnostic" and meta.get("preparation_imputation") == "none":
                continue
            self.graph_choice.addItem(stem.replace("_", " ").title(), available.get(stem))
        self.graph_choice.blockSignals(False); self._show_graph()
        self.status.setText("No features met the selected significance criteria." if tables["significant_results"].empty else "Historical results loaded.")

    def _filter_results(self, *_):
        if not self.statistics: return
        frame = self.statistics["tables"]["all_results"]
        choice = self.result_filter.currentText()
        column = {"Tested": ("tested", "TRUE"), "FDR significant": ("fdr_significant", "TRUE"),
            "Combined significant": ("combined_significant", "TRUE"), "Higher in Condition A": ("effect_direction", "higher_in_condition_A"),
            "Higher in Condition B": ("effect_direction", "higher_in_condition_B"),
            "Not significant": ("result_class", "Not significant"), "Not tested": ("tested", "FALSE")}.get(choice)
        if column: frame = frame[frame[column[0]] == column[1]]
        effect = self.effect_filter.currentIndex()
        if effect: frame = frame[frame["effect_size_pass"] == ("TRUE" if effect == 1 else "FALSE")]
        query = self.search.text().strip()
        if query:
            fields = [name for name in ("feature_id", "display_identifier", "gene_symbol", "uniprot_accession") if name in frame]
            mask = pd.Series(False, index=frame.index)
            for name in fields: mask |= frame[name].str.contains(query, case=False, regex=False, na=False)
            frame = frame[mask]
        sort = self.sort.currentData()
        if sort:
            frame = frame.copy()
            frame["_sort"] = pd.to_numeric(frame["effect" if sort == "abs_effect" else sort], errors="coerce")
            if sort == "abs_effect": frame["_sort"] = frame["_sort"].abs()
            frame = frame.sort_values("_sort", ascending=sort in ("adj.P.Val", "P.Value"), na_position="last", kind="stable")
        self._visible_results = frame.drop(columns=["_sort"], errors="ignore")
        self._fill(self.results_table, self._visible_results, RESULT_FIELDS)
        self._current_table_key = "all_results"

    def _show_feature(self):
        if not self.statistics or not self.results_table.selectedItems(): return
        item = self.results_table.selectedItems()[0]
        feature_id = item.data(Qt.ItemDataRole.UserRole)
        frame = self.statistics["tables"]["all_results"]
        matched = frame[frame.feature_id == feature_id]
        if matched.empty: return
        row = matched.iloc[0]
        lines = [f"{key.replace('_', ' ').title()}: {value or 'NA'}" for key, value in row.items()]
        if row.get("imputed_cell_count", "0") not in ("0", "", "NA"):
            lines += ["This feature contains imputed quantitative values from the selected Preparation run.",
                "The differential model does not propagate imputation uncertainty."]
        self.feature_detail.setPlainText("\n".join(lines))

    def _set_fit(self, fit):
        self._fit = fit; self._show_graph()

    def _show_graph(self, *_):
        path = self.graph_choice.currentData()
        if not path or not Path(path).is_file():
            self.graph_label.setPixmap(QPixmap()); self.graph_label.setText("Graph not available for this run."); return
        pixmap = QPixmap(str(path))
        if self._fit and not pixmap.isNull():
            pixmap = pixmap.scaled(self.graph_scroll.viewport().size(),
                Qt.AspectRatioMode.KeepAspectRatio, Qt.TransformationMode.SmoothTransformation)
        self.graph_label.setPixmap(pixmap)

    def _current_frame(self):
        index = self.result_tabs.currentIndex()
        return {0: self._frames.get("summary"), 1: getattr(self, "_visible_results", self._frames.get("all_results")),
            2: self._frames.get("significant_results"), 3: self._frames.get("qualitative_candidates"),
            4: self._frames.get("model_audit")}.get(index)

    def _save_path(self, title, name, filter):
        return QFileDialog.getSaveFileName(self, title, name, filter)[0]

    def _export_table(self):
        frame = self._current_frame()
        if frame is None: return
        path = self._save_path("Export current table", "differential_table.csv", "CSV files (*.csv)")
        if path: frame.to_csv(path, index=False)

    def _export_workbook(self):
        if not self.statistics: return
        source = self.statistics["workbook"]
        path = self._save_path("Export workbook", source.name, "Excel workbook (*.xlsx)")
        if path: shutil.copy2(source, path)

    def _export_graph(self):
        source = self.graph_choice.currentData()
        if not source: return
        source = Path(source)
        path = self._save_path("Export graph", source.name, "Images and PDF (*.png *.pdf)")
        if path:
            wanted = source.with_suffix(Path(path).suffix.lower())
            if wanted.is_file(): shutil.copy2(wanted, path)

    def _export_feature(self):
        if not self.statistics or not self.results_table.selectedItems(): return
        feature_id = self.results_table.selectedItems()[0].data(Qt.ItemDataRole.UserRole)
        path = self._save_path("Export selected feature details", f"{feature_id}.csv", "CSV files (*.csv)")
        if not path: return
        frame = self.statistics["tables"]["all_results"]
        selected = frame[frame.feature_id == feature_id]
        audit = self.statistics["tables"]["model_audit"]
        audit = audit[audit.feature_id == feature_id]
        for name, value in audit.iloc[0].items() if not audit.empty else ():
            selected = selected.assign(**{f"audit_{name}": value})
        selected.to_csv(path, index=False)

    def _explore_selected(self):
        if not self.statistics or not self.results_table.selectedItems():
            return
        feature_id = self.results_table.selectedItems()[0].data(Qt.ItemDataRole.UserRole)
        self.explore_requested.emit("differential", self.statistics["metadata"]["run_id"], feature_id)
    def _open_results(self):
        if self.statistics: QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.statistics["run_root"])))

    def _open_preparation(self):
        if self.preparation: QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.preparation["run_root"])))
