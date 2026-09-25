"""Experiment Comparison setup, frozen results, exports and history."""
from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd
from PySide6.QtCore import Qt, QSortFilterProxyModel, Signal
from PySide6.QtGui import QPixmap, QStandardItem, QStandardItemModel
from PySide6.QtWidgets import (QCheckBox, QComboBox, QFileDialog, QFormLayout, QHBoxLayout,
    QInputDialog, QLabel, QLineEdit, QMessageBox, QPushButton, QTabWidget, QTableView,
    QTextEdit, QVBoxLayout, QWidget)

from ..core import experiment_comparison as comparison
from ..core.differential_analysis import list_runs as list_differential_runs
from ..core.importer import read_table, xlsx_sheets


class FilterProxy(QSortFilterProxyModel):
    def filterAcceptsRow(self, row, parent):
        needle = self.filterRegularExpression().pattern().casefold()
        if not needle:
            return True
        model = self.sourceModel()
        return any(needle in str(model.index(row, column, parent).data() or "").casefold()
            for column in range(model.columnCount()))


class ResultTable(QWidget):
    def __init__(self, name: str):
        super().__init__()
        self.name = name
        self.frame = pd.DataFrame()
        self.search = QLineEdit()
        self.search.setPlaceholderText("Search rows (case-insensitive)")
        self.model = QStandardItemModel()
        self.proxy = FilterProxy()
        self.proxy.setSourceModel(self.model)
        self.table = QTableView()
        self.table.setModel(self.proxy)
        self.table.setSortingEnabled(True)
        self.table.setSelectionBehavior(QTableView.SelectRows)
        self.detail = QTextEdit()
        self.detail.setReadOnly(True)
        self.csv = QPushButton("Export CSV...")
        self.xlsx = QPushButton("Export XLSX...")
        bar = QHBoxLayout()
        bar.addWidget(self.csv)
        bar.addWidget(self.xlsx)
        layout = QVBoxLayout(self)
        layout.addWidget(self.search)
        layout.addWidget(self.table, 1)
        layout.addWidget(QLabel("Selected-row detail"))
        layout.addWidget(self.detail)
        layout.addLayout(bar)
        self.search.textChanged.connect(self.proxy.setFilterFixedString)
        self.table.selectionModel().selectionChanged.connect(self._detail)
        self.csv.clicked.connect(lambda: self._export("csv"))
        self.xlsx.clicked.connect(lambda: self._export("xlsx"))

    def set_frame(self, frame: pd.DataFrame):
        self.frame = frame.copy()
        self.model.clear()
        self.model.setHorizontalHeaderLabels([str(column) for column in frame.columns])
        for _, row in frame.iterrows():
            self.model.appendRow([QStandardItem("" if pd.isna(value) else str(value))
                for value in row])
        self.detail.clear()

    def selected_record(self) -> dict | None:
        selected = self.table.selectionModel().selectedRows()
        if not selected:
            return None
        source = self.proxy.mapToSource(selected[0]).row()
        return self.frame.iloc[source].to_dict() if source < len(self.frame) else None

    def _detail(self):
        record = self.selected_record()
        self.detail.setPlainText("\n".join(f"{key}: {value}" for key, value in record.items())
            if record else "")

    def _export(self, extension: str):
        path, _ = QFileDialog.getSaveFileName(self, f"Export {self.name}",
            f"{self.name.replace(' ', '_').lower()}.{extension}",
            "CSV (*.csv)" if extension == "csv" else "Excel workbook (*.xlsx)")
        if not path:
            return
        if extension == "csv":
            self.frame.to_csv(path, index=False)
        else:
            self.frame.to_excel(path, index=False)


class SourceBox(QWidget):
    config_changed = Signal()

    def __init__(self, side: str):
        super().__init__()
        self.side = side
        self.project = None
        self.path: Path | None = None
        self.sheet: str | None = None
        self.source_kind = "external"
        self.preview = pd.DataFrame()
        self.alias = QLineEdit(f"Dataset {side}")
        self.location = QLabel("No table selected.")
        self.location.setWordWrap(True)
        self.browse = QPushButton("Select file...")
        self.current = QPushButton("Use current project dataset")
        self.frozen = QPushButton("Select frozen run table...")
        self.column = QComboBox()
        self.kind = QComboBox()
        for label, value in (("UniProt", "uniprot"), ("Gene Symbol", "gene_symbol"),
                ("NCBI Gene / Entrez", "entrez"), ("Ensembl Gene", "ensembl_gene"),
                ("Ensembl Protein", "ensembl_protein"), ("RefSeq", "refseq_protein"),
                ("Original ID", "original"), ("Unknown", "unknown")):
            self.kind.addItem(label, value)
        self.inference = QLabel("Suggestions require review and confirmation.")
        self.confirmed = QCheckBox("I confirm this identifier column and type")
        form = QFormLayout(self)
        form.addRow("Alias", self.alias)
        form.addRow("Source", self.location)
        buttons = QHBoxLayout()
        buttons.addWidget(self.browse)
        buttons.addWidget(self.current)
        buttons.addWidget(self.frozen)
        form.addRow(buttons)
        form.addRow("Identifier column", self.column)
        form.addRow("Identifier type", self.kind)
        form.addRow("Suggested configuration", self.inference)
        form.addRow(self.confirmed)
        self.browse.clicked.connect(self._browse)
        self.current.clicked.connect(self._current)
        self.frozen.clicked.connect(self._frozen)
        self.column.currentIndexChanged.connect(lambda: self.confirmed.setChecked(False))
        self.kind.currentIndexChanged.connect(lambda: self.confirmed.setChecked(False))
        self.kind.currentIndexChanged.connect(lambda _index: self.config_changed.emit())

    def set_project(self, project):
        self.project = project
        self.current.setEnabled(bool(project and project.config.get("input", {}).get("processed_file")))
        self.frozen.setEnabled(bool(project))
        self.path = None
        self.location.setText("No table selected.")
        self.confirmed.setChecked(False)

    def _browse(self):
        path, _ = QFileDialog.getOpenFileName(self, f"Select dataset {self.side}", "",
            "Tables (*.csv *.tsv *.xlsx)")
        if path:
            self.set_source(Path(path))

    def _current(self):
        if self.project:
            relative = self.project.config.get("input", {}).get("processed_file")
            if relative:
                self.set_source(Path(self.project.root) / relative, source_kind="project")

    def _frozen(self):
        if not self.project:
            return
        choices = {}
        for run_id in comparison.list_runs(self.project):
            try:
                loaded = comparison.load_run(self.project, run_id)
                for side in ("a", "b"):
                    path = loaded["run_root"] / f"inputs/dataset_{side}" / Path(loaded["config"][side]["path"]).name
                    choices[f"Comparison {run_id} / dataset {side.upper()}"] = (path, f"comparison:{run_id}:{side}")
            except Exception:
                continue
        for run_id in list_differential_runs(self.project):
            path = Path(self.project.root) / "analyses/Differential/Statistics/runs" / run_id / "results/differential_results_all.csv"
            if path.is_file():
                choices[f"Differential {run_id} / all results"] = (path, f"differential:{run_id}")
        if not choices:
            QMessageBox.information(self, "No frozen tables", "No suitable frozen run table is available.")
            return
        label, accepted = QInputDialog.getItem(self, "Select frozen run table",
            f"Dataset {self.side}", list(choices), 0, False)
        if accepted:
            path, kind = choices[label]
            self.set_source(path, source_kind=kind)
    def set_source(self, path: Path, *, source_kind="external", sheet: str | None = None):
        selected = sheet
        if path.suffix.lower() == ".xlsx" and selected is None:
            sheets = xlsx_sheets(path)
            if len(sheets) > 1:
                selected, ok = QInputDialog.getItem(self, "Choose worksheet",
                    f"Dataset {self.side} worksheet", sheets, 0, False)
                if not ok:
                    return
        frame, selected = comparison._read_source(path, selected)
        self.path, self.sheet, self.source_kind, self.preview = path, selected, source_kind, frame
        self.location.setText(str(path) + (f" | sheet: {selected}" if selected else ""))
        self.column.clear()
        self.column.addItems([str(column) for column in frame.columns])
        guesses = comparison.suggestions(frame)
        preferred = next((name for name, item in guesses.items() if item.get("primary_identifier")), None)
        if preferred:
            self.column.setCurrentText(preferred)
            kind = guesses[preferred].get("identifier_type")
            index = self.kind.findData(kind)
            if index >= 0:
                self.kind.setCurrentIndex(index)
        quantitative = [name for name, item in guesses.items() if item.get("role") == "quantification"]
        self.inference.setText(f"Suggested: {preferred or 'none'} / "
            f"{guesses.get(preferred, {}).get('identifier_type', 'unknown')} | "
            f"Quantification: {', '.join(quantitative) or 'none'}")
        self.confirmed.setChecked(False)
        self.config_changed.emit()

    def spec(self) -> comparison.TableSpec:
        if self.path is None or not self.confirmed.isChecked():
            raise comparison.ComparisonError(f"Confirm dataset {self.side} source, column and type.")
        return comparison.TableSpec(str(self.path), self.column.currentText(),
            self.kind.currentData(), self.sheet, self.alias.text().strip() or f"Dataset {self.side}",
            self.source_kind)


class ExperimentComparisonPage(QWidget):
    biological_context_requested = Signal(str, str, str, object)

    def __init__(self):
        super().__init__()
        self.project = None
        self.loaded = None
        self.source_a, self.source_b = SourceBox("A"), SourceBox("B")
        self.comparison_type = QComboBox()
        for label, value in (("Original identifier", "original"), ("UniProt", "uniprot"),
                ("Gene Symbol", "gene_symbol"), ("NCBI Gene", "ncbi_gene"),
                ("Ensembl Gene", "ensembl_gene"), ("Ensembl Protein", "ensembl_protein"),
                ("RefSeq", "refseq")):
            self.comparison_type.addItem(label, value)
        self.mapping_path: Path | None = None
        self.mapping_label = QLabel("No mapping selected. Different ID types will require mapping.")
        self.mapping_button = QPushButton("Select persisted Mapping catalog...")
        self.differential_a, self.differential_b = QComboBox(), QComboBox()
        self.comparable = QCheckBox("These contrasts represent comparable biological directions")
        self.match_a = QCheckBox("Differential run A corresponds to dataset A")
        self.match_b = QCheckBox("Differential run B corresponds to dataset B")
        self.contrasts = QLabel("Select two persisted Differential runs to review their contrasts.")
        self.contrasts.setWordWrap(True)
        self.compare_button = QPushButton("Compare")
        self.status = QLabel("Select two independent datasets.")
        self.status.setWordWrap(True)
        setup = QWidget()
        form = QFormLayout(setup)
        form.addRow("Dataset A", self.source_a)
        form.addRow("Dataset B", self.source_b)
        form.addRow("Comparison identity", self.comparison_type)
        form.addRow("Mapping", self.mapping_label)
        form.addRow(self.mapping_button)
        form.addRow("Differential run A", self.differential_a)
        form.addRow("Differential run B", self.differential_b)
        form.addRow("Contrasts", self.contrasts)
        form.addRow(self.match_a)
        form.addRow(self.match_b)
        form.addRow(self.comparable)
        form.addRow(self.compare_button)
        form.addRow(self.status)
        self.tabs = QTabWidget()
        self.tabs.addTab(setup, "Setup")
        self.summary = QLabel("No comparison loaded.")
        self.summary.setWordWrap(True)
        self.tabs.addTab(self.summary, "Summary")
        self.views = {}
        for title, key in (("A only", "a_only"), ("B only", "b_only"),
                ("Shared", "shared"), ("Ambiguous / Unmapped", "ambiguous"),
                ("Differential behavior", "differential_comparison")):
            view = ResultTable(title)
            self.views[key] = view
            self.tabs.addTab(view, title)
        self.graph = QLabel("No graph available.")
        self.graph_choice = QComboBox()
        graph_page = QWidget()
        graph_layout = QVBoxLayout(graph_page)
        graph_layout.addWidget(self.graph_choice)
        graph_layout.addWidget(self.graph, 1)
        self.tabs.addTab(graph_page, "Graphs")
        self.graph_choice.currentIndexChanged.connect(self._show_graph)
        exports = QWidget()
        export_layout = QVBoxLayout(exports)
        self.export_workbook = QPushButton("Export consolidated workbook...")
        export_layout.addWidget(self.export_workbook)
        export_layout.addStretch()
        self.tabs.addTab(exports, "Exports")
        history = QWidget()
        history_layout = QVBoxLayout(history)
        self.history = QComboBox()
        self.load_history = QPushButton("Load selected comparison")
        history_layout.addWidget(self.history)
        history_layout.addWidget(self.load_history)
        history_layout.addStretch()
        self.tabs.addTab(history, "History")
        layout = QVBoxLayout(self)
        layout.addWidget(self.tabs)
        self.context_button = QPushButton("Biological context...")
        self.context_button.setEnabled(False)
        layout.addWidget(self.context_button)
        self.context_button.clicked.connect(self._open_context)
        self.mapping_button.clicked.connect(self._select_mapping)
        self.source_a.config_changed.connect(self._refresh_identity_options)
        self.source_b.config_changed.connect(self._refresh_identity_options)
        self.compare_button.clicked.connect(self._compare)
        self.export_workbook.clicked.connect(self._export_workbook)
        self.load_history.clicked.connect(self._load_history)
        self.differential_a.currentIndexChanged.connect(self._contrast_text)
        self.differential_b.currentIndexChanged.connect(self._contrast_text)
        self.comparable.setChecked(False)

    def set_project(self, project):
        self.project = project
        self.loaded = None
        self.mapping_path = None
        self.mapping_label.setText("No mapping selected. Different ID types will require mapping.")
        self.source_a.set_project(project)
        self.source_b.set_project(project)
        for combo in (self.differential_a, self.differential_b):
            combo.clear()
            combo.addItem("None", None)
            if project:
                for run_id in list_differential_runs(project):
                    combo.addItem(run_id, run_id)
        self.comparable.setChecked(False)
        self.match_a.setChecked(False)
        self.match_b.setChecked(False)
        self.history.clear()
        if project:
            for run_id in comparison.list_runs(project):
                self.history.addItem(run_id, run_id)
        self._refresh_identity_options()
        self.status.setText("Select two independent datasets." if project else "Open a project.")

    def _refresh_identity_options(self):
        choices = (("Original identifier", "original"), ("UniProt", "uniprot"),
            ("Gene Symbol", "gene_symbol"), ("NCBI Gene", "ncbi_gene"),
            ("Ensembl Gene", "ensembl_gene"), ("Ensembl Protein", "ensembl_protein"),
            ("RefSeq", "refseq"))
        left = comparison.canonical_type(self.source_a.kind.currentData() or "")
        right = comparison.canonical_type(self.source_b.kind.currentData() or "")
        columns = set()
        if self.mapping_path:
            try:
                columns = set(pd.read_csv(self.mapping_path, nrows=0).columns)
            except (OSError, ValueError):
                columns = set()
        def resolvable(source, target):
            return source == target or bool(columns and
                comparison.IDENTITY_COLUMNS.get(source) in columns and
                comparison.IDENTITY_COLUMNS.get(target) in columns)
        previous = self.comparison_type.currentData()
        self.comparison_type.clear()
        if self.source_a.path is not None and self.source_b.path is not None:
            for label, target in choices:
                if target == "original":
                    allowed = left == right
                else:
                    allowed = resolvable(left, target) and resolvable(right, target)
                if allowed:
                    self.comparison_type.addItem(label, target)
        index = self.comparison_type.findData(previous)
        if index >= 0:
            self.comparison_type.setCurrentIndex(index)
    def _select_mapping(self):
        if not self.project:
            return
        default = str(Path(self.project.root) / "mapping/tables/protein_catalog.csv")
        path, _ = QFileDialog.getOpenFileName(self, "Select persisted Mapping catalog",
            default, "CSV (*.csv)")
        if path:
            self.mapping_path = Path(path)
            self.mapping_label.setText(path)
            self._refresh_identity_options()

    def _contrast_text(self):
        if not self.project:
            return
        from ..core.differential_analysis import load_run
        lines = []
        for side, combo in (("A", self.differential_a), ("B", self.differential_b)):
            if combo.currentData():
                try:
                    meta = load_run(self.project, combo.currentData())["metadata"]
                    lines.append(f"Run {side}: {meta.get('condition_A', '?')} − "
                        f"{meta.get('condition_B', '?')} | {meta.get('comparison_direction', 'Condition A - Condition B')}")
                except Exception as error:
                    lines.append(f"Run {side}: unavailable ({error})")
        self.contrasts.setText("\n".join(lines) if lines else "No Differential runs selected.")
        self.comparable.setChecked(False)
        self.match_a.setChecked(False)
        self.match_b.setChecked(False)

    def _compare(self):
        if not self.project:
            return
        try:
            spec = comparison.ComparisonSpec(self.source_a.spec(), self.source_b.spec(),
                self.comparison_type.currentData(), str(self.mapping_path) if self.mapping_path else None,
                self.differential_a.currentData(), self.differential_b.currentData(),
                self.comparable.isChecked(), self.match_a.isChecked(), self.match_b.isChecked())
            root = comparison.compare(self.project, spec)
            self._display(comparison.load_run(self.project, root.name))
            self.history.insertItem(0, root.name, root.name)
            self.status.setText(f"Comparison ready: {root.name}")
        except Exception as error:
            self.status.setText(str(error))
            QMessageBox.warning(self, "Comparison unavailable", str(error))

    def _display(self, loaded):
        self.loaded = loaded
        summary = loaded["manifest"]["summary"]
        aliases = loaded["config"]
        lines = [f"{aliases['a']['alias']} vs {aliases['b']['alias']}",
            *(f"{name}: {value}" for name, value in summary.items())]
        differential = loaded["manifest"].get("differential") or {}
        if differential:
            lines.append(f"Comparable contrasts confirmed: {bool(aliases.get('comparable_contrasts'))}")
            for side in ("A", "B"):
                metadata = differential.get(side) or {}
                lines.append(f"Differential {side}: {aliases.get('differential_' + side.lower())} | "
                    f"{metadata.get('condition_A', '?')} - {metadata.get('condition_B', '?')} | "
                    f"{metadata.get('comparison_direction', 'Condition A - Condition B')}")
        self.summary.setText("\n".join(lines))
        for key, view in self.views.items():
            if key == "ambiguous":
                parts = []
                for name, side in (("ambiguous", "Ambiguous"), ("unmapped_a", "Unmapped A"),
                                   ("unmapped_b", "Unmapped B")):
                    frame = loaded["tables"].get(name, pd.DataFrame()).copy()
                    if not frame.empty:
                        frame.insert(0, "category", side)
                        parts.append(frame)
                view.set_frame(pd.concat(parts, ignore_index=True) if parts else pd.DataFrame())
            else:
                view.set_frame(loaded["tables"].get(key, pd.DataFrame()))
        self.graph_choice.clear()
        for name, filename in (("Membership counts", "membership.png"),
                ("Differential patterns", "differential_patterns.png")):
            path = loaded["run_root"] / "graphs" / filename
            if path.is_file():
                self.graph_choice.addItem(name, str(path))
        self._show_graph()
        self.context_button.setEnabled(any(str(loaded["config"][side].get("source_kind", "")).startswith("differential:") or (loaded["config"].get(f"differential_{side}") and loaded["config"].get(f"differential_match_{side}")) for side in ("a", "b")))
        self.tabs.setCurrentIndex(1)

    def _show_graph(self):
        path = self.graph_choice.currentData()
        self.graph.setPixmap(QPixmap(path) if path else QPixmap())
    def _open_context(self):
        if not self.loaded:
            return
        widget = self.tabs.currentWidget()
        if not isinstance(widget, ResultTable):
            QMessageBox.information(self, "Choose a protein", "Select a comparison table row first.")
            return
        record = widget.selected_record()
        if not record or not record.get("comparison_entity"):
            QMessageBox.information(self, "Choose a protein", "Select a comparison entity first.")
            return
        options = {}
        for side in ("a", "b"):
            source = self.loaded["config"][side]
            kind = str(source.get("source_kind", ""))
            present = str(record.get(f"{side.upper()}_present", "True")).lower() != "false"
            if not present:
                continue
            run_id = kind.split(":", 1)[1] if kind.startswith("differential:") else (
                self.loaded["config"].get(f"differential_{side}")
                if self.loaded["config"].get(f"differential_match_{side}") else None)
            if run_id:
                options[f"Dataset {side.upper()} - {source.get('alias', '')}"] = run_id
        if not options:
            QMessageBox.information(self, "Biological context unavailable",
                "The selected entity has no safely identified Differential source side.")
            return
        label, accepted = QInputDialog.getItem(self, "Choose context side",
            "Open Biological Context from", list(options), 0, False)
        if accepted:
            self.biological_context_requested.emit("differential", options[label],
                str(record["comparison_entity"]), None)
    def _load_history(self):
        if not self.project or not self.history.currentData():
            return
        try:
            self._display(comparison.load_run(self.project, self.history.currentData()))
            self.status.setText(f"Historical comparison loaded: {self.history.currentData()}")
        except Exception as error:
            QMessageBox.warning(self, "Historical comparison unavailable", str(error))

    def _export_workbook(self):
        if not self.loaded:
            return
        source = self.loaded["run_root"] / "experiment_comparison.xlsx"
        target, _ = QFileDialog.getSaveFileName(self, "Export consolidated workbook",
            source.name, "Excel workbook (*.xlsx)")
        if target:
            shutil.copy2(source, target)
