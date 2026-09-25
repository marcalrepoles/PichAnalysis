"""Functional comparison embedded in Experiment Comparison."""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import pandas as pd
from PySide6.QtCore import Signal
from PySide6.QtGui import QPixmap
from PySide6.QtWidgets import (QComboBox, QFileDialog, QHBoxLayout, QLabel,
    QMessageBox, QPushButton, QTabWidget, QVBoxLayout, QWidget)

from ..core import functional_comparison as functional


LABELS = {"go": "GO", "kegg": "KEGG", "reactome": "Reactome",
    "mitocarta": "MitoCarta", "mtdna": "mtDNA Evidence", "domains": "InterPro / Pfam",
    "complexes": "Complexes", "string": "STRING"}


class ModulePanel(QWidget):
    context_requested = Signal(str, str, str, object)
    changed = Signal()

    def __init__(self, module: str):
        super().__init__()
        # Import at construction time to reuse 20B's sortable/searchable table
        # without introducing an import cycle with its containing page.
        from .experiment_comparison_page import ResultTable
        self.module = module
        self.project = None
        self.comparison_run_id = None
        self.loaded = None
        self.run_a, self.run_b = QComboBox(), QComboBox()
        self.compare_button = QPushButton("Compare selected runs")
        self.history = QComboBox()
        self.load_button = QPushButton("Load functional comparison")
        self.warning = QLabel(f"No {LABELS[module]} comparison configured.")
        self.warning.setWordWrap(True)
        self.unit = QComboBox()
        self.views = {}
        result_tabs = QTabWidget()
        self.result_tabs = result_tabs
        for title, key in (("Master", "master"), ("A only", "a_only"),
                ("B only", "b_only"), ("Shared", "shared"), ("Enrichment", "enrichment"), ("Members", "members")):
            view = ResultTable(f"{LABELS[module]} {title}")
            result_tabs.addTab(view, title)
            self.views[key] = view
        self.graph = QLabel()
        self.export_workbook = QPushButton("Export functional workbook...")
        self.context_a = QPushButton("Open A biological context")
        self.context_b = QPushButton("Open B biological context")
        self.context_a.setEnabled(False)
        self.context_b.setEnabled(False)
        controls = QHBoxLayout()
        for label, combo in (("Run A", self.run_a), ("Run B", self.run_b)):
            controls.addWidget(QLabel(label))
            controls.addWidget(combo)
        controls.addWidget(self.compare_button)
        history_bar = QHBoxLayout()
        history_bar.addWidget(QLabel("History"))
        history_bar.addWidget(self.history)
        history_bar.addWidget(self.load_button)
        unit_bar = QHBoxLayout()
        unit_bar.addWidget(QLabel("Comparable unit"))
        unit_bar.addWidget(self.unit)
        unit_bar.addWidget(self.context_a)
        unit_bar.addWidget(self.context_b)
        layout = QVBoxLayout(self)
        layout.addLayout(controls)
        layout.addLayout(history_bar)
        layout.addWidget(self.warning)
        layout.addLayout(unit_bar)
        layout.addWidget(result_tabs, 1)
        layout.addWidget(self.graph)
        layout.addWidget(self.export_workbook)
        self.compare_button.clicked.connect(self._compare)
        self.load_button.clicked.connect(self._load_selected)
        self.unit.currentIndexChanged.connect(self._show_unit)
        self.export_workbook.clicked.connect(self._export_workbook)
        self.context_a.clicked.connect(lambda: self._open_context("A"))
        self.context_b.clicked.connect(lambda: self._open_context("B"))

    def set_comparison(self, project, comparison_run_id: str | None):
        self.project = project
        self.comparison_run_id = comparison_run_id
        self.loaded = None
        self.unit.clear()
        self.run_a.clear()
        self.run_b.clear()
        self.history.clear()
        self.run_a.addItem("Select run A", None)
        self.run_b.addItem("Select run B", None)
        if project and comparison_run_id:
            for run_id in functional.list_source_runs(project, self.module):
                self.run_a.addItem(run_id, run_id)
                self.run_b.addItem(run_id, run_id)
            for run_id in functional.list_functional_runs(project, comparison_run_id, self.module):
                self.history.addItem(run_id, run_id)
            selected_path = Path(project.root) / "analyses/experiment_comparison/runs" / comparison_run_id / \
                "functional/selected.json"
            if selected_path.is_file():
                try:
                    selected = json.loads(selected_path.read_text(encoding="utf-8")).get(self.module)
                    if selected:
                        self._display(functional.load_functional(project, comparison_run_id,
                            self.module, selected))
                        index = self.history.findData(selected)
                        if index >= 0:
                            self.history.setCurrentIndex(index)
                except (OSError, ValueError, KeyError, functional.FunctionalComparisonError):
                    self.warning.setText("Saved functional selection is incomplete or unavailable.")
        if self.loaded is None:
            self.warning.setText(f"No {LABELS[self.module]} comparison configured.")
            for view in self.views.values():
                view.set_frame(pd.DataFrame())
        self.compare_button.setEnabled(bool(project and comparison_run_id))
        self.context_a.setEnabled(bool(self.loaded and self.module in {"kegg", "reactome", "string"}))
        self.context_b.setEnabled(bool(self.loaded and self.module in {"kegg", "reactome", "string"}))

    def _remember(self, run_id: str):
        root = Path(self.project.root) / "analyses/experiment_comparison/runs" / \
            self.comparison_run_id / "functional"
        root.mkdir(parents=True, exist_ok=True)
        path = root / "selected.json"
        selected = json.loads(path.read_text(encoding="utf-8")) if path.is_file() else {}
        selected[self.module] = run_id
        temporary = path.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(selected, indent=2) + "\n", encoding="utf-8")
        temporary.replace(path)

    def _compare(self):
        if not self.project or not self.comparison_run_id:
            return
        if not self.run_a.currentData() or not self.run_b.currentData():
            self.warning.setText("Select both persisted runs explicitly.")
            return
        try:
            root = functional.compare(self.project, functional.FunctionalSpec(
                self.comparison_run_id, self.module,
                self.run_a.currentData(), self.run_b.currentData()))
            self.history.insertItem(0, root.name, root.name)
            self._remember(root.name)
            self._display(functional.load_functional(self.project, self.comparison_run_id,
                self.module, root.name))
            self.changed.emit()
        except Exception as error:
            self.warning.setText(str(error))
            QMessageBox.warning(self, "Functional comparison unavailable", str(error))

    def _load_selected(self):
        if not self.project or not self.history.currentData():
            return
        try:
            self._display(functional.load_functional(self.project, self.comparison_run_id,
                self.module, self.history.currentData()))
            self._remember(self.history.currentData())
            self.changed.emit()
        except Exception as error:
            self.warning.setText(str(error))
            QMessageBox.warning(self, "Functional history unavailable", str(error))

    def _display(self, result):
        self.loaded = result
        config = result["config"]
        self.warning.setText(f"A: {config['run_a']} | snapshot: {config.get('snapshot_A') or 'not recorded'}\n"
            f"B: {config['run_b']} | snapshot: {config.get('snapshot_B') or 'not recorded'}\n"
            f"Parameters A: {config['parameters']['A']}\nParameters B: {config['parameters']['B']}\n"
            + "\n".join(config.get("warnings", [])))
        self.run_a.setCurrentIndex(self.run_a.findData(config["run_a"]))
        self.run_b.setCurrentIndex(self.run_b.findData(config["run_b"]))
        self.unit.clear()
        for unit in result["manifest"]["units"]:
            self.unit.addItem(unit, unit)
        self._show_unit()
        self.context_a.setEnabled(self.module in {"kegg", "reactome", "string"})
        self.context_b.setEnabled(self.module in {"kegg", "reactome", "string"})

    def _show_unit(self):
        if not self.loaded or not self.unit.currentData():
            return
        unit = self.unit.currentData()
        tables = self.loaded["tables"][unit]
        for name, view in self.views.items():
            view.set_frame(tables[name])
        slug = "_".join(unit.lower().split())
        path = self.loaded["root"] / slug / "counts.png"
        self.graph.setPixmap(QPixmap(str(path)) if path.is_file() else QPixmap())

    def _selected_entity(self):
        view = self.result_tabs.currentWidget()
        if view not in (self.views[key] for key in ("master", "shared", "a_only", "b_only")):
            return None
        record = view.selected_record()
        return record if record and record.get("entity_id") else None

    def _open_context(self, side: str):
        if not self.loaded:
            return
        record = self._selected_entity()
        if not record:
            QMessageBox.information(self, "Choose a result", "Select a result row first.")
            return
        if str(record.get(f"present_{side}", "False")).lower() != "true":
            QMessageBox.information(self, "Not present", f"This result is not present in dataset {side}.")
            return
        config = self.loaded["config"]
        self.context_requested.emit(self.module, config[f"run_{side.lower()}"],
            str(record["entity_id"]), None)

    def _export_workbook(self):
        if not self.loaded:
            return
        source = self.loaded["root"] / "functional_comparison.xlsx"
        target, _ = QFileDialog.getSaveFileName(self, "Export functional workbook",
            source.name, "Excel workbook (*.xlsx)")
        if target:
            shutil.copy2(source, target)


class FunctionalComparisonWidget(QWidget):
    context_requested = Signal(str, str, str, object)

    def __init__(self):
        super().__init__()
        self.project = None
        self.comparison_run_id = None
        self.tabs = QTabWidget()
        self.overview = QLabel("No functional comparisons configured.")
        self.overview.setWordWrap(True)
        self.tabs.addTab(self.overview, "Overview")
        self.panels = {}
        for module in functional.MODULES:
            panel = ModulePanel(module)
            panel.context_requested.connect(self.context_requested)
            panel.changed.connect(self._overview)
            self.panels[module] = panel
            self.tabs.addTab(panel, LABELS[module])
        layout = QVBoxLayout(self)
        layout.addWidget(self.tabs)

    def set_comparison(self, project, comparison_run_id: str | None):
        self.project = project
        self.comparison_run_id = comparison_run_id
        for panel in self.panels.values():
            panel.set_comparison(project, comparison_run_id)
        self._overview()

    def _overview(self):
        rows = []
        for module, panel in self.panels.items():
            if panel.loaded:
                counts = panel.loaded["manifest"]["counts"]
                parts = []
                for unit, values in counts.items():
                    frame = panel.loaded["tables"][unit]["master"]
                    enriched = {label: int((frame.comparison_class == label).sum()) for label in
                        ("Enriched only in A", "Enriched only in B", "Enriched in both")}
                    parts.append(f"{unit} — A only {values['A only']}, shared {values['Shared']}, B only {values['B only']}"
                        + (f"; enriched only A {enriched['Enriched only in A']}, only B {enriched['Enriched only in B']}, both {enriched['Enriched in both']}"
                            if any(enriched.values()) else ""))
                rows.append(f"{LABELS[module]}: " + "; ".join(parts))
            else:
                rows.append(f"No {LABELS[module]} comparison configured.")
        self.overview.setText("\n".join(rows))
