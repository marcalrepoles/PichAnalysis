from __future__ import annotations

from collections import Counter
from pathlib import Path

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QPixmap
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QDoubleSpinBox, QFormLayout, QHBoxLayout, QLabel,
    QListWidget, QListWidgetItem, QPlainTextEdit, QProgressBar, QPushButton,
    QTableWidget, QTableWidgetItem, QTabWidget, QVBoxLayout, QWidget,
)

from ..core.presence_analysis import (
    PresenceOutputs, presence_readiness, quantification_types, quantitative_columns,
    replicate_counts, suggested_minimum,
)
from ..core.project import Project


class PresencePage(QWidget):
    run_requested = Signal(dict)
    export_table_requested = Signal()
    export_workbook_requested = Signal()
    export_graph_requested = Signal(str)
    open_graphs_requested = Signal()

    def __init__(self) -> None:
        super().__init__()
        self.project: Project | None = None
        self.outputs: PresenceOutputs | None = None
        self.quant_type = QComboBox()
        self.conditions = QListWidget()
        self.conditions.setSelectionMode(QListWidget.SelectionMode.MultiSelection)
        self.zero_missing = QCheckBox("Yes — zero means absence")
        self.zero_missing.setChecked(True)
        self.threshold = QDoubleSpinBox(); self.threshold.setRange(-1e15,1e15); self.threshold.setDecimals(6)
        self.rule_mode = QComboBox(); self.rule_mode.addItem("Minimum number of replicates","count"); self.rule_mode.addItem("Minimum fraction of replicates","fraction")
        self.rule_value = QDoubleSpinBox(); self.rule_value.setRange(.01,1e6); self.rule_value.setDecimals(3)
        self.predominant = QCheckBox("Calculate predominant exploratory category"); self.predominant.setChecked(True)
        self.column_table = QTableWidget(0,4); self.column_table.setHorizontalHeaderLabels(["Column","Condition","Replicate","Type"])
        self.method = QPlainTextEdit(); self.method.setReadOnly(True); self.method.setMaximumHeight(120)
        self.readiness = QLabel("Open a project."); self.readiness.setWordWrap(True)
        self.run_button = QPushButton("Run analysis"); self.run_button.setEnabled(False)
        self.progress = QProgressBar(); self.progress.setRange(0,0); self.progress.hide()
        form = QFormLayout(); form.addRow("Quantification used",self.quant_type); form.addRow("Conditions",self.conditions)
        form.addRow("Zero as absence",self.zero_missing); form.addRow("Detection threshold",self.threshold)
        form.addRow("Criterion",self.rule_mode); form.addRow("Minimum value",self.rule_value)
        form.addRow("Exploratory",self.predominant)
        config = QWidget(); config_layout=QVBoxLayout(config); config_layout.addLayout(form); config_layout.addWidget(self.column_table)
        config_layout.addWidget(QLabel("Method")); config_layout.addWidget(self.method); config_layout.addWidget(self.readiness)
        config_layout.addWidget(self.run_button); config_layout.addWidget(self.progress)
        self.summary=QLabel("No results available."); self.summary.setWordWrap(True)
        self.filter=QComboBox(); self.filter.addItem("All")
        self.table=QTableWidget()
        self.detail=QPlainTextEdit(); self.detail.setReadOnly(True)
        results=QWidget(); results_layout=QVBoxLayout(results); results_layout.addWidget(self.summary); results_layout.addWidget(self.filter)
        results_layout.addWidget(self.table,1); results_layout.addWidget(self.detail)
        actions=QHBoxLayout(); self.export_table=QPushButton("Export current table..."); self.export_workbook=QPushButton("Export workbook...")
        actions.addWidget(self.export_table); actions.addWidget(self.export_workbook); results_layout.addLayout(actions)
        self.graph_picker=QComboBox(); self.graph_image=QLabel("No graph available."); self.graph_image.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.graph_image.setMinimumHeight(300); self.graph_image.setScaledContents(False)
        self.export_graph=QPushButton("Export graph..."); self.open_graphs=QPushButton("Open graphs folder")
        graphs=QWidget(); graph_layout=QVBoxLayout(graphs); graph_layout.addWidget(self.graph_picker); graph_layout.addWidget(self.graph_image,1)
        graph_actions=QHBoxLayout(); graph_actions.addWidget(self.export_graph); graph_actions.addWidget(self.open_graphs); graph_layout.addLayout(graph_actions)
        tabs=QTabWidget(); tabs.addTab(config,"Configuration"); tabs.addTab(results,"Results"); tabs.addTab(graphs,"Graphs")
        layout=QVBoxLayout(self); layout.addWidget(tabs)
        self.quant_type.currentIndexChanged.connect(self._refresh_configuration)
        self.conditions.itemSelectionChanged.connect(self._refresh_method)
        self.rule_mode.currentIndexChanged.connect(self._rule_changed)
        for control in (self.zero_missing,self.predominant): control.toggled.connect(self._refresh_method)
        for control in (self.threshold,self.rule_value): control.valueChanged.connect(self._refresh_method)
        self.run_button.clicked.connect(self._emit_run)
        self.filter.currentTextChanged.connect(self._fill_table)
        self.table.currentCellChanged.connect(self._show_detail)
        self.graph_picker.currentIndexChanged.connect(self._show_graph)
        self.export_table.clicked.connect(self.export_table_requested); self.export_workbook.clicked.connect(self.export_workbook_requested)
        self.export_graph.clicked.connect(lambda: self.export_graph_requested.emit(str(self.graph_picker.currentData() or "")))
        self.open_graphs.clicked.connect(self.open_graphs_requested)

    def set_project(self, project: Project | None) -> None:
        self.project=project; self.quant_type.clear()
        if project: self.quant_type.addItems(quantification_types(project))
        self._refresh_configuration()

    def _refresh_configuration(self) -> None:
        self.conditions.clear(); self.column_table.setRowCount(0)
        if not self.project or not self.quant_type.currentText(): self.run_button.setEnabled(False); return
        kind=self.quant_type.currentText(); columns=[x for x in quantitative_columns(self.project) if x["quantification_type"]==kind]
        for condition,count in replicate_counts(self.project,kind).items():
            item=QListWidgetItem(f"{condition} — {count} replicate(s)"); item.setData(Qt.ItemDataRole.UserRole,condition); item.setSelected(True); self.conditions.addItem(item)
        self.column_table.setRowCount(len(columns))
        for row,item in enumerate(columns):
            for col,key in enumerate(("column","condition","replicate","quantification_type")): self.column_table.setItem(row,col,QTableWidgetItem(item[key]))
        self.rule_value.setValue(float(suggested_minimum(self.project,kind))); self._refresh_method()

    def _rule_changed(self) -> None:
        fraction=self.rule_mode.currentData()=="fraction"; self.rule_value.setMaximum(1 if fraction else 1e6)
        if fraction and self.rule_value.value()>1: self.rule_value.setValue(.5)
        self._refresh_method()

    def selected_conditions(self) -> list[str]:
        return [item.data(Qt.ItemDataRole.UserRole) for item in self.conditions.selectedItems()]

    def _refresh_method(self) -> None:
        if not self.project: return
        kind=self.quant_type.currentText(); state=presence_readiness(self.project,kind or None); selected=self.selected_conditions()
        counts=replicate_counts(self.project,kind) if kind else {}; description="; ".join(f"{c}: {counts[c]} replicate(s)" for c in selected)
        rule=(f"at least {self.rule_value.value():g} replicate(s)" if self.rule_mode.currentData()=="count" else f"minimum fraction {self.rule_value.value():g}")
        zero="Zero will be treated as absence." if self.zero_missing.isChecked() else "Zero may be considered detected at a zero threshold."
        self.method.setPlainText(f"{description}\nDetected when the value meets the threshold {self.threshold.value():g}. {zero}\nReproducible when detected in {rule}. Specific requires complete absence in the opposite condition.")
        self.readiness.setText(state.reason); self.run_button.setEnabled(state.ready and bool(selected))

    def _emit_run(self) -> None:
        self.run_requested.emit({"quantification_type":self.quant_type.currentText(),"selected_conditions":self.selected_conditions(),
            "zero_is_missing":self.zero_missing.isChecked(),"threshold":self.threshold.value(),"rule_mode":self.rule_mode.currentData(),
            "rule_value":self.rule_value.value(),"predominant":self.predominant.isChecked()})

    def set_running(self,running:bool) -> None: self.progress.setVisible(running); self.run_button.setEnabled(not running)

    def show_outputs(self,outputs:PresenceOutputs) -> None:
        self.outputs=outputs; counts=outputs.metadata.get("classification_counts",{})
        self.summary.setText(f"Rows analyzed: {outputs.metadata.get('total_entities',0)} | "+" | ".join(f"{k}: {v}" for k,v in counts.items()))
        self.filter.clear(); self.filter.addItem("All"); self.filter.addItems(sorted(counts)); self._fill_table()
        self.graph_picker.clear()
        for path in outputs.graphs: self.graph_picker.addItem(path.stem,str(path))

    def _fill_table(self) -> None:
        if not self.outputs:return
        frame=self.outputs.classification; chosen=self.filter.currentText()
        if chosen and chosen!="All": frame=frame[frame["classification"]==chosen]
        self.table.setRowCount(len(frame)); self.table.setColumnCount(len(frame.columns)); self.table.setHorizontalHeaderLabels(list(frame.columns))
        for r,row in enumerate(frame.itertuples(index=False,name=None)):
            for c,value in enumerate(row): self.table.setItem(r,c,QTableWidgetItem("" if value is None else str(value)))

    def _show_detail(self,row:int,*_args) -> None:
        if row<0:return
        headers=[self.table.horizontalHeaderItem(i).text() for i in range(self.table.columnCount())]
        wanted=("original_id","uniprot_accession","gene_symbol","protein_name","classification","ncbi_summary")
        self.detail.setPlainText("\n\n".join(f"{x}: {self.table.item(row,headers.index(x)).text()}" for x in wanted if x in headers and self.table.item(row,headers.index(x))))

    def _show_graph(self) -> None:
        path=self.graph_picker.currentData()
        if path:
            pixmap=QPixmap(str(path)); self.graph_image.setPixmap(pixmap.scaled(self.graph_image.size(),Qt.AspectRatioMode.KeepAspectRatio,Qt.TransformationMode.SmoothTransformation))
