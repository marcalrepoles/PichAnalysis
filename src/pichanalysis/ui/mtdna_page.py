"""Offline mtDNA Evidence analysis UI; all scientific values come from persisted R runs."""
from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd
from PySide6.QtCore import QThread, Qt, QUrl, Signal
from PySide6.QtGui import QDesktopServices, QPixmap
from PySide6.QtWidgets import (QAbstractItemView, QComboBox, QDoubleSpinBox, QFileDialog,
    QFormLayout, QGroupBox, QHBoxLayout, QLabel, QLineEdit, QMessageBox,
    QProgressBar, QPushButton, QSpinBox, QTabWidget, QTableWidget, QTableWidgetItem,
    QVBoxLayout, QWidget)

from ..core.mtdna_analysis import (MtdnaParameters, list_mtdna_runs, map_mtdna_experiment,
    read_mtdna_outputs, run_mtdna_analysis)


class MtdnaAnalysisWorker(QThread):
    succeeded = Signal(object)
    failed = Signal(str)
    target_outside = Signal(dict)

    def __init__(self, project, manager, runtime, run_id, parameters, runner=run_mtdna_analysis):
        super().__init__()
        self.project, self.manager, self.runtime = project, manager, runtime
        self.run_id, self.parameters, self.runner = run_id, parameters, runner

    def run(self):
        from ..core.mtdna_analysis import MtdnaTargetOutsideBackgroundError
        try: self.succeeded.emit(self.runner(self.project, self.manager, self.runtime,
                                              run_id=self.run_id, parameters=self.parameters))
        except MtdnaTargetOutsideBackgroundError as error: self.target_outside.emit(error.details)
        except Exception as error: self.failed.emit(str(error))


class MtdnaPage(QWidget):
    run_requested = Signal(dict)
    open_database_requested = Signal()
    TABLES = (("summary","Summary"),("entity_evidence","Entities"),
        ("category_frequency","Categories"),("enrichment_all","Enrichment"),
        ("enrichment_excluded","Excluded from enrichment"),("source_frequency","Sources"),
        ("source_agreement","Source agreement"),("source_count","Source count"),
        ("evidence_records","Evidence"),("comparison","Comparison"),
        ("mapping","Mapping"),("unmapped","Unmapped"),("ambiguous","Ambiguous"))

    def __init__(self, manager):
        super().__init__()
        self.manager, self.project, self.outputs = manager, None, None
        self._frames = {}; self.mapping_frame = pd.DataFrame(); self._running = False
        self.database = QLabel(); self.database.setWordWrap(True)
        self.open_database = QPushButton("Open Database Manager")
        self.target = QComboBox(); self.background = QComboBox()
        self.manual = QTableWidget(); self.manual.setColumnCount(6)
        self.manual.setHorizontalHeaderLabels(["Gene","NCBI Gene ID","Entity key","UniProt","Original identifiers","mtDNA evidence"])
        self.manual.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.manual.setSelectionMode(QAbstractItemView.SelectionMode.MultiSelection)
        self.background_manual = QTableWidget(); self.background_manual.setColumnCount(3)
        self.background_manual.setHorizontalHeaderLabels(["Gene","Entity key","Original identifier"])
        self.background_manual.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.background_manual.setSelectionMode(QAbstractItemView.SelectionMode.MultiSelection)
        self.presence_help = QLabel(); self.presence_help.setWordWrap(True)
        self.background_help = QLabel("Default: all uniquely resolved gene-level entities in the experiment, including entities without mtDNA evidence.")
        self.background_help.setWordWrap(True)
        defaults = MtdnaParameters()
        self.advanced = QGroupBox("Advanced options")
        self.fdr = QDoubleSpinBox(); self.fdr.setRange(0,1); self.fdr.setDecimals(4); self.fdr.setValue(defaults.fdr_cutoff)
        self.minimum_overlap = QSpinBox(); self.minimum_overlap.setRange(1,1000000); self.minimum_overlap.setValue(defaults.minimum_overlap)
        self.top_n = QSpinBox(); self.top_n.setRange(1,10000); self.top_n.setValue(defaults.top_n)
        advanced_form = QFormLayout(self.advanced)
        advanced_form.addRow("FDR threshold",self.fdr)
        advanced_form.addRow("Minimum overlap",self.minimum_overlap)
        advanced_form.addRow("Top N in plots",self.top_n)
        self.readiness = QLabel(); self.readiness.setWordWrap(True)
        self.run_button = QPushButton("5. Run analysis")
        self.progress = QProgressBar(); self.progress.setRange(0,0); self.progress.hide()
        self.running_text = QLabel()
        config = QWidget(); config_layout = QVBoxLayout(config)
        for widget in (QLabel("1. Database"),self.database,self.open_database,
            QLabel("2. Target set"),self.target,QLabel("Manual target selection"),self.manual,
            self.presence_help,QLabel("3. Background"),self.background,self.background_manual,
            self.background_help,QLabel("4. Analysis settings"),self.advanced,
            QLabel("Evidence related to mtDNA does not by itself establish direct mtDNA binding."),
            self.readiness,self.run_button,self.progress,self.running_text): config_layout.addWidget(widget)
        self.tabs = QTabWidget(); self.tabs.addTab(config,"Configuration")
        self.tables = {key:QTableWidget() for key,_ in self.TABLES}
        self.search = QLineEdit(); self.search.setPlaceholderText("Search gene, Gene ID, UniProt, category, GO term, MitoPathway or source")
        self.enrichment_filter = QComboBox(); self.enrichment_filter.addItem("All tested",False)
        self.enrichment_filter.addItem("Significant only",True)
        self.enrichment_status = QLabel()
        self.entity_detail = QTableWidget(); self.entity_title = QLabel("Select an entity to inspect independent evidence records.")
        self.category_entities = QTableWidget(); self.category_status = QLabel("Select a category to see target entities.")
        self.comparison_disclaimer = QLabel("This comparison is descriptive and does not by itself establish recruitment, loss, or differential association with mtDNA.")
        self.comparison_disclaimer.setWordWrap(True)
        for key,label in self.TABLES:
            table = self.tables[key]
            if key == "entity_evidence": widget = self._page(self.search,table,self.entity_title,self.entity_detail)
            elif key == "category_frequency": widget = self._page(table,self.category_status,self.category_entities)
            elif key == "enrichment_all": widget = self._page(self.enrichment_filter,self.enrichment_status,table)
            elif key == "comparison": widget = self._page(self.comparison_disclaimer,table)
            else: widget = table
            self.tabs.addTab(widget,label)
        self.graph_choice = QComboBox(); self.graph_preview = QLabel("No graphs are available for this run.")
        self.graph_preview.setAlignment(Qt.AlignmentFlag.AlignCenter); self.graph_preview.setMinimumHeight(280)
        self.tabs.addTab(self._page(self.graph_choice,self.graph_preview),"Graphs")
        self.history = QComboBox(); self.history_status = QLabel(); self.history_status.setWordWrap(True)
        self.tabs.addTab(self._page(self.history,self.history_status),"History")
        self.export_table = QPushButton("Export current table...")
        self.export_workbook = QPushButton("Export workbook...")
        self.export_graph = QPushButton("Export graph...")
        self.export_evidence = QPushButton("Export selected entity evidence...")
        self.open_results = QPushButton("Open mtDNA results folder")
        actions = QHBoxLayout()
        for button in (self.export_table,self.export_workbook,self.export_graph,self.export_evidence,self.open_results):
            actions.addWidget(button)
        outer = QVBoxLayout(self); outer.addWidget(QLabel("mtDNA Evidence — Homo sapiens"))
        outer.addWidget(self.tabs,1); outer.addLayout(actions)
        self.open_database.clicked.connect(self.open_database_requested)
        self.run_button.clicked.connect(lambda:self.run_requested.emit(self.parameters()))
        self.target.currentIndexChanged.connect(self._selection_changed)
        self.background.currentIndexChanged.connect(self._selection_changed)
        self.manual.itemSelectionChanged.connect(self._selection_changed)
        self.background_manual.itemSelectionChanged.connect(self._selection_changed)
        self.history.currentIndexChanged.connect(self._load_history)
        self.graph_choice.currentIndexChanged.connect(self._show_graph)
        self.search.textChanged.connect(self._filter_entities)
        self.enrichment_filter.currentIndexChanged.connect(self._filter_enrichment)
        self.tables["entity_evidence"].itemSelectionChanged.connect(self._entity_selected)
        self.tables["category_frequency"].itemSelectionChanged.connect(self._category_selected)
        self.export_table.clicked.connect(self._export_current)
        self.export_workbook.clicked.connect(self._export_workbook)
        self.export_graph.clicked.connect(self._export_graph)
        self.export_evidence.clicked.connect(self._export_evidence)
        self.open_results.clicked.connect(self._open_results)

    @staticmethod
    def _page(*widgets):
        page=QWidget();layout=QVBoxLayout(page)
        for widget in widgets:layout.addWidget(widget)
        return page

    @staticmethod
    def _fill(table, frame):
        table.clear();table.setRowCount(len(frame));table.setColumnCount(len(frame.columns))
        table.setHorizontalHeaderLabels([str(name) for name in frame.columns])
        for i, values in enumerate(frame.itertuples(index=False,name=None)):
            for j, value in enumerate(values):table.setItem(i,j,QTableWidgetItem("" if pd.isna(value) else str(value)))

    def set_project(self, project):
        self.project = project;self.outputs=None;self._frames.clear()
        self._populate_targets();self.refresh_readiness();self._refresh_history()

    def _populate_targets(self):
        self.mapping_frame=pd.DataFrame();classes=[]
        path=self.project.root/"analyses/presence_absence/tables/classification.csv" if self.project else None
        if path and path.is_file():
            try:classes=sorted(pd.read_csv(path).classification.dropna().astype(str).unique())
            except Exception:classes=[]
        for combo in (self.target,self.background):
            combo.blockSignals(True);combo.clear()
            combo.addItem("All uniquely resolved entities in the experiment" if combo is self.background else "All mapped entities","All mapped entities")
            for label in dict.fromkeys(["Reproducibly detected","Shared",*[x for x in classes if x.endswith("-specific")],"Sporadic"]):
                combo.addItem(label,label);combo.model().item(combo.count()-1).setEnabled(bool(classes))
            combo.addItem("Manual selection","Manual selection");combo.blockSignals(False)
        self.presence_help.setText("Presence/Absence-derived targets use existing classifications." if classes else
            "Run Presence / absence to enable derived targets.")
        if self.project and self.manager.mtdna_evidence.is_ready():
            try:self.mapping_frame=map_mtdna_experiment(self.project,self.manager.mtdna_evidence)
            except Exception:self.mapping_frame=pd.DataFrame()
        if not self.mapping_frame.empty:
            frame=self.mapping_frame.rename(columns={"gene_symbol":"Gene","ncbi_gene_id":"NCBI Gene ID",
                "entity_key":"Entity key","original_uniprot":"UniProt",
                "original_identifier":"Original identifiers","evidence_status":"mtDNA evidence"})
            self._fill(self.manual,frame[["Gene","NCBI Gene ID","Entity key","UniProt","Original identifiers","mtDNA evidence"]])
            self._fill(self.background_manual,frame[["Gene","Entity key","Original identifiers"]])
        else:
            self.manual.setRowCount(0);self.background_manual.setRowCount(0)

    def refresh_readiness(self):
        db=self.manager.mtdna_evidence;snapshot=db.active_snapshot();manifest=db.manifest(snapshot) if snapshot else {}
        counts=manifest.get("counts",{})
        self.database.setText("\n".join(("mtDNA Evidence — Homo sapiens",
            f"Status: {'Ready' if snapshot else 'Not installed'}",f"Snapshot: {snapshot.name if snapshot else '—'}",
            f"Evidence entities: {counts.get('entities','—')}",f"Evidence records: {counts.get('evidence_records','—')}",
            f"MitoCarta snapshot: {manifest.get('mitocarta_snapshot_id','—')}",
            f"Gene Ontology snapshot: {manifest.get('go_snapshot_id','—')}",
            f"NCBI reference: {manifest.get('ncbi_accession','—')}")))
        ready=bool(snapshot and self.project and not self.mapping_frame.empty and not self._running)
        self.readiness.setText("Ready for offline analysis." if ready else "mtDNA Evidence database is not ready." if not snapshot else
            "An experimental mapping catalog is required.")
        self.run_button.setEnabled(ready)

    def parameters(self):
        selected=lambda table:tuple(sorted({int(self.mapping_frame.iloc[index.row()].source_row)
            for index in table.selectionModel().selectedRows()})) if not self.mapping_frame.empty else ()
        return dict(target_selection=self.target.currentData(),background_selection=self.background.currentData(),
            manual_rows=selected(self.manual),background_manual_rows=selected(self.background_manual),
            minimum_overlap=self.minimum_overlap.value(),fdr_cutoff=self.fdr.value(),top_n=self.top_n.value())

    def _selection_changed(self):
        self.manual.setVisible(self.target.currentData()=="Manual selection")
        self.background_manual.setVisible(self.background.currentData()=="Manual selection")

    def set_running(self,running):
        self._running=running;self.progress.setVisible(running)
        self.running_text.setText("Running offline mtDNA Evidence analysis..." if running else "")
        self.refresh_readiness()

    def show_outputs(self,outputs):
        self.outputs=outputs;self._frames=outputs["tables"]
        for key,table in self.tables.items():self._fill(table,self._frames[key])
        self._filter_enrichment();self._filter_entities()
        self.enrichment_status.setText("No mtDNA evidence categories met the selected FDR threshold." if
            self._frames["enrichment_significant"].empty else "Category enrichment is relative to the remaining experimental background.")
        self.graph_choice.blockSignals(True);self.graph_choice.clear()
        for path in outputs["plots"]:self.graph_choice.addItem(path.stem,path)
        self.graph_choice.blockSignals(False);self._show_graph()
        self._refresh_history(select=outputs["metadata"]["run_id"])

    def _refresh_history(self,select=None):
        self.history.blockSignals(True);self.history.clear()
        if self.project:
            for run_id in list_mtdna_runs(self.project):
                try:
                    metadata=read_mtdna_outputs(self.project,run_id)["metadata"]
                    label=f"{metadata.get('created_at','')} | {metadata.get('target_definition','')} | {metadata.get('background_definition','')} | {metadata.get('snapshot_id','')} | overlap {metadata.get('minimum_overlap','')} | FDR {metadata.get('fdr_cutoff','')}"
                except Exception:label=f"{run_id} — incomplete"
                self.history.addItem(label,run_id)
        if select:
            index=self.history.findData(select)
            if index>=0:self.history.setCurrentIndex(index)
        self.history.blockSignals(False)

    def _load_history(self):
        if not self.project:return
        run_id=self.history.currentData()
        if not run_id:return
        try:self.show_outputs(read_mtdna_outputs(self.project,run_id));self.history_status.setText(f"Loaded historical run {run_id} offline.")
        except Exception as error:self.history_status.setText(f"Incomplete historical run: {error}")

    def _filter_entities(self):
        if "entity_evidence" not in self._frames:return
        frame=self._frames["entity_evidence"];term=self.search.text().strip()
        if term:
            evidence=self._frames.get("evidence_records",pd.DataFrame())
            matching=set()
            if not evidence.empty:
                mask=evidence.astype(str).apply(lambda column:column.str.contains(term,case=False,regex=False)).any(axis=1)
                matching=set(evidence.loc[mask,"entity_key"].astype(str))
            own=frame.astype(str).apply(lambda column:column.str.contains(term,case=False,regex=False)).any(axis=1)
            frame=frame[own|frame.entity_key.astype(str).isin(matching)]
        self._fill(self.tables["entity_evidence"],frame)

    def _filter_enrichment(self):
        if "enrichment_all" not in self._frames:return
        key="enrichment_significant" if self.enrichment_filter.currentData() else "enrichment_all"
        self._fill(self.tables["enrichment_all"],self._frames[key])

    @staticmethod
    def _selected_value(table,column):
        rows=table.selectionModel().selectedRows()
        headers=[table.horizontalHeaderItem(i).text() for i in range(table.columnCount())]
        if not rows or column not in headers:return ""
        item=table.item(rows[0].row(),headers.index(column))
        return item.text() if item else ""

    def _entity_selected(self):
        key=self._selected_value(self.tables["entity_evidence"],"entity_key")
        if not key or "evidence_records" not in self._frames:return
        frame=self._frames["evidence_records"]
        self._fill(self.entity_detail,frame[frame.entity_key.astype(str)==key])
        self.entity_title.setText(f"Evidence related to mtDNA for {key}; independent source records, not a binding score.")

    def _category_selected(self):
        category=self._selected_value(self.tables["category_frequency"],"category")
        if not category or "category_to_entities" not in self._frames:return
        frame=self._frames["category_to_entities"]
        self._fill(self.category_entities,frame[frame.category.astype(str)==category])
        self.category_status.setText(f"Target entities in {category}")

    def _show_graph(self):
        path=self.graph_choice.currentData()
        if path and Path(path).is_file():
            pixmap=QPixmap(str(path));self.graph_preview.setPixmap(pixmap.scaled(1100,650,
                Qt.AspectRatioMode.KeepAspectRatio,Qt.TransformationMode.SmoothTransformation))
        else:self.graph_preview.setText("No graphs are available for this run.")

    def _export_current(self):
        if not self.outputs:return
        widget=self.tabs.currentWidget();key=next((key for key,table in self.tables.items() if widget is table or table.parent() is widget),None)
        if not key or key not in self._frames:return
        path,_=QFileDialog.getSaveFileName(self,"Export table","mtDNA_table.csv","CSV (*.csv)")
        if path:self._frames[key].to_csv(path,index=False)

    def _export_workbook(self):
        if not self.outputs or not self.outputs["workbook"]:return
        path,_=QFileDialog.getSaveFileName(self,"Export workbook","mtDNA_analysis.xlsx","Excel (*.xlsx)")
        if path:shutil.copy2(self.outputs["workbook"],path)

    def _export_graph(self):
        source=self.graph_choice.currentData()
        if not source:return
        path,_=QFileDialog.getSaveFileName(self,"Export graph",Path(source).name,"PNG (*.png)")
        if path:shutil.copy2(source,path)

    def _export_evidence(self):
        if not self.outputs:return
        key=self._selected_value(self.tables["entity_evidence"],"entity_key")
        if not key:return
        path,_=QFileDialog.getSaveFileName(self,"Export entity evidence",f"{key.replace(':','_')}_evidence.csv","CSV (*.csv)")
        if path:
            frame=self._frames["evidence_records"]
            frame[frame.entity_key.astype(str)==key].to_csv(path,index=False)

    def _open_results(self):
        if self.project:QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.project.root/"analyses/mtDNA")))
