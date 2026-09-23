"""Complex Portal analysis UI: scientific values are read only from persisted runs."""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import pandas as pd
from PySide6.QtCore import QThread, Qt, QUrl, Signal
from PySide6.QtGui import QDesktopServices, QPixmap
from PySide6.QtWidgets import (QAbstractItemView, QComboBox, QDoubleSpinBox, QFileDialog,
    QFormLayout, QGroupBox, QHBoxLayout, QLabel, QLineEdit, QListWidget,
    QListWidgetItem, QMessageBox, QProgressBar, QPushButton, QSpinBox, QTabWidget,
    QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget)

from ..core.complex_analysis import (ComplexParameters, ComplexTargetOutsideBackgroundError,
    list_complex_runs, map_complex_experiment, read_complex_outputs, run_complex_analysis)

COVERAGE_LABELS = {
    "complete_protein_component_coverage": "Complete protein-component coverage",
    "partial_protein_component_coverage": "Partial protein-component coverage",
    "no_detected_protein_components": "No detected protein components",
    "not_applicable_no_protein_components": "Not applicable — no protein components",
}


class ComplexAnalysisWorker(QThread):
    succeeded = Signal(object)
    failed = Signal(str)
    target_outside = Signal(dict)

    def __init__(self, project, manager, runtime, run_id, parameters, runner=run_complex_analysis):
        super().__init__()
        self.project, self.manager, self.runtime = project, manager, runtime
        self.run_id, self.parameters, self.runner = run_id, parameters, runner

    def run(self):
        try:
            self.succeeded.emit(self.runner(self.project, self.manager, self.runtime,
                                            run_id=self.run_id, parameters=self.parameters))
        except ComplexTargetOutsideBackgroundError as error:
            self.target_outside.emit(error.details)
        except Exception as error:
            self.failed.emit(str(error))


class ComplexPage(QWidget):
    run_requested = Signal(dict)
    open_database_requested = Signal()
    TABLES = (("summary","Summary"),("complex_coverage","Complex coverage"),
              ("component_groups","Component groups"),("alternative_component_groups","Alternative component groups"),
              ("direct_participants","Direct participants"),("expanded_protein_components","Expanded protein components"),
              ("complex_frequency","Frequency"),("complex_enrichment_all","Enrichment"),
              ("complex_enrichment_excluded","Excluded from enrichment"),("protein_to_complexes","Proteins"),
              ("complex_to_proteins","Complex to proteins"),("nonprotein_participants","Non-protein participants"),
              ("nested_complexes","Nested complexes"),("mapping","Mapping"),("unmapped","Unmapped"),("ambiguous","Ambiguous"))

    def __init__(self, manager):
        super().__init__()
        self.manager, self.project, self.outputs = manager, None, None
        self.mapping_frame = pd.DataFrame(); self._frames = {}; self._running = False
        defaults = ComplexParameters()
        self.database = QLabel(); self.database.setWordWrap(True)
        self.readiness = QLabel(); self.readiness.setWordWrap(True)
        self.open_database = QPushButton("Open Database Manager")
        self.open_database_folder = QPushButton("Open Complex Portal database folder")
        self.target = QComboBox(); self.background = QComboBox()
        self.manual = QTableWidget(); self.manual.setColumnCount(5)
        self.manual.setHorizontalHeaderLabels(["Gene","Original UniProt","Canonical UniProt","Protein name","Original ID"])
        self.manual.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.manual.setSelectionMode(QAbstractItemView.SelectionMode.MultiSelection)
        self.background_manual = QListWidget(); self.background_manual.setSelectionMode(QAbstractItemView.SelectionMode.MultiSelection)
        self.presence_help = QLabel(); self.presence_help.setWordWrap(True)
        self.target_count = QLabel("Uniquely resolved canonical proteins: 0")
        self.advanced = QGroupBox("Advanced options"); self.advanced.setCheckable(True); self.advanced.setChecked(False)
        self.fdr = QDoubleSpinBox(); self.fdr.setRange(0,1); self.fdr.setDecimals(4); self.fdr.setValue(defaults.fdr_cutoff)
        self.minimum_overlap = QSpinBox(); self.minimum_overlap.setRange(1,1000000); self.minimum_overlap.setValue(defaults.minimum_overlap)
        self.top_n = QSpinBox(); self.top_n.setRange(1,10000); self.top_n.setValue(defaults.top_n)
        form = QFormLayout(self.advanced)
        form.addRow("FDR threshold",self.fdr); form.addRow("Minimum overlap",self.minimum_overlap); form.addRow("Top N in plots",self.top_n)
        self.run_button = QPushButton("5. Run analysis")
        self.progress = QProgressBar(); self.progress.setRange(0,0); self.progress.hide()
        self.running_text = QLabel()
        config = QWidget(); layout = QVBoxLayout(config)
        for widget in (QLabel("1. Database"),self.database,self.open_database,self.open_database_folder,
            self._help("This analysis uses manually curated human complexes from Complex Portal. Predicted hu.MAP/MuSIC complexes are not included."),
            QLabel("2. Target set"),self.target,
            QLabel("Manual target selection"),self.manual,
            self.target_count,self.presence_help,QLabel("3. Background"),self.background,self.background_manual,
            self._help("The background represents experimentally eligible proteins, including proteins without curated Complex Portal membership."),
            QLabel("4. Analysis settings"),self.advanced,
            self._help("Protein-component coverage describes which curated protein component groups were detected. It does not demonstrate that the complex is physically assembled."),
            self._help("Alternative component proteins belong to one component group for coverage. Detecting multiple alternatives does not increase the number of covered component groups."),
            self._help("Frequency counts unique target proteins annotated as possible members of the complex. It is distinct from component-group coverage."),
            self._help("Enrichment tests whether proteins annotated as possible components of a curated complex are over-represented in the selected target relative to the remaining experimental background."),
            self._help("Curated stoichiometry is descriptive metadata and is not interpreted as proteomic abundance."),
            self.readiness,self.run_button,self.progress,self.running_text): layout.addWidget(widget)
        self.tabs = QTabWidget(); self.tabs.addTab(config,"Configuration")
        self.tables = {key:QTableWidget() for key,_ in self.TABLES}
        self.search = QLineEdit(); self.search.setPlaceholderText("Search complex ID, name, alias, UniProt or gene")
        self.coverage_filter = QComboBox(); self.coverage_filter.addItem("All","")
        for code,label in COVERAGE_LABELS.items(): self.coverage_filter.addItem(label,code)
        self.enrichment_filter = QComboBox(); self.enrichment_filter.addItem("All tested",False); self.enrichment_filter.addItem("Significant only",True)
        self.enrichment_status = QLabel()
        self.protein_memberships = QTableWidget()
        self.protein_status = QLabel("Select a protein to see all curated complex memberships.")
        for key,label in self.TABLES:
            widget = self.tables[key]
            if key=="complex_coverage": widget=self._page(self.search,self.coverage_filter,widget)
            elif key=="complex_enrichment_all": widget=self._page(self.enrichment_filter,self.enrichment_status,widget)
            elif key=="protein_to_complexes": widget=self._page(widget,self.protein_status,self.protein_memberships)
            self.tabs.addTab(widget,label)
        self.detail_title = QLabel("Select a complex in a results table."); self.detail_title.setWordWrap(True)
        self.coverage_progress = QProgressBar(); self.coverage_progress.setFormat("Protein-component coverage: no complex selected")
        self.detail_tables = {key:QTableWidget() for key in ("component_groups","complex_to_proteins","direct_participants",
            "expanded_protein_components","nonprotein_participants","nested_complexes")}
        details = QTabWidget()
        for key,label in (("component_groups","Component groups"),("complex_to_proteins","Possible / detected proteins"),
            ("direct_participants","Direct participants"),("expanded_protein_components","Expanded protein components"),
            ("nonprotein_participants","Non-protein participants"),("nested_complexes","Nested complexes")):
            details.addTab(self.detail_tables[key],label)
        self.tabs.addTab(self._page(self.detail_title,self.coverage_progress,details),"Complex details")
        self.graph_choice = QComboBox(); self.graph_preview = QLabel("No plots are available for this run.")
        self.graph_preview.setAlignment(Qt.AlignmentFlag.AlignCenter); self.graph_preview.setMinimumHeight(280)
        self.tabs.addTab(self._page(self.graph_choice,self.graph_preview),"Graphs")
        self.history = QComboBox(); self.history_status = QLabel(); self.history_status.setWordWrap(True)
        self.tabs.addTab(self._page(self.history,self.history_status),"History")
        self.export_table = QPushButton("Export current table..."); self.export_workbook = QPushButton("Export workbook...")
        self.export_graph = QPushButton("Export graph..."); self.export_details = QPushButton("Export complex details...")
        self.open_results = QPushButton("Open Complexes results folder")
        actions = QHBoxLayout()
        for button in (self.export_table,self.export_workbook,self.export_graph,self.export_details,self.open_results):actions.addWidget(button)
        outer=QVBoxLayout(self); outer.addWidget(QLabel("Complexes — Homo sapiens")); outer.addWidget(self.tabs,1); outer.addLayout(actions)
        self.open_database.clicked.connect(self.open_database_requested)
        self.open_database_folder.clicked.connect(self._open_database_folder)
        self.target.currentIndexChanged.connect(self._selection_changed)
        self.background.currentIndexChanged.connect(self._selection_changed)
        self.manual.itemSelectionChanged.connect(self._selection_changed)
        self.run_button.clicked.connect(lambda:self.run_requested.emit(self.parameters()))
        self.coverage_filter.currentIndexChanged.connect(self._filter_coverage)
        self.enrichment_filter.currentIndexChanged.connect(self._filter_enrichment)
        self.search.textChanged.connect(self._filter_coverage)
        self.history.currentIndexChanged.connect(self._load_history)
        self.graph_choice.currentIndexChanged.connect(self._show_graph)
        for table in self.tables.values():table.itemSelectionChanged.connect(lambda t=table:self._table_selected(t))
        self.export_table.clicked.connect(self._export_current); self.export_workbook.clicked.connect(self._export_workbook)
        self.export_graph.clicked.connect(self._export_graph); self.export_details.clicked.connect(self._export_details)
        self.open_results.clicked.connect(self._open_results)

    @staticmethod
    def _help(text):
        widget=QLabel(text); widget.setWordWrap(True); return widget

    @staticmethod
    def _page(*widgets):
        page=QWidget(); layout=QVBoxLayout(page)
        for widget in widgets:layout.addWidget(widget)
        return page

    def set_project(self,project):
        self.project=project; self._clear_results(); self._populate_targets(); self.refresh_readiness(); self._refresh_history()

    def _clear_results(self):
        self.outputs=None; self._frames.clear()
        for table in (*self.tables.values(),*self.detail_tables.values(),self.protein_memberships):
            table.clear(); table.setRowCount(0); table.setColumnCount(0)
        self.graph_choice.clear(); self.graph_preview.setPixmap(QPixmap()); self.graph_preview.setText("No plots are available for this run.")
        self.detail_title.setText("Select a complex in a results table.")
        self.protein_status.setText("Select a protein to see all curated complex memberships.")
        self.coverage_progress.setFormat("Protein-component coverage: no complex selected")

    def _populate_targets(self):
        self.mapping_frame=pd.DataFrame()
        path=self.project.root/"analyses/presence_absence/tables/classification.csv" if self.project else None
        classes=[]
        if path and path.is_file():
            try:classes=sorted(pd.read_csv(path).classification.dropna().astype(str).unique())
            except Exception:classes=[]
        for combo in (self.target,self.background):
            combo.blockSignals(True); combo.clear(); combo.addItem("All uniquely resolved proteins in the experiment" if combo is self.background else "All mapped proteins","All mapped proteins")
            for label in dict.fromkeys(["Reproducibly detected","Shared",*[x for x in classes if x.endswith("-specific")],"Sporadic"]):
                combo.addItem(label,label); combo.model().item(combo.count()-1).setEnabled(bool(classes))
            combo.addItem("Manual selection","Manual selection"); combo.blockSignals(False)
        self.presence_help.setText("Presence/Absence-derived selections use existing results." if classes else "Run Presence / absence to enable derived selections.")
        self.manual.setRowCount(0); self.background_manual.clear()
        names = {}
        catalog = self.project.root/"mapping/tables/protein_catalog.csv" if self.project else None
        if catalog and catalog.is_file():
            try:
                source = pd.read_csv(catalog,dtype=str).fillna("")
                if "protein_name" in source:
                    names = dict(zip(source.source_row.astype(str),source.protein_name.astype(str)))
            except Exception: pass
        if self.project and self.manager.complex_portal.is_ready():
            try:self.mapping_frame=map_complex_experiment(self.project,self.manager.complex_portal)
            except Exception:pass
        for row in self.mapping_frame.to_dict("records"):
            values=[row["gene_symbol"],row["original_uniprot"],row["canonical_uniprot"],
                    names.get(str(row["source_row"]),""),row["original_id"]]
            index=self.manual.rowCount(); self.manual.insertRow(index)
            for column,value in enumerate(values):
                cell=QTableWidgetItem(str(value))
                cell.setData(Qt.ItemDataRole.UserRole,int(row["source_row"]))
                self.manual.setItem(index,column,cell)
            if row["isoform_normalized"]:
                self.manual.item(index,1).setToolTip(
                    f"Original: {row['original_uniprot']}\nComplex Portal lookup: {row['canonical_uniprot']}")
            value=" | ".join(str(x) for x in values)
            item=QListWidgetItem(value); item.setData(Qt.ItemDataRole.UserRole,int(row["source_row"]))
            self.background_manual.addItem(item)
        self._selection_changed()

    def _selection_changed(self):
        self.manual.setEnabled(self.target.currentData()=="Manual selection")
        self.background_manual.setEnabled(self.background.currentData()=="Manual selection")
        count=self.mapping_frame.loc[self.mapping_frame.mapping_status=="mapped_unique","canonical_uniprot"].nunique() if not self.mapping_frame.empty else 0
        self.target_count.setText(f"Uniquely resolved canonical proteins: {count}")
        self.run_button.setEnabled(bool(count) and not self._running and self.readiness.text()=="Ready for Complex Portal analysis.")

    def refresh_readiness(self):
        db=self.manager.complex_portal; snap=db.active_snapshot(); m=db.manifest(snap) if snap else {}
        self.database.setText(f"Complex Portal — Homo sapiens\nStatus: {'Ready' if snap else 'Not installed'}\nScope: Manually curated complexes\nRelease: {m.get('complex_portal_release','—')}\nSnapshot: {snap.name if snap else '—'}\nCurated complexes: {m.get('complex_count','—')}\nUnique proteins: {m.get('unique_protein_count','—')}")
        if not snap:reason="Complex Portal database is not installed."
        elif not self.project:reason="Open a project to analyze complexes."
        elif str(self.project.config.get("organism_tax_id") or "")!="9606":reason="Complex Portal analysis currently supports Homo sapiens only."
        elif not (self.project.root/"mapping/tables/protein_catalog.csv").is_file():reason="Run protein mapping before Complex Portal analysis."
        else:reason="Ready for Complex Portal analysis."
        self.readiness.setText(reason); self._selection_changed()

    def parameters(self):
        return dict(target_selection=self.target.currentData(),background_selection=self.background.currentData(),
            manual_rows=tuple(int(self.manual.item(row,0).data(Qt.ItemDataRole.UserRole)) for row in sorted({item.row() for item in self.manual.selectedItems()})),
            background_manual_rows=tuple(int(item.data(Qt.ItemDataRole.UserRole)) for item in self.background_manual.selectedItems()),
            minimum_overlap=self.minimum_overlap.value(),fdr_cutoff=self.fdr.value(),top_n=self.top_n.value())

    def set_running(self,running):
        self._running=running; self.progress.setVisible(running)
        self.running_text.setText("Running complex analysis..." if running else "")
        self._selection_changed()

    def _fill(self,table,frame):
        self._frames[table]=frame.copy(); table.clear(); table.setRowCount(len(frame)); table.setColumnCount(len(frame.columns))
        table.setHorizontalHeaderLabels([str(name).replace("_"," ").title() for name in frame.columns])
        for i,values in enumerate(frame.itertuples(index=False,name=None)):
            for j,value in enumerate(values):
                if pd.isna(value):value=""
                elif str(value) in COVERAGE_LABELS:value=COVERAGE_LABELS[str(value)]
                elif str(frame.columns[j])=="group_covered":value="Covered" if str(value).lower() in ("true","1") else "Not detected"
                elif str(frame.columns[j]) in ("possible_uniprot_accessions","detected_uniprot_accessions","possible_accessions","detected_accessions"):value=str(value).replace(";"," | ")
                elif str(value)=="0" and "stoichiometry" in str(frame.columns[j]).lower():value="Unknown"
                table.setItem(i,j,QTableWidgetItem(str(value)))

    def show_outputs(self,outputs):
        self._clear_results(); self.outputs=outputs; tables=outputs["tables"]; run=outputs["run_root"]
        for key,widget in self.tables.items():
            frame=tables.get(key)
            if frame is None:
                path=run/"inputs"/f"{key}.csv"; frame=pd.read_csv(path) if path.is_file() else pd.DataFrame()
            self._fill(widget,frame)
        self._filter_coverage(); self._filter_enrichment()
        if tables["complex_enrichment_significant"].empty:
            self.enrichment_status.setText("No complexes met the selected FDR threshold.")
        else:self.enrichment_status.setText(f"{len(tables['complex_enrichment_significant'])} complexes met the selected FDR threshold.")
        self.graph_choice.blockSignals(True); self.graph_choice.clear()
        for path in sorted(outputs["plots"]):
            if path.suffix.lower()==".png":self.graph_choice.addItem(path.stem.replace("_"," ").title(),{"png":path,"pdf":path.with_suffix(".pdf")})
        self.graph_choice.blockSignals(False); self._show_graph()
        m=outputs["metadata"]
        self.history_status.setText(f"Loaded run {m.get('run_id')} | Snapshot {m.get('snapshot_id')} | Release {m.get('complex_portal_release')}")
        self._refresh_history(m.get("run_id"))

    def _filter_coverage(self):
        table=self.tables["complex_coverage"]; frame=self._frames.get(table)
        if frame is None:return
        code=self.coverage_filter.currentData(); query=self.search.text().strip().casefold()
        details=self.outputs["tables"].get("complex_details",pd.DataFrame()) if self.outputs else pd.DataFrame()
        proteins=self.outputs["tables"].get("protein_to_complexes",pd.DataFrame()) if self.outputs else pd.DataFrame()
        for i,row in frame.reset_index(drop=True).iterrows():
            visible=not code or str(row.get("coverage_class",""))==code
            if query:
                values=[str(row.get(k,"")) for k in ("complex_id","complex_name","expected_possible_proteins","detected_member_proteins")]
                if "complex_id" in details:
                    selected=details[details.complex_id.astype(str)==str(row.get("complex_id"))]
                    if not selected.empty:values += [str(selected.iloc[0].get(k,"")) for k in ("aliases","recommended_name")]
                if "complex_id" in proteins:
                    selected=proteins[proteins.complex_id.astype(str)==str(row.get("complex_id"))]
                    values += selected.get("gene_symbol",pd.Series(dtype=str)).astype(str).tolist()
                visible=visible and any(query in value.casefold() for value in values)
            table.setRowHidden(i,not visible)

    def _filter_enrichment(self):
        table=self.tables["complex_enrichment_all"]; frame=self._frames.get(table)
        if frame is None:return
        significant=self.outputs["tables"].get("complex_enrichment_significant",pd.DataFrame()) if self.outputs else pd.DataFrame()
        ids=set(significant.complex_id.astype(str)) if "complex_id" in significant else set()
        for i,row in frame.reset_index(drop=True).iterrows():table.setRowHidden(i,bool(self.enrichment_filter.currentData()) and str(row.get("complex_id")) not in ids)

    def _table_selected(self,table):
        frame=self._frames.get(table); row=table.currentRow()
        if frame is not None and 0<=row<len(frame):
            if table is self.tables["protein_to_complexes"] and "canonical_uniprot" in frame:
                protein=str(frame.iloc[row]["canonical_uniprot"])
                related=frame[frame.canonical_uniprot.astype(str)==protein]
                self.protein_status.setText(f"{protein}: {len(related)} curated complex memberships. Multiple memberships are not mapping ambiguity.")
                self._fill(self.protein_memberships,related)
            if "complex_id" in frame:self.show_complex_details(str(frame.iloc[row]["complex_id"]))

    def show_complex_details(self,complex_id):
        if not self.outputs:return
        tables=self.outputs["tables"]; run=self.outputs["run_root"]
        coverage=tables["complex_coverage"]; selected=coverage[coverage.complex_id.astype(str)==complex_id]
        if selected.empty:return
        record=selected.iloc[0]; total=int(record.total_protein_component_groups); covered=int(record.covered_component_groups)
        self.coverage_progress.setMaximum(max(total,1)); self.coverage_progress.setValue(covered)
        self.coverage_progress.setFormat(f"Protein-component coverage: {covered} / {total} groups")
        detail=tables.get("complex_details",pd.DataFrame())
        detail=detail[detail.complex_id.astype(str)==complex_id] if "complex_id" in detail else pd.DataFrame()
        value=detail.iloc[0] if not detail.empty else record
        lines=[f"Complex Portal ID: {complex_id}",f"Recommended name: {value.get('recommended_name',record.get('complex_name',''))}",
            f"Description: {value.get('description','')}",f"Aliases: {value.get('aliases','')}",
            f"Coverage class: {COVERAGE_LABELS.get(str(record.coverage_class),record.coverage_class)}",
            f"Protein-component coverage: {covered} / {total} groups",f"Confidence term: {value.get('confidence_name','')}",
            f"Evidence: {value.get('experimental_evidence','')}",f"Source: {value.get('source','')}"]
        if bool(record.get("has_nonprotein_participants")) and total and covered==total:
            lines.append("All curated protein component groups were detected. Non-protein components are not assessed by this proteomic analysis.")
        self.detail_title.setText("\n".join(lines))
        for key,widget in self.detail_tables.items():
            frame=tables.get(key)
            if frame is None:
                path=run/"inputs"/f"{key}.csv"; frame=pd.read_csv(path) if path.is_file() else pd.DataFrame()
            column="parent_complex_id" if key=="nested_complexes" else "complex_id"
            self._fill(widget,frame[frame[column].astype(str)==complex_id] if column in frame else pd.DataFrame())

    def _show_graph(self):
        data=self.graph_choice.currentData()
        if not data or not data["png"].is_file():self.graph_preview.setPixmap(QPixmap());self.graph_preview.setText("No plots are available for this run.");return
        self.graph_preview.setPixmap(QPixmap(str(data["png"])).scaled(900,580,Qt.AspectRatioMode.KeepAspectRatio,Qt.TransformationMode.SmoothTransformation))

    def _refresh_history(self,selected=None):
        self.history.blockSignals(True); self.history.clear()
        if self.project:
            for run_id in list_complex_runs(self.project):
                try:
                    m=json.loads((self.project.root/"analyses/Complexes/runs"/run_id/"metadata.json").read_text(encoding="utf-8"))
                    label=f"{m.get('created_at','—')} | Target {m.get('target_definition','—')} | Background {m.get('background_definition','—')} | Snapshot {m.get('snapshot_id','—')} | Release {m.get('complex_portal_release','—')} | Overlap ≥{m.get('minimum_overlap','—')} | FDR {m.get('fdr_cutoff','—')}"
                except Exception:label=f"{run_id} | Incomplete metadata"
                self.history.addItem(label,run_id)
        self.history.setCurrentIndex(self.history.findData(selected) if selected else -1);self.history.blockSignals(False)

    def _load_history(self,index):
        if not self.project or index<0:return
        try:self.show_outputs(read_complex_outputs(self.project,self.history.itemData(index)))
        except Exception as error:self._clear_results();self.history_status.setText(f"Historical run is incomplete: {error}")

    def _current_table(self):
        current=self.tabs.currentWidget()
        if isinstance(current,QTableWidget):return current
        tables=current.findChildren(QTableWidget) if current else []
        return tables[0] if tables else None

    def _export_current(self):
        frame=self._frames.get(self._current_table())
        if frame is None:return QMessageBox.information(self,"Export table","The current panel does not contain an exportable table.")
        path,_=QFileDialog.getSaveFileName(self,"Export current table","complex_table.csv","CSV (*.csv)")
        if path:frame.to_csv(path,index=False)

    def _export_workbook(self):
        if not self.outputs:return
        path,_=QFileDialog.getSaveFileName(self,"Export workbook","Complex_analysis.xlsx","Excel workbook (*.xlsx)")
        if path:shutil.copy2(self.outputs["workbook"],path)

    def _export_graph(self):
        data=self.graph_choice.currentData()
        if not data:return
        path,selected=QFileDialog.getSaveFileName(self,"Export graph","complex_plot.png","PNG (*.png);;PDF (*.pdf)")
        if path:
            source=data["pdf" if "PDF" in selected or Path(path).suffix.lower()==".pdf" else "png"]
            if source.is_file():shutil.copy2(source,path)
            else:QMessageBox.information(self,"Export graph","This run has no graph in the selected format.")

    def _export_details(self):
        groups=self._frames.get(self.detail_tables["component_groups"])
        if groups is None or groups.empty:
            return QMessageBox.information(self,"Export complex details","Select a complex first.")
        sections=[]
        for key,table in self.detail_tables.items():
            frame=self._frames.get(table)
            if frame is not None and not frame.empty:
                section=frame.copy()
                section.insert(0,"section",key.replace("_"," ").title())
                sections.append(section)
        path,_=QFileDialog.getSaveFileName(self,"Export complex details","complex_details.csv","CSV (*.csv)")
        if path:pd.concat(sections,ignore_index=True,sort=False).to_csv(path,index=False)

    def _open_results(self):
        if self.outputs:QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.outputs["run_root"])))

    def _open_database_folder(self):
        snapshot=self.manager.complex_portal.active_snapshot()
        if snapshot:QDesktopServices.openUrl(QUrl.fromLocalFile(str(snapshot)))
