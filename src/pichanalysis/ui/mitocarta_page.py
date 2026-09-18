from __future__ import annotations

import json
from pathlib import Path
from typing import Callable

import pandas as pd
from PySide6.QtCore import QThread,Qt,QUrl,Signal
from PySide6.QtGui import QDesktopServices,QPixmap
from PySide6.QtWidgets import (QAbstractItemView,QComboBox,QDoubleSpinBox,QFileDialog,QFormLayout,QGroupBox,QHBoxLayout,QLabel,QListWidget,QListWidgetItem,QMessageBox,QProgressBar,QPushButton,QSpinBox,QTableWidget,QTableWidgetItem,QTabWidget,QVBoxLayout,QWidget)

from ..core.database_manager import DatabaseManager
from ..core.mitocarta_analysis import (MissingMitoCartaOutputError,MitoCartaOutputs,MitoCartaParameters,MitoCartaTargetOutsideBackgroundError,export_mitocarta_artifact,list_mitocarta_runs,mitocarta_compartment_genes,mitocarta_gene_pathways,mitocarta_gene_subcompartments,mitocarta_pathway_genes,mitocarta_readiness,mitocarta_target_genes,read_mitocarta_outputs,run_mitocarta_analysis)
from ..core.project import Project
from ..core.r_runtime import RRuntime


class MitoCartaAnalysisWorker(QThread):
    succeeded=Signal(object);target_outside_background=Signal(dict);failed=Signal(str)
    def __init__(self,project:Project,manager:DatabaseManager,runtime:RRuntime,run_id:str,parameters:MitoCartaParameters,runner:Callable=run_mitocarta_analysis):super().__init__();self.project=project;self.manager=manager;self.runtime=runtime;self.run_id=run_id;self.parameters=parameters;self.runner=runner
    def run(self):
        try:self.succeeded.emit(self.runner(self.project,self.manager,self.runtime,run_id=self.run_id,parameters=self.parameters))
        except MitoCartaTargetOutsideBackgroundError as error:self.target_outside_background.emit(error.details)
        except Exception as error:self.failed.emit(str(error))


class MitoCartaPage(QWidget):
    run_requested=Signal(dict);open_database_requested=Signal()
    def __init__(self,manager:DatabaseManager):
        super().__init__();self.manager=manager;self.project=None;self.outputs=None;self._current_frames={}
        self.database=QLabel();self.database.setWordWrap(True);self.open_database=QPushButton("Open Database Manager")
        self.target=QComboBox();self.background=QComboBox();self.background_help=QLabel("Background represents the genes that could have been observed in this experiment.");self.background_help.setWordWrap(True)
        self.manual=QListWidget();self.manual.setSelectionMode(QAbstractItemView.SelectionMode.MultiSelection);self.manual.setMinimumHeight(100);self.manual_help=QLabel("Manual selection preserves the stable experimental source row; display labels are not scientific keys.");self.manual_help.setWordWrap(True);self.presence_help=QLabel();self.presence_help.setWordWrap(True)
        self.advanced=QGroupBox("Advanced options");self.advanced.setCheckable(True);self.advanced.setChecked(False);defaults=MitoCartaParameters();self.fdr=QDoubleSpinBox();self.fdr.setRange(.0001,1);self.fdr.setDecimals(4);self.fdr.setValue(defaults.fdr_cutoff);self.minimum=QSpinBox();self.minimum.setRange(1,100000);self.minimum.setValue(defaults.minimum_overlap);self.top_n=QSpinBox();self.top_n.setRange(1,1000);self.top_n.setValue(defaults.top_n);advanced=QFormLayout(self.advanced);advanced.addRow("FDR threshold",self.fdr);advanced.addRow("Minimum overlap",self.minimum);advanced.addRow("Top N pathways in plots",self.top_n);self.advanced.toggled.connect(self._toggle_advanced);self._toggle_advanced(False)
        self.method_help=QLabel("MitoCarta membership: Shows which experimentally resolved genes belong to MitoCarta3.0.\nOverall enrichment: Tests whether MitoCarta genes are more frequent in the selected target than in the remaining experimental background.\nSub-compartment enrichment: Tests mitochondrial sub-compartments using only MitoCarta-positive genes in the experimental background.\nMitoPathway enrichment: Tests MitoPathways using only MitoCarta-positive genes in the experimental background.");self.method_help.setWordWrap(True)
        self.readiness=QLabel();self.readiness.setWordWrap(True);self.run_button=QPushButton("5. Run MitoCarta analysis");self.progress=QProgressBar();self.progress.setRange(0,0);self.progress.hide();self.running_text=QLabel()
        config=QWidget();form=QFormLayout();form.addRow("2. Target set",self.target);form.addRow("Manual selection",self.manual);form.addRow("3. Background",self.background);layout=QVBoxLayout(config);layout.addWidget(QLabel("1. Database"));layout.addWidget(self.database);layout.addWidget(self.open_database);layout.addLayout(form);layout.addWidget(self.manual_help);layout.addWidget(self.presence_help);layout.addWidget(self.background_help);layout.addWidget(QLabel("4. Analysis settings"));layout.addWidget(self.advanced);layout.addWidget(self.method_help);layout.addWidget(self.readiness);layout.addWidget(self.run_button);layout.addWidget(self.progress);layout.addWidget(self.running_text)
        self.summary=QTableWidget();self.membership=QTableWidget();self.membership_filter=QComboBox();self.membership_filter.addItem("All resolved genes","all");self.membership_filter.addItem("MitoCarta only","mito");self.membership_filter.addItem("Non-MitoCarta only","non");membership_page=self._page(self.membership_filter,self.membership)
        self.overall=QTableWidget();self.overall_note=QLabel("Reference represents background \\ target. Values are loaded from the persisted Fisher exact test.");overall_page=self._page(self.overall_note,self.overall)
        self.sub_frequency=QTableWidget();self.sub_choice=QComboBox();self.sub_genes=QTableWidget();sub_frequency_page=self._page(QLabel("Genes may belong to more than one sub-compartment."),self.sub_frequency,QLabel("Select sub-compartment"),self.sub_choice,self.sub_genes)
        self.sub_enrichment=QTableWidget();self.sub_filter=QComboBox();self.sub_filter.addItem("All tested","all");self.sub_filter.addItem("Significant only","significant");self.no_sub_significant=QLabel("No sub-compartments met the selected FDR threshold.");sub_enrichment_page=self._page(self.sub_filter,self.no_sub_significant,self.sub_enrichment)
        self.path_frequency=QTableWidget();self.path_choice=QComboBox();self.path_genes=QTableWidget();path_frequency_page=self._page(self.path_frequency,QLabel("Select MitoPathway"),self.path_choice,self.path_genes)
        self.path_enrichment=QTableWidget();self.path_filter=QComboBox();self.path_filter.addItem("All tested","all");self.path_filter.addItem("Significant only","significant");self.no_path_significant=QLabel("No MitoPathways met the selected FDR threshold.");self.path_excluded=QTableWidget();path_tabs=QTabWidget();path_tabs.addTab(self.path_enrichment,"Tested");path_tabs.addTab(self.path_excluded,"Excluded by minimum overlap");path_enrichment_page=self._page(self.path_filter,self.no_path_significant,path_tabs)
        self.gene_choice=QComboBox();self.gene_compartments=QTableWidget();self.gene_pathways=QTableWidget();gene_tabs=QTabWidget();gene_tabs.addTab(self.gene_compartments,"Sub-compartments");gene_tabs.addTab(self.gene_pathways,"MitoPathways");genes_page=self._page(QLabel("Select a MitoCarta gene from the target"),self.gene_choice,gene_tabs)
        self.mapping=QTableWidget();self.unmapped=QTableWidget();self.ambiguous=QTableWidget()
        self.graph_choice=QComboBox();self.graph_preview=QLabel("No plots are available for this run.");self.graph_preview.setAlignment(Qt.AlignmentFlag.AlignCenter);self.graph_preview.setMinimumHeight(300);graphs_page=self._page(self.graph_choice,self.graph_preview)
        self.history=QComboBox();self.history_status=QLabel();history_page=self._page(QLabel("Previous MitoCarta runs"),self.history,self.history_status)
        self.tabs=QTabWidget();
        for widget,name in ((config,"Configuration"),(self.summary,"Summary"),(membership_page,"Membership"),(overall_page,"Overall enrichment"),(sub_frequency_page,"Sub-compartment frequency"),(sub_enrichment_page,"Sub-compartment enrichment"),(path_frequency_page,"MitoPathway frequency"),(path_enrichment_page,"MitoPathway enrichment"),(genes_page,"Genes"),(self.mapping,"Mapping"),(self.unmapped,"Unmapped"),(self.ambiguous,"Ambiguous"),(graphs_page,"Graphs"),(history_page,"History")):self.tabs.addTab(widget,name)
        self.export_table=QPushButton("Export current table...");self.export_workbook=QPushButton("Export workbook...");self.export_graph=QPushButton("Export graph...");self.open_folder=QPushButton("Open MitoCarta results folder");actions=QHBoxLayout()
        for button in (self.export_table,self.export_workbook,self.export_graph,self.open_folder):button.setEnabled(False);actions.addWidget(button)
        outer=QVBoxLayout(self);outer.addWidget(QLabel("MitoCarta3.0 — Homo sapiens"));outer.addWidget(self.tabs,1);outer.addLayout(actions)
        self.open_database.clicked.connect(self.open_database_requested);self.run_button.clicked.connect(self._emit_run);self.target.currentIndexChanged.connect(self._selection_changed);self.background.currentIndexChanged.connect(self._selection_changed);self.membership_filter.currentIndexChanged.connect(self._show_membership);self.sub_filter.currentIndexChanged.connect(self._show_sub_enrichment);self.path_filter.currentIndexChanged.connect(self._show_path_enrichment);self.gene_choice.currentIndexChanged.connect(self._show_gene_navigation);self.sub_choice.currentIndexChanged.connect(self._show_compartment_genes);self.path_choice.currentIndexChanged.connect(self._show_pathway_genes);self.graph_choice.currentIndexChanged.connect(self._show_graph);self.history.currentIndexChanged.connect(self._load_history);self.export_table.clicked.connect(self._export_current_table);self.export_workbook.clicked.connect(self._export_workbook);self.export_graph.clicked.connect(self._export_graph);self.open_folder.clicked.connect(self._open_results_folder)

    @staticmethod
    def _page(*widgets):
        page=QWidget();layout=QVBoxLayout(page)
        for widget in widgets:layout.addWidget(widget)
        return page
    def _toggle_advanced(self,expanded):
        for child in self.advanced.findChildren(QWidget):child.setVisible(expanded)
    def set_project(self,project):
        self.project=project;self.outputs=None;self._populate_sets();self._populate_manual();self.refresh_readiness();self._refresh_history()
        if project and (project.root/"analyses/MitoCarta/latest_metadata.json").is_file():
            try:self.show_outputs(read_mitocarta_outputs(project))
            except MissingMitoCartaOutputError as error:self.history_status.setText(f"Latest MitoCarta run is incomplete: {error}")
    def refresh_readiness(self):
        snapshot=self.manager.mitocarta.active_snapshot();manifest=self.manager.mitocarta.manifest(snapshot) if snapshot else {};ready,reason=mitocarta_readiness(self.project,self.manager);self.database.setText(f"Status: {'Ready' if snapshot else 'Not installed'}\nVersion: {manifest.get('version','—')}\nSnapshot: {snapshot.name if snapshot else '—'}");self.readiness.setText("MitoCarta3.0 is not installed." if not snapshot else reason);self.run_button.setEnabled(ready)
    def _populate_sets(self):
        self.target.clear();self.background.clear();self.target.addItem("All mapped entities","All mapped entities");self.background.addItem("All gene-resolved entities in the experiment","All mapped entities")
        path=self.project.root/"analyses/presence_absence/tables/classification.csv" if self.project else Path();classes=[]
        if path.is_file():
            try:classes=sorted(pd.read_csv(path).classification.dropna().astype(str).unique())
            except (OSError,ValueError,AttributeError):classes=[]
        dependent=["Reproducibly detected","Shared",*[value for value in classes if value.endswith("-specific")],"Sporadic"]
        for label in dict.fromkeys(dependent):
            self.target.addItem(label,label);self.background.addItem(label,label)
            if not classes:self.target.model().item(self.target.count()-1).setEnabled(False);self.background.model().item(self.background.count()-1).setEnabled(False)
        self.target.addItem("Manual selection","Manual selection");self.background.addItem("Manual selection","Manual selection");self.presence_help.setText("Presence/Absence-derived sets use existing results." if classes else "Run Presence / absence to enable reproducibility, shared, condition-specific, and sporadic sets.")
    def _populate_manual(self):
        self.manual.clear();path=self.project.root/"mapping/tables/protein_catalog.csv" if self.project else Path()
        if not path.is_file():return
        try:frame=pd.read_csv(path).drop_duplicates("source_row")
        except (OSError,ValueError,KeyError):return
        for row in frame.itertuples(index=False):
            source=int(row.source_row);item=QListWidgetItem(f"{getattr(row,'gene_symbol','')} | NCBI {getattr(row,'ncbi_gene_id','')} | UniProt {getattr(row,'uniprot_accession','')} | {getattr(row,'original_id','')}");item.setData(Qt.ItemDataRole.UserRole,source);self.manual.addItem(item)
        self._selection_changed()
    def _selection_changed(self):manual=self.target.currentData()=="Manual selection" or self.background.currentData()=="Manual selection";self.manual.setEnabled(manual);self.manual_help.setEnabled(manual)
    def parameters(self):return {"target_selection":self.target.currentData(),"background_selection":self.background.currentData(),"manual_rows":tuple(int(item.data(Qt.ItemDataRole.UserRole)) for item in self.manual.selectedItems()),"minimum_overlap":self.minimum.value(),"fdr_cutoff":self.fdr.value(),"top_n":self.top_n.value()}
    def _emit_run(self):
        if (self.target.currentData()=="Manual selection" or self.background.currentData()=="Manual selection") and not self.manual.selectedItems():self.readiness.setText("Select at least one experimental entity for Manual selection.");return
        self.run_requested.emit(self.parameters())
    def set_running(self,running):self.progress.setVisible(running);self.running_text.setText("Running MitoCarta analysis..." if running else "");self.run_button.setEnabled(not running and mitocarta_readiness(self.project,self.manager)[0])
    def _fill(self,table,frame):
        self._current_frames[table]=frame.copy();table.clear();table.setRowCount(len(frame));table.setColumnCount(len(frame.columns));table.setHorizontalHeaderLabels([str(column) for column in frame.columns])
        for row,values in enumerate(frame.itertuples(index=False,name=None)):
            for column,value in enumerate(values):table.setItem(row,column,QTableWidgetItem("" if pd.isna(value) else str(value)))
    def _clear_results(self):
        self.outputs=None;self._current_frames.clear()
        for table in (self.summary,self.membership,self.overall,self.sub_frequency,self.sub_enrichment,self.sub_genes,self.path_frequency,self.path_enrichment,self.path_excluded,self.path_genes,self.gene_compartments,self.gene_pathways,self.mapping,self.unmapped,self.ambiguous):self._fill(table,pd.DataFrame())
        self.gene_choice.clear();self.sub_choice.clear();self.path_choice.clear();self.graph_choice.clear();self.graph_preview.setText("No plots are available for this run.")
    def show_outputs(self,outputs):
        self.outputs=outputs;self._fill(self.summary,outputs.summary);self._fill(self.overall,outputs.overall);self._fill(self.sub_frequency,outputs.subcompartment_frequency);self._fill(self.path_frequency,outputs.pathway_frequency);self._fill(self.mapping,outputs.mapping);self._fill(self.unmapped,outputs.unmapped);self._fill(self.ambiguous,outputs.ambiguous);self._fill(self.path_excluded,outputs.pathway_excluded);self._show_membership();self._show_sub_enrichment();self._show_path_enrichment();self.no_sub_significant.setVisible(outputs.subcompartment_significant.empty);self.no_path_significant.setVisible(outputs.pathway_significant.empty)
        target_genes=mitocarta_target_genes(outputs);self.gene_choice.blockSignals(True);self.gene_choice.clear()
        for row in target_genes.itertuples(index=False):self.gene_choice.addItem(f"{row.gene_symbol} — {row.canonical_gene_key}",str(row.canonical_gene_key))
        self.gene_choice.blockSignals(False);self.sub_choice.clear();[self.sub_choice.addItem(str(value),str(value)) for value in outputs.subcompartment_frequency.get("Subcompartment",pd.Series(dtype=str))];self.path_choice.clear();[self.path_choice.addItem(str(value),str(value)) for value in outputs.pathway_frequency.get("MitoPathway",pd.Series(dtype=str))]
        self.graph_choice.blockSignals(True);self.graph_choice.clear();pngs={path.stem:path for path in outputs.graphs if path.suffix.lower()==".png"};pdfs={path.stem:path for path in outputs.graphs if path.suffix.lower()==".pdf"}
        for stem,path in sorted(pngs.items()):self.graph_choice.addItem(stem.replace("_"," ").title(),{"png":path,"pdf":pdfs.get(stem)})
        self.graph_choice.blockSignals(False);self._show_gene_navigation();self._show_compartment_genes();self._show_pathway_genes();self._show_graph();self.history_status.setText(f"Loaded MitoCarta run: {outputs.metadata.get('run_id','unknown')}");[button.setEnabled(True) for button in (self.export_table,self.export_workbook,self.open_folder)];self.export_graph.setEnabled(self.graph_choice.count()>0);self._refresh_history(outputs.metadata.get("run_id"))
    def _show_membership(self):
        if not self.outputs:return self._fill(self.membership,pd.DataFrame())
        frame=self.outputs.membership.copy();flag=frame.is_mitocarta.astype(str).str.lower().isin({"true","1","t"});mode=self.membership_filter.currentData();self._fill(self.membership,frame if mode=="all" else frame[flag if mode=="mito" else ~flag])
    def _show_sub_enrichment(self):self._fill(self.sub_enrichment,pd.DataFrame() if not self.outputs else (self.outputs.subcompartment_significant if self.sub_filter.currentData()=="significant" else self.outputs.subcompartment_enrichment))
    def _show_path_enrichment(self):self._fill(self.path_enrichment,pd.DataFrame() if not self.outputs else (self.outputs.pathway_significant if self.path_filter.currentData()=="significant" else self.outputs.pathway_enrichment))
    def _show_gene_navigation(self):
        key=self.gene_choice.currentData();self._fill(self.gene_compartments,mitocarta_gene_subcompartments(self.outputs,key) if self.outputs and key else pd.DataFrame());self._fill(self.gene_pathways,mitocarta_gene_pathways(self.outputs,key) if self.outputs and key else pd.DataFrame())
    def _show_compartment_genes(self):self._fill(self.sub_genes,mitocarta_compartment_genes(self.outputs,self.sub_choice.currentData()) if self.outputs and self.sub_choice.currentData() else pd.DataFrame())
    def _show_pathway_genes(self):self._fill(self.path_genes,mitocarta_pathway_genes(self.outputs,self.path_choice.currentData()) if self.outputs and self.path_choice.currentData() else pd.DataFrame())
    def _show_graph(self):
        data=self.graph_choice.currentData()
        if not data or not Path(data.get("png","")).is_file():self.graph_preview.clear();self.graph_preview.setText("No plots are available for this run.");return
        self.graph_preview.setPixmap(QPixmap(str(data["png"])).scaled(900,600,Qt.AspectRatioMode.KeepAspectRatio,Qt.TransformationMode.SmoothTransformation))
    def _refresh_history(self,selected=None):
        self.history.blockSignals(True);self.history.clear()
        if self.project:
            for run_id in list_mitocarta_runs(self.project):
                path=self.project.root/"analyses/MitoCarta/runs"/run_id/"metadata.json";label=run_id
                try:
                    meta=json.loads(path.read_text(encoding="utf-8"));label=f"{run_id} | {meta.get('target_definition','Unknown target')} | {meta.get('background_definition','Unknown background')} | MitoCarta {meta.get('mitocarta_version','Unknown')} | Snapshot {meta.get('snapshot_id','Unknown')}"
                except (OSError,ValueError):label=f"{run_id} | Incomplete metadata"
                self.history.addItem(label,run_id)
        self.history.setCurrentIndex(self.history.findData(selected) if selected else -1);self.history.blockSignals(False)
    def _load_history(self,index):
        if not self.project or index<0:return
        try:self.show_outputs(read_mitocarta_outputs(self.project,self.history.itemData(index)))
        except MissingMitoCartaOutputError as error:self._clear_results();self.history_status.setText(f"Historical run is incomplete: {error}")
    def _export_current_table(self):
        widget=self.tabs.currentWidget();table=widget if isinstance(widget,QTableWidget) else {2:self.membership,3:self.overall,4:self.sub_frequency,5:self.sub_enrichment,6:self.path_frequency,7:self.path_enrichment}.get(self.tabs.currentIndex());frame=self._current_frames.get(table) if table else None
        if frame is None:return QMessageBox.information(self,"Export table","The current panel does not contain an exportable table.")
        filename,_=QFileDialog.getSaveFileName(self,"Export current table","mitocarta_table.csv","CSV (*.csv)")
        if filename:frame.to_csv(filename,index=False)
    def _export_workbook(self):
        if not self.outputs:return
        filename,_=QFileDialog.getSaveFileName(self,"Export workbook","MitoCarta_analysis.xlsx","Excel workbook (*.xlsx)")
        if filename:export_mitocarta_artifact(self.outputs.workbook,Path(filename))
    def _export_graph(self):
        data=self.graph_choice.currentData()
        if not data:return
        filename,selected=QFileDialog.getSaveFileName(self,"Export graph","mitocarta_plot.png","PNG (*.png);;PDF (*.pdf)")
        if filename:
            suffix=".pdf" if "PDF" in selected or Path(filename).suffix.lower()==".pdf" else ".png";source=data.get(suffix[1:]);destination=Path(filename) if Path(filename).suffix else Path(filename).with_suffix(suffix)
            if source:export_mitocarta_artifact(source,destination)
    def _open_results_folder(self):
        if self.outputs and not QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.outputs.run_root))):QMessageBox.warning(self,"MitoCarta results","Could not open the MitoCarta results folder.")
