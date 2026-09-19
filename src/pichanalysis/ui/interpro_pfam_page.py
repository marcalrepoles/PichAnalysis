from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Callable

import pandas as pd
from PySide6.QtCore import QThread, Qt, QUrl, Signal
from PySide6.QtGui import QDesktopServices, QPixmap
from PySide6.QtWidgets import (QAbstractItemView,QComboBox,QDoubleSpinBox,QFileDialog,QFormLayout,QGroupBox,QHBoxLayout,QLabel,QListWidget,QListWidgetItem,QMessageBox,QProgressBar,QPushButton,QSpinBox,QTabWidget,QTableWidget,QTableWidgetItem,QVBoxLayout,QWidget)

from ..core.database_manager import DatabaseManager
from ..core.interpro_database import InterProCancelled
from ..core.interpro_pfam_analysis import (InterProPfamOutputs,InterProPfamParameters,InterProPfamTargetOutsideBackgroundError,architectures_for_protein,interpro_for_protein,interpro_locations_for_protein,interpro_pfam_experiment_accessions,list_interpro_pfam_runs,pfam_for_protein,pfam_locations_for_protein,proteins_for_interpro,proteins_for_pfam,read_interpro_pfam_outputs,run_interpro_pfam_analysis,select_interpro_pfam_set,_resolved_accessions)
from ..core.project import Project
from ..core.r_runtime import RRuntime


class InterProPfamAnalysisWorker(QThread):
    succeeded=Signal(object);target_outside_background=Signal(dict);failed=Signal(str)
    def __init__(self,project,manager,runtime,run_id,parameters,runner:Callable=run_interpro_pfam_analysis):super().__init__();self.project=project;self.manager=manager;self.runtime=runtime;self.run_id=run_id;self.parameters=parameters;self.runner=runner
    def run(self):
        try:self.succeeded.emit(self.runner(self.project,self.manager,self.runtime,run_id=self.run_id,parameters=self.parameters))
        except InterProPfamTargetOutsideBackgroundError as error:self.target_outside_background.emit(error.details)
        except Exception as error:self.failed.emit(str(error))


class AnnotationSetBuildWorker(QThread):
    progress=Signal(dict);succeeded=Signal(str);failed=Signal(str);canceled=Signal(str)
    def __init__(self,database,project_path,accessions):super().__init__();import threading;self.database=database;self.project_path=project_path;self.accessions=tuple(accessions);self._cancel=threading.Event()
    def run(self):
        try:self.succeeded.emit(str(self.database.build_annotation_set(self.project_path,self.accessions,progress=self.progress.emit,cancel_requested=self._cancel.is_set)))
        except InterProCancelled as error:self.canceled.emit(str(error))
        except Exception as error:self.failed.emit(str(error))
    def cancel(self):self._cancel.set()


class InterProPfamPage(QWidget):
    run_requested=Signal(dict);open_database_requested=Signal()
    def __init__(self,manager:DatabaseManager):
        super().__init__();self.manager=manager;self.project=None;self.outputs=None;self.build_worker=None;self._frames={}
        self.metadata=QLabel();self.metadata.setWordWrap(True);self.open_database=QPushButton("Open Database Manager")
        self.annotation_set=QComboBox();self.set_status=QLabel();self.set_status.setWordWrap(True);self.build_set=QPushButton("Build annotation set");self.cancel_build=QPushButton("Cancel");self.cancel_build.setEnabled(False);self.build_progress=QProgressBar();self.build_progress.setRange(0,1);self.build_text=QLabel("No annotation acquisition is running.");self.open_set=QPushButton("Open annotation set folder")
        self.target=QComboBox();self.background=QComboBox();self.background_help=QLabel("Background represents the proteins that could have been observed in this experiment.");self.manual=QListWidget();self.manual.setSelectionMode(QAbstractItemView.SelectionMode.MultiSelection);self.presence_help=QLabel()
        defaults=InterProPfamParameters();self.advanced=QGroupBox("Advanced options");self.advanced.setCheckable(True);self.advanced.setChecked(False);self.fdr=QDoubleSpinBox();self.fdr.setRange(.0001,1);self.fdr.setDecimals(4);self.fdr.setValue(defaults.fdr_cutoff);self.minimum=QSpinBox();self.minimum.setRange(1,100000);self.minimum.setValue(defaults.minimum_overlap);self.top_n=QSpinBox();self.top_n.setRange(1,1000);self.top_n.setValue(defaults.top_n);advanced=QFormLayout(self.advanced);advanced.addRow("FDR threshold",self.fdr);advanced.addRow("Minimum overlap",self.minimum);advanced.addRow("Top N features in plots",self.top_n);self.advanced.toggled.connect(self._toggle_advanced);self._toggle_advanced(False)
        help_text=QLabel("Frequency: Shows how many proteins in the selected target contain each annotation.\nEnrichment: Tests whether an annotation is more frequent in the selected target than in the remaining experimental background.\nLocations: A protein may contain the same domain more than once. Repeated locations do not increase membership counts.\nArchitecture: Shows the ordered arrangement of structural annotations along each protein.");help_text.setWordWrap(True)
        self.readiness=QLabel();self.readiness.setWordWrap(True);self.run_button=QPushButton("5. Run InterPro/Pfam analysis");self.progress=QProgressBar();self.progress.setRange(0,0);self.progress.hide();self.running_text=QLabel()
        config=QWidget();layout=QVBoxLayout(config);layout.addWidget(QLabel("1. Annotation data"));layout.addWidget(self.metadata);layout.addWidget(self.open_database);layout.addWidget(self.annotation_set);layout.addWidget(self.set_status);row=QHBoxLayout();[row.addWidget(x) for x in (self.build_set,self.cancel_build,self.open_set)];layout.addLayout(row);layout.addWidget(self.build_progress);layout.addWidget(self.build_text);form=QFormLayout();form.addRow("2. Target set",self.target);form.addRow("Manual selection",self.manual);form.addRow("3. Background",self.background);layout.addLayout(form);layout.addWidget(self.presence_help);layout.addWidget(self.background_help);layout.addWidget(QLabel("4. Analysis settings"));layout.addWidget(self.advanced);layout.addWidget(help_text);layout.addWidget(self.readiness);layout.addWidget(self.run_button);layout.addWidget(self.progress);layout.addWidget(self.running_text)
        self.summary=QTableWidget();self.coverage=QTableWidget();self.interpro_frequency=QTableWidget();self.interpro_type=QComboBox();self.interpro_enrichment=QTableWidget();self.interpro_enrichment_filter=QComboBox();self.interpro_enrichment_filter.addItems(["All tested","Significant only"]);self.interpro_excluded=QTableWidget();self.no_interpro=QLabel("No InterPro entries met the selected FDR threshold.")
        self.pfam_frequency=QTableWidget();self.pfam_enrichment=QTableWidget();self.pfam_enrichment_filter=QComboBox();self.pfam_enrichment_filter.addItems(["All tested","Significant only"]);self.pfam_excluded=QTableWidget();self.no_pfam=QLabel("No Pfam entries met the selected FDR threshold.")
        self.protein=QComboBox();self.protein_interpro=QTableWidget();self.protein_pfam=QTableWidget();self.protein_interpro_locations=QTableWidget();self.protein_pfam_locations=QTableWidget();protein_tabs=QTabWidget();[protein_tabs.addTab(table,name) for table,name in ((self.protein_interpro,"InterPro annotations"),(self.protein_interpro_locations,"InterPro locations"),(self.protein_pfam,"Pfam annotations"),(self.protein_pfam_locations,"Pfam locations"))];proteins_page=self._page(self.protein,protein_tabs)
        self.interpro_feature=QComboBox();self.interpro_proteins=QTableWidget();self.pfam_feature=QComboBox();self.pfam_proteins=QTableWidget();feature_page=self._page(QLabel("InterPro feature → proteins"),self.interpro_feature,self.interpro_proteins,QLabel("Pfam feature → proteins"),self.pfam_feature,self.pfam_proteins)
        self.repeated=QTableWidget();repeated_page=self._page(QLabel("Repeated occurrences are shown separately as locations but count once for frequency and enrichment."),self.repeated)
        self.pfam_architecture=QTableWidget();self.interpro_architecture=QTableWidget();self.pfam_arch_frequency=QTableWidget();self.interpro_arch_frequency=QTableWidget();self.integration=QTableWidget();arch_tabs=QTabWidget();[arch_tabs.addTab(table,name) for table,name in ((self.pfam_architecture,"Pfam architecture"),(self.interpro_architecture,"InterPro structural architecture"),(self.pfam_arch_frequency,"Pfam architecture frequency"),(self.interpro_arch_frequency,"InterPro architecture frequency"),(self.integration,"Pfam-InterPro integration"))];architecture_page=self._page(QLabel("Structural architecture includes InterPro domain and repeat entries. Overlapping annotations remain visible."),arch_tabs)
        self.mapping=QTableWidget();self.unmapped=QTableWidget();self.ambiguous=QTableWidget();self.graph_choice=QComboBox();self.graph_preview=QLabel("No plots are available for this run.");self.graph_preview.setAlignment(Qt.AlignmentFlag.AlignCenter);self.graph_preview.setMinimumHeight(300);graphs=self._page(self.graph_choice,self.graph_preview);self.history=QComboBox();self.history_status=QLabel();history_page=self._page(self.history,self.history_status)
        ip_freq=self._page(self.interpro_type,self.interpro_frequency);ip_enr_tabs=QTabWidget();ip_enr_tabs.addTab(self.interpro_enrichment,"Tested");ip_enr_tabs.addTab(self.interpro_excluded,"Excluded by minimum overlap");ip_enr=self._page(self.interpro_enrichment_filter,self.no_interpro,ip_enr_tabs);ip_tabs=QTabWidget();ip_tabs.addTab(ip_freq,"Frequency");ip_tabs.addTab(ip_enr,"Enrichment")
        pf_enr_tabs=QTabWidget();pf_enr_tabs.addTab(self.pfam_enrichment,"Tested");pf_enr_tabs.addTab(self.pfam_excluded,"Excluded by minimum overlap");pf_tabs=QTabWidget();pf_tabs.addTab(self.pfam_frequency,"Frequency");pf_tabs.addTab(self._page(self.pfam_enrichment_filter,self.no_pfam,pf_enr_tabs),"Enrichment")
        self.tabs=QTabWidget();
        for widget,name in ((config,"Configuration"),(self.summary,"Summary"),(self.coverage,"Annotation coverage"),(ip_tabs,"InterPro"),(pf_tabs,"Pfam"),(proteins_page,"Proteins"),(feature_page,"Feature navigation"),(architecture_page,"Architecture"),(repeated_page,"Repeated features"),(self.mapping,"Mapping"),(self.unmapped,"Unmapped"),(self.ambiguous,"Ambiguous"),(graphs,"Graphs"),(history_page,"History")):self.tabs.addTab(widget,name)
        self.export_table=QPushButton("Export current table...");self.export_workbook=QPushButton("Export workbook...");self.export_graph=QPushButton("Export graph...");self.open_results=QPushButton("Open InterPro/Pfam results folder");actions=QHBoxLayout();[actions.addWidget(x) for x in (self.export_table,self.export_workbook,self.export_graph,self.open_results)];outer=QVBoxLayout(self);outer.addWidget(QLabel("InterPro / Pfam"));outer.addWidget(self.tabs,1);outer.addLayout(actions)
        self.open_database.clicked.connect(self.open_database_requested);self.build_set.clicked.connect(self._build_annotation_set);self.cancel_build.clicked.connect(self._cancel_build);self.open_set.clicked.connect(self._open_set_folder);self.run_button.clicked.connect(lambda:self.run_requested.emit(self.parameters()));self.annotation_set.currentIndexChanged.connect(self.refresh_readiness);self.target.currentIndexChanged.connect(self._selection_changed);self.background.currentIndexChanged.connect(self.refresh_readiness);self.interpro_type.currentIndexChanged.connect(self._show_interpro_frequency);self.interpro_enrichment_filter.currentIndexChanged.connect(self._show_interpro_enrichment);self.pfam_enrichment_filter.currentIndexChanged.connect(self._show_pfam_enrichment);self.protein.currentIndexChanged.connect(self._show_protein);self.interpro_feature.currentIndexChanged.connect(self._show_feature);self.pfam_feature.currentIndexChanged.connect(self._show_feature);self.graph_choice.currentIndexChanged.connect(self._show_graph);self.history.currentIndexChanged.connect(self._load_history);self.export_table.clicked.connect(self._export_current);self.export_workbook.clicked.connect(self._export_workbook);self.export_graph.clicked.connect(self._export_graph);self.open_results.clicked.connect(self._open_results)
    @staticmethod
    def _page(*widgets):page=QWidget();layout=QVBoxLayout(page);[layout.addWidget(x) for x in widgets];return page
    def _toggle_advanced(self,value):
        for child in self.advanced.findChildren(QWidget):child.setVisible(value)
    def set_project(self,project):self.project=project;self.outputs=None;self._populate_sets();self._populate_targets();self._populate_manual();self.refresh_readiness();self._refresh_history()
    def _populate_sets(self):
        self.annotation_set.clear()
        if self.project:
            for sid in self.manager.interpro.list_annotation_sets(self.project.root):
                if self.manager.interpro.is_annotation_set_ready(self.project.root,sid):self.annotation_set.addItem(sid,sid)
    def _populate_targets(self):
        self.target.clear();self.background.clear();self.target.addItem("All mapped proteins","All mapped proteins");self.background.addItem("All annotation-set-covered proteins in the experiment","All mapped proteins");path=self.project.root/"analyses/presence_absence/tables/classification.csv" if self.project else Path();classes=[]
        if path.is_file():
            try:classes=sorted(pd.read_csv(path).classification.dropna().astype(str).unique())
            except Exception:classes=[]
        for label in dict.fromkeys(["Reproducibly detected","Shared",*[x for x in classes if x.endswith("-specific")],"Sporadic"]):
            self.target.addItem(label,label);self.background.addItem(label,label);enabled=bool(classes);self.target.model().item(self.target.count()-1).setEnabled(enabled);self.background.model().item(self.background.count()-1).setEnabled(enabled)
        self.target.addItem("Manual selection","Manual selection");self.background.addItem("Manual selection","Manual selection");self.presence_help.setText("Presence/Absence-derived sets use existing results." if classes else "Run Presence / absence to enable derived protein sets.")
    def _populate_manual(self):
        self.manual.clear();path=self.project.root/"mapping/tables/protein_catalog.csv" if self.project else Path()
        if not path.is_file():return
        frame=pd.read_csv(path).drop_duplicates("source_row")
        for row in frame.itertuples(index=False):item=QListWidgetItem(f"{getattr(row,'uniprot_accession','')} | {getattr(row,'gene_symbol','')} | {getattr(row,'protein_name','')} | {getattr(row,'original_id','')}");item.setData(Qt.ItemDataRole.UserRole,int(row.source_row));self.manual.addItem(item)
    def _selection_changed(self):self.manual.setEnabled(self.target.currentData()=="Manual selection" or self.background.currentData()=="Manual selection");self.refresh_readiness()
    def refresh_readiness(self):
        manifest=self.manager.interpro.metadata_manifest();self.metadata.setText(f"InterPro release: {manifest.get('interpro_release','Unknown')}\nPfam release: {manifest.get('pfam_release','Unknown')}\nMetadata status: {'Ready' if manifest else 'Not installed'}")
        if not manifest:self.readiness.setText("InterPro/Pfam metadata is not installed.");self.run_button.setEnabled(False);return
        sid=self.annotation_set.currentData()
        if not self.project or not sid:self.set_status.setText("Annotation set status: Unavailable");self.readiness.setText("Build or select a Ready annotation set.");self.run_button.setEnabled(False);return
        if self.background.currentData() is None:self.run_button.setEnabled(False);return
        annotation=self.manager.interpro.load_annotation_set(self.project.root,sid);required=set(_resolved_accessions(select_interpro_pfam_set(self.project,self.background.currentData())));covered=set(annotation.proteins.requested_accession.astype(str));missing=required-covered;self.set_status.setText(f"Annotation set: {sid}\nStatus: Ready\nProteins: {len(annotation.proteins)}\nInterPro release: {annotation.manifest.get('interpro_release')}\nPfam release: {annotation.manifest.get('pfam_release')}")
        self.readiness.setText("The selected annotation set does not cover all required proteins." if missing else "Ready for InterPro/Pfam analysis.");self.run_button.setEnabled(not missing)
    def parameters(self):return {"annotation_set_id":self.annotation_set.currentData(),"target_selection":self.target.currentData(),"background_selection":self.background.currentData(),"manual_rows":tuple(int(x.data(Qt.ItemDataRole.UserRole)) for x in self.manual.selectedItems()),"minimum_overlap":self.minimum.value(),"fdr_cutoff":self.fdr.value(),"top_n":self.top_n.value()}
    def _build_annotation_set(self):
        if not self.project:return
        accessions=interpro_pfam_experiment_accessions(self.project);answer=QMessageBox.question(self,"Build annotation set",f"Proteins to annotate: {len(accessions)}\n\nBuild a new immutable annotation set using the official InterPro service?")
        if answer!=QMessageBox.StandardButton.Yes:return
        self.build_worker=AnnotationSetBuildWorker(self.manager.interpro,self.project.root,accessions);self.build_worker.progress.connect(self._build_progress);self.build_worker.succeeded.connect(self._build_finished);self.build_worker.failed.connect(self._build_failed);self.build_worker.canceled.connect(self._build_canceled);self.build_worker.finished.connect(self.build_worker.deleteLater);self.build_worker.start();self.build_set.setEnabled(False);self.cancel_build.setEnabled(True)
    def _build_progress(self,p):self.build_progress.setRange(0,max(1,int(p.get("total",1))));self.build_progress.setValue(int(p.get("done",0)));self.build_text.setText(str(p.get("message","")))
    def _build_finished(self,path):self.build_worker=None;self.build_set.setEnabled(True);self.cancel_build.setEnabled(False);self.build_text.setText("Annotation set is Ready.");self._populate_sets();self.annotation_set.setCurrentIndex(self.annotation_set.findData(Path(path).name));self.refresh_readiness()
    def _build_failed(self,message):self.build_worker=None;self.build_set.setEnabled(True);self.cancel_build.setEnabled(False);self.build_text.setText("Acquisition failed.");QMessageBox.warning(self,"Annotation acquisition failed",message)
    def _build_canceled(self,message):self.build_worker=None;self.build_set.setEnabled(True);self.cancel_build.setEnabled(False);self.build_text.setText(message)
    def _cancel_build(self):
        if self.build_worker:self.build_text.setText("Cancellation requested.");self.build_worker.cancel()
    def is_running(self):return bool(self.build_worker and self.build_worker.isRunning())
    def set_running(self,running):self.progress.setVisible(running);self.running_text.setText("Running InterPro/Pfam analysis..." if running else "");self.run_button.setEnabled(not running and "Ready for" in self.readiness.text())
    def _fill(self,table,frame):
        self._frames[table]=frame.copy();table.clear();table.setRowCount(len(frame));table.setColumnCount(len(frame.columns));table.setHorizontalHeaderLabels([str(x).replace("_"," ") for x in frame.columns])
        for r,values in enumerate(frame.itertuples(index=False,name=None)):
            for c,value in enumerate(values):table.setItem(r,c,QTableWidgetItem("" if pd.isna(value) else str(value)))
    def show_outputs(self,o:InterProPfamOutputs):
        self.outputs=o;self._fill(self.summary,o.summary)
        try:coverage=self.manager.interpro.load_annotation_set(self.project.root,o.metadata.get("annotation_set_id")).proteins
        except Exception:coverage=o.coverage
        self._fill(self.coverage,coverage);self._fill(self.interpro_excluded,o.interpro_excluded);self._fill(self.pfam_frequency,o.pfam_frequency);self._fill(self.pfam_excluded,o.pfam_excluded);self._fill(self.repeated,o.repeated_features);self._fill(self.pfam_architecture,o.pfam_architecture);self._fill(self.interpro_architecture,o.interpro_architecture);self._fill(self.pfam_arch_frequency,o.pfam_architecture_frequency);self._fill(self.interpro_arch_frequency,o.interpro_architecture_frequency);self._fill(self.integration,o.integration);self._fill(self.mapping,o.mapping);self._fill(self.unmapped,o.unmapped);self._fill(self.ambiguous,o.ambiguous)
        self.interpro_type.blockSignals(True);self.interpro_type.clear();self.interpro_type.addItem("All entry types",None);[self.interpro_type.addItem(str(x).replace("_"," ").title(),str(x)) for x in sorted(o.interpro_frequency.Entry_type.dropna().astype(str).unique())];self.interpro_type.blockSignals(False);self._show_interpro_frequency();self._show_interpro_enrichment();self._show_pfam_enrichment();self.no_interpro.setVisible(o.interpro_significant.empty);self.no_pfam.setVisible(o.pfam_significant.empty)
        proteins=sorted(set(o.mapping.loc[o.mapping.mapping_status=="mapped_unique","uniprot_accession"].dropna().astype(str)));self.protein.clear();[self.protein.addItem(x,f"UP:{x}") for x in proteins];self.interpro_feature.clear();[self.interpro_feature.addItem(x,x) for x in o.interpro_frequency.InterPro_ID.astype(str)];self.pfam_feature.clear();[self.pfam_feature.addItem(x,x) for x in o.pfam_frequency.Pfam_ID.astype(str)];self._show_protein();self._show_feature()
        self.graph_choice.clear();pngs={p.stem:p for p in o.plots if p.suffix.lower()==".png"};pdfs={p.stem:p for p in o.plots if p.suffix.lower()==".pdf"};[self.graph_choice.addItem(stem.replace("_"," ").title(),{"png":path,"pdf":pdfs.get(stem)}) for stem,path in sorted(pngs.items())];self._show_graph();self.history_status.setText(f"Loaded run {o.metadata.get('run_id')} using annotation set {o.metadata.get('annotation_set_id')}");self._refresh_history(o.metadata.get("run_id"))
    def _show_interpro_frequency(self):
        if not self.outputs:return
        frame=self.outputs.interpro_frequency;kind=self.interpro_type.currentData();self._fill(self.interpro_frequency,frame if kind is None else frame[frame.Entry_type.astype(str)==kind])
    def _show_interpro_enrichment(self):self._fill(self.interpro_enrichment,pd.DataFrame() if not self.outputs else (self.outputs.interpro_significant if self.interpro_enrichment_filter.currentText()=="Significant only" else self.outputs.interpro_enrichment))
    def _show_pfam_enrichment(self):self._fill(self.pfam_enrichment,pd.DataFrame() if not self.outputs else (self.outputs.pfam_significant if self.pfam_enrichment_filter.currentText()=="Significant only" else self.outputs.pfam_enrichment))
    def _show_protein(self):
        if not self.outputs:return
        accession=str(self.protein.currentData() or "").removeprefix("UP:");self._fill(self.protein_interpro,interpro_for_protein(self.outputs,accession));self._fill(self.protein_pfam,pfam_for_protein(self.outputs,accession));self._fill(self.protein_interpro_locations,interpro_locations_for_protein(self.outputs,accession));self._fill(self.protein_pfam_locations,pfam_locations_for_protein(self.outputs,accession))
    def _show_feature(self):
        if not self.outputs:return
        self._fill(self.interpro_proteins,proteins_for_interpro(self.outputs,self.interpro_feature.currentData()));self._fill(self.pfam_proteins,proteins_for_pfam(self.outputs,self.pfam_feature.currentData()))
    def _show_graph(self):
        data=self.graph_choice.currentData()
        if not data or not Path(data["png"]).is_file():self.graph_preview.setText("No plots are available for this run.");return
        self.graph_preview.setPixmap(QPixmap(str(data["png"])).scaled(900,600,Qt.AspectRatioMode.KeepAspectRatio,Qt.TransformationMode.SmoothTransformation))
    def _refresh_history(self,selected=None):
        self.history.blockSignals(True);self.history.clear()
        if self.project:
            for rid in list_interpro_pfam_runs(self.project):
                try:m=json.loads((self.project.root/"analyses/InterPro_Pfam/runs"/rid/"metadata.json").read_text());label=f"{rid} | {m.get('target_definition')} | {m.get('background_definition')} | Set {m.get('annotation_set_id')} | InterPro {m.get('interpro_release')} | Pfam {m.get('pfam_release')}"
                except Exception:label=f"{rid} | Incomplete metadata"
                self.history.addItem(label,rid)
        self.history.setCurrentIndex(self.history.findData(selected) if selected else -1);self.history.blockSignals(False)
    def _load_history(self,index):
        if not self.project or index<0:return
        try:self.show_outputs(read_interpro_pfam_outputs(self.project,self.history.itemData(index)))
        except Exception as error:self.outputs=None;self.history_status.setText(f"Historical run is incomplete: {error}")
    def _export_current(self):
        tables=[x for x in self.tabs.currentWidget().findChildren(QTableWidget)] if not isinstance(self.tabs.currentWidget(),QTableWidget) else [self.tabs.currentWidget()];table=tables[0] if tables else None;frame=self._frames.get(table)
        if frame is None:return QMessageBox.information(self,"Export table","The current panel does not contain an exportable table.")
        path,_=QFileDialog.getSaveFileName(self,"Export current table","interpro_pfam_table.csv","CSV (*.csv)");path and frame.to_csv(path,index=False)
    def _export_workbook(self):
        if not self.outputs:return
        path,_=QFileDialog.getSaveFileName(self,"Export workbook","InterPro_Pfam_analysis.xlsx","Excel workbook (*.xlsx)");path and shutil.copy2(self.outputs.workbook,path)
    def _export_graph(self):
        data=self.graph_choice.currentData()
        if not data:return
        path,selected=QFileDialog.getSaveFileName(self,"Export graph","interpro_pfam_plot.png","PNG (*.png);;PDF (*.pdf)")
        if path:
            source=data.get("pdf" if "PDF" in selected or Path(path).suffix.lower()==".pdf" else "png");source and shutil.copy2(source,path)
    def _open_results(self):
        if self.outputs:QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.outputs.run_root)))
    def _open_set_folder(self):
        sid=self.annotation_set.currentData()
        if self.project and sid:QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.project.root/"annotations/InterPro/sets"/sid)))
