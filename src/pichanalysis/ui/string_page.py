from __future__ import annotations

import json
import shutil
from pathlib import Path

import pandas as pd
from PySide6.QtCore import QThread, Qt, QUrl, Signal
from PySide6.QtGui import QDesktopServices, QPixmap
from PySide6.QtWidgets import (
    QAbstractItemView, QComboBox, QFileDialog, QFormLayout, QGroupBox, QHBoxLayout,
    QLabel, QListWidget, QListWidgetItem, QMessageBox, QProgressBar, QPushButton,
    QSpinBox, QTabWidget, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget,
)

from ..core.database_manager import DatabaseManager
from ..core.string_analysis import (
    StringExpansionLimitError, StringParameters, StringTargetOutsideBackgroundError,
    list_string_runs, map_experiment, read_string_outputs, run_string_analysis,
    string_edge_details, string_node_neighbors, string_node_support,
)
from .string_viewer import StringNetworkViewer


class StringAnalysisWorker(QThread):
    succeeded = Signal(object); failed = Signal(str); expansion_limit = Signal(dict); target_outside = Signal(dict)
    def __init__(self, project, manager, runtime, run_id, parameters, runner=run_string_analysis):
        super().__init__();self.project=project;self.manager=manager;self.runtime=runtime
        self.run_id=run_id;self.parameters=parameters;self.runner=runner
    def run(self):
        try:self.succeeded.emit(self.runner(self.project,self.manager,self.runtime,run_id=self.run_id,parameters=self.parameters))
        except StringExpansionLimitError as error:self.expansion_limit.emit(error.details)
        except StringTargetOutsideBackgroundError as error:self.target_outside.emit(error.details)
        except Exception as error:self.failed.emit(str(error))


class StringPage(QWidget):
    run_requested = Signal(dict); open_database_requested = Signal()
    def __init__(self, manager:DatabaseManager):
        super().__init__();self.manager=manager;self.project=None;self.outputs=None;self.mapping_frame=pd.DataFrame();self._frames={};self._running=False
        defaults=StringParameters()
        self.database=QLabel();self.database.setWordWrap(True)
        self.open_database=QPushButton("Open Database Manager")
        self.readiness=QLabel();self.readiness.setWordWrap(True)
        self.target=QComboBox();self.manual=QListWidget();self.manual.setSelectionMode(QAbstractItemView.SelectionMode.MultiSelection)
        self.manual.setMinimumHeight(130);self.seed_count=QLabel("Selected seeds: 0");self.presence_help=QLabel()
        self.network_type=QComboBox();self.network_type.addItem("Functional association","functional");self.network_type.addItem("Physical","physical")
        self.network_help=QLabel("Functional association includes functional relationships and should not be interpreted automatically as direct physical interaction.");self.network_help.setWordWrap(True)
        self.threshold_preset=QComboBox()
        for label,value in (("Exploratory — 150","Exploratory"),("Functional — 400","Functional"),("Robust — 700","Robust"),("Custom","Custom")):self.threshold_preset.addItem(label,value)
        self.threshold_preset.setCurrentIndex(1);self.custom_threshold=QSpinBox();self.custom_threshold.setRange(0,1000);self.custom_threshold.setValue(defaults.threshold);self.custom_threshold.setEnabled(False)
        self.threshold_value=QLabel("Confidence threshold: 400")
        self.max_hop=QComboBox();self.max_hop.addItem("0 — Seeds only",0);self.max_hop.addItem("1 — Direct neighbors",1);self.max_hop.addItem("2 — Second-degree neighbors",2);self.max_hop.setCurrentIndex(1)
        self.selection_mode=QComboBox();self.selection_mode.addItem("Union","union");self.selection_mode.addItem("Strict common direct neighbors","strict_common")
        self.selection_help=QLabel();self.selection_help.setWordWrap(True)
        self.advanced=QGroupBox("Advanced options");self.advanced.setCheckable(True);self.advanced.setChecked(False)
        self.minimum_hub_degree=QSpinBox();self.minimum_hub_degree.setRange(0,100000);self.minimum_hub_degree.setValue(defaults.minimum_hub_degree)
        self.max_external_nodes=QSpinBox();self.max_external_nodes.setRange(0,10000000);self.max_external_nodes.setValue(defaults.max_external_nodes)
        self.top_n=QSpinBox();self.top_n.setRange(1,1000);self.top_n.setValue(defaults.top_n)
        advanced_form=QFormLayout(self.advanced);advanced_form.addRow("Minimum hub degree",self.minimum_hub_degree);advanced_form.addRow("Maximum external nodes",self.max_external_nodes);advanced_form.addRow("Top N metrics in plots",self.top_n)
        self.guard_help=QLabel("Maximum external nodes is an operational safety limit. It does not rank or truncate the network.");self.guard_help.setWordWrap(True)
        self.run_button=QPushButton("6. Run STRING analysis");self.progress=QProgressBar();self.progress.setRange(0,0);self.progress.hide();self.running_text=QLabel()
        config=QWidget();configuration=QVBoxLayout(config)
        configuration.addWidget(QLabel("1. Database"));configuration.addWidget(self.database);configuration.addWidget(self.open_database)
        form=QFormLayout();form.addRow("2. Seeds / target",self.target);form.addRow("Manual seed selection",self.manual)
        form.addRow("3. Network",self.network_type);form.addRow("Confidence threshold",self.threshold_preset);form.addRow("Custom threshold",self.custom_threshold)
        form.addRow("4. Maximum hop",self.max_hop);form.addRow("Direct-neighbor selection",self.selection_mode)
        configuration.addLayout(form)
        for widget in (self.seed_count,self.presence_help,self.network_help,self.threshold_value,self.selection_help,QLabel("5. Analysis settings"),self.advanced,self.guard_help,self.readiness,self.run_button,self.progress,self.running_text):configuration.addWidget(widget)
        self.tables={name:QTableWidget() for name in ("summary","seeds","internal_nodes","internal_edges","degree1_nodes","degree2_nodes","common_direct_neighbors","seed_support","node_metrics","hubs","components","edge_evidence","evidence_summary","mapping","unmapped","ambiguous","network_summary","expanded_nodes","expanded_edges")}
        self.viewer=StringNetworkViewer();self.node_details=QLabel("Select a node in the viewer.");self.node_details.setWordWrap(True)
        self.edge_details=QLabel("Select an edge in the viewer.");self.edge_details.setWordWrap(True)
        self.neighbor_table=QTableWidget();self.support_table=QTableWidget()
        network_tabs=QTabWidget()
        for key,label in (("expanded_nodes","Expanded nodes"),("expanded_edges","Expanded edges"),("internal_nodes","Internal nodes"),("internal_edges","Internal edges"),("degree1_nodes","Degree 1"),("degree2_nodes","Degree 2")):
            network_tabs.addTab(self.tables[key],label)
        self.common_status=QLabel();common_page=self._page(self.common_status,self.tables["common_direct_neighbors"])
        self.hub_status=QLabel();hubs_page=self._page(self.hub_status,self.tables["hubs"])
        self.metrics_scope=QLabel();metrics_page=self._page(self.metrics_scope,self.tables["node_metrics"])
        evidence_tabs=QTabWidget();evidence_tabs.addTab(self.tables["edge_evidence"],"Edge evidence");evidence_tabs.addTab(self.tables["evidence_summary"],"Evidence summary")
        details=QTabWidget();details.addTab(self._page(self.node_details,self.neighbor_table,self.support_table),"Node details");details.addTab(self._page(self.edge_details),"Edge details")
        viewer_page=self._page(self.viewer,details)
        self.graph_choice=QComboBox();self.graph_preview=QLabel("No plots are available for this run.");self.graph_preview.setAlignment(Qt.AlignmentFlag.AlignCenter);self.graph_preview.setMinimumHeight(300)
        graphs_page=self._page(self.graph_choice,self.graph_preview)
        self.history=QComboBox();self.history_status=QLabel();self.history_status.setWordWrap(True);history_page=self._page(self.history,self.history_status)
        self.tabs=QTabWidget()
        for widget,label in ((config,"Configuration"),(viewer_page,"Network"),(self.tables["summary"],"Summary"),(self.tables["seeds"],"Seeds"),(network_tabs,"Internal / expansion"),(common_page,"Common neighbors"),(self.tables["seed_support"],"Seed support"),(metrics_page,"Metrics"),(hubs_page,"Hubs"),(self.tables["components"],"Components"),(evidence_tabs,"Evidence"),(self.tables["mapping"],"Mapping"),(self.tables["unmapped"],"Unmapped"),(self.tables["ambiguous"],"Ambiguous"),(graphs_page,"Graphs"),(history_page,"History")):self.tabs.addTab(widget,label)
        self.export_table=QPushButton("Export current table...");self.export_workbook=QPushButton("Export workbook...");self.export_graph=QPushButton("Export graph...");self.export_network_data=QPushButton("Export network data...");self.open_results=QPushButton("Open STRING results folder")
        actions=QHBoxLayout()
        for widget in (self.export_table,self.export_workbook,self.export_graph,self.export_network_data,self.open_results):actions.addWidget(widget)
        outer=QVBoxLayout(self);outer.addWidget(QLabel("STRING — Homo sapiens"));outer.addWidget(self.tabs,1);outer.addLayout(actions)
        self.open_database.clicked.connect(self.open_database_requested)
        self.target.currentIndexChanged.connect(self._selection_changed);self.manual.itemSelectionChanged.connect(self._selection_changed)
        self.threshold_preset.currentIndexChanged.connect(self._threshold_changed);self.custom_threshold.valueChanged.connect(self._threshold_changed)
        self.selection_mode.currentIndexChanged.connect(self._selection_changed)
        self.run_button.clicked.connect(lambda:self.run_requested.emit(self.parameters()))
        self.history.currentIndexChanged.connect(self._load_history)
        self.graph_choice.currentIndexChanged.connect(self._show_graph)
        self.viewer.node_selected.connect(self.show_node_details);self.viewer.edge_selected.connect(self.show_edge_details)
        self.export_table.clicked.connect(self._export_current);self.export_workbook.clicked.connect(self._export_workbook)
        self.export_graph.clicked.connect(self._export_graph);self.export_network_data.clicked.connect(self._export_network_data)
        self.open_results.clicked.connect(self._open_results)
        self._threshold_changed();self._selection_changed()

    @staticmethod
    def _page(*widgets):
        page=QWidget();layout=QVBoxLayout(page)
        for widget in widgets:layout.addWidget(widget)
        return page

    def set_project(self,project):
        self.project=project;self._clear_results();self._populate_targets();self._populate_manual();self.refresh_readiness();self._refresh_history()

    def _clear_results(self):
        self.outputs = None
        self.viewer.outputs = None
        self.viewer.title.setText("No STRING run loaded.")
        self.viewer.legend.clear()
        self.viewer.render()
        self._frames.clear()
        for widget in self.tables.values():
            widget.clear()
            widget.setRowCount(0)
            widget.setColumnCount(0)
        for widget in (self.neighbor_table, self.support_table):
            widget.clear()
            widget.setRowCount(0)
            widget.setColumnCount(0)
        self.graph_choice.clear()
        self.graph_preview.setPixmap(QPixmap())
        self.graph_preview.setText("No plots are available for this run.")
        self.node_details.setText("Select a node in the viewer.")
        self.edge_details.setText("Select an edge in the viewer.")
    def _populate_targets(self):
        self.target.blockSignals(True);self.target.clear();self.target.addItem("All mapped proteins","All mapped proteins")
        path=self.project.root/"analyses/presence_absence/tables/classification.csv" if self.project else None
        classes=[]
        if path and path.is_file():
            try:classes=sorted(pd.read_csv(path).classification.dropna().astype(str).unique())
            except Exception:classes=[]
        for label in dict.fromkeys(["Reproducibly detected","Shared",*[x for x in classes if x.endswith("-specific")],"Sporadic"]):
            self.target.addItem(label,label);self.target.model().item(self.target.count()-1).setEnabled(bool(classes))
        self.target.addItem("Manual selection","Manual selection");self.target.blockSignals(False)
        self.presence_help.setText("Presence/Absence-derived seed sets use existing results." if classes else "Run Presence / absence to enable derived seed sets.")

    def _populate_manual(self):
        self.manual.clear();self.mapping_frame=pd.DataFrame()
        if not self.project or not self.manager.string.is_ready():return
        try:self.mapping_frame=map_experiment(self.project,self.manager.string)
        except Exception:return
        for row in self.mapping_frame.drop_duplicates("source_row").to_dict("records"):
            item=QListWidgetItem(f"{row.get('gene_symbol','')} | {row.get('uniprot_accession','')} | {row.get('string_protein_id','')} | {row.get('preferred_name','')} | {row.get('original_id','')}")
            item.setData(Qt.ItemDataRole.UserRole,int(row["source_row"]));self.manual.addItem(item)

    def _selection_changed(self):
        self.manual.setEnabled(self.target.currentData()=="Manual selection")
        self._update_seed_count();self._mode_help()
    def _update_seed_count(self):
        if self.mapping_frame.empty:self.seed_count.setText("Selected seeds: 0");return
        if self.target.currentData()=="Manual selection":rows={str(item.data(Qt.ItemDataRole.UserRole)) for item in self.manual.selectedItems()}
        elif self.target.currentData()=="All mapped proteins":rows=set(self.mapping_frame.source_row.astype(str))
        else:
            path=self.project.root/"analyses/presence_absence/tables/classification.csv" if self.project else None
            if path and path.is_file():
                frame=pd.read_csv(path,dtype=str).fillna("");label=self.target.currentData()
                rows=set(frame.loc[frame.classification!="Not reproducibly detected" if label=="Reproducibly detected" else frame.classification==label,"source_row"])
            else:rows=set()
        selected=self.mapping_frame[self.mapping_frame.source_row.astype(str).isin(rows) & (self.mapping_frame.mapping_status=="mapped_unique")]
        count=selected.string_protein_id.nunique();self.seed_count.setText(f"Selected seeds: {count}")
        self._seed_count_value=count
        self.run_button.setEnabled(not self._running and bool(count) and "Ready for" in self.readiness.text())
    def _mode_help(self):
        if self.selection_mode.currentData()=="strict_common":
            suffix=" With one seed, this equals that seed's direct neighbors." if getattr(self,"_seed_count_value",0)==1 else ""
            self.selection_help.setText("Includes only proteins connected directly to every selected seed at the selected confidence threshold."+suffix)
        else:self.selection_help.setText("Includes proteins connected directly to at least one selected seed.")
    def _threshold_changed(self):
        custom=self.threshold_preset.currentData()=="Custom";self.custom_threshold.setEnabled(custom)
        value=self.custom_threshold.value() if custom else {"Exploratory":150,"Functional":400,"Robust":700}.get(self.threshold_preset.currentData(),400)
        self.threshold_value.setText(f"Confidence threshold: {value} (PichAnalysis preset)" if not custom else f"Confidence threshold: {value} (custom)")

    def refresh_readiness(self):
        snap=self.manager.string.active_snapshot();m=self.manager.string.manifest(snap) if snap else {}
        self.database.setText(f"Status: {'Ready' if snap else 'Not installed'}\nSTRING version: {m.get('string_version','—')}\nSnapshot: {snap.name if snap else '—'}\nFunctional network: {m.get('functional_edge_count','—')} edges\nPhysical network: {m.get('physical_edge_count','—')} edges")
        if not snap:reason="STRING database is not installed."
        elif not self.project:reason="Open a project to analyze STRING networks."
        elif str(self.project.config.get("organism_tax_id") or "")!="9606":reason="STRING analysis currently supports Homo sapiens only."
        elif not (self.project.root/"mapping/tables/protein_catalog.csv").is_file():reason="Run protein mapping before STRING analysis."
        else:reason="Ready for STRING analysis."
        self.readiness.setText(reason);self._update_seed_count()

    def parameters(self):
        return dict(target_selection=self.target.currentData(),manual_rows=tuple(int(x.data(Qt.ItemDataRole.UserRole)) for x in self.manual.selectedItems()),network_type=self.network_type.currentData(),threshold_preset=self.threshold_preset.currentData(),combined_score_threshold=self.custom_threshold.value() if self.threshold_preset.currentData()=="Custom" else None,max_hop=self.max_hop.currentData(),degree1_selection_mode=self.selection_mode.currentData(),minimum_hub_degree=self.minimum_hub_degree.value(),max_external_nodes=self.max_external_nodes.value(),top_n=self.top_n.value())
    def set_running(self,running):
        self._running=running;self.progress.setVisible(running);self.running_text.setText("Running STRING analysis..." if running else "");self._update_seed_count()

    def _fill(self,table,frame):
        self._frames[table]=frame.copy();table.clear();table.setRowCount(len(frame));table.setColumnCount(len(frame.columns))
        table.setHorizontalHeaderLabels([str(name).replace("_"," ").title() for name in frame.columns])
        for row,values in enumerate(frame.itertuples(index=False,name=None)):
            for col,value in enumerate(values):table.setItem(row,col,QTableWidgetItem("NA" if pd.isna(value) else str(value)))
    def show_outputs(self,outputs):
        self.outputs=outputs;tables=outputs["tables"];metadata=outputs["metadata"]
        for name,widget in self.tables.items():self._fill(widget,tables[name])
        metrics=tables["node_metrics"]
        self.metrics_scope.setText(f"Metrics graph scope: {metadata.get('metrics_graph_scope','—')}. Hop level and network degree are different quantities.")
        self.hub_status.setText(f"Hub definition: network degree ≥ {metadata.get('minimum_hub_degree','—')}" if not tables["hubs"].empty else f"Hub definition: network degree ≥ {metadata.get('minimum_hub_degree','—')}. No nodes met the current hub definition.")
        self.common_status.setText("No proteins were directly connected to all selected seeds at this threshold." if metadata.get("degree1_selection_mode")=="strict_common" and tables["common_direct_neighbors"].empty else "Common direct neighbors are connected to every selected seed at this threshold.")
        self.viewer.load_outputs(outputs)
        self.graph_choice.blockSignals(True);self.graph_choice.clear()
        for path in sorted(outputs["plots"]):
            if path.suffix.lower()==".png":self.graph_choice.addItem(path.stem.replace("_"," ").title(),{"png":path,"pdf":path.with_suffix(".pdf") if path.with_suffix(".pdf").is_file() else None})
        self.graph_choice.blockSignals(False);self._show_graph()
        self.node_details.setText("Select a node in the viewer.");self.edge_details.setText("Select an edge in the viewer.")
        self.history_status.setText(f"Loaded run {metadata.get('run_id')} | Snapshot {metadata.get('snapshot_id')} | STRING {metadata.get('string_version')}")
        self._refresh_history(metadata.get("run_id"))
    def show_node_details(self,string_id):
        if not self.outputs:return
        tables=self.outputs["tables"];nodes=tables["expanded_nodes"];metrics=tables["node_metrics"]
        n=nodes[nodes.string_protein_id.astype(str)==str(string_id)];m=metrics[metrics.string_protein_id.astype(str)==str(string_id)]
        if n.empty:return
        row=n.iloc[0];metric=m.iloc[0] if not m.empty else pd.Series(dtype=object)
        fields=(("Gene",row.get("gene_symbol")),("UniProt",row.get("uniprot_accession")),("STRING ID",string_id),("Preferred name",row.get("preferred_name")),("STRING annotation",row.get("annotation")),("Hop level",row.get("hop_level")),("Seed / external","Seed" if row.get("hop_level")==0 else "External"),("Network degree",metric.get("network_degree")),("Betweenness",metric.get("betweenness")),("Closeness",metric.get("closeness")),("Clustering coefficient",metric.get("local_clustering_coefficient")),("Component",metric.get("component_id")),("Seed support count",row.get("seed_support_count")))
        text="\n".join(f"{label}: {'NA' if pd.isna(value) else value}" for label,value in fields)
        if row.get("hop_level")==2:text+="\nHop-2 support indicates recorded reachability, not necessarily a direct seed-to-node edge."
        self.node_details.setText(text);self._fill(self.neighbor_table,string_node_neighbors(self.outputs,string_id));self._fill(self.support_table,string_node_support(self.outputs,string_id))
    def show_edge_details(self,a,b):
        if not self.outputs:return
        frame=string_edge_details(self.outputs,a,b)
        if frame.empty:return
        row=frame.iloc[0];lines=[f"Protein A: {a}",f"Protein B: {b}",f"Network type: {self.outputs['metadata'].get('network_type')}",f"STRING confidence score: {row.get('combined_score')}"]
        lines += [f"{name}: {'NA' if pd.isna(value) else value}" for name,value in row.items() if name not in {"protein_a","protein_b","combined_score"}]
        self.edge_details.setText("\n".join(lines))

    def _show_graph(self):
        data=self.graph_choice.currentData()
        if not data or not Path(data["png"]).is_file():self.graph_preview.setText("No plots are available for this run.");self.graph_preview.setPixmap(QPixmap());return
        self.graph_preview.setPixmap(QPixmap(str(data["png"])).scaled(900,600,Qt.AspectRatioMode.KeepAspectRatio,Qt.TransformationMode.SmoothTransformation))
    def _refresh_history(self,selected=None):
        self.history.blockSignals(True);self.history.clear()
        if self.project:
            for run_id in list_string_runs(self.project):
                try:
                    m=json.loads((self.project.root/"analyses/STRING/runs"/run_id/"metadata.json").read_text(encoding="utf-8"))
                    label=f"{m.get('created_at','—')} | {m.get('target_definition','—')} | {m.get('network_type','—')} | ≥{m.get('combined_score_threshold','—')} ({m.get('threshold_preset','—')}) | Hop {m.get('max_hop','—')} | {m.get('degree1_selection_mode','—')} | Snapshot {m.get('snapshot_id','—')} | STRING {m.get('string_version','—')}"
                except Exception:label=f"{run_id} | Incomplete metadata"
                self.history.addItem(label,run_id)
        self.history.setCurrentIndex(self.history.findData(selected) if selected else -1);self.history.blockSignals(False)
    def _load_history(self,index):
        if not self.project or index<0:return
        try:self.show_outputs(read_string_outputs(self.project,self.history.itemData(index)))
        except Exception as error:self._clear_results();self.history_status.setText(f"Historical run is incomplete: {error}")

    def _current_table(self):
        current=self.tabs.currentWidget()
        if isinstance(current,QTableWidget):return current
        if isinstance(current,QTabWidget):
            child=current.currentWidget()
            if isinstance(child,QTableWidget):return child
        tables=current.findChildren(QTableWidget) if current else []
        return tables[0] if tables else None
    def _export_current(self):
        table=self._current_table();frame=self._frames.get(table)
        if frame is None:return QMessageBox.information(self,"Export table","The current panel does not contain an exportable table.")
        path,_=QFileDialog.getSaveFileName(self,"Export current table","string_table.csv","CSV (*.csv)")
        if path:frame.to_csv(path,index=False)
    def _export_workbook(self):
        if not self.outputs:return
        path,_=QFileDialog.getSaveFileName(self,"Export workbook","STRING_analysis.xlsx","Excel workbook (*.xlsx)")
        if path:shutil.copy2(self.outputs["workbook"],path)
    def _export_graph(self):
        data=self.graph_choice.currentData()
        if not data:return
        path,selected=QFileDialog.getSaveFileName(self,"Export graph","string_plot.png","PNG (*.png);;PDF (*.pdf)")
        if path:
            source=data.get("pdf" if "PDF" in selected or Path(path).suffix.lower()==".pdf" else "png")
            if source:shutil.copy2(source,path)
    def _export_network_data(self):
        if not self.outputs:return
        folder=QFileDialog.getExistingDirectory(self,"Export full network data")
        if folder:
            run=self.outputs["run_root"]
            shutil.copy2(run/"networks/expanded_nodes.csv",Path(folder)/"string_nodes.csv")
            shutil.copy2(run/"networks/expanded_edges.csv",Path(folder)/"string_edges.csv")
    def _open_results(self):
        if self.outputs:QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.outputs["run_root"])))
