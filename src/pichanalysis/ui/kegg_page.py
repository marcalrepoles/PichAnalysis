from __future__ import annotations
from pathlib import Path
from PySide6.QtCore import Signal,Qt
from PySide6.QtGui import QPixmap
import pandas as pd
from PySide6.QtWidgets import QComboBox,QDoubleSpinBox,QFileDialog,QFormLayout,QHBoxLayout,QLabel,QLineEdit,QMessageBox,QProgressBar,QPushButton,QScrollArea,QSpinBox,QTableWidget,QTableWidgetItem,QTabWidget,QVBoxLayout,QWidget
from ..core.database_manager import DatabaseManager
from ..core.go_analysis import available_sets
from ..core.kegg_analysis import KEGGOutputs,kegg_readiness,read_kegg_outputs
from ..core.project import Project
from ..core.pathway_viewer import render_overlay,export_overlay,pathway_proteins,protein_pathways
from ..core.kegg_analysis import list_kegg_runs
class KEGGPage(QWidget):
 run_requested=Signal(dict);open_database_requested=Signal();export_requested=Signal(str)
 def __init__(self,manager:DatabaseManager):
  super().__init__();self.manager=manager;self.project=None;self.outputs=None
  self.database=QLabel();self.source=QLabel("Source: Local KEGG snapshot");self.open_database=QPushButton("Open Database Manager")
  self.target=QComboBox();self.background=QComboBox();self.policy=QComboBox();self.policy.addItem("Unambiguous mappings only","unique-only");self.policy.addItem("Include all mapped candidates — exploratory","all-candidates")
  self.manual=QLineEdit();self.manual.setPlaceholderText("Comma-separated source row numbers");self.fdr=QDoubleSpinBox();self.fdr.setRange(.000001,1);self.fdr.setValue(.05);self.pvalue=QDoubleSpinBox();self.pvalue.setRange(.000001,1);self.pvalue.setValue(1);self.minimum=QSpinBox();self.minimum.setRange(1,10000);self.minimum.setValue(3);self.top=QSpinBox();self.top.setRange(1,200);self.top.setValue(20)
  form=QFormLayout();form.addRow("1. Select protein set",self.target);form.addRow("Manual selection",self.manual);form.addRow("2. Select background",self.background);form.addRow("3. Select mapping policy",self.policy);form.addRow("FDR cutoff",self.fdr);form.addRow("P-value cutoff",self.pvalue);form.addRow("Minimum gene count",self.minimum);form.addRow("Top N",self.top)
  self.help=QLabel("Pathway frequency shows how many genes from the selected protein set belong to each KEGG pathway. Pathway enrichment tests whether a pathway contains more genes from the selected set than expected from the selected background. The background represents the genes that could have been selected in this experiment. FDR controls multiple-testing errors.");self.help.setWordWrap(True)
  self.ready=QLabel();self.run=QPushButton("5. Run KEGG analysis");self.progress=QProgressBar();self.progress.setRange(0,0);self.progress.hide();config=QWidget();cl=QVBoxLayout(config)
  for w in (self.database,self.source,self.open_database):cl.addWidget(w)
  cl.addLayout(form);cl.addWidget(self.help);cl.addWidget(self.ready);cl.addWidget(self.run);cl.addWidget(self.progress)
  self.summary=QTableWidget();self.frequency=QTableWidget();self.enrichment=QTableWidget();self.pathway=QComboBox();self.image=QLabel("Pathway images are not installed for the active KEGG snapshot.");self.image.setAlignment(Qt.AlignmentFlag.AlignCenter);self.proteins=QTableWidget();self.scale=1.0;self.overlay=None
  self.fit=QPushButton("Fit to Window");self.actual=QPushButton("Actual Size");self.zoom_in=QPushButton("Zoom In");self.zoom_out=QPushButton("Zoom Out");self.export=QPushButton("Export Highlighted Pathway...");controls=QHBoxLayout()
  for b in (self.fit,self.actual,self.zoom_in,self.zoom_out,self.export):controls.addWidget(b)
  scroll=QScrollArea();scroll.setWidget(self.image);scroll.setWidgetResizable(False);self.scroll=scroll
  viewer=QWidget();vl=QVBoxLayout(viewer);vl.addWidget(self.pathway);vl.addLayout(controls);vl.addWidget(scroll,1);vl.addWidget(QLabel("KGML is required to highlight experimental genes on this pathway."));vl.addWidget(self.proteins)
  self.protein=QComboBox();self.protein_paths=QTableWidget();self.open_path=QPushButton("Open in Pathway Viewer");protein_view=QWidget();pl=QVBoxLayout(protein_view);pl.addWidget(QLabel("Select protein/gene"));pl.addWidget(self.protein);pl.addWidget(self.protein_paths);pl.addWidget(self.open_path)
  self.history=QComboBox();history=QWidget();hl=QVBoxLayout(history);hl.addWidget(QLabel("Previous KEGG runs"));hl.addWidget(self.history);hl.addStretch()
  tabs=QTabWidget();self.tabs=tabs;tabs.addTab(config,"Configuration");tabs.addTab(self.summary,"Summary");tabs.addTab(self.frequency,"Pathway Frequency");tabs.addTab(self.enrichment,"Pathway Enrichment");tabs.addTab(viewer,"Pathway Viewer");tabs.addTab(protein_view,"Proteins");tabs.addTab(history,"History");layout=QVBoxLayout(self);layout.addWidget(tabs)
  self.run.clicked.connect(self._emit);self.open_database.clicked.connect(self.open_database_requested);self.pathway.currentIndexChanged.connect(self._show_pathway);self.fit.clicked.connect(self.fit_to_window);self.actual.clicked.connect(lambda:self._set_scale(1));self.zoom_in.clicked.connect(lambda:self._set_scale(self.scale*1.25));self.zoom_out.clicked.connect(lambda:self._set_scale(self.scale/1.25));self.export.clicked.connect(self._export);self.protein.currentIndexChanged.connect(self._show_protein_paths);self.open_path.clicked.connect(self._open_selected_path);self.history.currentTextChanged.connect(self._load_history)
 def set_project(self,project:Project|None):
  self.project=project;self.target.clear();self.background.clear()
  if project:
   for key,label in available_sets(project).items():self.target.addItem(label,key);self.background.addItem(label,key)
   i=self.background.findData("mapped");self.background.setCurrentIndex(max(i,0))
   self.background.insertItem(0,"All KEGG-mapped genes in the experiment","all_experiment");self.background.setCurrentIndex(0);self.history.addItems(list_kegg_runs(project))
  state=kegg_readiness(project,self.manager);ok,msg=self.manager.pathway_analysis_compatibility();snap=self.manager.active_snapshot();self.database.setText(f"KEGG database\nStatus: {msg}\nSnapshot: {snap.name if snap else '—'}")
  self.ready.setText(state.reason);self.run.setEnabled(state.ready)
 def _emit(self):
  try:rows=[int(x.strip()) for x in self.manual.text().split(",") if x.strip()]
  except ValueError:self.ready.setText("Manual selection must contain row numbers.");return
  self.run_requested.emit({"target_selection":self.target.currentData(),"background_selection":self.background.currentData(),"manual_rows":rows,"mapping_policy":self.policy.currentData(),"fdr_cutoff":self.fdr.value(),"p_cutoff":self.pvalue.value(),"min_count":self.minimum.value(),"top_n":self.top.value()})
 def set_running(self,running):self.progress.setVisible(running);self.run.setEnabled(not running)
 def _fill(self,table,frame):
  table.setRowCount(len(frame));table.setColumnCount(len(frame.columns));table.setHorizontalHeaderLabels(list(frame.columns))
  for r,row in enumerate(frame.itertuples(index=False,name=None)):
   for c,v in enumerate(row):table.setItem(r,c,QTableWidgetItem("" if str(v)=="nan" else str(v)))
 def show_outputs(self,o:KEGGOutputs):
  self.outputs=o;self._fill(self.summary,o.summary);self._fill(self.frequency,o.frequency);self._fill(self.enrichment,o.enrichment);self.pathway.clear()
  for row in o.frequency.itertuples():self.pathway.addItem(f"{row.pathway_id} — {getattr(row,'name','')}",row.pathway_id)
  self.protein.clear()
  for row in o.mapping[o.mapping.included_in_analysis.astype(bool)].itertuples():self.protein.addItem(f"{row.gene_symbol} — {row.kegg_gene_id}",str(row.kegg_gene_id))
 def _show_pathway(self):
  pid=self.pathway.currentData();snap=self.manager.active_snapshot()
  if not pid or not snap:return
  image=snap/"images"/f"{pid}.png"
  if image.is_file():
   nodes_path=self.outputs.root/"pathways/pathway_nodes.csv" if self.outputs else None
   if nodes_path and nodes_path.is_file():
    try:
     nodes=pd.read_csv(nodes_path);hits=nodes[nodes["pathway_id"]==pid];self.overlay=render_overlay(image,hits);self._set_scale(self.scale)
    except (OSError,ValueError,KeyError):pass
  if self.outputs:
   membership_path=self.outputs.root/"mapping/pathway_membership.csv"
   if membership_path.is_file():
    members=pd.read_csv(membership_path);members=members[members["pathway_id"]==pid]
    wanted=[x for x in ("gene_symbol","kegg_gene_id","uniprot_accession","protein_name","original_id","uniprot_function","ncbi_summary") if x in members.columns];self._fill(self.proteins,members[wanted])
 def _set_scale(self,value):
  self.scale=max(.1,min(8.0,float(value)))
  if self.overlay is not None:self.image.setPixmap(QPixmap.fromImage(self.overlay).scaled(self.overlay.width()*self.scale,self.overlay.height()*self.scale,Qt.AspectRatioMode.KeepAspectRatio,Qt.TransformationMode.SmoothTransformation));self.image.adjustSize()
 def fit_to_window(self):
  if self.overlay is not None:self._set_scale(min(self.scroll.viewport().width()/self.overlay.width(),self.scroll.viewport().height()/self.overlay.height()))
 def _export(self):
  if self.overlay is None:return
  pid=self.pathway.currentData();internal=self.project.root/"analyses/KEGG/pathways"/pid/"highlighted.png";export_overlay(self.overlay,internal);name,_=QFileDialog.getSaveFileName(self,"Export Highlighted Pathway",f"{pid}_highlighted.png","PNG (*.png)")
  if name:
   try:export_overlay(self.overlay,Path(name))
   except OSError as e:QMessageBox.warning(self,"Export failed",str(e))
 def _show_protein_paths(self):
  if not self.project or not self.protein.currentData():return
  p=self.outputs.root/"mapping/pathway_membership.csv" if self.outputs else Path()
  if p.is_file():self._fill(self.protein_paths,protein_pathways(pd.read_csv(p),self.protein.currentData()))
 def _open_selected_path(self):
  row=self.protein_paths.currentRow()
  if row>=0:
   pid=self.protein_paths.item(row,0).text();i=self.pathway.findData(pid)
   if i>=0:self.pathway.setCurrentIndex(i);self.tabs.setCurrentIndex(4)
 def _load_history(self,run_id):
  if not self.project or not run_id:return
  try:self.show_outputs(read_kegg_outputs(self.project,run_id))
  except RuntimeError as error:self.ready.setText(f"Historical run is incomplete: {error}")
