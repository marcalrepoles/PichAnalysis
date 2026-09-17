import json
from pathlib import Path
import pandas as pd
from PySide6.QtGui import QColor,QImage
from PySide6.QtWidgets import QApplication
from pichanalysis.core.kegg_analysis import read_kegg_outputs
from pichanalysis.core.pathway_viewer import export_overlay
from pichanalysis.core.project import create_project
from pichanalysis.ui.kegg_page import KEGGPage
from pichanalysis.ui.main_window import MainWindow

def app():return QApplication.instance() or QApplication([])
class Manager:
 def active_snapshot(self):return None
 def pathway_analysis_compatibility(self):return False,"No active snapshot"
def test_zoom_limits_actual_fit_and_export_resolution(tmp_path):
 app();page=KEGGPage(Manager());page.overlay=QImage(400,200,QImage.Format.Format_RGB32);page.overlay.fill(QColor("white"));page.scroll.resize(202,102)
 page._set_scale(1);page.zoom_in.click();assert page.scale==1.25;page.zoom_out.click();assert page.scale==1
 page._set_scale(99);assert page.scale==8;page._set_scale(.001);assert page.scale==.1;page.actual.click();assert page.scale==1
 page.scroll.viewport().resize(200,100);page.fit_to_window();assert abs(page.scale-.5)<.03
 page._set_scale(4);out=export_overlay(page.overlay,tmp_path/"x.png");saved=QImage(str(out));assert (saved.width(),saved.height())==(400,200)
def _write_run(root,name,value):
 run=root/"analyses/KEGG/runs"/name
 for d in ("mapping","frequency","enrichment","graphs"): (run/d).mkdir(parents=True,exist_ok=True)
 pd.DataFrame({"included_in_analysis":[True],"gene_symbol":[value],"kegg_gene_id":["1"]}).to_csv(run/"mapping/kegg_mapping.csv",index=False)
 pd.DataFrame({"pathway_id":["hsa00010"],"name":[value],"gene_count":[1]}).to_csv(run/"frequency/kegg_pathway_frequency.csv",index=False)
 pd.DataFrame({"pathway_id":["hsa00010"],"name":[value]}).to_csv(run/"enrichment/kegg_pathway_all.csv",index=False)
 pd.DataFrame({"metric":["run"],"value":[value]}).to_csv(run/"summary.csv",index=False);(run/"metadata.json").write_text(json.dumps({"run_id":name}),encoding="utf-8")
def test_history_a_b_a_and_incomplete(tmp_path):
 project=create_project(tmp_path,"History");_write_run(project.root,"A","alpha");_write_run(project.root,"B","beta")
 assert read_kegg_outputs(project,"A").summary.value.iloc[0]=="alpha";assert read_kegg_outputs(project,"B").summary.value.iloc[0]=="beta";assert read_kegg_outputs(project,"A").summary.value.iloc[0]=="alpha"
 (project.root/"analyses/KEGG/runs/broken").mkdir();
 try:read_kegg_outputs(project,"broken");assert False
 except RuntimeError:pass
def test_target_background_cancel_and_continue(monkeypatch):
 app();window=MainWindow();window.kegg_worker=object();calls=[];monkeypatch.setattr(window,"_run_kegg",lambda p,allow=False:calls.append((p,allow)));monkeypatch.setattr(window.analyses_page.kegg_page,"set_running",lambda x:None)
 monkeypatch.setattr(window,"_confirm_target_adjustment",lambda n:False);window._kegg_failed("2 target KEGG gene(s) are outside","",{"x":1});assert calls==[]
 window.kegg_worker=object();monkeypatch.setattr(window,"_confirm_target_adjustment",lambda n:True);window._kegg_failed("2 target KEGG gene(s) are outside","",{"x":1});assert calls==[({"x":1},True)]
