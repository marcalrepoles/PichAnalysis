from pathlib import Path
import pandas as pd
from PySide6.QtGui import QColor,QImage
from pichanalysis.core.pathway_viewer import export_overlay,pathway_proteins,protein_pathways,render_overlay
def test_overlay_multiple_hits_and_export(tmp_path):
 source=tmp_path/"base.png";image=QImage(100,100,QImage.Format.Format_RGB32);image.fill(QColor("white"));assert image.save(str(source))
 nodes=pd.DataFrame([{"x":20,"y":20,"width":10,"height":10,"is_target_hit":True},{"x":70,"y":70,"width":12,"height":12,"is_target_hit":True},{"x":50,"y":50,"width":10,"height":10,"is_target_hit":False}])
 rendered=render_overlay(source,nodes);assert rendered.size()==image.size();assert rendered.pixelColor(15,15)!=QColor("white")
 destination=export_overlay(rendered,tmp_path/"out"/"highlighted.png");assert destination.is_file() and source.read_bytes()!=destination.read_bytes()
def test_bidirectional_navigation_helpers():
 frame=pd.DataFrame({"pathway_id":["hsa00010","hsa00020","hsa00010"],"name":["A","B","A"],"kegg_gene_id":["1","1","2"]})
 assert len(pathway_proteins(frame,"hsa00010"))==2
 assert set(protein_pathways(frame,"1").pathway_id)=={"hsa00010","hsa00020"}
