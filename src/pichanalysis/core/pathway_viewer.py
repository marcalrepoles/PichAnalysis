from __future__ import annotations
from pathlib import Path
import pandas as pd
from PySide6.QtCore import QRectF
from PySide6.QtGui import QColor,QImage,QPainter,QPen

def render_overlay(source:Path,nodes:pd.DataFrame,color:str="#ff1744",width:int=4)->QImage:
 image=QImage(str(source))
 if image.isNull():raise ValueError("Pathway image could not be loaded.")
 result=image.copy();p=QPainter(result);p.setPen(QPen(QColor(color),width))
 for row in nodes.itertuples():
  if bool(row.is_target_hit):p.drawRect(QRectF(float(row.x)-float(row.width)/2,float(row.y)-float(row.height)/2,float(row.width),float(row.height)))
 p.end();return result
def export_overlay(image:QImage,destination:Path)->Path:
 destination.parent.mkdir(parents=True,exist_ok=True)
 if not image.save(str(destination),"PNG"):raise OSError(f"Could not save highlighted pathway to {destination}")
 return destination
def pathway_proteins(membership:pd.DataFrame,pathway_id:str)->pd.DataFrame:return membership[membership.pathway_id==pathway_id].copy()
def protein_pathways(membership:pd.DataFrame,kegg_gene_id:str)->pd.DataFrame:
 columns=[x for x in ("pathway_id","name") if x in membership.columns]
 return membership[membership.kegg_gene_id.astype(str)==str(kegg_gene_id)][columns].drop_duplicates()
