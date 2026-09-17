"""Manual online KEGG integration smoke. Not collected by pytest."""
from __future__ import annotations
import hashlib,tempfile
from pathlib import Path
import pandas as pd
from pichanalysis.core.database_manager import DatabaseManager
from pichanalysis.core.kegg_analysis import prepare_kegg_arguments,read_kegg_outputs
from pichanalysis.core.organism import set_organism
from pichanalysis.core.pathway_viewer import export_overlay,render_overlay
from pichanalysis.core.project import create_project
from pichanalysis.core.r_runtime import RRuntime

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
 root=Path(tempfile.mkdtemp(prefix="pichanalysis-kegg-real-e2e-"));manager=DatabaseManager(root/"db");snapshot=manager.download(["kgml","images"],pathway_subset=["hsa00010","hsa00020"])
 links=pd.read_csv(snapshot/"tables/gene_to_pathway.tsv",sep="\t",dtype=str);conv=pd.read_csv(snapshot/"tables/kegg_to_ncbi_geneid.tsv",sep="\t",dtype=str)
 chosen=links[links.pathway_id.isin(["hsa00010","hsa00020"])].merge(conv,left_on="gene_id",right_on="kegg_gene_id").drop_duplicates("kegg_gene_id").head(8);assert len(chosen)>=3
 project=create_project(root,"Smoke");set_organism(project,"Homo sapiens","9606");catalog=pd.DataFrame({"source_row":range(1,len(chosen)+1),"original_id":chosen.ncbi_gene_id,"ncbi_gene_id":chosen.ncbi_gene_id,"uniprot_accession":"","gene_symbol":"","protein_name":"","mapping_status":"mapped"});(project.root/"mapping/tables").mkdir(parents=True,exist_ok=True);catalog.to_csv(project.root/"mapping/tables/protein_catalog.csv",index=False)
 args=prepare_kegg_arguments(project,manager,target_selection="all_experiment",background_selection="all_experiment",manual_rows=[],mapping_policy="unique-only",fdr_cutoff=.05,p_cutoff=1,min_count=1,top_n=20,run_id="real_smoke")
 result=RRuntime().run(Path(__file__).parents[1]/"r_scripts/04_kegg_analysis.R",*args);assert result.returncode==0,result.stderr
 outputs=read_kegg_outputs(project);assert len(outputs.mapping)>0 and len(outputs.frequency)>0 and len(outputs.enrichment)>0
 membership=pd.read_csv(project.root/"analyses/KEGG/mapping/pathway_membership.csv");nodes=pd.read_csv(project.root/"analyses/KEGG/pathways/pathway_nodes.csv");assert len(membership)>0 and nodes.is_target_hit.astype(bool).any()
 pid=str(nodes[nodes.is_target_hit.astype(bool)].pathway_id.iloc[0]);original=snapshot/"images"/f"{pid}.png";before=sha(original);image=render_overlay(original,nodes[nodes.pathway_id==pid]);exported=export_overlay(image,root/"external_highlighted.png");assert before==sha(original) and exported.is_file() and image.width()>0
 assert (project.root/"analyses/KEGG/KEGG_analysis.xlsx").is_file();assert outputs.metadata["snapshot_id"]==snapshot.name
 print(f"REAL_KEGG_SMOKE_OK root={root} snapshot={snapshot.name} pathways=hsa00010,hsa00020 mapped={len(outputs.mapping)} membership={len(membership)} nodes={len(nodes)} png_sha256={before} exported_sha256={sha(exported)}")
if __name__=="__main__":main()
