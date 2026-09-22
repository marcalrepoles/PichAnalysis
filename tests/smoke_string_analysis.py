"""Run the controlled offline STRING scientific smoke explicitly."""
import tempfile
from pathlib import Path

from test_string_analysis import scientific_fixture, IDS
from pichanalysis.core.r_runtime import RRuntime
from pichanalysis.core.string_analysis import StringParameters, run_string_analysis

def main():
    root=Path(tempfile.mkdtemp(prefix="pichanalysis-string-analysis-"))
    project,manager=scientific_fixture(root);manager.string.provider=None
    outputs={}
    for preset in ("Exploratory","Functional","Robust"):
        outputs[preset]=run_string_analysis(project,manager,RRuntime(),run_id=preset.lower(),parameters=StringParameters(max_hop=2,threshold_preset=preset),timeout=120)
    edges=[len(outputs[p]["tables"]["expanded_edges"]) for p in ("Exploratory","Functional","Robust")]
    neighbors=[len(outputs[p]["tables"]["degree1_nodes"])+len(outputs[p]["tables"]["degree2_nodes"]) for p in ("Exploratory","Functional","Robust")]
    assert edges[2]<=edges[1]<=edges[0] and neighbors[2]<=neighbors[1]<=neighbors[0]
    strict=run_string_analysis(project,manager,RRuntime(),run_id="strict",parameters=StringParameters(max_hop=2,degree1_selection_mode="strict_common"),timeout=120)
    assert strict["tables"]["common_direct_neighbors"].string_protein_id.tolist()==[IDS["X"]]
    assert len(strict["tables"]["hubs"])>=1 and len(strict["plots"])==12
    print({"temporary_path":str(root),"edge_counts_150_400_700":edges,"external_node_counts_150_400_700":neighbors,"strict_common":[IDS["X"]],"runs":4})
if __name__=="__main__":main()
