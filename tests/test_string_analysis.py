import gzip
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from pichanalysis.core.database_manager import DatabaseManager
from pichanalysis.core.r_runtime import RRuntime
from pichanalysis.core.string_analysis import (
    StringParameters, StringExpansionLimitError, StringTargetOutsideBackgroundError,
    StringPresenceRequiredError, map_experiment, prepare_string_run, read_string_outputs,
    run_string_analysis, list_string_runs,
)

IDS={name:f"9606.ENSP{name}" for name in "ABCXYZW"}
def _gz(path,lines):
    with gzip.open(path,"wt",encoding="utf-8") as handle:handle.write("\n".join(lines)+"\n")

def scientific_fixture(tmp_path):
    source=tmp_path/"raw";source.mkdir()
    _gz(source/"9606.protein.info.vfixture.txt.gz",["#string_protein_id\tpreferred_name\tprotein_size\tannotation"]+[f"{id}\t{name}\t100\t{name} protein" for name,id in IDS.items()])
    aliases=["#string_protein_id\talias\tsource"]+[f"{IDS[name]}\tP{name}\tUniProt_AC" for name in "ABC"]
    aliases += [f"{IDS['A']}\tDUP\tUniProt_AC",f"{IDS['B']}\tDUP\tUniProt_AC"]
    _gz(source/"9606.protein.aliases.vfixture.txt.gz",aliases)
    header="protein1 protein2 experiments experiments_transferred textmining combined_score"
    edge=lambda a,b,score:f"{IDS[a]} {IDS[b]} {score-100} 7 23 {score}"
    functional=[edge("A","B",900),edge("A","X",800),edge("B","X",600),edge("C","X",750),edge("A","Y",500),edge("B","Y",500),edge("A","Z",500),edge("C","Z",500),edge("X","W",450),edge("Z","W",450)]
    _gz(source/"9606.protein.links.full.vfixture.txt.gz",[header,*functional])
    _gz(source/"9606.protein.physical.links.full.vfixture.txt.gz",[header,edge("A","X",850)])
    manager=DatabaseManager(tmp_path/"db");manager.string.install_from_directory(source)
    project=SimpleNamespace(root=tmp_path/"project",config={"organism_tax_id":"9606"})
    catalog=project.root/"mapping/tables/protein_catalog.csv";catalog.parent.mkdir(parents=True)
    pd.DataFrame([dict(source_row=i,original_id=f"row{i}",uniprot_accession=accession,gene_symbol=gene) for i,(accession,gene) in enumerate([("PA","A"),("PB","B"),("PC","C"),("DUP","ambiguous"),("NONE","unmapped"),("PA","A")],1)]).to_csv(catalog,index=False)
    return project,manager

def test_mapping_and_strict_thresholds(tmp_path):
    project,manager=scientific_fixture(tmp_path);mapping=map_experiment(project,manager.string)
    assert set(mapping.mapping_status)=={"mapped_unique","ambiguous","unmapped"}
    assert len(set(mapping.loc[mapping.mapping_status=="mapped_unique","string_protein_id"]))==3
    run,_=prepare_string_run(project,manager,"strict400",StringParameters(max_hop=2,degree1_selection_mode="strict_common"))
    common=pd.read_csv(run/"networks/common_direct_neighbors.csv")
    assert common.string_protein_id.tolist()==[IDS["X"]]
    d1=pd.read_csv(run/"networks/degree1_nodes.csv");assert d1.string_protein_id.tolist()==[IDS["X"]]
    d2=pd.read_csv(run/"networks/degree2_nodes.csv");assert set(d2.string_protein_id)=={IDS["W"]}
    support=pd.read_csv(run/"networks/seed_support.csv");assert len(support[support.string_protein_id==IDS["X"]])==3
    run,_=prepare_string_run(project,manager,"strict700",StringParameters(max_hop=1,degree1_selection_mode="strict_common",threshold_preset="Robust"))
    assert pd.read_csv(run/"networks/degree1_nodes.csv").empty
    assert pd.read_csv(run/"networks/common_direct_neighbors.csv").empty

def test_shortest_hop_internal_and_guard(tmp_path):
    project,manager=scientific_fixture(tmp_path)
    run,_=prepare_string_run(project,manager,"union",StringParameters(max_hop=2))
    internal=pd.read_csv(run/"networks/internal_edges.csv");assert len(internal)==1
    internal_nodes=pd.read_csv(run/"networks/internal_nodes.csv");assert set(internal_nodes.string_protein_id)=={IDS[x] for x in "ABC"}
    d1=pd.read_csv(run/"networks/degree1_nodes.csv");d2=pd.read_csv(run/"networks/degree2_nodes.csv")
    assert set(d1.string_protein_id)=={IDS[x] for x in "XYZ"}
    assert set(d2.string_protein_id)=={IDS["W"]}
    assert not set(d1.string_protein_id)&set(d2.string_protein_id)
    with pytest.raises(StringExpansionLimitError) as err:prepare_string_run(project,manager,"guard",StringParameters(max_hop=2,max_external_nodes=2))
    assert err.value.details["observed"]==3

def test_persistence_r_outputs_and_snapshots(tmp_path):
    project,manager=scientific_fixture(tmp_path);runtime=RRuntime()
    if not runtime.available:pytest.skip("Rscript unavailable")
    a=run_string_analysis(project,manager,runtime,run_id="a",parameters=StringParameters(max_hop=2),timeout=120)
    assert a["tables"]["network_summary"].iloc[0].expanded_node_count==7
    assert a["tables"]["internal_nodes"].network_degree.tolist().count(0)==1
    assert len(a["plots"])==12 and a["workbook"].is_file()
    assert a["metadata"]["snapshot_id"]==manager.string.active_snapshot().name
    assert (project.root/"analyses/STRING/mapping/string_mapping.csv").is_file()
    assert (project.root/"analyses/STRING/STRING_analysis.xlsx").is_file()
    assert (project.root/"scripts/runs/a_string/R_session_info.txt").is_file()
    b=run_string_analysis(project,manager,runtime,run_id="b",parameters=StringParameters(max_hop=1,network_type="physical"),timeout=120)
    assert len(b["tables"]["expanded_edges"])==1
    assert len(read_string_outputs(project,"a")["tables"]["expanded_edges"])>1
    assert list_string_runs(project)==["a","b"]

def test_errors(tmp_path):
    project,manager=scientific_fixture(tmp_path)
    with pytest.raises(StringPresenceRequiredError):prepare_string_run(project,manager,"presence",StringParameters(target_selection="Shared"))
    with pytest.raises(StringTargetOutsideBackgroundError):prepare_string_run(project,manager,"outside",StringParameters(background_selection="Manual selection",background_manual_rows=(1,)))


def test_zero_edge_single_seed_and_history_snapshot(tmp_path):
    project,manager=scientific_fixture(tmp_path);runtime=RRuntime()
    if not runtime.available:pytest.skip("Rscript unavailable")
    first=manager.string.active_snapshot().name
    isolated=run_string_analysis(project,manager,runtime,run_id="isolated",parameters=StringParameters(target_selection="Manual selection",manual_rows=(3,),max_hop=0,network_type="physical",degree1_selection_mode="strict_common"),timeout=120)
    assert isolated["tables"]["network_summary"].iloc[0].isolated_seed_count==1
    assert isolated["tables"]["hubs"].empty
    assert isolated["tables"]["node_metrics"].iloc[0].network_degree==0
    second=manager.string.install_from_directory(tmp_path/"raw")
    assert second.name!=first
    newer=run_string_analysis(project,manager,runtime,run_id="newer",parameters=StringParameters(max_hop=0),timeout=120)
    assert newer["metadata"]["snapshot_id"]==second.name
    assert read_string_outputs(project,"isolated")["metadata"]["snapshot_id"]==first