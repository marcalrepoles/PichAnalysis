import gzip, json, sqlite3
from pathlib import Path
import pytest

from pichanalysis.core.database_registry import DatabaseState
from pichanalysis.core.databases.string import StringDownloadCancelled, StringProvider
from pichanalysis.core.string_database import StringDatabase

EVIDENCE="neighborhood neighborhood_transferred fusion cooccurence homology coexpression coexpression_transferred experiments experiments_transferred database database_transferred textmining textmining_transferred combined_score"

def gz(path,text):
    with gzip.open(path,"wt",encoding="utf-8",newline="") as f:f.write(text)

def fixture(root,conflict=False,many=0):
    root.mkdir(parents=True,exist_ok=True); ids=[f"9606.ENSP{i:05d}" for i in range(1,5)]
    gz(root/"9606.protein.info.vfixture.txt.gz","#string_protein_id\tpreferred_name\tprotein_size\tannotation\n"+"\n".join(f"{x}\tG{i}\t100\tProtein {i}" for i,x in enumerate(ids,1))+"\n")
    gz(root/"9606.protein.aliases.vfixture.txt.gz","#string_protein_id\talias\tsource\n"+f"{ids[0]}\tP00001\tUniProt_AC\n{ids[1]}\tP00002-2\tUniProt_AC\n{ids[0]}\tAMB\tUniProt_AC\n{ids[1]}\tAMB\tUniProt_AC\n{ids[2]}\tG3\tEnsembl\n")
    values=lambda score:"0 1 2 3 4 5 6 7 8 9 10 11 12 "+str(score)
    rows=[f"{ids[0]} {ids[1]} {values(900)}",f"{ids[1]} {ids[0]} {values(700 if conflict else 900)}",f"{ids[0]} {ids[2]} {values(500)}",f"{ids[0]} {ids[3]} {values(200)}",f"{ids[0]} {ids[0]} {values(999)}"]
    for i in range(many):rows.append(f"9606.X{i:05d} 9606.Y{i:05d} {values(i%1001)}")
    header="protein1 protein2 "+EVIDENCE+"\n";gz(root/"9606.protein.links.full.vfixture.txt.gz",header+"\n".join(rows)+"\n")
    gz(root/"9606.protein.physical.links.full.vfixture.txt.gz",header+f"{ids[1]} {ids[2]} {values(800)}\n")
    return ids

def test_provider_version_and_names(monkeypatch):
    provider=StringProvider();monkeypatch.setattr(provider,"_read",lambda request:json.dumps([{"string_version":"12.0","string_stable_address":"https://version-12-0.string-db.org"}]).encode())
    assert provider.detect_version()["string_version"]=="12.0"
    assert provider.filenames("12.0")[0]=="9606.protein.links.full.v12.0.txt.gz"

def test_import_mapping_queries_and_integrity(tmp_path):
    ids=fixture(tmp_path/"source");db=StringDatabase(tmp_path/"db");snap=db.install_from_directory(tmp_path/"source")
    assert db.active_snapshot()==snap and db.validate_snapshot()
    manifest=db.manifest();assert manifest["functional_edge_count"]==3 and manifest["physical_edge_count"]==1 and manifest["self_loops"]["functional"]==1
    assert len(manifest["files"])==4 and all(x["sha256"] and x["headers"] for x in manifest["files"])
    assert db.lookup_uniprot("P00001")["mapping_status"]=="mapped_unique"
    assert db.lookup_uniprot("AMB")["mapping_status"]=="ambiguous"
    assert db.lookup_uniprot("NONE")["mapping_status"]=="unmapped"
    assert db.lookup_uniprot("P00002-2")["candidate_string_ids"]==[ids[1]]
    assert db.lookup_alias("G3")["candidate_sources"]==["Ensembl"]
    assert db.lookup_alias("G1")["mapping_status"]=="unmapped"
    assert db.get_protein(ids[0])["preferred_name"]=="G1"
    assert db.get_edge(ids[1],ids[0])["combined_score"]==900
    assert [len(db.get_neighbors(ids[0],min_combined_score=x)) for x in (0,400,700)]==[3,2,1]
    assert not db.get_neighbors(ids[0],"physical") and db.get_neighbors(ids[1],"physical")
    assert len(db.get_internal_edges(ids[:3]))==2
    assert db.get_neighbors_for_many(ids[:2])[ids[0]]
    with sqlite3.connect(db.database_path()) as c:
        columns=[r[1] for r in c.execute("PRAGMA table_info(functional_edges)")]
        assert "experiments" in columns and "experiments_transferred" in columns

def test_conflicting_reciprocal_edge_fails(tmp_path):
    fixture(tmp_path/"source",True);db=StringDatabase(tmp_path/"db")
    with pytest.raises(ValueError,match="Conflicting reciprocal"):db.install_from_directory(tmp_path/"source")
    assert db.active_snapshot() is None and db.state()==DatabaseState.ERROR

def test_safe_update_failure_and_cancel(tmp_path):
    fixture(tmp_path/"a");db=StringDatabase(tmp_path/"db");active=db.install_from_directory(tmp_path/"a")
    fixture(tmp_path/"bad",True)
    with pytest.raises(ValueError):db.install_from_directory(tmp_path/"bad")
    assert db.active_snapshot()==active
    fixture(tmp_path/"cancel")
    with pytest.raises(StringDownloadCancelled):db.install_from_directory(tmp_path/"cancel",cancel_requested=lambda:True)
    assert db.active_snapshot()==active

def test_incremental_import_thousands_of_edges(tmp_path):
    fixture(tmp_path/"source",many=2500);db=StringDatabase(tmp_path/"db");db.install_from_directory(tmp_path/"source")
    assert db.manifest()["functional_edge_count"]==2503
