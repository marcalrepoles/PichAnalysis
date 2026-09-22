"""Explicit real-data STRING smoke; excluded from normal pytest collection."""
from __future__ import annotations
import argparse, json, tempfile
from pathlib import Path
from pichanalysis.core.string_database import StringDatabase

def main(destination=None):
    root=Path(destination) if destination else Path(tempfile.mkdtemp(prefix="pichanalysis-string-"))
    db=StringDatabase(root);snap=db.download_and_build(progress=lambda p:print(p,flush=True));m=db.manifest(snap)
    # Deterministically use the first auditable UniProt alias found in the local database.
    import sqlite3
    with sqlite3.connect(db.database_path()) as c:
        alias,string_id=c.execute("SELECT alias,string_protein_id FROM aliases WHERE source='UniProt_AC' ORDER BY alias,string_protein_id LIMIT 1").fetchone()
        node=c.execute("SELECT protein_a FROM functional_edges ORDER BY protein_a LIMIT 1").fetchone()[0]
    report={"temporary_path":str(root),"snapshot_id":snap.name,**m,"mapping":db.lookup_uniprot(alias),"protein":db.get_protein(string_id),"functional_neighbors":len(db.get_neighbors(node,"functional")),"physical_neighbors":len(db.get_neighbors(node,"physical")),"threshold_0":len(db.get_neighbors(node,"functional",0)),"threshold_700":len(db.get_neighbors(node,"functional",700)),"internal_edges":len(db.get_internal_edges([node,*[e["protein_b"] if e["protein_a"]==node else e["protein_a"] for e in db.get_neighbors(node)[:5]]]))}
    assert report["functional_neighbors"]>0 and report["threshold_700"]<=report["threshold_0"]
    db.provider=None;assert db.lookup_uniprot(alias) and db.get_protein(string_id) and db.get_neighbors(node)
    print(json.dumps(report,indent=2))
if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--destination");main(parser.parse_args().destination)

