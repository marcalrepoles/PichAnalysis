"""Managed, immutable human Gene Ontology ontology and annotation snapshots."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import shutil
import sqlite3
import urllib.request
import uuid
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from .database_registry import DatabaseState

ONTOLOGY_URL = "https://current.geneontology.org/ontology/go-basic.obo"
ANNOTATION_URL = "https://current.geneontology.org/annotations/gaf/HUMAN-uniprot.gaf.gz"
GO_FILES = {"go-basic.obo": ONTOLOGY_URL, "HUMAN-uniprot.gaf.gz": ANNOTATION_URL}
GAF_COLUMNS = ("database", "db_object_id", "gene_symbol", "qualifier", "go_id", "reference",
    "evidence_code", "with_from", "aspect", "db_object_name", "db_object_synonym",
    "db_object_type", "taxon", "date", "assigned_by", "annotation_extension", "gene_product_form_id")


class GeneOntologyError(RuntimeError): pass
class GeneOntologyCancelled(GeneOntologyError): pass


def _now(): return datetime.now(timezone.utc).isoformat()


def _sha(path):
    digest=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b""):digest.update(chunk)
    return digest.hexdigest()


def parse_obo(path):
    """Parse official OBO term stanzas, preserving is_a versus part_of."""
    header={};terms=[];relations=[];term=None
    def flush():
        if not term or not term.get("go_id"):return
        terms.append({"go_id":term["go_id"],"name":term.get("name",""),
            "namespace":term.get("namespace",""),"is_obsolete":term.get("is_obsolete","false")})
        relations.extend(term.get("relations",[]))
    with Path(path).open(encoding="utf-8-sig") as stream:
        for raw in stream:
            line=raw.strip()
            if line=="[Term]":flush();term={"relations":[]};continue
            if line.startswith("["):
                flush();term=None;continue
            if not line or line.startswith("!") or ": " not in line:continue
            key,value=line.split(": ",1)
            if term is None:
                if key in {"data-version","date","format-version"}:header[key]=value
                continue
            if key=="id":term["go_id"]=value
            elif key in {"name","namespace","is_obsolete"}:term[key]=value
            elif key=="is_a":term["relations"].append({"child_go_id":term.get("go_id",""),
                "parent_go_id":value.split()[0],"relation":"is_a"})
            elif key=="relationship" and value.startswith("part_of "):
                term["relations"].append({"child_go_id":term.get("go_id",""),
                    "parent_go_id":value.split()[1],"relation":"part_of"})
    flush()
    if not terms or not relations or any(not row["child_go_id"] for row in relations):
        raise GeneOntologyError("GO ontology has no valid terms or relations.")
    return pd.DataFrame(terms),pd.DataFrame(relations).drop_duplicates(),header


def _gaf_stream(path):
    return gzip.open(path,"rt",encoding="utf-8") if str(path).endswith(".gz") else Path(path).open(encoding="utf-8")


def parse_gaf(path,output_csv=None,cancel_requested=None):
    """Stream GAF 2.2, retaining audit fields and requiring human taxonomy."""
    cancel_requested=cancel_requested or (lambda:False)
    metadata={};rows=[];count=0;not_count=0;excluded_taxa={}
    writer=None;output=None
    if output_csv:
        output=Path(output_csv).open("w",encoding="utf-8",newline="")
        writer=csv.DictWriter(output,fieldnames=GAF_COLUMNS);writer.writeheader()
    try:
        with _gaf_stream(path) as stream:
            for line in stream:
                if count%10000==0 and cancel_requested():raise GeneOntologyCancelled("Gene Ontology processing was canceled.")
                if line.startswith("!"):
                    value=line[1:].strip()
                    if ":" in value:
                        key,rest=value.split(":",1)
                        metadata[key.strip().lower().replace(" ","_").replace("-","_")]=rest.strip()
                    continue
                if not line.strip():continue
                parts=line.rstrip("\r\n").split("\t")
                if len(parts)!=len(GAF_COLUMNS):raise GeneOntologyError(f"Expected 17 GAF columns, found {len(parts)}.")
                row=dict(zip(GAF_COLUMNS,parts))
                if row["taxon"]!="taxon:9606":
                    excluded_taxa[row["taxon"]]=excluded_taxa.get(row["taxon"],0)+1
                    continue
                if not row["go_id"].startswith("GO:"):raise GeneOntologyError("Human GAF contains an invalid GO ID.")
                if "NOT" in row["qualifier"].split("|"):not_count+=1
                if writer:writer.writerow(row)
                else:rows.append(row)
                count+=1
    finally:
        if output:output.close()
    if count==0:raise GeneOntologyError("Human GAF contains no annotations.")
    metadata.update(annotation_count=count,not_count=not_count,excluded_taxa=excluded_taxa)
    return (pd.DataFrame(rows,columns=GAF_COLUMNS) if not writer else None),metadata


class GeneOntologyProvider:
    def download(self,url,destination,progress=None,cancel_requested=None):
        cancel_requested=cancel_requested or (lambda:False)
        destination=Path(destination);destination.parent.mkdir(parents=True,exist_ok=True)
        part=destination.with_name(destination.name+".part")
        request=urllib.request.Request(url,headers={"User-Agent":"PichAnalysis/GO-snapshot"})
        with urllib.request.urlopen(request,timeout=90) as response,part.open("wb") as target:
            headers={"Last-Modified":response.headers.get("Last-Modified",""),"ETag":response.headers.get("ETag","")}
            total=int(response.headers.get("Content-Length") or 0);done=0
            while chunk:=response.read(1024*1024):
                if cancel_requested():raise GeneOntologyCancelled("Gene Ontology download was canceled.")
                target.write(chunk);done+=len(chunk)
                if progress:progress({"message":f"Downloading {destination.name}","bytes_downloaded":done,"bytes_total":total})
        os.replace(part,destination)
        return {"source_url":url,"retrieved_at":_now(),"last_modified":headers["Last-Modified"],
            "etag":headers["ETag"],"size":destination.stat().st_size,"sha256":_sha(destination)}


class GeneOntologyDatabase:
    def __init__(self,root,provider=None):
        self.root=Path(root)/"gene_ontology"/"human"
        self.snapshots=self.root/"snapshots";self.active_pointer=self.root/"active_snapshot.json"
        self.provider=provider or GeneOntologyProvider()

    def active_snapshot(self):
        try:
            sid=json.loads(self.active_pointer.read_text(encoding="utf-8"))["snapshot_id"]
            snap=self.snapshots/sid
            return snap if snap.is_dir() and self.manifest(snap).get("status")==DatabaseState.READY else None
        except (OSError,ValueError,KeyError):return None

    def is_ready(self):return self.active_snapshot() is not None

    def manifest(self,snapshot=None):
        try:return json.loads(((snapshot or self.active_snapshot())/"manifest.json").read_text(encoding="utf-8"))
        except (OSError,TypeError,ValueError):return {}

    def database_path(self):
        snap=self.active_snapshot()
        return snap/"gene_ontology.sqlite" if snap else None

    @staticmethod
    def _atomic_json(path,payload):
        path.parent.mkdir(parents=True,exist_ok=True)
        part=path.with_name(path.name+".part");part.write_text(json.dumps(payload,indent=2),encoding="utf-8")
        os.replace(part,path)

    def _new_snapshot(self):
        sid=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")+"-"+uuid.uuid4().hex[:8]
        snap=self.snapshots/sid
        (snap/"raw").mkdir(parents=True);(snap/"tables").mkdir()
        return snap

    def build_from_directory(self,source,progress=None,cancel_requested=None):
        """Fixture/local import; source files copied byte-for-byte."""
        source=Path(source)
        return self._build(lambda name,dest:shutil.copy2(source/name,dest),progress,cancel_requested,local=True)

    def download_and_build(self,progress=None,cancel_requested=None):
        return self._build(lambda name,dest:self.provider.download(GO_FILES[name],dest,progress,cancel_requested),
                           progress,cancel_requested,local=False)

    def _build(self,acquire,progress,cancel_requested,local):
        cancel_requested=cancel_requested or (lambda:False)
        emit=lambda message:progress({"message":message}) if progress else None
        snap=self._new_snapshot()
        try:
            info={}
            for name,url in GO_FILES.items():
                if cancel_requested():raise GeneOntologyCancelled("Gene Ontology processing was canceled.")
                emit(f"Retrieving {name}...")
                path=snap/"raw"/name
                response=acquire(name,path)
                info[name]=response if isinstance(response,dict) else {"source_url":url if not local else str(Path(response)),
                    "retrieved_at":_now(),"last_modified":"","etag":"","size":path.stat().st_size,"sha256":_sha(path)}
            emit("Parsing GO ontology...")
            terms,relations,obo_meta=parse_obo(snap/"raw"/"go-basic.obo")
            terms.to_csv(snap/"tables"/"go_terms.csv",index=False)
            relations.to_csv(snap/"tables"/"go_relations.csv",index=False)
            if cancel_requested():raise GeneOntologyCancelled("Gene Ontology processing was canceled.")
            emit("Parsing human GO annotations...")
            _,gaf_meta=parse_gaf(snap/"raw"/"HUMAN-uniprot.gaf.gz",snap/"tables"/"go_annotations.csv",cancel_requested)
            emit("Indexing Gene Ontology database...")
            sqlite=snap/"gene_ontology.sqlite"
            with sqlite3.connect(sqlite) as conn:
                terms.to_sql("terms",conn,index=False,if_exists="replace")
                relations.to_sql("relations",conn,index=False,if_exists="replace")
                for chunk in pd.read_csv(snap/"tables"/"go_annotations.csv",dtype=str,chunksize=25000,keep_default_na=False):
                    if cancel_requested():raise GeneOntologyCancelled("Gene Ontology processing was canceled.")
                    chunk.to_sql("annotations",conn,index=False,if_exists="append")
                pd.DataFrame([{"key":"schema_version","value":"1"}]).to_sql("metadata",conn,index=False)
                for table,column in (("terms","go_id"),("relations","child_go_id"),("relations","parent_go_id"),
                    ("annotations","go_id"),("annotations","gene_symbol"),("annotations","db_object_id"),("annotations","taxon")):
                    conn.execute(f'CREATE INDEX "idx_{table}_{column}" ON "{table}" ("{column}")')
            manifest={"database":"Gene Ontology","snapshot_id":snap.name,"status":DatabaseState.INCOMPLETE,
                "created_at":_now(),"organism":"Homo sapiens","tax_id":"9606",
                "ontology_version":obo_meta.get("data-version",""),"ontology_date":obo_meta.get("date",""),
                "annotation_version":gaf_meta.get("gaf_version",gaf_meta.get("generated_by","")),
                "annotation_date":gaf_meta.get("date_generated",""),"sources":info,
                "counts":{"terms":len(terms),"relations":len(relations),"annotations":gaf_meta["annotation_count"],
                          "not_annotations":gaf_meta["not_count"],"excluded_nonhuman_annotations":sum(gaf_meta["excluded_taxa"].values())},
                "excluded_taxa":gaf_meta["excluded_taxa"],
                "normalized_hashes":{name:_sha(snap/"tables"/name) for name in ("go_terms.csv","go_relations.csv","go_annotations.csv")},
                "sqlite_sha256":_sha(sqlite)}
            self._atomic_json(snap/"manifest.json",manifest)
            emit("Validating Gene Ontology snapshot...")
            self.validate_snapshot(snap,allow_incomplete=True)
            if cancel_requested():raise GeneOntologyCancelled("Gene Ontology processing was canceled.")
            manifest.update(status=DatabaseState.READY,completed_at=_now())
            self._atomic_json(snap/"manifest.json",manifest)
            self._atomic_json(self.active_pointer,{"snapshot_id":snap.name,"activated_at":_now()})
            return snap
        except Exception as error:
            self._atomic_json(snap/"manifest.json",{"database":"Gene Ontology","snapshot_id":snap.name,
                "status":DatabaseState.CANCELLED if isinstance(error,GeneOntologyCancelled) else DatabaseState.ERROR,
                "last_error":str(error)})
            raise

    def validate_snapshot(self,snapshot=None,allow_incomplete=False):
        snap=Path(snapshot or self.active_snapshot() or "")
        m=self.manifest(snap)
        if m.get("status") not in ({DatabaseState.INCOMPLETE,DatabaseState.READY} if allow_incomplete else {DatabaseState.READY}):
            raise GeneOntologyError("Gene Ontology snapshot is not Ready.")
        if m.get("tax_id")!="9606":raise GeneOntologyError("Gene Ontology snapshot is not human.")
        for name in GO_FILES:
            path=snap/"raw"/name;info=m.get("sources",{}).get(name,{})
            if not path.is_file() or _sha(path)!=info.get("sha256") or path.stat().st_size!=info.get("size"):
                raise GeneOntologyError(f"Gene Ontology raw source is invalid: {name}")
        for name,expected in m.get("normalized_hashes",{}).items():
            path=snap/"tables"/name
            if not path.is_file() or _sha(path)!=expected:raise GeneOntologyError(f"Gene Ontology table is invalid: {name}")
        sqlite=snap/"gene_ontology.sqlite"
        if not sqlite.is_file() or _sha(sqlite)!=m.get("sqlite_sha256"):
            raise GeneOntologyError("Gene Ontology SQLite hash mismatch.")
        with sqlite3.connect(sqlite) as conn:
            tables={row[0] for row in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
            if not {"terms","relations","annotations","metadata"}.issubset(tables):
                raise GeneOntologyError("Gene Ontology SQLite tables are missing.")
            for name in ("terms","relations","annotations"):
                if conn.execute(f'SELECT COUNT(*) FROM "{name}"').fetchone()[0]!=m["counts"][name]:
                    raise GeneOntologyError(f"Gene Ontology row count mismatch: {name}")
            if conn.execute("SELECT COUNT(*) FROM annotations WHERE taxon!='taxon:9606'").fetchone()[0]:
                raise GeneOntologyError("Gene Ontology contains non-human annotations.")
        return True

    def _query(self,sql,params=(),snapshot=None):
        path=(Path(snapshot)/"gene_ontology.sqlite") if snapshot else self.database_path()
        if not path or not path.is_file():raise GeneOntologyError("Gene Ontology snapshot is not installed.")
        with sqlite3.connect(path) as conn:return pd.read_sql_query(sql,conn,params=params)

    def get_term(self,go_id,snapshot=None):return self._query("SELECT * FROM terms WHERE go_id=?",(go_id,),snapshot)

    def get_descendants(self,go_id,relations=("is_a","part_of"),snapshot=None):
        allowed=tuple(value for value in relations if value in {"is_a","part_of"})
        if not allowed:return pd.DataFrame(columns=("go_id","distance"))
        edges=self._query("SELECT child_go_id,parent_go_id,relation FROM relations",snapshot=snapshot)
        children={}
        for row in edges.itertuples(index=False):
            if row.relation in allowed:children.setdefault(row.parent_go_id,set()).add(row.child_go_id)
        seen={go_id:0};queue=[go_id]
        for parent in queue:
            for child in sorted(children.get(parent,())):
                if child not in seen:seen[child]=seen[parent]+1;queue.append(child)
        return pd.DataFrame(({"go_id":key,"distance":distance} for key,distance in seen.items()))

    def get_annotations_for_terms(self,go_ids,snapshot=None):
        ids=sorted(set(go_ids))
        if not ids:return pd.DataFrame(columns=GAF_COLUMNS)
        frames=[]
        for offset in range(0,len(ids),800):
            part=ids[offset:offset+800];holders=",".join("?" for _ in part)
            frames.append(self._query(f"SELECT * FROM annotations WHERE go_id IN ({holders})",part,snapshot))
        return pd.concat(frames,ignore_index=True) if frames else pd.DataFrame(columns=GAF_COLUMNS)
