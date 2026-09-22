from __future__ import annotations

import csv, gzip, hashlib, json, os, shutil, sqlite3
from datetime import datetime, timezone
from pathlib import Path
from .database_registry import DatabaseState
from .databases.string import StringDownloadCancelled, StringProvider

REQUIRED_TABLES={"proteins","aliases","functional_edges","physical_edges","metadata"}
def _now():return datetime.now(timezone.utc).isoformat()
def _sha(path):
 h=hashlib.sha256()
 with Path(path).open("rb") as f:
  while chunk:=f.read(1024*1024):h.update(chunk)
 return h.hexdigest()
def _atomic_json(path,payload):
 path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);part=path.with_name(path.name+".part");part.write_text(json.dumps(payload,indent=2,ensure_ascii=False),encoding="utf-8");os.replace(part,path)
def _q(name):return '"'+str(name).replace('"','""')+'"'

class StringDatabase:
 def __init__(self,root,provider=None):self.root=Path(root)/"string"/"human";self.snapshots=self.root/"snapshots";self.active_pointer=self.root/"active_snapshot.json";self.provider=provider or StringProvider()
 def create_staging_snapshot(self,version=None,stable_address=None):
  meta={"string_version":version,"string_stable_address":stable_address} if version else self.provider.detect_version();sid=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ");snap=self.snapshots/sid
  for folder in ("raw","tables"):(snap/folder).mkdir(parents=True,exist_ok=True)
  manifest={"database":"STRING","organism":"Homo sapiens","tax_id":"9606","string_version":meta["string_version"],"stable_address":meta["string_stable_address"],"snapshot_id":sid,"status":DatabaseState.INCOMPLETE,"created_at":_now(),"license":"CC BY 4.0","files":[],"edge_directionality":"undirected","duplicate_pair_policy":"canonical pair; identical AB/BA merged; conflicts fail","self_loop_policy":"excluded and counted","score_storage":"official raw integer scores; no recomputation"};_atomic_json(snap/"manifest.json",manifest);return snap
 def manifest(self,snapshot=None):
  target=snapshot or self.active_snapshot()
  try:return json.loads((target/"manifest.json").read_text(encoding="utf-8")) if target else {}
  except (OSError,ValueError):return {}
 def _write(self,snap,m):m["updated_at"]=_now();_atomic_json(snap/"manifest.json",m)
 def mark_downloading(self,snap):m=self.manifest(snap);m["status"]=DatabaseState.DOWNLOADING;self._write(snap,m)
 def mark_incomplete(self,snap):m=self.manifest(snap);m["status"]=DatabaseState.CANCELLED;self._write(snap,m)
 def mark_error(self,snap,error):m=self.manifest(snap);m.update(status=DatabaseState.ERROR,last_error=str(error));self._write(snap,m)
 def state(self):
  if self.active_snapshot():return DatabaseState.READY
  manifests=sorted(self.snapshots.glob("*/manifest.json"),reverse=True) if self.snapshots.exists() else []
  return DatabaseState(self.manifest(manifests[0].parent).get("status",DatabaseState.ERROR)) if manifests else DatabaseState.NOT_INSTALLED
 def active_snapshot(self):
  try:
   snap=self.snapshots/json.loads(self.active_pointer.read_text(encoding="utf-8"))["snapshot_id"]
   return snap if snap.is_dir() and self.manifest(snap).get("status")==DatabaseState.READY else None
  except (OSError,ValueError,KeyError):return None
 def is_ready(self):return self.active_snapshot() is not None
 def version(self):return self.manifest().get("string_version")
 def database_path(self):
  snap=self.active_snapshot();return snap/"string_network.sqlite" if snap else None
 def source_hashes(self):return {x["filename"]:x["sha256"] for x in self.manifest().get("files",[])}
 def download_and_build(self,progress=None,cancel_requested=None):
  snap=self.create_staging_snapshot();self.mark_downloading(snap);m=self.manifest(snap)
  try:
   files=[]
   for i,name in enumerate(self.provider.filenames(m["string_version"]),1):
    if cancel_requested and cancel_requested():raise StringDownloadCancelled("STRING download was canceled.")
    files.append(self.provider.download(m["stable_address"],name,snap/"raw"/name,lambda n,d,t:self._emit(progress,"Downloading",n,i,4,d,t),cancel_requested))
   m["files"]=files;self._write(snap,m);return self.finalize_snapshot(snap,progress,cancel_requested)
  except StringDownloadCancelled:self.mark_incomplete(snap);raise
  except Exception as exc:self.mark_error(snap,exc);raise
 def install_from_directory(self,source,version="fixture",stable_address="fixture://official",progress=None,cancel_requested=None):
  snap=self.create_staging_snapshot(version,stable_address);self.mark_downloading(snap);m=self.manifest(snap);files=[]
  try:
   for src in sorted(Path(source).glob("*.txt.gz")):
    dst=snap/"raw"/src.name;shutil.copyfile(src,dst);files.append({"filename":src.name,"source":str(src),"size":dst.stat().st_size,"sha256":_sha(dst)})
   m["files"]=files;self._write(snap,m);return self.finalize_snapshot(snap,progress,cancel_requested)
  except StringDownloadCancelled:self.mark_incomplete(snap);raise
  except Exception as exc:self.mark_error(snap,exc);raise
 def finalize_snapshot(self,snap,progress=None,cancel_requested=None):
  m=self.manifest(snap);raw=snap/"raw";names={p.name:p for p in raw.glob("*.txt.gz")}
  def choose(token):
   hits=[p for n,p in names.items() if token in n]
   if len(hits)!=1:raise ValueError(f"Expected exactly one STRING {token} file.")
   return hits[0]
  functional=choose("protein.links.full");physical=choose("protein.physical.links.full");info=choose("protein.info");aliases=choose("protein.aliases")
  headers={p.name:self._header(p) for p in (functional,physical,info,aliases)}
  for item in m["files"]:item["headers"]=headers.get(item["filename"],[])
  db=snap/"string_network.sqlite";conn=sqlite3.connect(db)
  try:
   conn.execute("PRAGMA journal_mode=OFF");conn.execute("PRAGMA synchronous=OFF")
   self._emit(progress,"Parsing","Parsing protein information...");pc=self._parse_proteins(conn,info,snap/"tables/proteins.csv",cancel_requested)
   self._emit(progress,"Parsing","Parsing aliases...");ac=self._parse_aliases(conn,aliases,snap/"tables/aliases.csv",cancel_requested)
   self._emit(progress,"Indexing","Indexing functional network...");fc,fs=self._parse_edges(conn,functional,"functional_edges",cancel_requested)
   self._emit(progress,"Indexing","Indexing physical network...");hc,hs=self._parse_edges(conn,physical,"physical_edges",cancel_requested)
   alias_sources=[r[0] for r in conn.execute("SELECT DISTINCT source FROM aliases ORDER BY source")]
   uniprot_alias_sources=[source for source in alias_sources if source == "UniProt_AC"]
   conn.execute("CREATE TABLE metadata(key TEXT PRIMARY KEY,value TEXT NOT NULL)")
   for key,value in (("string_version",m["string_version"]),("tax_id","9606")):conn.execute("INSERT INTO metadata VALUES(?,?)",(key,str(value)))
   conn.executescript("CREATE INDEX idx_alias ON aliases(alias);CREATE INDEX idx_alias_protein ON aliases(string_protein_id);CREATE INDEX idx_protein_id ON proteins(string_protein_id);CREATE INDEX idx_fe_a ON functional_edges(protein_a);CREATE INDEX idx_fe_b ON functional_edges(protein_b);CREATE INDEX idx_fe_score ON functional_edges(combined_score);CREATE INDEX idx_pe_a ON physical_edges(protein_a);CREATE INDEX idx_pe_b ON physical_edges(protein_b);CREATE INDEX idx_pe_score ON physical_edges(combined_score);");conn.commit()
  finally:conn.close()
  m.update(status=DatabaseState.READY,completed_at=_now(),retrieved_at=_now(),protein_count=pc,alias_count=ac,functional_edge_count=fc,physical_edge_count=hc,self_loops={"functional":fs,"physical":hs},sqlite_sha256=_sha(db),sqlite_size=db.stat().st_size,alias_sources=alias_sources,uniprot_alias_sources=uniprot_alias_sources);self._write(snap,m);self.validate_snapshot(snap);_atomic_json(self.active_pointer,{"snapshot_id":snap.name,"activated_at":_now()});return snap
 @staticmethod
 def _header(path):
  with gzip.open(path,"rt",encoding="utf-8",newline="") as f:
   line=f.readline().rstrip("\r\n");return line.split("\t") if "\t" in line else line.split()
 @staticmethod
 def _rows(path):
  f=gzip.open(path,"rt",encoding="utf-8",newline="");first=f.readline().rstrip("\r\n");tab="\t" in first;header=first.split("\t") if tab else first.split()
  def iterator():
   try:
    for line in f:
     values=line.rstrip("\r\n").split("\t") if tab else line.split()
     if values:yield dict(zip(header,values))
   finally:f.close()
  return header,iterator()
 def _parse_proteins(self,conn,path,csv_path,cancel):
  h,rows=self._rows(path);idcol="#string_protein_id" if "#string_protein_id" in h else "string_protein_id" if "string_protein_id" in h else None
  if not idcol:raise ValueError("Invalid STRING protein.info header.")
  cols=[idcol,*[x for x in ("preferred_name","protein_size","annotation") if x in h]];conn.execute("CREATE TABLE proteins("+",".join(_q("string_protein_id" if x==idcol else x)+" TEXT" for x in cols)+", PRIMARY KEY(string_protein_id))");count=0
  with csv_path.open("w",encoding="utf-8",newline="") as out:
   w=csv.writer(out);w.writerow(["string_protein_id",*[x for x in cols if x!=idcol]])
   for row in rows:
    if cancel and cancel():raise StringDownloadCancelled("STRING indexing was canceled.")
    vals=[row.get(x,"") for x in cols];conn.execute("INSERT OR REPLACE INTO proteins VALUES("+",".join("?"*len(vals))+")",vals);w.writerow(vals);count+=1
  return count
 def _parse_aliases(self,conn,path,csv_path,cancel):
  h,rows=self._rows(path);idcol="#string_protein_id" if "#string_protein_id" in h else "string_protein_id"
  if not {idcol,"alias","source"}.issubset(h):raise ValueError("Invalid STRING protein.aliases header.")
  conn.execute("CREATE TABLE aliases(string_protein_id TEXT,alias TEXT,source TEXT)");count=0
  with csv_path.open("w",encoding="utf-8",newline="") as out:
   w=csv.writer(out);w.writerow(["string_protein_id","alias","source"])
   for row in rows:
    if cancel and cancel():raise StringDownloadCancelled("STRING indexing was canceled.")
    vals=(row[idcol],row["alias"],row["source"]);conn.execute("INSERT INTO aliases VALUES(?,?,?)",vals);w.writerow(vals);count+=1
  return count
 def _parse_edges(self,conn,path,table,cancel):
  h,rows=self._rows(path);a="#protein1" if "#protein1" in h else "protein1";b="protein2"
  if a not in h or b not in h or "combined_score" not in h:raise ValueError(f"Invalid STRING edge header in {path.name}.")
  scores=[x for x in h if x not in (a,b)];conn.execute(f"CREATE TABLE {table}(protein_a TEXT,protein_b TEXT,"+",".join(_q(x)+" INTEGER" for x in scores)+",PRIMARY KEY(protein_a,protein_b))");insert=f"INSERT INTO {table} VALUES("+",".join("?"*(2+len(scores)))+")";count=loops=0
  for i,row in enumerate(rows):
   if cancel and i%1000==0 and cancel():raise StringDownloadCancelled("STRING indexing was canceled.")
   x,y=row[a],row[b]
   if x==y:loops+=1;continue
   pair=tuple(sorted((x,y)));vals=tuple(int(row[s]) for s in scores)
   try:conn.execute(insert,(*pair,*vals));count+=1
   except sqlite3.IntegrityError:
    existing=conn.execute(f"SELECT "+",".join(_q(s) for s in scores)+f" FROM {table} WHERE protein_a=? AND protein_b=?",pair).fetchone()
    if tuple(existing)!=vals:raise ValueError(f"Conflicting reciprocal STRING edge for {pair[0]} and {pair[1]}.")
   if count and count%100000==0:conn.commit()
  return count,loops
 def validate_snapshot(self,snapshot=None):
  snap=snapshot or self.active_snapshot();m=self.manifest(snap)
  if not snap or m.get("status")!=DatabaseState.READY:raise ValueError("STRING snapshot is not Ready.")
  if len(m.get("files",[]))!=4:raise ValueError("STRING snapshot does not contain four raw files.")
  for item in m["files"]:
   p=snap/"raw"/item["filename"]
   if not p.is_file() or _sha(p)!=item["sha256"]:raise ValueError(f"STRING raw file hash mismatch: {item['filename']}")
  db=snap/"string_network.sqlite"
  if _sha(db)!=m.get("sqlite_sha256"):raise ValueError("STRING SQLite hash mismatch.")
  with sqlite3.connect(db) as c:
   tables={r[0] for r in c.execute("SELECT name FROM sqlite_master WHERE type='table'")}
   if not REQUIRED_TABLES.issubset(tables):raise ValueError("STRING SQLite is missing required tables.")
   indexes={r[0] for r in c.execute("SELECT name FROM sqlite_master WHERE type='index'")}
   required_indexes={"idx_alias","idx_alias_protein","idx_protein_id","idx_fe_a","idx_fe_b","idx_fe_score","idx_pe_a","idx_pe_b","idx_pe_score"}
   if not required_indexes.issubset(indexes):raise ValueError("STRING SQLite is missing required indexes.")
   for table,key in (("proteins","protein_count"),("aliases","alias_count"),("functional_edges","functional_edge_count"),("physical_edges","physical_edge_count")):
    if c.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]!=m[key]:raise ValueError(f"STRING {table} count mismatch.")
  return True
 def _connect(self):
  path=self.database_path()
  if not path:raise RuntimeError("STRING database is not installed.")
  c=sqlite3.connect(path);c.row_factory=sqlite3.Row;return c
 def lookup_alias(self,identifier,input_type="alias",sources=None):
  with self._connect() as c:rows=c.execute("SELECT string_protein_id,source FROM aliases WHERE alias=? ORDER BY string_protein_id,source",(identifier,)).fetchall()
  if sources is not None:rows=[r for r in rows if r["source"] in sources]
  ids=sorted({r["string_protein_id"] for r in rows});return {"input_identifier":identifier,"input_type":input_type,"candidate_string_ids":ids,"candidate_sources":sorted({r["source"] for r in rows}),"mapping_status":"unmapped" if not ids else "mapped_unique" if len(ids)==1 else "ambiguous"}
 def lookup_uniprot(self,accession):
  sources=self.manifest().get("uniprot_alias_sources",["UniProt_AC"])
  return self.lookup_alias(accession,"uniprot",sources)
 def get_protein(self,string_id):
  with self._connect() as c:r=c.execute("SELECT * FROM proteins WHERE string_protein_id=?",(string_id,)).fetchone();return dict(r) if r else None
 def get_edge(self,a,b,network_type="functional"):
  table=self._table(network_type);a,b=sorted((a,b))
  with self._connect() as c:r=c.execute(f"SELECT * FROM {table} WHERE protein_a=? AND protein_b=?",(a,b)).fetchone();return dict(r) if r else None
 def get_neighbors(self,string_id,network_type="functional",min_combined_score=0):
  table=self._table(network_type)
  with self._connect() as c:rows=c.execute(f"SELECT * FROM {table} WHERE combined_score>=? AND (protein_a=? OR protein_b=?) ORDER BY combined_score DESC,protein_a,protein_b",(int(min_combined_score),string_id,string_id)).fetchall();return [dict(r) for r in rows]
 def get_neighbors_for_many(self,string_ids,network_type="functional",min_combined_score=0):return {x:self.get_neighbors(x,network_type,min_combined_score) for x in tuple(dict.fromkeys(string_ids))}
 def get_internal_edges(self,string_ids,network_type="functional",min_combined_score=0):
  ids=tuple(dict.fromkeys(string_ids))
  if not ids:return []
  table=self._table(network_type);marks=",".join("?"*len(ids))
  with self._connect() as c:rows=c.execute(f"SELECT * FROM {table} WHERE combined_score>=? AND protein_a IN ({marks}) AND protein_b IN ({marks}) ORDER BY protein_a,protein_b",(int(min_combined_score),*ids,*ids)).fetchall();return [dict(r) for r in rows]
 @staticmethod
 def _table(kind):
  if kind not in {"functional","physical"}:raise ValueError("network_type must be 'functional' or 'physical'.")
  return kind+"_edges"
 @staticmethod
 def _emit(callback,stage,message,file_index=0,file_total=0,done=0,total=0):
  if callback:callback({"stage":stage,"message":message,"filename":message,"file_index":file_index,"file_total":file_total,"bytes_downloaded":done,"bytes_total":total})
