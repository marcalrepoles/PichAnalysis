from __future__ import annotations
import hashlib,json,os,shutil
from datetime import datetime,timezone
from pathlib import Path
from .database_registry import DatabaseState
from .databases.reactome import BASE_URL,CORE_FILES,DATA_LICENSE

def _now():return datetime.now(timezone.utc).isoformat()
def _sha(path:Path):
 h=hashlib.sha256()
 with path.open("rb") as f:
  while chunk:=f.read(1024*1024):h.update(chunk)
 return h.hexdigest()
def _rows(path:Path):return [x.split("\t") for x in path.read_text(encoding="utf-8").splitlines() if x.strip()]
def validate_reactome_core(raw:Path)->None:
 for name in CORE_FILES:
  path=raw/name
  if not path.is_file() or not path.stat().st_size:raise ValueError(f"Missing or empty Reactome file: {name}")
  rows=_rows(path)
  if name=="ReactomePathways.txt" and not any(len(r)>=3 and r[0].startswith("R-HSA-") and r[1] and "Homo sapiens" in r[2] for r in rows):raise ValueError(f"Invalid Reactome structure: {name}")
  if name=="ReactomePathwaysRelation.txt" and not any(len(r)>=2 and r[0].startswith("R-") and r[1].startswith("R-") for r in rows):raise ValueError(f"Invalid Reactome structure: {name}")
  if name in {"UniProt2Reactome.txt","NCBI2Reactome.txt"} and not any(len(r)>=2 and r[0] and r[1].startswith("R-HSA-") for r in rows):raise ValueError(f"Invalid Reactome structure: {name}")
  if name=="humanPathwaysWithDiagrams.txt" and not any(any(v.startswith("R-HSA-") for v in r) for r in rows):raise ValueError(f"Invalid Reactome structure: {name}")
  if name=="pathway2summation.txt" and not any(len(r)>=2 and r[0].startswith("R-HSA-") and r[1] for r in rows):raise ValueError(f"Invalid Reactome structure: {name}")

class ReactomeDatabase:
 def __init__(self,root:Path):
  self.root=Path(root)/"reactome"/"human";self.snapshots=self.root/"snapshots";self.active_pointer=self.root/"active_snapshot.json"
 def create_staging_snapshot(self)->Path:
  sid=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ");snap=self.snapshots/sid;(snap/"raw").mkdir(parents=True);self._write_manifest(snap,{"database":"Reactome","organism_name":"Homo sapiens","tax_id":"9606","reactome_species_prefix":"R-HSA-","snapshot_id":sid,"release_version":"unknown","status":DatabaseState.INCOMPLETE,"download_started_at":_now(),"download_completed_at":None,"source":BASE_URL,"data_license":DATA_LICENSE,"files":[]});return snap
 def install_from_directory(self,source:Path)->Path:
  snap=self.create_staging_snapshot();raw=snap/"raw"
  try:
   for name in CORE_FILES:shutil.copy2(Path(source)/name,raw/name)
   return self.finalize_snapshot(snap)
  except Exception as e:
   if self.manifest(snap).get("status")!=DatabaseState.ERROR:self.mark_error(snap,str(e))
   raise
 def finalize_snapshot(self,snap:Path)->Path:
  try:
   raw=Path(snap)/"raw";validate_reactome_core(raw);files=[{"filename":n,"size":(raw/n).stat().st_size,"sha256":_sha(raw/n),"source_url":f"{BASE_URL}/{n}"} for n in CORE_FILES]
   completed=_now();m=self.manifest(snap);m.update(status=DatabaseState.READY,download_completed_at=completed,retrieved_at=completed,files=files);self._write_manifest(snap,m);self._atomic_json(self.active_pointer,{"snapshot_id":Path(snap).name,"activated_at":completed});return Path(snap)
  except Exception as e:
   self.mark_error(Path(snap),str(e));raise
 def mark_error(self,snap:Path,message:str)->None:
  m=self.manifest(snap);m.update(status=DatabaseState.ERROR,last_error=message);self._write_manifest(snap,m)
 def mark_downloading(self,snap:Path)->None:
  m=self.manifest(snap);m.update(status=DatabaseState.DOWNLOADING);self._write_manifest(snap,m)
 def mark_incomplete(self,snap:Path)->None:
  m=self.manifest(snap);m.update(status=DatabaseState.INCOMPLETE);self._write_manifest(snap,m)
 def state(self):
  if self.is_core_ready():return DatabaseState.READY
  manifests=sorted(self.snapshots.glob("*/manifest.json"),reverse=True) if self.snapshots.exists() else []
  if not manifests:return DatabaseState.NOT_INSTALLED
  try:
   state=DatabaseState(json.loads(manifests[0].read_text(encoding="utf-8")).get("status",DatabaseState.ERROR))
   return DatabaseState.INCOMPLETE if state==DatabaseState.CANCELLED else state
  except (OSError,ValueError,json.JSONDecodeError):return DatabaseState.ERROR
 def active_snapshot(self):
  try:
   snap=self.snapshots/json.loads(self.active_pointer.read_text(encoding="utf-8"))["snapshot_id"]
   return snap if self.manifest(snap).get("status")==DatabaseState.READY else None
  except (OSError,KeyError,json.JSONDecodeError):return None
 def is_core_ready(self):return self.active_snapshot() is not None
 def manifest(self,snapshot:Path|None=None):
  try:return json.loads(((snapshot or self.active_snapshot())/"manifest.json").read_text(encoding="utf-8"))
  except (OSError,TypeError,json.JSONDecodeError):return {}
 def core_paths(self):
  snap=self.active_snapshot();return {n:snap/"raw"/n for n in CORE_FILES} if snap else {}
 def source_hashes(self):return {x["filename"]:x["sha256"] for x in self.manifest().get("files",[])}
 @staticmethod
 def _atomic_json(path,payload):
  path.parent.mkdir(parents=True,exist_ok=True);part=path.with_name(path.name+".part");part.write_text(json.dumps(payload,indent=2),encoding="utf-8");os.replace(part,path)
 def _write_manifest(self,snap,payload):self._atomic_json(snap/"manifest.json",payload)
