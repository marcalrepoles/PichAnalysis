from __future__ import annotations
import hashlib,os,time,urllib.error,urllib.request,tarfile,re,shutil
from pathlib import PurePosixPath
from pathlib import Path
BASE_URL="https://reactome.org/download/current"
CORE_FILES=("ReactomePathways.txt","ReactomePathwaysRelation.txt","UniProt2Reactome.txt","NCBI2Reactome.txt","humanPathwaysWithDiagrams.txt","pathway2summation.txt")
DATA_LICENSE="CC0";DIAGRAM_LICENSE="CC BY 4.0"
DIAGRAM_ARCHIVE="diagrams.png.tgz";DIAGRAM_URL=f"{BASE_URL}/{DIAGRAM_ARCHIVE}"
class ReactomeDownloadCancelled(RuntimeError):
 pass
class ReactomeProvider:
 def __init__(self,retries=3,timeout=60):self.retries=retries;self.timeout=timeout
 def download(self,name:str,destination:Path,progress=None,cancel_requested=None)->dict:
  error=None
  for attempt in range(self.retries+1):
   part=destination.with_name(destination.name+".part")
   try:
    if cancel_requested and cancel_requested():raise ReactomeDownloadCancelled("Reactome Core Data download was canceled.")
    req=urllib.request.Request(f"{BASE_URL}/{name}",headers={"User-Agent":"PichAnalysis/0.1"});digest=hashlib.sha256();size=0
    with urllib.request.urlopen(req,timeout=self.timeout) as response,part.open("wb") as out:
     total=int(response.headers.get("Content-Length",0));
     while True:
      if cancel_requested and cancel_requested():raise ReactomeDownloadCancelled("Reactome Core Data download was canceled.")
      chunk=response.read(1024*1024)
      if not chunk:break
      out.write(chunk);digest.update(chunk);size+=len(chunk);progress and progress(name,size,total)
    os.replace(part,destination);return {"name":name,"source":f"{BASE_URL}/{name}","file_size":size,"sha256":digest.hexdigest()}
   except ReactomeDownloadCancelled:
    part.unlink(missing_ok=True);raise
   except (OSError,urllib.error.URLError) as exc:
    error=exc;part.unlink(missing_ok=True)
    if attempt<self.retries:time.sleep(min(2**attempt,4))
  raise RuntimeError(f"Reactome download failed: {error}")
def human_rows(path:Path,id_column:int=1):
 rows=[]
 for line in path.read_text(encoding="utf-8").splitlines():
  fields=line.split("\t")
  if len(fields)>id_column and fields[id_column].startswith("R-HSA-"):rows.append(fields)
 return rows
def validate_core(root:Path)->None:
 for name in CORE_FILES:
  path=root/name
  if not path.is_file() or path.stat().st_size==0:raise ValueError(f"Missing or empty Reactome file: {name}")
def safe_extract_tar(archive:Path,destination:Path,progress=None,cancel_requested=None)->None:
 destination.mkdir(parents=True,exist_ok=True);base=destination.resolve()
 with tarfile.open(archive,"r:*") as bundle:
  members=bundle.getmembers()
  for index,member in enumerate(members,1):
   if cancel_requested and cancel_requested():raise ReactomeDownloadCancelled("Reactome diagram installation was canceled.")
   raw=member.name.replace("\\","/");parts=PurePosixPath(raw).parts;target=(destination/Path(*parts)).resolve()
   if not raw or raw.startswith("/") or re.match(r"^[A-Za-z]:",raw) or ".." in parts or (target!=base and base not in target.parents):raise ValueError(f"Unsafe archive member: {member.name}")
   if member.issym() or member.islnk() or not (member.isfile() or member.isdir()):raise ValueError(f"Unsafe archive member type: {member.name}")
   if member.isdir():target.mkdir(parents=True,exist_ok=True)
   else:
    target.parent.mkdir(parents=True,exist_ok=True);source=bundle.extractfile(member)
    if source is None:raise ValueError(f"Unreadable archive member: {member.name}")
    with source,target.open("wb") as output:shutil.copyfileobj(source,output,1024*1024)
   if progress:progress(index,len(members),member.name)
