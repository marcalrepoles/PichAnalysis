from __future__ import annotations

import hashlib
import json
import os
import shutil
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Iterable

import pandas as pd

from .database_registry import DatabaseState

API_BASE = "https://www.ebi.ac.uk/interpro/api/"
LAYERS = ("interpro", "pfam")
TABLES = (
    "proteins.csv", "protein_interpro_membership.csv", "interpro_locations.csv",
    "protein_pfam_membership.csv", "pfam_locations.csv", "interpro_entries.csv",
    "pfam_entries.csv", "interpro_hierarchy.csv",
)


class InterProError(RuntimeError): pass
class InterProCancelled(InterProError): pass
class InterProRequestError(InterProError): pass
class InvalidAnnotationSet(InterProError): pass


def _now() -> str: return datetime.now(timezone.utc).isoformat()
def _sha_bytes(data: bytes) -> str: return hashlib.sha256(data).hexdigest()
def _sha_file(path: Path) -> str:
    digest=hashlib.sha256()
    with path.open("rb") as stream:
        while chunk:=stream.read(1024*1024):digest.update(chunk)
    return digest.hexdigest()
def _safe_accession(accession: str) -> str:
    return urllib.parse.quote(accession, safe="").replace("%", "_")
def _atomic_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True,exist_ok=True);part=path.with_name(path.name+".part")
    part.write_text(json.dumps(value,indent=2,ensure_ascii=False),encoding="utf-8");os.replace(part,path)


class InterProProvider:
    def __init__(self, *, base_url: str=API_BASE, retries: int=3, timeout: int=60,
                 backoff: float=.5, opener=None):
        self.base_url=base_url.rstrip("/")+"/";self.retries=retries;self.timeout=timeout
        self.backoff=backoff;self.opener=opener or urllib.request.urlopen

    def fetch(self, url: str, cancel_requested: Callable[[],bool]|None=None) -> tuple[dict,int,bytes]:
        error=None
        for attempt in range(self.retries+1):
            if cancel_requested and cancel_requested():raise InterProCancelled("InterPro/Pfam acquisition was canceled.")
            try:
                request=urllib.request.Request(url,headers={"Accept":"application/json","User-Agent":"PichAnalysis/0.1"})
                with self.opener(request,timeout=self.timeout) as response:
                    raw=response.read();return json.loads(raw.decode("utf-8")),int(getattr(response,"status",200)),raw
            except urllib.error.HTTPError as exc:
                if exc.code==404:return {"count":0,"next":None,"results":[]},404,b"{\"count\":0,\"next\":null,\"results\":[]}"
                error=exc
                if exc.code not in {429,500,502,503,504}:break
                retry_after=exc.headers.get("Retry-After") if exc.headers else None
                delay=float(retry_after) if retry_after and retry_after.isdigit() else self.backoff*(2**attempt)
            except (OSError,urllib.error.URLError,json.JSONDecodeError) as exc:
                error=exc;delay=self.backoff*(2**attempt)
            if attempt<self.retries:time.sleep(min(delay,8))
        raise InterProRequestError(f"Official InterPro API request failed after {self.retries+1} attempt(s): {error}")

    def metadata(self, cancel_requested=None) -> tuple[dict,bytes]:
        payload,_status,raw=self.fetch(self.base_url,cancel_requested);return payload,raw

    def annotations(self, accession: str, layer: str, cancel_requested=None) -> tuple[dict,int,bytes,str]:
        if layer not in LAYERS:raise ValueError(f"Unsupported InterPro layer: {layer}")
        url=f"{self.base_url}entry/{layer}/protein/uniprot/{urllib.parse.quote(accession,safe='-')}/?page_size=200"
        combined={"count":0,"next":None,"previous":None,"results":[]};raw_pages=[];status=200
        while url:
            payload,status,raw=self.fetch(url,cancel_requested);raw_pages.append(json.loads(raw.decode("utf-8")))
            combined["results"].extend(payload.get("results",[]));combined["count"]=payload.get("count",len(combined["results"]));url=payload.get("next")
            if cancel_requested and cancel_requested():raise InterProCancelled("InterPro/Pfam acquisition was canceled.")
        encoded=json.dumps({"pages":raw_pages,"combined":combined},ensure_ascii=False,separators=(",",":")).encode("utf-8")
        return combined,status,encoded,f"{self.base_url}entry/{layer}/protein/uniprot/{urllib.parse.quote(accession,safe='-')}/"


@dataclass(frozen=True)
class AnnotationSet:
    root: Path
    manifest: dict
    proteins: pd.DataFrame
    interpro_membership: pd.DataFrame
    interpro_locations: pd.DataFrame
    pfam_membership: pd.DataFrame
    pfam_locations: pd.DataFrame
    interpro_entries: pd.DataFrame
    pfam_entries: pd.DataFrame
    interpro_hierarchy: pd.DataFrame


class InterProDatabase:
    def __init__(self, root: Path, provider: InterProProvider|None=None):
        self.root=Path(root)/"interpro";self.metadata_root=self.root/"metadata";self.cache_root=self.root/"cache";self.provider=provider or InterProProvider()

    @property
    def active_metadata_pointer(self):return self.metadata_root/"active.json"
    def metadata_snapshot(self):
        try:
            sid=json.loads(self.active_metadata_pointer.read_text(encoding="utf-8"))["snapshot_id"];path=self.metadata_root/"snapshots"/sid
            return path if path.is_dir() else None
        except (OSError,KeyError,json.JSONDecodeError):return None
    def metadata_manifest(self):
        snap=self.metadata_snapshot()
        try:return json.loads((snap/"manifest.json").read_text(encoding="utf-8")) if snap else {}
        except (OSError,json.JSONDecodeError):return {}
    def state(self):
        if self.metadata_snapshot():return DatabaseState.READY
        attempts=sorted((self.metadata_root/"staging").glob("*/manifest.json"),reverse=True) if (self.metadata_root/"staging").exists() else []
        if not attempts:return DatabaseState.NOT_INSTALLED
        try:return DatabaseState(json.loads(attempts[0].read_text(encoding="utf-8")).get("status",DatabaseState.ERROR))
        except (OSError,ValueError,json.JSONDecodeError):return DatabaseState.ERROR
    def cached_protein_count(self):
        return len({path.parent.name for path in self.cache_root.glob("*/*/*.json") if not path.name.endswith(".meta.json")}) if self.cache_root.exists() else 0

    def download_metadata(self, cancel_requested=None) -> Path:
        sid=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ");staging=self.metadata_root/"staging"/sid;final=self.metadata_root/"snapshots"/sid
        staging.mkdir(parents=True);manifest={"database":"InterPro / Pfam","snapshot_id":sid,"status":str(DatabaseState.DOWNLOADING),"api_base":self.provider.base_url,"created_at":_now()};_atomic_json(staging/"manifest.json",manifest)
        try:
            payload,raw=self.provider.metadata(cancel_requested)
            if cancel_requested and cancel_requested():raise InterProCancelled("InterPro/Pfam metadata download was canceled.")
            (staging/"raw").mkdir();(staging/"tables").mkdir();(staging/"raw/api_root.json").write_bytes(raw)
            databases=payload.get("databases",{});rows=[{"database":key,"name":value.get("name"),"version":value.get("version"),"release_date":value.get("releaseDate"),"type":value.get("type")} for key,value in databases.items()]
            pd.DataFrame(rows).to_csv(staging/"tables/database_versions.csv",index=False)
            interpro=databases.get("interpro",{});pfam=databases.get("pfam",{});context="interpro-{}__pfam-{}".format(interpro.get("version") or "unknown",pfam.get("version") or "unknown")
            manifest.update(status=str(DatabaseState.READY),completed_at=_now(),interpro_release=interpro.get("version") or "unknown",pfam_release=pfam.get("version") or "unknown",release_context=context,files={"raw/api_root.json":_sha_file(staging/"raw/api_root.json"),"tables/database_versions.csv":_sha_file(staging/"tables/database_versions.csv")})
            _atomic_json(staging/"manifest.json",manifest);final.parent.mkdir(parents=True,exist_ok=True);os.replace(staging,final);_atomic_json(self.active_metadata_pointer,{"snapshot_id":sid,"activated_at":_now()});return final
        except Exception as exc:
            manifest.update(status=str(DatabaseState.INCOMPLETE if isinstance(exc,InterProCancelled) else DatabaseState.ERROR),last_error=str(exc));_atomic_json(staging/"manifest.json",manifest);raise

    def _cache_paths(self, context: str, accession: str, layer: str):
        base=self.cache_root/context/_safe_accession(accession);return base/f"{layer}.json",base/f"{layer}.meta.json"
    def _acquire(self, accession, layer, context, releases, cancel_requested=None):
        raw_path,meta_path=self._cache_paths(context,accession,layer)
        if raw_path.is_file() and meta_path.is_file():
            meta=json.loads(meta_path.read_text(encoding="utf-8"))
            if meta.get("release_context")==context and meta.get("response_sha256")==_sha_file(raw_path):return json.loads(raw_path.read_text(encoding="utf-8"))["combined"],meta,raw_path
        payload,status,raw,url=self.provider.annotations(accession,layer,cancel_requested);raw_path.parent.mkdir(parents=True,exist_ok=True);part=raw_path.with_name(raw_path.name+".part");part.write_bytes(raw);os.replace(part,raw_path)
        meta={"requested_accession":accession,"retrieved_at":_now(),"source_url":url,"http_status":status,"response_sha256":_sha_bytes(raw),"release_context":context,"interpro_release":releases[0],"pfam_release":releases[1]};_atomic_json(meta_path,meta);return payload,meta,raw_path

    def build_annotation_set(self, project_path: Path, accessions: Iterable[str], *, annotation_set_id: str|None=None, progress=None, cancel_requested=None) -> Path:
        metadata=self.metadata_manifest()
        if not metadata:raise InterProError("InterPro/Pfam metadata is not installed.")
        proteins=tuple(dict.fromkeys(str(value).strip() for value in accessions if str(value).strip()))
        if not proteins:raise InterProError("At least one UniProt accession is required.")
        set_id=annotation_set_id or datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ");base=Path(project_path)/"annotations"/"InterPro";staging=base/"staging"/set_id;final=base/"sets"/set_id
        if staging.exists() or final.exists():raise InterProError(f"Annotation set already exists: {set_id}")
        for folder in ("raw/interpro","raw/pfam","tables","metadata"): (staging/folder).mkdir(parents=True,exist_ok=True)
        context=metadata.get("release_context","unknown");releases=(metadata.get("interpro_release","unknown"),metadata.get("pfam_release","unknown"));manifest={"database":"InterPro / Pfam","organism_context":"protein accessions supplied explicitly","annotation_set_id":set_id,"status":str(DatabaseState.DOWNLOADING),"created_at":_now(),"interpro_release":releases[0],"pfam_release":releases[1],"api_base":self.provider.base_url,"metadata_snapshot":self.metadata_snapshot().name,"requested_protein_count":len(proteins)};_atomic_json(staging/"manifest.json",manifest)
        protein_rows=[];results={};hash_rows=[]
        try:
            total=len(proteins)*2;done=0
            for accession in proteins:
                results[accession]={};statuses={};counts={}
                for layer in LAYERS:
                    if cancel_requested and cancel_requested():raise InterProCancelled("InterPro/Pfam annotation set build was canceled.")
                    progress and progress({"done":done,"total":total,"message":f"Fetching {layer.title()} annotations for {accession}","accession":accession,"layer":layer})
                    payload,meta,cache_path=self._acquire(accession,layer,context,releases,cancel_requested);destination=staging/"raw"/layer/f"{_safe_accession(accession)}.json";shutil.copy2(cache_path,destination)
                    hash_rows.append({"requested_accession":accession,"layer":layer,"relative_path":str(destination.relative_to(staging)).replace("\\","/"),"sha256":_sha_file(destination),"source_url":meta["source_url"],"http_status":meta["http_status"],"retrieved_at":meta["retrieved_at"]})
                    results[accession][layer]=payload;counts[layer]=len(payload.get("results",[]));statuses[layer]="no_matches" if counts[layer]==0 else "annotated";done+=1
                    progress and progress({"done":done,"total":total,"message":f"Protein {proteins.index(accession)+1} of {len(proteins)}: {layer.title()} complete","accession":accession,"layer":layer})
                annotation_status="annotated" if any(counts.values()) else "no_matches";protein_rows.append({"requested_accession":accession,"interpro_query_status":statuses["interpro"],"pfam_query_status":statuses["pfam"],"interpro_entry_count":counts["interpro"],"pfam_entry_count":counts["pfam"],"annotation_status":annotation_status})
            frames=self._normalize(results);frames["proteins"]=pd.DataFrame(protein_rows);frames["hashes"]=pd.DataFrame(hash_rows)
            filenames={"proteins":"proteins.csv","interpro_membership":"protein_interpro_membership.csv","interpro_locations":"interpro_locations.csv","pfam_membership":"protein_pfam_membership.csv","pfam_locations":"pfam_locations.csv","interpro_entries":"interpro_entries.csv","pfam_entries":"pfam_entries.csv","interpro_hierarchy":"interpro_hierarchy.csv"}
            for key,name in filenames.items():frames[key].to_csv(staging/"tables"/name,index=False)
            frames["hashes"].to_csv(staging/"raw_response_hashes.csv",index=False);shutil.copy2(self.metadata_snapshot()/"manifest.json",staging/"metadata/metadata_manifest.json")
            table_hashes={name:_sha_file(staging/"tables"/name) for name in filenames.values()};counts=pd.DataFrame(protein_rows).annotation_status.value_counts().to_dict()
            manifest.update(status=str(DatabaseState.READY),completed_at=_now(),annotated_protein_count=int(counts.get("annotated",0)),no_match_protein_count=int(counts.get("no_matches",0)),failed_protein_count=0,raw_response_hash_index_sha256=_sha_file(staging/"raw_response_hashes.csv"),metadata_manifest_sha256=_sha_file(staging/"metadata/metadata_manifest.json"),table_hashes=table_hashes)
            _atomic_json(staging/"manifest.json",manifest);self._validate_root(staging);final.parent.mkdir(parents=True,exist_ok=True);os.replace(staging,final);_atomic_json(base/"active.json",{"annotation_set_id":set_id,"activated_at":_now()});return final
        except Exception as exc:
            manifest.update(status=str(DatabaseState.INCOMPLETE if isinstance(exc,InterProCancelled) else DatabaseState.ERROR),last_error=str(exc),failed_protein_count=max(0,len(proteins)-len(protein_rows)));_atomic_json(staging/"manifest.json",manifest);raise

    @staticmethod
    def _normalize(results):
        im=[];il=[];pm=[];pl=[];ie={};pe={};hier=[]
        for accession,layers in results.items():
            for layer,payload in layers.items():
                for result in payload.get("results",[]):
                    metadata=result.get("metadata") or {};entry=str(metadata.get("accession") or "");name=metadata.get("name");etype=metadata.get("type");integrated=metadata.get("integrated");parent=metadata.get("parent")
                    if not entry:continue
                    if layer=="interpro":
                        im.append({"uniprot_accession":accession,"interpro_id":entry,"entry_name":name,"entry_type":etype});ie[entry]={"interpro_id":entry,"name":name,"type":etype,"parent":parent}
                        if parent:hier.append({"parent_interpro_id":parent,"child_interpro_id":entry})
                    else:pm.append({"uniprot_accession":accession,"pfam_id":entry,"pfam_name":name,"integrated_interpro_id":integrated,"entry_type":etype});pe[entry]={"pfam_id":entry,"name":name,"integrated_interpro_id":integrated,"type":etype}
                    for protein in result.get("proteins") or []:
                        for location_index,location in enumerate(protein.get("entry_protein_locations") or [],1):
                            for fragment_index,fragment in enumerate(location.get("fragments") or [],1):
                                row={"uniprot_accession":accession,"start":fragment.get("start"),"end":fragment.get("end"),"fragment":fragment.get("dc-status"),"location_index":location_index,"fragment_index":fragment_index,"representative":location.get("representative"),"model":location.get("model"),"score":location.get("score")}
                                if layer=="interpro":row["interpro_id"]=entry;il.append(row)
                                else:row["pfam_id"]=entry;pl.append(row)
        def frame(rows,columns):return pd.DataFrame(rows,columns=columns).drop_duplicates().reset_index(drop=True)
        return {"interpro_membership":frame(im,["uniprot_accession","interpro_id","entry_name","entry_type"]),"interpro_locations":frame(il,["uniprot_accession","interpro_id","start","end","fragment","location_index","fragment_index","representative","model","score"]),"pfam_membership":frame(pm,["uniprot_accession","pfam_id","pfam_name","integrated_interpro_id","entry_type"]),"pfam_locations":frame(pl,["uniprot_accession","pfam_id","start","end","fragment","location_index","fragment_index","representative","model","score"]),"interpro_entries":frame(list(ie.values()),["interpro_id","name","type","parent"]),"pfam_entries":frame(list(pe.values()),["pfam_id","name","integrated_interpro_id","type"]),"interpro_hierarchy":frame(hier,["parent_interpro_id","child_interpro_id"])}

    def list_annotation_sets(self, project_path: Path):
        root=Path(project_path)/"annotations/InterPro/sets";return sorted(path.name for path in root.iterdir() if path.is_dir()) if root.is_dir() else []
    def _set_root(self,project_path,set_id=None):
        base=Path(project_path)/"annotations/InterPro"
        if set_id is None:
            try:set_id=json.loads((base/"active.json").read_text(encoding="utf-8"))["annotation_set_id"]
            except (OSError,KeyError,json.JSONDecodeError):return None
        return base/"sets"/set_id
    def _validate_root(self,root):
        try:manifest=json.loads((Path(root)/"manifest.json").read_text(encoding="utf-8"))
        except (OSError,json.JSONDecodeError) as exc:raise InvalidAnnotationSet(f"Invalid annotation set manifest: {exc}") from exc
        if manifest.get("status")!=str(DatabaseState.READY) or manifest.get("failed_protein_count")!=0:raise InvalidAnnotationSet("Annotation set is not scientifically Ready.")
        for name in TABLES:
            path=Path(root)/"tables"/name
            if not path.is_file() or manifest.get("table_hashes",{}).get(name)!=_sha_file(path):raise InvalidAnnotationSet(f"Missing or invalid annotation table: {name}")
        index=Path(root)/"raw_response_hashes.csv"
        if not index.is_file() or manifest.get("raw_response_hash_index_sha256")!=_sha_file(index):raise InvalidAnnotationSet("Raw response hash index is missing or invalid.")
        metadata=Path(root)/"metadata/metadata_manifest.json"
        if not metadata.is_file() or manifest.get("metadata_manifest_sha256")!=_sha_file(metadata):raise InvalidAnnotationSet("Metadata provenance is missing or invalid.")
        for row in pd.read_csv(index).itertuples(index=False):
            raw=Path(root)/str(row.relative_path)
            if not raw.is_file() or _sha_file(raw)!=str(row.sha256):raise InvalidAnnotationSet(f"Raw response is missing or invalid: {row.relative_path}")
        return manifest
    def validate_annotation_set(self,project_path,set_id=None):
        root=self._set_root(project_path,set_id)
        try:self._validate_root(root);return True
        except (InvalidAnnotationSet,TypeError):return False
    def is_annotation_set_ready(self,project_path,set_id=None):return self.validate_annotation_set(project_path,set_id)
    def load_annotation_set(self,project_path,set_id=None):
        root=self._set_root(project_path,set_id)
        if root is None:raise InvalidAnnotationSet("No active InterPro annotation set.")
        manifest=self._validate_root(root);frames={name:pd.read_csv(root/"tables"/filename) for name,filename in {"proteins":"proteins.csv","interpro_membership":"protein_interpro_membership.csv","interpro_locations":"interpro_locations.csv","pfam_membership":"protein_pfam_membership.csv","pfam_locations":"pfam_locations.csv","interpro_entries":"interpro_entries.csv","pfam_entries":"pfam_entries.csv","interpro_hierarchy":"interpro_hierarchy.csv"}.items()}
        return AnnotationSet(root,manifest,**frames)
    def find_interpro_memberships(self,project_path,accession,set_id=None):return self.load_annotation_set(project_path,set_id).interpro_membership.query("uniprot_accession == @accession").copy()
    def find_pfam_memberships(self,project_path,accession,set_id=None):return self.load_annotation_set(project_path,set_id).pfam_membership.query("uniprot_accession == @accession").copy()
    def find_interpro_locations(self,project_path,accession,set_id=None):return self.load_annotation_set(project_path,set_id).interpro_locations.query("uniprot_accession == @accession").copy()
    def find_pfam_locations(self,project_path,accession,set_id=None):return self.load_annotation_set(project_path,set_id).pfam_locations.query("uniprot_accession == @accession").copy()
