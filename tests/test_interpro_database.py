from __future__ import annotations

import json
import urllib.error
from pathlib import Path

import pandas as pd
import pytest

from pichanalysis.core.database_registry import DatabaseState
from pichanalysis.core.interpro_database import (
    InterProCancelled, InterProDatabase, InterProProvider, InvalidAnnotationSet,
)


def metadata_payload(interpro="110.0",pfam="38.2"):
    return {"databases":{"interpro":{"name":"InterPro","version":interpro,"releaseDate":"2026-09-03T00:00:00Z","type":"entry"},"pfam":{"name":"Pfam","version":pfam,"releaseDate":"2026-03-26T00:00:00Z","type":"entry"}}}

def hit(entry,name,etype,locations,integrated=None):
    return {"metadata":{"accession":entry,"name":name,"source_database":"pfam" if entry.startswith("PF") else "interpro","type":etype,"integrated":integrated,"parent":None},"proteins":[{"accession":"p1","entry_protein_locations":[{"fragments":[{"start":start,"end":end,"dc-status":"CONTINUOUS"}],"representative":False,"model":entry,"score":score} for start,end,score in locations]}]}

class FakeProvider:
    base_url="https://www.ebi.ac.uk/interpro/api/"
    def __init__(self,releases=("110.0","38.2"),fail=None):self.releases=releases;self.fail=fail;self.calls=[]
    def metadata(self,cancel_requested=None):
        if cancel_requested and cancel_requested():raise InterProCancelled("cancelled")
        payload=metadata_payload(*self.releases);return payload,json.dumps(payload).encode()
    def annotations(self,accession,layer,cancel_requested=None):
        self.calls.append((accession,layer))
        if cancel_requested and cancel_requested():raise InterProCancelled("cancelled")
        if accession==self.fail:raise RuntimeError("planned HTTP error")
        if accession=="EMPTY" or "-2" in accession:results=[]
        elif layer=="interpro":results=[hit("IPR000001","Example family","family",[(10,50,1.0),(100,140,.5)]),hit("IPR000002","Overlapping domain","domain",[(40,80,.2)])]
        else:results=[hit("PF00001","Repeated Pfam","domain",[(10,50,1.0),(100,140,.5)],"IPR000001"),hit("PF99999","Unintegrated Pfam","repeat",[(40,80,.2)],None)]
        payload={"count":len(results),"next":None,"results":results};raw=json.dumps({"pages":[payload],"combined":payload}).encode();return payload,200,raw,f"{self.base_url}{layer}/{accession}"

def ready_database(tmp_path,provider=None):
    database=InterProDatabase(tmp_path/"global",provider or FakeProvider());database.download_metadata();return database

def test_metadata_cache_normalization_repeats_overlaps_and_offline_loading(tmp_path):
    provider=FakeProvider();database=ready_database(tmp_path,provider);project=tmp_path/"project"
    root=database.build_annotation_set(project,["P1","EMPTY","P12345-2"],annotation_set_id="A")
    loaded=database.load_annotation_set(project,"A");assert root==loaded.root and loaded.manifest["status"]==DatabaseState.READY
    assert len(database.find_interpro_memberships(project,"P1","A"))==2 and len(database.find_pfam_memberships(project,"P1","A"))==2
    assert len(database.find_interpro_locations(project,"P1","A").query("interpro_id == 'IPR000001'"))==2
    assert len(database.find_pfam_locations(project,"P1","A").query("pfam_id == 'PF00001'"))==2
    pfam=database.find_pfam_memberships(project,"P1","A");assert pd.isna(pfam.loc[pfam.pfam_id=="PF99999","integrated_interpro_id"]).all()
    locations=database.find_pfam_locations(project,"P1","A");assert ((locations.start==10)&(locations.end==50)).any() and ((locations.start==40)&(locations.end==80)).any()
    proteins=loaded.proteins.set_index("requested_accession");assert proteins.loc["EMPTY","annotation_status"]=="no_matches" and proteins.loc["P12345-2","annotation_status"]=="no_matches"
    assert loaded.interpro_entries.set_index("interpro_id").loc["IPR000001","type"]=="family" and loaded.pfam_entries.set_index("pfam_id").loc["PF99999","type"]=="repeat"
    calls=len(provider.calls);database.build_annotation_set(project,["P1"],annotation_set_id="B");assert len(provider.calls)==calls
    database.provider=object();assert len(database.find_pfam_locations(project,"P1","A"))==3
    assert database.list_annotation_sets(project)==["A","B"] and database.validate_annotation_set(project,"A")

def test_release_unknown_mismatch_failure_cancel_and_immutability(tmp_path):
    provider=FakeProvider((None,None));database=ready_database(tmp_path,provider);assert database.metadata_manifest()["release_context"]=="interpro-unknown__pfam-unknown"
    project=tmp_path/"project";first=database.build_annotation_set(project,["P1"],annotation_set_id="A");original=(first/"manifest.json").read_bytes()
    with pytest.raises(Exception):database.build_annotation_set(project,["P1"],annotation_set_id="A")
    assert (first/"manifest.json").read_bytes()==original
    changed=FakeProvider(("111.0","39.0"));database.provider=changed;database.download_metadata();database.build_annotation_set(project,["P1"],annotation_set_id="B");assert changed.calls==[("P1","interpro"),("P1","pfam")]
    failed=ready_database(tmp_path/"failure",FakeProvider(fail="BAD"))
    with pytest.raises(RuntimeError,match="planned HTTP error"):failed.build_annotation_set(tmp_path/"failed-project",["BAD"],annotation_set_id="failed")
    assert failed.list_annotation_sets(tmp_path/"failed-project")==[]
    cancelled=ready_database(tmp_path/"cancel",FakeProvider())
    with pytest.raises(InterProCancelled):cancelled.build_annotation_set(tmp_path/"cancel-project",["P1"],annotation_set_id="cancelled",cancel_requested=lambda:True)
    assert cancelled.list_annotation_sets(tmp_path/"cancel-project")==[]

def test_hash_validation_and_a_b_a(tmp_path):
    database=ready_database(tmp_path);project=tmp_path/"project";database.build_annotation_set(project,["P1"],annotation_set_id="A");database.build_annotation_set(project,["EMPTY"],annotation_set_id="B")
    a=database.load_annotation_set(project,"A");b=database.load_annotation_set(project,"B");again=database.load_annotation_set(project,"A");pd.testing.assert_frame_equal(a.pfam_locations,again.pfam_locations);assert len(a.pfam_membership)>len(b.pfam_membership)
    (a.root/"tables/proteins.csv").write_text("tampered",encoding="utf-8");assert not database.is_annotation_set_ready(project,"A")
    with pytest.raises(InvalidAnnotationSet):database.load_annotation_set(project,"A")

def test_provider_retry_404_and_pagination(monkeypatch):
    attempts=[]
    class Response:
        status=200
        def __init__(self,payload):self.payload=payload
        def __enter__(self):return self
        def __exit__(self,*a):pass
        def read(self):return json.dumps(self.payload).encode()
    def opener(request,timeout):
        attempts.append(request.full_url)
        if len(attempts)==1:raise urllib.error.HTTPError(request.full_url,429,"rate",{},None)
        return Response({"count":0,"next":None,"results":[]})
    monkeypatch.setattr("time.sleep",lambda *_:None);provider=InterProProvider(retries=1,opener=opener);payload,status,_=provider.fetch("https://example.invalid")
    assert status==200 and payload["count"]==0 and len(attempts)==2
    def missing(request,timeout):raise urllib.error.HTTPError(request.full_url,404,"missing",{},None)
    payload,status,_=InterProProvider(opener=missing).fetch("https://example.invalid");assert status==404 and payload["results"]==[]
