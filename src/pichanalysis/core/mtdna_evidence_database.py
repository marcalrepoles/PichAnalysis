"""Immutable, source-specific human mtDNA evidence snapshots.

MitoCarta and GO are read from existing local snapshots. Only the NCBI RefSeq
record may be fetched; all scientific claims remain individual evidence rows.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import sqlite3
import urllib.request
import uuid
import xml.etree.ElementTree as ET
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from .database_registry import DatabaseState

NCBI_ACCESSION = "NC_012920.1"
NCBI_URL = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            "?db=nuccore&id=NC_012920.1&rettype=gb&retmode=xml")
GO_SEEDS = {
    "GO:0042645": ("mitochondrial_nucleoid", "physical_localization"),
    "GO:0000262": ("mitochondrial_chromosome", "physical_localization"),
    "GO:0032042": ("mtDNA_metabolism", "general_mtDNA_metabolism"),
    "GO:0006264": ("mtDNA_replication", "replication"),
    "GO:0043504": ("mtDNA_repair", "repair"),
    "GO:0090139": ("mtDNA_packaging", "packaging"),
    "GO:0006390": ("mitochondrial_transcription", "transcription"),
    "GO:1901858": ("mtDNA_metabolism_regulation", "regulation"),
    "GO:1903108": ("mitochondrial_transcription_regulation", "regulation"),
}
CATEGORIES = (
    ("mtDNA_encoded", "mtDNA encoded", "genomic_origin", "Protein encoded by the mitochondrial genome; no binding or maintenance inference."),
    ("mtDNA_maintenance", "mtDNA maintenance", "maintenance", "Curated mtDNA maintenance pathway."),
    ("mtDNA_replication", "mtDNA replication", "replication", "Mitochondrial DNA replication evidence."),
    ("mitochondrial_nucleoid", "Mitochondrial nucleoid", "physical_localization", "Specific nucleoid evidence."),
    ("mitochondrial_chromosome", "Mitochondrial chromosome", "physical_localization", "Specific mitochondrial chromosome evidence."),
    ("mtDNA_repair", "mtDNA repair", "repair", "Mitochondrial DNA repair evidence."),
    ("mtDNA_packaging", "mtDNA packaging", "packaging", "Mitochondrial chromosome packaging evidence."),
    ("mtDNA_modification", "mtDNA modification", "modification", "Mitochondrial DNA modification evidence."),
    ("mtDNA_stability_decay", "mtDNA stability and decay", "stability_decay", "Mitochondrial DNA stability or decay evidence."),
    ("mtDNA_metabolism", "Mitochondrial DNA metabolism", "general_mtDNA_metabolism", "Mitochondrial DNA metabolic process evidence."),
    ("mitochondrial_transcription", "Mitochondrial transcription", "transcription", "Transcription from the mitochondrial genome."),
    ("mtDNA_metabolism_regulation", "mtDNA metabolism regulation", "regulation", "Regulation, not direct execution, of mtDNA metabolism."),
    ("mitochondrial_transcription_regulation", "Mitochondrial transcription regulation", "regulation", "Regulation, not direct execution, of mitochondrial transcription."),
)
EVIDENCE_COLUMNS = ("entity_key","ncbi_gene_id","gene_symbol","uniprot","source","source_snapshot",
    "evidence_class","evidence_category","source_identifier","source_name","source_detail",
    "evidence_code","reference","qualifier","assigned_by","seed_go_id","annotated_go_id",
    "relation_path","source_path","protein_id","coordinates")
TABLES = ("entities","evidence_records","category_membership","categories","mitocarta_evidence",
          "go_evidence","mtdna_encoded","source_agreement","go_exclusions","metadata")


class MtdnaEvidenceError(RuntimeError): pass
class MissingMtdnaDependency(MtdnaEvidenceError): pass
class MtdnaBuildCancelled(MtdnaEvidenceError): pass


def _now(): return datetime.now(timezone.utc).isoformat()


def _sha(path):
    digest=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b""):digest.update(chunk)
    return digest.hexdigest()


def _text(value): return "" if pd.isna(value) else str(value).strip()


def _gene_id(value):
    value=_text(value)
    return value[:-2] if value.endswith(".0") and value[:-2].isdigit() else value


def _uniprots(value):
    return sorted({token.strip() for token in re.split(r"[;,|\s]+",_text(value)) if token.strip() and token.strip().lower()!="nan"})


def parse_ncbi_refseq(payload: bytes):
    """Parse official EFetch GenBank XML, deriving protein products from CDS."""
    try:root=ET.fromstring(payload)
    except ET.ParseError as error:raise MtdnaEvidenceError(f"Invalid NCBI GenBank XML: {error}") from error
    seq=root.find(".//GBSeq")
    if seq is None:raise MtdnaEvidenceError("NCBI GenBank XML contains no GBSeq.")
    accession=seq.findtext("GBSeq_accession-version","")
    organism=seq.findtext("GBSeq_organism","")
    if accession!=NCBI_ACCESSION or organism!="Homo sapiens":
        raise MtdnaEvidenceError(f"Unexpected NCBI reference: {accession} / {organism}.")
    taxon=set(); products=[]
    for feature in seq.findall(".//GBFeature"):
        qualifiers=defaultdict(list)
        for item in feature.findall(".//GBQualifier"):
            qualifiers[item.findtext("GBQualifier_name","")].append(item.findtext("GBQualifier_value","") or "")
        if feature.findtext("GBFeature_key")=="source":
            taxon.update(value for value in qualifiers["db_xref"] if value.startswith("taxon:"))
        if feature.findtext("GBFeature_key")!="CDS":continue
        ids=[value.split(":",1)[1] for value in qualifiers["db_xref"] if value.startswith("GeneID:")]
        gene=(qualifiers["gene"] or [""])[0]
        if not gene or not ids:raise MtdnaEvidenceError("NCBI mitochondrial CDS lacks gene or GeneID.")
        products.append({"gene_symbol":gene,"ncbi_gene_id":ids[0],"protein_id":(qualifiers["protein_id"] or [""])[0],
            "product":(qualifiers["product"] or [""])[0],"coordinates":feature.findtext("GBFeature_location","") or ""})
    if taxon!={"taxon:9606"} or not products:
        raise MtdnaEvidenceError("NCBI mitochondrial record requires taxon:9606 and protein CDS features.")
    return products


def mitocarta_category(path):
    levels=[part.strip().casefold() for part in str(path).split(">")]
    if len(levels)<2 or levels[0]!="mitochondrial central dogma":return None
    if levels[1]=="mtrna metabolism":
        return "mitochondrial_transcription" if len(levels)>=3 and levels[2]=="transcription" else None
    if levels[1]!="mtdna maintenance":return None
    suffix=" ".join(levels[2:])
    for word,category in (("replication","mtDNA_replication"),("nucleoid","mitochondrial_nucleoid"),
                          ("repair","mtDNA_repair"),("modification","mtDNA_modification"),
                          ("stability","mtDNA_stability_decay"),("decay","mtDNA_stability_decay")):
        if word in suffix:return category
    return "mtDNA_maintenance"


def go_ancestor_paths(term, relations):
    """Only is_a and part_of edges; keep the actual relation chain."""
    queue=[(term,())]; seen={term}
    while queue:
        current,path=queue.pop(0)
        yield current,path
        for parent,relation in relations.get(current,()):
            if relation in {"is_a","part_of"} and parent not in seen:
                seen.add(parent);queue.append((parent,(*path,f"{relation}:{parent}")))


def _read_source_table(snapshot,name):
    path=Path(snapshot)/"tables"/name
    if not path.is_file():raise MissingMtdnaDependency(f"Required local source table is missing: {path}")
    return pd.read_csv(path,dtype=str).fillna("")


class MtdnaEvidenceDatabase:
    def __init__(self,root,provider=None):
        self.database_root=Path(root)
        self.root=self.database_root/"mtdna_evidence"/"human"
        self.snapshots=self.root/"snapshots"
        self.active_pointer=self.root/"active_snapshot.json"
        self.provider=provider

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
        return snap/"mtdna_evidence.sqlite" if snap else None

    def dependencies(self,mitocarta=None,go_snapshot=None):
        from .mitocarta_database import MitoCartaDatabase
        mitocarta=mitocarta or MitoCartaDatabase(self.database_root)
        mito=mitocarta.active_snapshot()
        go=Path(go_snapshot) if go_snapshot else self._active_go_snapshot()
        return {"mitocarta_ready":bool(mito),"mitocarta_snapshot":mito,
            "go_ready":bool(go and (go/"manifest.json").is_file() and self._go_manifest(go).get("status")==DatabaseState.READY and self._go_manifest(go).get("tax_id")=="9606"),
            "go_snapshot":go}

    def _active_go_snapshot(self):
        pointer=self.database_root/"gene_ontology"/"human"/"active_snapshot.json"
        try:return self.database_root/"gene_ontology"/"human"/"snapshots"/json.loads(pointer.read_text(encoding="utf-8"))["snapshot_id"]
        except (OSError,ValueError,KeyError):return None

    @staticmethod
    def _go_manifest(snap):
        try:return json.loads((Path(snap)/"manifest.json").read_text(encoding="utf-8"))
        except (OSError,ValueError):return {}

    def _fetch_ncbi(self,cancel_requested):
        if cancel_requested():raise MtdnaBuildCancelled("mtDNA evidence build was canceled.")
        if self.provider:return self.provider.fetch(NCBI_URL)
        request=urllib.request.Request(NCBI_URL,headers={"User-Agent":"PichAnalysis/mtDNA-evidence (NCBI E-utilities)"})
        with urllib.request.urlopen(request,timeout=45) as response:
            data=response.read(5_000_001)
        if len(data)>5_000_000:raise MtdnaEvidenceError("NCBI reference exceeds expected size.")
        return data

    def build(self,mitocarta=None,go_snapshot=None,ncbi_xml=None,progress=None,cancel_requested=None):
        """Build in staging, validate, then atomically promote and activate."""
        from .mitocarta_database import MitoCartaDatabase
        mitocarta=mitocarta or MitoCartaDatabase(self.database_root)
        cancel_requested=cancel_requested or (lambda:False)
        emit=lambda message:progress({"message":message}) if progress else None
        emit("Checking source databases...")
        deps=self.dependencies(mitocarta,go_snapshot)
        if not deps["mitocarta_ready"]:raise MissingMtdnaDependency("MitoCarta database is required.")
        if not deps["go_ready"]:raise MissingMtdnaDependency("Gene Ontology database is required.")
        mito,go=deps["mitocarta_snapshot"],deps["go_snapshot"]
        mm=mitocarta.manifest(mito); gm=self._go_manifest(go)
        for name in ("mitocarta_genes.csv","mitopathway_membership.csv","mitopathway_hierarchy.csv"):
            if not (mito/"tables"/name).is_file():raise MissingMtdnaDependency(f"MitoCarta source table missing: {name}")
        for name in ("go_terms.csv","go_relations.csv","go_annotations.csv"):
            if not (go/"tables"/name).is_file():raise MissingMtdnaDependency(f"Gene Ontology source table missing: {name}")
        sid=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")+"-"+uuid.uuid4().hex[:8]
        stage=self.snapshots/sid;final=stage
        (stage/"raw").mkdir(parents=True);(stage/"tables").mkdir()
        try:
            if cancel_requested():raise MtdnaBuildCancelled("mtDNA evidence build was canceled.")
            emit("Retrieving NCBI mitochondrial reference...")
            payload=Path(ncbi_xml).read_bytes() if ncbi_xml else self._fetch_ncbi(cancel_requested)
            products=parse_ncbi_refseq(payload)
            raw=stage/"raw"/"NC_012920.1.xml";raw.write_bytes(payload)
            emit("Reading MitoCarta mtDNA pathways...")
            mito_rows=self._mitocarta_rows(mito,mito.name)
            if cancel_requested():raise MtdnaBuildCancelled("mtDNA evidence build was canceled.")
            emit("Reading Gene Ontology mtDNA annotations...")
            go_rows,excluded=self._go_rows(go,go.name)
            if cancel_requested():raise MtdnaBuildCancelled("mtDNA evidence build was canceled.")
            emit("Reconciling identifiers...")
            ncbi_rows=[]
            for item in products:
                ncbi_rows.append(self._record(ncbi_gene_id=item["ncbi_gene_id"],gene_symbol=item["gene_symbol"],
                    source="NCBI mtDNA genome",source_snapshot=NCBI_ACCESSION,evidence_class="genomic_origin",
                    evidence_category="mtDNA_encoded",source_identifier=NCBI_ACCESSION,
                    source_name=item["product"],source_detail=item["coordinates"],protein_id=item["protein_id"],
                    coordinates=item["coordinates"]))
            records=self._reconcile(mito_rows+go_rows+ncbi_rows)
            if not records:raise MtdnaEvidenceError("No mtDNA evidence records were derived.")
            emit("Building evidence index...")
            frames=self._frames(records,excluded)
            for name,frame in frames.items():frame.to_csv(stage/"tables"/f"{name}.csv",index=False)
            mito_hashes={item["filename"]:item["sha256"] for item in mm.get("files",[]) if "filename" in item and "sha256" in item}
            go_hashes={name:_sha(go/"tables"/name) for name in ("go_terms.csv","go_relations.csv","go_annotations.csv")}
            manifest={"database":"mtDNA evidence","snapshot_id":sid,"status":DatabaseState.INCOMPLETE,
                "created_at":_now(),"organism":"Homo sapiens","tax_id":"9606",
                "mitocarta_snapshot_id":mito.name,"mitocarta_snapshot_path":str(mito),"mitocarta_source_hashes":mito_hashes,
                "mitocarta_manifest_sha256":_sha(mito/"manifest.json"),
                "go_snapshot_id":go.name,"go_snapshot_path":str(go),"go_source_hashes":go_hashes,
                "go_manifest_sha256":_sha(go/"manifest.json"),
                "go_source_manifest":gm,"ncbi_accession":NCBI_ACCESSION,"ncbi_source_url":NCBI_URL,
                "ncbi_raw_size":len(payload),"ncbi_raw_sha256":_sha(raw),
                "counts":{**{name:len(frame) for name,frame in frames.items()},
                    "mitocarta_supported_genes":frames["mitocarta_evidence"].entity_key.nunique(),
                    "go_supported_genes":frames["go_evidence"].entity_key.nunique(),
                    "mtdna_encoded_genes":frames["mtdna_encoded"].entity_key.nunique()},
                "category_counts":frames["category_membership"].category.value_counts().to_dict()}
            emit("Indexing database...")
            self._write_sqlite(stage,frames,manifest)
            manifest["sqlite_sha256"]=_sha(stage/"mtdna_evidence.sqlite")
            (stage/"manifest.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
            emit("Validating snapshot...")
            self.validate_snapshot(stage,allow_incomplete=True)
            if cancel_requested():raise MtdnaBuildCancelled("mtDNA evidence build was canceled.")
            manifest["status"]=DatabaseState.READY;manifest["completed_at"]=_now()
            (stage/"manifest.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
            self._atomic_json(self.active_pointer,{"snapshot_id":sid,"activated_at":_now()})
            return final
        except Exception as error:
            if stage.is_dir():
                (stage/"manifest.json").write_text(json.dumps({"status":DatabaseState.CANCELLED if isinstance(error,MtdnaBuildCancelled) else DatabaseState.ERROR,
                    "snapshot_id":sid,"last_error":str(error)},indent=2),encoding="utf-8")
            raise

    @staticmethod
    def _record(**fields):
        return {name:_text(fields.get(name,"")) for name in EVIDENCE_COLUMNS}

    def _mitocarta_rows(self,snapshot,sid):
        genes=_read_source_table(snapshot,"mitocarta_genes.csv")
        members=_read_source_table(snapshot,"mitopathway_membership.csv")
        hierarchy=_read_source_table(snapshot,"mitopathway_hierarchy.csv")
        paths=dict(zip(hierarchy.mitopathway,hierarchy.pathway_full))
        by_symbol=defaultdict(list)
        for row in genes.to_dict("records"):
            if _text(row.get("is_mitocarta")).lower() in ("true","1"):by_symbol[_text(row.get("gene_symbol"))].append(row)
        result=[]
        for row in members.to_dict("records"):
            path=paths.get(row["mitopathway"],"");category=mitocarta_category(path)
            if not category:continue
            for gene in by_symbol.get(row["gene_symbol"],[]):
                result.append(self._record(ncbi_gene_id=_gene_id(gene.get("ncbi_gene_id")),gene_symbol=row["gene_symbol"],
                    uniprot=gene.get("uniprot_accession"),source="MitoCarta",source_snapshot=sid,
                    evidence_class=dict((c[0],c[2]) for c in CATEGORIES)[category],evidence_category=category,
                    source_identifier=row["mitopathway"],source_name=path,source_path=path,
                    source_detail=f"Evidence: {_text(gene.get('evidence'))}; Subcompartment: {_text(gene.get('sub_compartment_raw'))}"))
        return result

    def _go_rows(self,snapshot,sid):
        terms=_read_source_table(snapshot,"go_terms.csv")
        relations_frame=_read_source_table(snapshot,"go_relations.csv")
        annotations=_read_source_table(snapshot,"go_annotations.csv")
        names=dict(zip(terms.go_id,terms.name))
        relations=defaultdict(list)
        for row in relations_frame.to_dict("records"):
            relations[row["child_go_id"]].append((row["parent_go_id"],row["relation"]))
        records=[];excluded=[]
        for annotation in annotations.to_dict("records"):
            qualifier=_text(annotation.get("qualifier"))
            annotated=_text(annotation.get("go_id"))
            paths={ancestor:path for ancestor,path in go_ancestor_paths(annotated,relations)}
            if "NOT" in re.split(r"[|,\s]+",qualifier):
                if any(seed in paths for seed in GO_SEEDS):excluded.append(annotation)
                continue
            for seed,(category,evidence_class) in GO_SEEDS.items():
                if seed not in paths:continue
                # Regulation and direct execution are distinct claims.
                if evidence_class!="regulation" and any(paths.get(reg) is not None for reg in ("GO:1901858","GO:1903108")):
                    continue
                records.append(self._record(ncbi_gene_id=_gene_id(annotation.get("ncbi_gene_id")),
                    gene_symbol=annotation.get("gene_symbol"),uniprot=annotation.get("uniprot") or (annotation.get("db_object_id") if annotation.get("database")=="UniProtKB" else ""),
                    source="Gene Ontology",source_snapshot=sid,evidence_class=evidence_class,
                    evidence_category=category,source_identifier=annotated,source_name=names.get(annotated,""),
                    source_detail=" > ".join(paths[seed]),evidence_code=annotation.get("evidence_code"),
                    reference=annotation.get("reference"),qualifier=qualifier,assigned_by=annotation.get("assigned_by"),
                    seed_go_id=seed,annotated_go_id=annotated,relation_path=" > ".join(paths[seed])))
        return records,excluded

    @staticmethod
    def _reconcile(records):
        symbol_ids=defaultdict(set)
        for row in records:
            if row["gene_symbol"] and row["ncbi_gene_id"]:symbol_ids[row["gene_symbol"].upper()].add(row["ncbi_gene_id"])
        for row in records:
            identifier=row["ncbi_gene_id"];symbol=row["gene_symbol"].upper()
            if not identifier and len(symbol_ids[symbol])==1:identifier=next(iter(symbol_ids[symbol]))
            row["ncbi_gene_id"]=identifier
            row["entity_key"]=f"NCBI:{identifier}" if identifier else f"SYMBOL:{symbol}"
        return records

    @staticmethod
    def _frames(records,excluded):
        evidence=pd.DataFrame(records,columns=EVIDENCE_COLUMNS).drop_duplicates()
        categories=pd.DataFrame(CATEGORIES,columns=("category","display_name","evidence_class","description"))
        membership=evidence[["entity_key","evidence_category"]].drop_duplicates().rename(columns={"evidence_category":"category"})
        entities=[];agreement=[]
        for key,group in evidence.groupby("entity_key",sort=True):
            sources=set(group.source);cats=sorted(set(group.evidence_category));uniprots=sorted({u for value in group.uniprot for u in _uniprots(value)})
            entity={"entity_key":key,"ncbi_gene_id":next((v for v in group.ncbi_gene_id if v),""),
                "gene_symbol":next((v for v in group.gene_symbol if v),""),"description":next((v for v in group.source_name if v),""),
                "uniprot_accessions":";".join(uniprots),"evidence_source_count":len(sources),
                "evidence_record_count":len(group),"category_count":len(cats),
                "has_mtdna_related_evidence":True,"is_mtdna_encoded":"mtDNA_encoded" in cats,
                "symbol_fallback":key.startswith("SYMBOL:")}
            entities.append(entity)
            agreement.append({"entity_key":key,"mitocarta_evidence":"MitoCarta" in sources,
                "go_evidence":"Gene Ontology" in sources,"mtdna_encoded":"NCBI mtDNA genome" in sources,
                "source_count":len(sources),"categories":";".join(cats)})
        return {"entities":pd.DataFrame(entities),"evidence_records":evidence,"category_membership":membership,
            "categories":categories,"mitocarta_evidence":evidence[evidence.source=="MitoCarta"],
            "go_evidence":evidence[evidence.source=="Gene Ontology"],
            "mtdna_encoded":evidence[evidence.source=="NCBI mtDNA genome"],
            "source_agreement":pd.DataFrame(agreement),"go_exclusions":pd.DataFrame(excluded,columns=("ncbi_gene_id","gene_symbol","uniprot","go_id","evidence_code","reference","qualifier","assigned_by")),
            "metadata":pd.DataFrame([{"key":"schema_version","value":"1"}])}

    @staticmethod
    def _write_sqlite(snapshot,frames,manifest):
        path=snapshot/"mtdna_evidence.sqlite"
        with sqlite3.connect(path) as conn:
            for name,frame in frames.items():frame.to_sql(name,conn,index=False,if_exists="replace")
            for table,column in (("entities","entity_key"),("entities","ncbi_gene_id"),("entities","gene_symbol"),
                ("category_membership","entity_key"),("category_membership","category"),
                ("evidence_records","entity_key"),("evidence_records","uniprot"),
                ("evidence_records","source"),("evidence_records","source_identifier"),
                ("evidence_records","source_path")):
                conn.execute(f'CREATE INDEX "idx_{table}_{column}" ON "{table}" ("{column}")')

    def validate_snapshot(self,snapshot=None,allow_incomplete=False):
        snap=Path(snapshot or self.active_snapshot() or "")
        manifest=self.manifest(snap)
        if manifest.get("status") not in ({DatabaseState.READY,DatabaseState.INCOMPLETE} if allow_incomplete else {DatabaseState.READY}):
            raise MtdnaEvidenceError("mtDNA evidence snapshot is not Ready.")
        if manifest.get("tax_id")!="9606" or manifest.get("ncbi_accession")!=NCBI_ACCESSION:
            raise MtdnaEvidenceError("Invalid mtDNA evidence organism or NCBI reference.")
        raw=snap/"raw"/"NC_012920.1.xml";sqlite=snap/"mtdna_evidence.sqlite"
        if not raw.is_file() or _sha(raw)!=manifest.get("ncbi_raw_sha256") or not sqlite.is_file() or _sha(sqlite)!=manifest.get("sqlite_sha256"):
            raise MtdnaEvidenceError("mtDNA evidence raw or SQLite hash mismatch.")
        parse_ncbi_refseq(raw.read_bytes())
        for key,path_key in (("mitocarta_snapshot_id","mitocarta_snapshot_path"),("go_snapshot_id","go_snapshot_path")):
            path=Path(manifest.get(path_key,""))
            if path.name!=manifest.get(key) or not path.is_dir():raise MtdnaEvidenceError(f"Missing or mismatched source snapshot: {key}")
        mito=Path(manifest["mitocarta_snapshot_path"])
        if _sha(mito/"manifest.json")!=manifest.get("mitocarta_manifest_sha256"):
            raise MtdnaEvidenceError("MitoCarta source manifest hash mismatch.")
        if json.loads((mito/"manifest.json").read_text(encoding="utf-8")).get("status")!=DatabaseState.READY:
            raise MtdnaEvidenceError("MitoCarta source snapshot is not Ready.")
        go=Path(manifest["go_snapshot_path"])
        if _sha(go/"manifest.json")!=manifest.get("go_manifest_sha256"):
            raise MtdnaEvidenceError("GO source manifest hash mismatch.")
        go_sources=self._go_manifest(go).get("sources",{})
        for name,info in go_sources.items():
            path=go/"raw"/name
            if not path.is_file() or _sha(path)!=info.get("sha256"):
                raise MtdnaEvidenceError("GO raw source hash mismatch.")
        if self._go_manifest(go).get("status")!=DatabaseState.READY:
            raise MtdnaEvidenceError("GO source snapshot is not Ready.")
        for name,expected in manifest.get("mitocarta_source_hashes",{}).items():
            path=mito/"raw"/name
            if not path.is_file() or _sha(path)!=expected:raise MtdnaEvidenceError("MitoCarta source hash mismatch.")
        for name,expected in manifest.get("go_source_hashes",{}).items():
            if _sha(Path(manifest["go_snapshot_path"])/"tables"/name)!=expected:raise MtdnaEvidenceError("GO source hash mismatch.")
        if not manifest.get("mitocarta_source_hashes") or not manifest.get("go_source_hashes"):
            raise MtdnaEvidenceError("Source hashes are missing.")
        with sqlite3.connect(sqlite) as conn:
            present={row[0] for row in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
            if not set(TABLES).issubset(present):raise MtdnaEvidenceError("mtDNA evidence SQLite tables are incomplete.")
            if conn.execute("SELECT COUNT(*) FROM evidence_records").fetchone()[0]<1 or conn.execute("SELECT COUNT(*) FROM categories").fetchone()[0]<1:
                raise MtdnaEvidenceError("mtDNA evidence index has no evidence or categories.")
        return True

    @staticmethod
    def _atomic_json(path,payload):
        path.parent.mkdir(parents=True,exist_ok=True);part=path.with_name(path.name+".part")
        part.write_text(json.dumps(payload,indent=2),encoding="utf-8");os.replace(part,path)

    def _query(self,sql,parameters=(),snapshot=None):
        path=(Path(snapshot)/"mtdna_evidence.sqlite") if snapshot else self.database_path()
        if not path or not path.is_file():raise MtdnaEvidenceError("mtDNA evidence database is not installed.")
        with sqlite3.connect(path) as conn:return pd.read_sql_query(sql,conn,params=parameters)

    def lookup_gene(self,value,snapshot=None):
        value=str(value).strip()
        return self._query("SELECT * FROM entities WHERE ncbi_gene_id=? OR UPPER(gene_symbol)=UPPER(?) OR entity_key=?",
                           (value,value,value),snapshot)

    def lookup_uniprot(self,accession,snapshot=None):
        accession=str(accession).strip()
        frame=self._query("SELECT DISTINCT e.* FROM entities e JOIN evidence_records r ON e.entity_key=r.entity_key WHERE (';'||REPLACE(r.uniprot, '|', ';')||';') LIKE ?",
                          (f"%;{accession};%",),snapshot)
        return {"status":"ambiguous" if len(frame)>1 else "unique" if len(frame)==1 else "unmapped","entities":frame}

    def get_evidence_for_entity(self,key,snapshot=None):return self._query("SELECT * FROM evidence_records WHERE entity_key=?",(key,),snapshot)
    def get_categories_for_entity(self,key,snapshot=None):return self._query("SELECT * FROM category_membership WHERE entity_key=?",(key,),snapshot)
    def get_entities_for_category(self,category,snapshot=None):return self._query("SELECT e.* FROM entities e JOIN category_membership c ON e.entity_key=c.entity_key WHERE c.category=?",(category,),snapshot)
    def get_mitocarta_evidence(self,key,snapshot=None):return self._query("SELECT * FROM mitocarta_evidence WHERE entity_key=?",(key,),snapshot)
    def get_go_evidence(self,key,snapshot=None):return self._query("SELECT * FROM go_evidence WHERE entity_key=?",(key,),snapshot)
    def is_mtdna_encoded(self,key,snapshot=None):return not self._query("SELECT entity_key FROM mtdna_encoded WHERE entity_key=? LIMIT 1",(key,),snapshot).empty
    def source_agreement(self,key=None,snapshot=None):return self._query("SELECT * FROM source_agreement"+(" WHERE entity_key=?" if key else ""),(key,) if key else (),snapshot)

    def search_entities(self,term,snapshot=None):
        needle=f"%{str(term).strip()}%"
        return self._query("SELECT DISTINCT e.* FROM entities e LEFT JOIN evidence_records r ON r.entity_key=e.entity_key "
            "WHERE UPPER(e.gene_symbol) LIKE UPPER(?) OR e.ncbi_gene_id LIKE ? OR e.uniprot_accessions LIKE ? "
            "OR UPPER(r.source_name) LIKE UPPER(?) OR UPPER(r.source_path) LIKE UPPER(?) "
            "OR UPPER(r.evidence_category) LIKE UPPER(?) OR r.source_identifier LIKE ?",
            (needle,)*7,snapshot)
