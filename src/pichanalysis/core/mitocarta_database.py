from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
from datetime import datetime, timezone
from tempfile import mkdtemp
from pathlib import Path

import pandas as pd

from .database_registry import DatabaseState
from .databases.mitocarta import BASE_URL, FILES, VERSION, validate_downloads

REQUIRED_COLUMNS = {"HumanGeneID", "Symbol", "Description", "MitoCarta3.0_List"}
COLUMN_MAP = {
    "Symbol": "gene_symbol", "Description": "gene_description", "HumanGeneID": "ncbi_gene_id",
    "UniProt": "uniprot_accession", "Synonyms": "synonyms", "MitoCarta2.0_Score": "maestro_score",
    "MitoCarta3.0_Evidence": "evidence", "MitoCarta3.0_SubMitoLocalization": "sub_compartment_raw",
    "MitoCarta3.0_MitoPathways": "mitopathways_raw", "Tissues": "tissues_raw",
}


def _now(): return datetime.now(timezone.utc).isoformat()


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(1024 * 1024): digest.update(chunk)
    return digest.hexdigest()


def inspect_workbook(path: Path) -> tuple[str, list[str], pd.DataFrame, list[str]]:
    engine = "openpyxl" if path.read_bytes()[:2] == b"PK" else "xlrd"
    try: excel = pd.ExcelFile(path, engine=engine)
    except Exception as exc: raise ValueError(f"Invalid MitoCarta workbook: {exc}") from exc
    candidates = []
    for sheet in excel.sheet_names:
        frame = pd.read_excel(path, sheet_name=sheet, engine=engine)
        columns = [str(column).strip() for column in frame.columns]
        frame.columns = columns
        if REQUIRED_COLUMNS.issubset(columns):
            positives = frame["MitoCarta3.0_List"].fillna("").astype(str).str.strip().str.casefold().eq("mitocarta3.0").sum()
            candidates.append((len(frame), positives < len(frame), sheet, columns, frame))
    if not candidates:
        raise ValueError("Invalid MitoCarta workbook: no sheet contains the required official headers: " + ", ".join(sorted(REQUIRED_COLUMNS)))
    _, _, sheet, columns, frame = max(candidates, key=lambda item: (item[1], item[0]))
    return sheet, columns, frame, list(excel.sheet_names)


def parse_workbook(path: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    sheet, columns, source, sheets = inspect_workbook(path)
    missing = [name for name in COLUMN_MAP if name not in source]
    selected = {old: new for old, new in COLUMN_MAP.items() if old in source}
    genes = source[list(selected)].rename(columns=selected).copy()
    genes["is_mitocarta"] = source["MitoCarta3.0_List"].fillna("").astype(str).str.strip().str.casefold().eq("mitocarta3.0")
    genes = genes[genes["gene_symbol"].notna() & genes["gene_symbol"].astype(str).str.strip().ne("")].copy()
    genes["gene_symbol"] = genes["gene_symbol"].astype(str).str.strip()
    memberships = []
    for row in genes.loc[genes.is_mitocarta, ["gene_symbol", "sub_compartment_raw"]].itertuples(index=False):
        if pd.isna(row.sub_compartment_raw): continue
        for compartment in re.split(r"\s*\|\s*", str(row.sub_compartment_raw).strip()):
            if compartment: memberships.append({"gene_symbol": row.gene_symbol, "subcompartment": compartment})
    subcompartments = pd.DataFrame(memberships, columns=["gene_symbol", "subcompartment"]).drop_duplicates()
    metadata = {"workbook_sheet_used": sheet, "workbook_sheets": sheets, "source_columns": [*selected, "MitoCarta3.0_List"], "workbook_columns": columns,
                "unused_documented_columns": missing, "total_workbook_rows": len(source),
                "mitocarta_gene_count": int(genes.is_mitocarta.sum())}
    return genes, subcompartments, metadata


def parse_gmx(path: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    rows = [line.rstrip("\r\n").split("\t") for line in path.read_text(encoding="utf-8-sig").splitlines()]
    if len(rows) < 3 or len(rows[0]) != len(rows[1]): raise ValueError("Invalid MitoPathways GMX structure.")
    width = len(rows[0]); membership, hierarchy = [], []
    for column in range(width):
        name = rows[0][column].strip(); description = rows[1][column].strip()
        if not name: raise ValueError(f"Invalid MitoPathways GMX: empty pathway name in column {column + 1}.")
        raw_levels = [part.strip() for part in (description.split(" > ") if " > " in description else description.split(".")) if part.strip()]
        full = " > ".join(level.replace("_", " ") for level in raw_levels) if raw_levels else name.replace("_", " ")
        hierarchy_row = {"mitopathway": name, "pathway_full": full}
        for index, level in enumerate(raw_levels, 1): hierarchy_row[f"pathway_level_{index}"] = level.replace("_", " ")
        hierarchy.append(hierarchy_row)
        for row in rows[2:]:
            gene = row[column].strip() if column < len(row) else ""
            if gene: membership.append({"mitopathway": name, "gene_symbol": gene})
    members = pd.DataFrame(membership, columns=["mitopathway", "gene_symbol"]).drop_duplicates().sort_values(["mitopathway", "gene_symbol"])
    hierarchy_frame = pd.DataFrame(hierarchy).sort_values("mitopathway")
    return members, hierarchy_frame, {"mitopathway_count": width, "pathway_membership_count": len(members)}


class MitoCartaDatabase:
    def __init__(self, root: Path):
        self.root = Path(root) / "mitocarta" / "human"; self.snapshots = self.root / "snapshots"; self.active_pointer = self.root / "active_snapshot.json"

    def create_staging_snapshot(self) -> Path:
        self.snapshots.mkdir(parents=True, exist_ok=True)
        snap = Path(mkdtemp(prefix=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ_"), dir=self.snapshots))
        sid = snap.name
        (snap / "raw").mkdir(); (snap / "tables").mkdir()
        self._write_manifest(snap, {"database":"MitoCarta", "version":VERSION, "organism_name":"Homo sapiens", "tax_id":"9606", "snapshot_id":sid, "status":DatabaseState.INCOMPLETE, "download_started_at":_now(), "download_completed_at":None, "retrieved_at":None, "source":BASE_URL, "files":[]})
        return snap

    def install_from_directory(self, source: Path) -> Path:
        snap = self.create_staging_snapshot()
        try:
            for name in FILES: shutil.copy2(Path(source) / name, snap / "raw" / name)
            return self.finalize_snapshot(snap)
        except Exception as exc:
            self.mark_error(snap, str(exc)); raise

    def finalize_snapshot(self, snap: Path, file_info: list[dict] | None = None) -> Path:
        try:
            raw, tables = Path(snap) / "raw", Path(snap) / "tables"; validate_downloads(raw)
            genes, compartments, workbook_meta = parse_workbook(raw / FILES[0]); membership, hierarchy, gmx_meta = parse_gmx(raw / FILES[1])
            genes.to_csv(tables / "mitocarta_genes.csv", index=False); compartments.to_csv(tables / "subcompartment_membership.csv", index=False)
            membership.to_csv(tables / "mitopathway_membership.csv", index=False); hierarchy.to_csv(tables / "mitopathway_hierarchy.csv", index=False)
            files = file_info or [{"filename":name,"source_url":f"{BASE_URL}/{name}","size":(raw/name).stat().st_size,"sha256":_sha(raw/name)} for name in FILES]
            completed = _now(); manifest = self.manifest(snap); manifest.update(status=DatabaseState.READY, download_completed_at=completed, retrieved_at=completed, files=files, **workbook_meta, **gmx_meta, subcompartment_membership_count=len(compartments)); self._write_manifest(snap, manifest)
            self._atomic_json(self.active_pointer, {"snapshot_id":Path(snap).name,"activated_at":completed}); return Path(snap)
        except Exception as exc:
            self.mark_error(Path(snap), str(exc)); raise

    def mark_downloading(self, snap: Path):
        manifest=self.manifest(snap);manifest["status"]=DatabaseState.DOWNLOADING;self._write_manifest(snap,manifest)
    def mark_incomplete(self, snap: Path):
        manifest=self.manifest(snap);manifest["status"]=DatabaseState.INCOMPLETE;self._write_manifest(snap,manifest)
    def mark_error(self, snap: Path, message: str):
        manifest=self.manifest(snap);manifest.update(status=DatabaseState.ERROR,last_error=message);self._write_manifest(snap,manifest)
    def active_snapshot(self):
        try:
            snap=self.snapshots/json.loads(self.active_pointer.read_text(encoding="utf-8"))["snapshot_id"]
            return snap if self.manifest(snap).get("status")==DatabaseState.READY else None
        except (OSError,KeyError,json.JSONDecodeError):return None
    def manifest(self,snapshot:Path|None=None):
        try:return json.loads(((snapshot or self.active_snapshot())/"manifest.json").read_text(encoding="utf-8"))
        except (OSError,TypeError,json.JSONDecodeError):return {}
    def state(self):
        if self.is_ready():return DatabaseState.READY
        manifests=sorted(self.snapshots.glob("*/manifest.json"),reverse=True) if self.snapshots.exists() else []
        if not manifests:return DatabaseState.NOT_INSTALLED
        try:return DatabaseState(json.loads(manifests[0].read_text(encoding="utf-8")).get("status",DatabaseState.ERROR))
        except (OSError,ValueError,json.JSONDecodeError):return DatabaseState.ERROR
    def is_ready(self):return self.active_snapshot() is not None
    def version(self):return self.manifest().get("version") if self.is_ready() else None
    def table_path(self,name:str):
        allowed={"genes":"mitocarta_genes.csv","pathways":"mitopathway_membership.csv","subcompartments":"subcompartment_membership.csv","hierarchy":"mitopathway_hierarchy.csv"};snap=self.active_snapshot();return snap/"tables"/allowed[name] if snap and name in allowed else None
    def genes_path(self):return self.table_path("genes")
    def pathway_membership_path(self):return self.table_path("pathways")
    def subcompartment_membership_path(self):return self.table_path("subcompartments")
    def source_hashes(self):return {item["filename"]:item["sha256"] for item in self.manifest().get("files",[])}
    def gene_count(self):return int(self.manifest().get("mitocarta_gene_count",0))
    def pathway_count(self):return int(self.manifest().get("mitopathway_count",0))
    @staticmethod
    def _atomic_json(path,payload):
        path.parent.mkdir(parents=True,exist_ok=True);part=path.with_name(path.name+".part");part.write_text(json.dumps(payload,indent=2),encoding="utf-8");os.replace(part,path)
    def _write_manifest(self,snap,payload):self._atomic_json(Path(snap)/"manifest.json",payload)
