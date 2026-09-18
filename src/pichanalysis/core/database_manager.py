from __future__ import annotations

import csv
import json
import logging
import os
import shutil
import threading
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Iterable

from platformdirs import user_data_path

from .database_registry import DatabaseState, KEGG_HSA
from .databases.kegg import (
    KEGGProvider, parse_gene_pathway_links, parse_genes, parse_pathways,
    validate_text, VALIDATORS, parse_conversion,
)
from .reactome_database import ReactomeDatabase
from .mitocarta_database import MitoCartaDatabase

PATHWAY_ANALYSIS_TABLES=("ncbi_geneid_to_kegg.tsv","kegg_to_ncbi_geneid.tsv","uniprot_to_kegg.tsv","kegg_to_uniprot.tsv")


ProgressCallback = Callable[[dict], None]


def default_database_root() -> Path:
    return Path(user_data_path("PichAnalysis", "PichAnalysis")) / "databases"


class DatabaseManager:
    def __init__(self, root: Path | None = None, provider: KEGGProvider | None = None):
        self.root = Path(root) if root else default_database_root()
        self.database_root = self.root / "kegg" / "hsa"
        self.snapshots_root = self.database_root / "snapshots"
        self.provider = provider or KEGGProvider()
        self.reactome = ReactomeDatabase(self.root)
        self.mitocarta = MitoCartaDatabase(self.root)
        self._cancel = threading.Event()
        self.logger = logging.getLogger("pichanalysis.database_manager")
        if not self.logger.handlers:
            log_dir = self.root.parent / "logs"
            log_dir.mkdir(parents=True, exist_ok=True)
            handler = logging.FileHandler(log_dir / "database_manager.log", encoding="utf-8")
            handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
            self.logger.addHandler(handler)
        self.logger.setLevel(logging.INFO)

    @property
    def active_pointer(self) -> Path:
        return self.database_root / "active.json"

    def cancel(self) -> None:
        self._cancel.set()

    def active_snapshot(self) -> Path | None:
        try:
            payload = json.loads(self.active_pointer.read_text(encoding="utf-8"))
            path = self.snapshots_root / payload["snapshot_id"]
            return path if path.is_dir() and self.manifest(path).get("status") == DatabaseState.READY else None
        except (OSError, KeyError, ValueError, json.JSONDecodeError):
            return None

    def manifest(self, snapshot: Path | None = None) -> dict:
        target = (snapshot or self.active_snapshot())
        if not target:
            return {}
        try:
            return json.loads((target / "manifest.json").read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return {}

    def state(self) -> DatabaseState:
        if self.active_snapshot():
            return DatabaseState.READY
        candidates = sorted(self.snapshots_root.glob("*/manifest.json"), reverse=True) if self.snapshots_root.exists() else []
        if not candidates:
            return DatabaseState.NOT_INSTALLED
        try:
            return DatabaseState(json.loads(candidates[0].read_text(encoding="utf-8")).get("status", DatabaseState.ERROR))
        except (OSError, ValueError, json.JSONDecodeError):
            return DatabaseState.ERROR

    def latest_incomplete_snapshot(self) -> Path | None:
        if not self.snapshots_root.exists():
            return None
        for manifest_path in sorted(self.snapshots_root.glob("*/manifest.json"), reverse=True):
            manifest = self.manifest(manifest_path.parent)
            if manifest.get("status") in {DatabaseState.INCOMPLETE, DatabaseState.CANCELLED, DatabaseState.ERROR}:
                return manifest_path.parent
        return None

    def paths(self) -> dict[str, Path | None]:
        snapshot = self.active_snapshot()
        return {
            "root": snapshot,
            "manifest": snapshot / "manifest.json" if snapshot else None,
            "tables": snapshot / "tables" if snapshot else None,
            "kgml": snapshot / "kgml" if snapshot else None,
            "images": snapshot / "images" if snapshot else None,
        }

    def pathway_analysis_compatibility(self) -> tuple[bool,str]:
        snapshot=self.active_snapshot()
        if snapshot is None:return False,"No active snapshot"
        if self.manifest(snapshot).get("status")!=DatabaseState.READY:return False,"Incomplete database"
        if not all((snapshot/"tables"/name).is_file() for name in PATHWAY_ANALYSIS_TABLES):
            return False,"Missing identifier mapping tables"
        return True,"Ready for pathway analysis"

    def download(self, components: Iterable[str], progress: ProgressCallback | None = None,
                 *, resume: bool = False, pathway_subset: Iterable[str] | None = None) -> Path:
        selected = tuple(dict.fromkeys(("core", *components)))
        snapshot = self.latest_incomplete_snapshot() if resume else None
        updating = self.active_snapshot() is not None and snapshot is None
        if snapshot is None:
            snapshot_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
            snapshot = self.snapshots_root / snapshot_id
        else:
            snapshot_id = snapshot.name
        snapshot.mkdir(parents=True, exist_ok=True)
        for folder in ("metadata", "tables", "entries", "kgml", "images"):
            (snapshot / folder).mkdir(exist_ok=True)
        self._cancel.clear()
        existing_manifest = self.manifest(snapshot)
        if resume and existing_manifest:
            selected = tuple(existing_manifest.get("components", selected))
        manifest = existing_manifest or {
            "schema_version": 1, "database_schema_version": 2, "database_id": "kegg", "display_name": KEGG_HSA.display_name,
            "organism": {"name": KEGG_HSA.organism_name, "code": "hsa", "taxonomy_id": "9606"},
            "snapshot_id": snapshot_id, "created_at": datetime.now(timezone.utc).isoformat(),
            "source_base_url": self.provider.base_url, "components": list(selected),
            "completed_items": [], "failed_items": [], "counts": {},
        }
        manifest["status"] = DatabaseState.UPDATING if updating else DatabaseState.DOWNLOADING
        manifest["updated_at"] = datetime.now(timezone.utc).isoformat()
        self._write_manifest(snapshot, manifest)
        self.logger.info("Starting snapshot %s with components %s", snapshot_id, selected)
        try:
            pathways, genes, links = self._download_core(snapshot, manifest, progress)
            wanted = list(pathway_subset) if pathway_subset is not None else [row["pathway_id"] for row in pathways]
            jobs = [(component, pathway_id) for component in selected if component != "core" for pathway_id in wanted]
            total = 5 + len(jobs)
            done = 5
            for component, pathway_id in jobs:
                if self._cancel.is_set():
                    manifest["status"] = DatabaseState.CANCELLED
                    self._write_manifest(snapshot, manifest)
                    self._emit(progress, manifest, done, total, f"Cancelled after {done} of {total} items")
                    return snapshot
                key = f"{component}:{pathway_id}"
                target = snapshot / component / self._filename(component, pathway_id)
                if key not in manifest["completed_items"] or not self._valid_existing(component, target):
                    try:
                        data = self.provider.fetch(self.provider.pathway_routes(pathway_id)[component])
                        VALIDATORS[component](data)
                        self._atomic_bytes(target, data)
                        if key not in manifest["completed_items"]:
                            manifest["completed_items"].append(key)
                        manifest["failed_items"] = [item for item in manifest["failed_items"] if item.get("item") != key]
                    except Exception as error:
                        manifest["failed_items"].append({"item": key, "error": str(error)})
                        self.logger.exception("Failed item %s", key)
                done += 1
                self._write_manifest(snapshot, manifest)
                self._emit(progress, manifest, done, total, f"Downloaded {done} of {total} items")
            manifest["counts"] = {"pathways": len(pathways), "genes": len(genes), "gene_pathway_links": len(links), "items": done}
            manifest["status"] = DatabaseState.INCOMPLETE if manifest["failed_items"] else DatabaseState.READY
            manifest["completed_at"] = datetime.now(timezone.utc).isoformat()
            self._write_manifest(snapshot, manifest)
            if manifest["status"] == DatabaseState.READY:
                self._atomic_json(self.active_pointer, {"snapshot_id": snapshot_id, "activated_at": manifest["completed_at"]})
            self._emit(progress, manifest, done, total, str(manifest["status"]))
            return snapshot
        except Exception as error:
            manifest["status"] = DatabaseState.ERROR
            manifest["last_error"] = str(error)
            self._write_manifest(snapshot, manifest)
            self.logger.exception("Snapshot %s failed", snapshot_id)
            raise

    def _download_core(self, snapshot: Path, manifest: dict, progress: ProgressCallback | None):
        definitions = (
            ("core:pathways", "list/pathway/hsa", "metadata/pathways.raw.tsv", parse_pathways),
            ("core:genes", "list/hsa", "metadata/genes.raw.tsv", parse_genes),
            ("core:links", "link/pathway/hsa", "metadata/gene_to_pathway.raw.tsv", parse_gene_pathway_links),
            ("core:ncbi", "conv/ncbi-geneid/hsa", "metadata/ncbi_geneid.raw.tsv", lambda x:parse_conversion(x,"ncbi-geneid:")),
            ("core:uniprot", "conv/uniprot/hsa", "metadata/uniprot.raw.tsv", lambda x:parse_conversion(x,"up:")),
        )
        parsed = []
        for index, (key, route, relative, parser) in enumerate(definitions, 1):
            raw_path = snapshot / relative
            if key not in manifest["completed_items"] or not raw_path.is_file():
                data = self.provider.fetch(route)
                validate_text(data)
                parser(data)
                self._atomic_bytes(raw_path, data)
                if key not in manifest["completed_items"]:
                    manifest["completed_items"].append(key)
            data = raw_path.read_bytes()
            parsed.append(parser(data))
            self._write_manifest(snapshot, manifest)
            self._emit(progress, manifest, index, 5, f"Downloaded core table {index} of 5")
            if self._cancel.is_set():
                manifest["status"] = DatabaseState.CANCELLED
                self._write_manifest(snapshot, manifest)
                return (*parsed, *([] for _ in range(5 - len(parsed))))[:3]
        pathways, genes, links, ncbi, uniprot = parsed
        self._write_tsv(snapshot / "tables/pathways.tsv", pathways, ("pathway_id", "name", "organism_code"))
        self._write_tsv(snapshot / "tables/genes.tsv", genes, ("gene_id", "kegg_gene_id", "symbol", "raw_description"))
        forward = [{"gene_id": gene, "pathway_id": pathway} for gene, pathway in links]
        inverse = [{"pathway_id": pathway, "gene_id": gene} for gene, pathway in links]
        self._write_tsv(snapshot / "tables/gene_to_pathway.tsv", forward, ("gene_id", "pathway_id"))
        self._write_tsv(snapshot / "tables/pathway_to_gene.tsv", inverse, ("pathway_id", "gene_id"))
        self._write_tsv(snapshot/"tables/ncbi_geneid_to_kegg.tsv",[{"ncbi_gene_id":x,"kegg_gene_id":k} for x,k in ncbi],("ncbi_gene_id","kegg_gene_id"))
        self._write_tsv(snapshot/"tables/kegg_to_ncbi_geneid.tsv",[{"kegg_gene_id":k,"ncbi_gene_id":x} for x,k in ncbi],("kegg_gene_id","ncbi_gene_id"))
        self._write_tsv(snapshot/"tables/uniprot_to_kegg.tsv",[{"uniprot_accession":x,"kegg_gene_id":k} for x,k in uniprot],("uniprot_accession","kegg_gene_id"))
        self._write_tsv(snapshot/"tables/kegg_to_uniprot.tsv",[{"kegg_gene_id":k,"uniprot_accession":x} for x,k in uniprot],("kegg_gene_id","uniprot_accession"))
        try:self._atomic_bytes(snapshot/"metadata/kegg_info.txt",self.provider.fetch("info/kegg"))
        except Exception as error:self.logger.warning("Could not download KEGG release information: %s",error)
        return pathways, genes, links

    @staticmethod
    def _filename(component: str, pathway_id: str) -> str:
        return pathway_id + {"entries": ".txt", "kgml": ".xml", "images": ".png"}[component]

    @staticmethod
    def _valid_existing(component: str, path: Path) -> bool:
        if not path.is_file():
            return False
        try:
            VALIDATORS[component](path.read_bytes())
            return True
        except (OSError, ValueError):
            return False

    @staticmethod
    def _atomic_bytes(path: Path, data: bytes) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        part = path.with_name(path.name + ".part")
        part.write_bytes(data)
        os.replace(part, path)

    @staticmethod
    def _atomic_json(path: Path, payload: dict) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        part = path.with_name(path.name + ".part")
        part.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
        os.replace(part, path)

    def _write_manifest(self, snapshot: Path, manifest: dict) -> None:
        manifest["updated_at"] = datetime.now(timezone.utc).isoformat()
        self._atomic_json(snapshot / "manifest.json", manifest)

    def _write_tsv(self, path: Path, rows: list[dict], fields: tuple[str, ...]) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        part = path.with_name(path.name + ".part")
        with part.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(part, path)

    @staticmethod
    def _emit(callback: ProgressCallback | None, manifest: dict, done: int, total: int, message: str) -> None:
        if callback:
            callback({"done": done, "total": total, "message": message, "status": str(manifest["status"]), "snapshot_id": manifest["snapshot_id"]})

    def disk_info(self) -> tuple[int, int, int]:
        anchor = self.root if self.root.exists() else self.root.parent
        anchor.mkdir(parents=True, exist_ok=True)
        return shutil.disk_usage(anchor)
