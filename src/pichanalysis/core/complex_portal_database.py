"""Immutable local snapshots of curated human Complex Portal ComplexTab.

This module stores source-derived facts and supports lookup only. It does not
infer complex presence, coverage, enrichment, or protein abundance.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import shutil
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

from .database_registry import DatabaseState
from .databases.complex_portal import (
    CITATION_URL, CURATED_HUMAN_URL, TERMS_URL, ComplexPortalCancelled,
    ComplexPortalProvider,
)

ACCESSION = re.compile(r"^(CPX-\d+)(?:\.(\d+))?$")
UNIPROT = re.compile(r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$")
TOKEN = re.compile(r"^(.+?)\((\d+)\)$")
REQUIRED = {
    "complex ac", "recommended name", "taxonomy identifier",
    "identifiers (and stoichiometry) of molecules in complex", "source",
}
TABLE_COLUMNS = {
    "complexes": ("complex_id", "complex_accession_raw", "complex_version", "recommended_name", "aliases", "taxonomy_id", "confidence_accession", "confidence_name", "experimental_evidence", "description", "properties", "assembly", "ligands", "disease", "agonist", "antagonist", "comment", "source", "cross_references", "go_annotations", "source_row_json"),
    "complex_aliases": ("complex_id", "alias"),
    "direct_participants": ("complex_id", "participant_accession", "participant_type", "stoichiometry_raw", "stoichiometry", "stoichiometry_known", "occurrence_index"),
    "expanded_protein_components": ("complex_id", "uniprot_accession", "participant_accession_raw", "participant_type", "stoichiometry_raw", "stoichiometry", "stoichiometry_known", "occurrence_index"),
    "protein_membership": ("complex_id", "uniprot_accession", "stoichiometry_raw", "stoichiometry", "stoichiometry_known"),
    "nonprotein_participants": ("complex_id", "participant_accession", "participant_type", "stoichiometry_raw", "stoichiometry", "stoichiometry_known", "occurrence_index"),
    "nested_complexes": ("parent_complex_id", "child_complex_id", "stoichiometry_raw", "stoichiometry", "stoichiometry_known", "occurrence_index"),
}
INDEXES = {
    "idx_complex_id": "CREATE INDEX idx_complex_id ON complexes(complex_id)",
    "idx_complex_name": "CREATE INDEX idx_complex_name ON complexes(recommended_name)",
    "idx_complex_alias": "CREATE INDEX idx_complex_alias ON complex_aliases(alias)",
    "idx_direct_accession": "CREATE INDEX idx_direct_accession ON direct_participants(participant_accession)",
    "idx_member_uniprot": "CREATE INDEX idx_member_uniprot ON protein_membership(uniprot_accession)",
    "idx_member_complex": "CREATE INDEX idx_member_complex ON protein_membership(complex_id)",
    "idx_nonprotein_complex": "CREATE INDEX idx_nonprotein_complex ON nonprotein_participants(complex_id)",
    "idx_nested_parent": "CREATE INDEX idx_nested_parent ON nested_complexes(parent_complex_id)",
}


def _now():
    return datetime.now(timezone.utc).isoformat()


def _sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _atomic_json(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(path.name + ".part")
    part.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    os.replace(part, path)


def _field(row, *names):
    by_name = {key.lstrip("#").strip().casefold(): value for key, value in row.items() if key}
    for name in names:
        value = by_name.get(name.casefold(), "")
        if value and value != "-":
            return value.strip()
    return ""


def _parts(value):
    return [part.strip() for part in value.split("|") if part.strip() and part.strip() != "-"]


def _stoich(token):
    match = TOKEN.fullmatch(token.strip())
    if not match:
        raise ValueError(f"Invalid ComplexTab participant or stoichiometry: {token!r}")
    accession, raw = match.groups()
    count = int(raw)
    return accession.strip(), raw, count if count else None, bool(count)


def _type(accession):
    upper = accession.upper()
    if ACCESSION.fullmatch(upper):
        return "complex"
    if upper.startswith("CHEBI:"):
        return "chemical"
    if upper.startswith("URS") or upper.startswith("RNACENTRAL:"):
        return "RNA"
    if UNIPROT.fullmatch(upper):
        return "protein"
    return "other/unknown"


def canonicalize_uniprot_for_complex_portal(accession: str) -> dict:
    """Normalize a valid isoform query only; never rewrite experimental data."""
    original = str(accession).strip()
    match = re.fullmatch(r"(.+)-([1-9]\d*)", original.upper())
    if match and UNIPROT.fullmatch(match.group(1)):
        return {"input_accession": original, "lookup_accession": match.group(1), "isoform_normalized": True}
    return {"input_accession": original, "lookup_accession": original, "isoform_normalized": False}


class ComplexPortalDatabase:
    def __init__(self, root, provider=None):
        self.root = Path(root) / "complex_portal" / "human"
        self.snapshots = self.root / "snapshots"
        self.active_pointer = self.root / "active_snapshot.json"
        self.provider = provider or ComplexPortalProvider()

    def create_staging_snapshot(self):
        snapshot_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        snapshot = self.snapshots / snapshot_id
        (snapshot / "raw").mkdir(parents=True)
        (snapshot / "tables").mkdir()
        self._write(snapshot, {
            "database": "Complex Portal", "organism": "Homo sapiens", "tax_id": "9606",
            "scope": "manually_curated", "snapshot_id": snapshot_id,
            "status": DatabaseState.INCOMPLETE, "created_at": _now(),
            "complex_portal_release": "unknown", "source_url": CURATED_HUMAN_URL,
            "license": "not explicitly stated in retrieved metadata",
            "license_or_terms_url": TERMS_URL, "citation_url": CITATION_URL,
        })
        return snapshot

    def _write(self, snapshot, manifest):
        manifest["updated_at"] = _now()
        _atomic_json(snapshot / "manifest.json", manifest)

    def manifest(self, snapshot=None):
        target = snapshot or self.active_snapshot()
        if target is None:
            return {}
        try:
            return json.loads((Path(target) / "manifest.json").read_text(encoding="utf-8"))
        except (OSError, ValueError):
            return {}

    def active_snapshot(self):
        try:
            snapshot = self.snapshots / json.loads(self.active_pointer.read_text(encoding="utf-8"))["snapshot_id"]
            return snapshot if snapshot.is_dir() and self.manifest(snapshot).get("status") == DatabaseState.READY else None
        except (OSError, ValueError, KeyError):
            return None

    def is_ready(self):
        return self.active_snapshot() is not None

    def state(self):
        if self.is_ready():
            return DatabaseState.READY
        manifests = sorted(self.snapshots.glob("*/manifest.json"), reverse=True) if self.snapshots.is_dir() else []
        return DatabaseState(self.manifest(manifests[0].parent).get("status", DatabaseState.ERROR)) if manifests else DatabaseState.NOT_INSTALLED

    def database_path(self):
        snapshot = self.active_snapshot()
        return snapshot / "complex_portal.sqlite" if snapshot else None

    def source_hashes(self):
        item = self.manifest().get("raw_file", {})
        return {item["filename"]: item["sha256"]} if item else {}

    def mark_downloading(self, snapshot):
        manifest = self.manifest(snapshot)
        if manifest.get("status") == DatabaseState.READY:
            raise ValueError("Ready Complex Portal snapshots are immutable.")
        manifest["status"] = DatabaseState.DOWNLOADING
        self._write(snapshot, manifest)

    def mark_incomplete(self, snapshot):
        manifest = self.manifest(snapshot)
        if manifest.get("status") == DatabaseState.READY:
            raise ValueError("Ready Complex Portal snapshots are immutable.")
        manifest["status"] = DatabaseState.CANCELLED
        self._write(snapshot, manifest)

    def mark_error(self, snapshot, error):
        manifest = self.manifest(snapshot)
        if manifest.get("status") == DatabaseState.READY:
            raise ValueError("Ready Complex Portal snapshots are immutable.")
        manifest.update(status=DatabaseState.ERROR, last_error=str(error))
        self._write(snapshot, manifest)

    @staticmethod
    def _emit(progress, stage, message, **extra):
        if progress:
            progress({"stage": stage, "message": message, **extra})

    def download_and_build(self, progress=None, cancel_requested=None):
        snapshot = self.create_staging_snapshot()
        self.mark_downloading(snapshot)
        try:
            self._emit(progress, "download", "Downloading Complex Portal data...")
            info = self.provider.download(snapshot / "raw" / "9606.tsv",
                progress=lambda done, total: self._emit(progress, "download", "Downloading Complex Portal data...", bytes_downloaded=done, bytes_total=total),
                cancel_requested=cancel_requested)
            return self.finalize_snapshot(snapshot, info, progress, cancel_requested)
        except ComplexPortalCancelled:
            self.mark_incomplete(snapshot)
            raise
        except Exception as error:
            self.mark_error(snapshot, error)
            raise

    def install_from_file(self, source, *, source_url="fixture://curated-human", progress=None, cancel_requested=None):
        """Offline fixture importer; same validation and activation path as the real provider."""
        source = Path(source)
        snapshot = self.create_staging_snapshot()
        self.mark_downloading(snapshot)
        try:
            if "predicted" in source.name.casefold() or "predicted" in source_url.casefold():
                raise ValueError("Predicted Complex Portal data cannot be activated as a curated snapshot.")
            if cancel_requested and cancel_requested():
                raise ComplexPortalCancelled("Complex Portal installation was canceled.")
            target = snapshot / "raw" / "9606.tsv"
            part = target.with_name(target.name + ".part")
            shutil.copyfile(source, part)
            os.replace(part, target)
            info = {"filename": target.name, "source_url": source_url, "resolved_url": source_url,
                    "last_modified": None, "etag": None, "content_length": target.stat().st_size,
                    "size": target.stat().st_size, "sha256": _sha(target)}
            return self.finalize_snapshot(snapshot, info, progress, cancel_requested)
        except ComplexPortalCancelled:
            self.mark_incomplete(snapshot)
            raise
        except Exception as error:
            self.mark_error(snapshot, error)
            raise

    def finalize_snapshot(self, snapshot, source_info, progress=None, cancel_requested=None):
        snapshot = Path(snapshot)
        manifest = self.manifest(snapshot)
        if manifest.get("status") != DatabaseState.DOWNLOADING:
            raise ValueError("Complex Portal snapshot is not in a downloadable state.")
        raw = snapshot / "raw" / "9606.tsv"
        if not raw.is_file() or raw.stat().st_size == 0:
            raise ValueError("Curated human ComplexTab is missing or empty.")
        if _sha(raw) != source_info.get("sha256"):
            raise ValueError("Complex Portal raw source hash mismatch.")
        if "predicted" in str(source_info.get("source_url", "")).casefold():
            raise ValueError("Predicted Complex Portal data cannot be activated as a curated snapshot.")
        self._emit(progress, "validate", "Validating source file...")
        with raw.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            headers = reader.fieldnames or []
            normalized = {header.lstrip("#").strip().casefold() for header in headers}
            if not REQUIRED.issubset(normalized) or not ({"expanded participant list", "expanded component list"} & normalized):
                raise ValueError("Unrecognized curated human ComplexTab header.")
            manifest.update(retrieved_at=_now(), source_url=source_info["source_url"],
                resolved_url=source_info.get("resolved_url"), http_last_modified=source_info.get("last_modified"),
                http_etag=source_info.get("etag"), http_content_length=source_info.get("content_length"),
                detected_headers=headers, raw_file={"filename": raw.name, "size": raw.stat().st_size,
                    "sha256": source_info["sha256"], "detected_headers": headers})
            rows = {name: [] for name in TABLE_COLUMNS}
            memberships = {}
            ids = set()
            self._emit(progress, "parse", "Parsing complexes...")
            for index, row in enumerate(reader, 1):
                if cancel_requested and index % 100 == 0 and cancel_requested():
                    raise ComplexPortalCancelled("Complex Portal parsing was canceled.")
                if None in row:
                    raise ValueError(f"ComplexTab row {index} has extra fields.")
                accession_raw = _field(row, "Complex ac")
                match = ACCESSION.fullmatch(accession_raw)
                if not match:
                    raise ValueError(f"Invalid Complex Portal accession at row {index}: {accession_raw!r}")
                complex_id, version = match.groups()
                if complex_id in ids:
                    raise ValueError(f"Duplicate Complex Portal accession: {complex_id}")
                ids.add(complex_id)
                taxonomy = _field(row, "Taxonomy identifier")
                if taxonomy != "9606":
                    raise ValueError(f"Curated human ComplexTab contains non-human tax ID {taxonomy!r} at {complex_id}.")
                evidence = _field(row, "Evidence Code", "Confidence")
                if evidence.startswith(("ECO:0007653", "ECO:0008004")):
                    raise ValueError("Predicted Complex Portal data cannot be activated as a curated snapshot.")
                evidence_match = re.fullmatch(r"(ECO:\d+)\((.*)\)", evidence)
                source = _field(row, "Source")
                if any(label in source.casefold() for label in ("hu.map", "humap", "music predicted")):
                    raise ValueError("Predicted Complex Portal data cannot be activated as a curated snapshot.")
                aliases = _field(row, "Aliases for complex")
                rows["complexes"].append({
                    "complex_id": complex_id, "complex_accession_raw": accession_raw, "complex_version": version or "",
                    "recommended_name": _field(row, "Recommended name"), "aliases": aliases, "taxonomy_id": taxonomy,
                    "confidence_accession": evidence_match.group(1) if evidence_match else evidence,
                    "confidence_name": evidence_match.group(2) if evidence_match else "",
                    "experimental_evidence": _field(row, "Experimental evidence"),
                    "description": _field(row, "Description"), "properties": _field(row, "Complex properties"),
                    "assembly": _field(row, "Complex assembly"), "ligands": _field(row, "Ligand", "Ligands"),
                    "disease": _field(row, "Disease"), "agonist": _field(row, "Agonist"),
                    "antagonist": _field(row, "Antagonist"), "comment": _field(row, "Comment", "Comments"),
                    "source": source, "cross_references": _field(row, "Cross references"),
                    "go_annotations": _field(row, "Go Annotations"),
                    "source_row_json": json.dumps(row, ensure_ascii=False),
                })
                for alias in dict.fromkeys(_parts(aliases)):
                    rows["complex_aliases"].append({"complex_id": complex_id, "alias": alias})
                direct = _field(row, "Identifiers (and stoichiometry) of molecules in complex")
                for occurrence, token in enumerate(_parts(direct), 1):
                    accession, raw_count, count, known = _stoich(token)
                    kind = _type(accession)
                    record = {"complex_id": complex_id, "participant_accession": accession,
                        "participant_type": kind, "stoichiometry_raw": raw_count,
                        "stoichiometry": count, "stoichiometry_known": int(known), "occurrence_index": occurrence}
                    rows["direct_participants"].append(record)
                    if kind == "complex":
                        rows["nested_complexes"].append({"parent_complex_id": complex_id,
                            "child_complex_id": ACCESSION.fullmatch(accession).group(1),
                            "stoichiometry_raw": raw_count, "stoichiometry": count,
                            "stoichiometry_known": int(known), "occurrence_index": occurrence})
                    elif kind != "protein":
                        rows["nonprotein_participants"].append(record.copy())
                expanded = _field(row, "Expanded participant list", "Expanded component list")
                for occurrence, token in enumerate(_parts(expanded), 1):
                    accession, raw_count, count, known = _stoich(token)
                    kind = "protein_set" if accession.startswith("[") and accession.endswith("]") else _type(accession)
                    record = {"complex_id": complex_id,
                        "uniprot_accession": accession if kind == "protein" else "",
                        "participant_accession_raw": accession, "participant_type": kind,
                        "stoichiometry_raw": raw_count, "stoichiometry": count,
                        "stoichiometry_known": int(known), "occurrence_index": occurrence}
                    rows["expanded_protein_components"].append(record)
                    if kind == "protein":
                        key = (complex_id, accession)
                        prior = memberships.get(key)
                        if prior is not None and prior["stoichiometry_raw"] != raw_count:
                            raise ValueError(f"Conflicting stoichiometry for {complex_id} / {accession}.")
                        memberships[key] = {key_: record[key_] for key_ in TABLE_COLUMNS["protein_membership"]}
            rows["protein_membership"] = list(memberships.values())
        if not rows["complexes"] or not rows["protein_membership"]:
            raise ValueError("Curated human ComplexTab has no complexes or protein memberships.")
        self._emit(progress, "parse", "Parsing participants...")
        self._emit(progress, "membership", "Building protein membership...")
        for name, columns in TABLE_COLUMNS.items():
            if cancel_requested and cancel_requested():
                raise ComplexPortalCancelled("Complex Portal processing was canceled.")
            path = snapshot / "tables" / f"{name}.csv"
            with path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=columns)
                writer.writeheader()
                writer.writerows(rows[name])
        self._emit(progress, "index", "Indexing database...")
        database = snapshot / "complex_portal.sqlite"
        with sqlite3.connect(database) as conn:
            for name, columns in TABLE_COLUMNS.items():
                if cancel_requested and cancel_requested():
                    raise ComplexPortalCancelled("Complex Portal indexing was canceled.")
                definitions = ",".join(f'"{column}" TEXT' for column in columns)
                conn.execute(f'CREATE TABLE "{name}" ({definitions})')
                placeholders = ",".join("?" for _ in columns)
                conn.executemany(f'INSERT INTO "{name}" VALUES ({placeholders})',
                    (["" if record[column] is None else str(record[column]) for column in columns] for record in rows[name]))
            conn.execute("CREATE TABLE metadata(key TEXT PRIMARY KEY, value TEXT NOT NULL)")
            conn.executemany("INSERT INTO metadata VALUES (?, ?)", ((key, value) for key, value in
                (("database", "Complex Portal"), ("tax_id", "9606"), ("scope", "manually_curated"))))
            for sql in INDEXES.values():
                if cancel_requested and cancel_requested():
                    raise ComplexPortalCancelled("Complex Portal indexing was canceled.")
                conn.execute(sql)
        counts = {
            "complex_count": len(rows["complexes"]),
            "direct_participant_count": len(rows["direct_participants"]),
            "protein_membership_count": len(rows["protein_membership"]),
            "unique_protein_count": len({r["uniprot_accession"] for r in rows["protein_membership"]}),
            "nonprotein_participant_count": len(rows["nonprotein_participants"]),
            "nested_complex_count": len(rows["nested_complexes"]),
            "known_stoichiometry_count": sum(r["stoichiometry_known"] for r in rows["direct_participants"]),
            "unknown_stoichiometry_count": sum(not r["stoichiometry_known"] for r in rows["direct_participants"]),
        }
        manifest.update(**counts, sqlite_size=database.stat().st_size, sqlite_sha256=_sha(database),
            completed_at=_now(), status=DatabaseState.DOWNLOADING)
        self._write(snapshot, manifest)
        self._emit(progress, "validate", "Validating snapshot...")
        self.validate_snapshot(snapshot, allow_staging=True)
        if cancel_requested and cancel_requested():
            raise ComplexPortalCancelled("Complex Portal activation was canceled.")
        manifest["status"] = DatabaseState.READY
        self._write(snapshot, manifest)
        _atomic_json(self.active_pointer, {"snapshot_id": snapshot.name, "activated_at": _now()})
        return snapshot

    def validate_snapshot(self, snapshot=None, *, allow_staging=False):
        snapshot = Path(snapshot) if snapshot else self.active_snapshot()
        manifest = self.manifest(snapshot)
        if not snapshot or manifest.get("status") not in ({DatabaseState.READY, DatabaseState.DOWNLOADING} if allow_staging else {DatabaseState.READY}):
            raise ValueError("Complex Portal snapshot is not Ready.")
        if manifest.get("database") != "Complex Portal" or manifest.get("tax_id") != "9606" or manifest.get("scope") != "manually_curated":
            raise ValueError("Complex Portal manifest scope is invalid.")
        raw_info = manifest.get("raw_file", {})
        raw = snapshot / "raw" / raw_info.get("filename", "")
        if not raw.is_file() or raw.stat().st_size != raw_info.get("size") or _sha(raw) != raw_info.get("sha256"):
            raise ValueError("Complex Portal raw file hash mismatch.")
        if not REQUIRED.issubset({header.lstrip("#").strip().casefold() for header in manifest.get("detected_headers", [])}):
            raise ValueError("Complex Portal detected headers are incomplete.")
        database = snapshot / "complex_portal.sqlite"
        if not database.is_file() or _sha(database) != manifest.get("sqlite_sha256"):
            raise ValueError("Complex Portal SQLite hash mismatch.")
        with sqlite3.connect(database) as conn:
            tables = {row[0] for row in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
            indexes = {row[0] for row in conn.execute("SELECT name FROM sqlite_master WHERE type='index'")}
            if not (set(TABLE_COLUMNS) | {"metadata"}).issubset(tables) or not set(INDEXES).issubset(indexes):
                raise ValueError("Complex Portal SQLite tables or indexes are incomplete.")
            for name, key in (("complexes", "complex_count"), ("direct_participants", "direct_participant_count"),
                ("protein_membership", "protein_membership_count"), ("nonprotein_participants", "nonprotein_participant_count"),
                ("nested_complexes", "nested_complex_count")):
                if conn.execute(f'SELECT COUNT(*) FROM "{name}"').fetchone()[0] != manifest.get(key):
                    raise ValueError(f"Complex Portal {name} count mismatch.")
            unique = conn.execute("SELECT COUNT(DISTINCT uniprot_accession) FROM protein_membership").fetchone()[0]
            known = conn.execute("SELECT COUNT(*) FROM direct_participants WHERE stoichiometry_known='1'").fetchone()[0]
            unknown = conn.execute("SELECT COUNT(*) FROM direct_participants WHERE stoichiometry_known='0'").fetchone()[0]
            if unique != manifest.get("unique_protein_count") or known != manifest.get("known_stoichiometry_count") or unknown != manifest.get("unknown_stoichiometry_count"):
                raise ValueError("Complex Portal protein or stoichiometry count mismatch.")
            metadata = dict(conn.execute("SELECT key, value FROM metadata"))
            if metadata != {"database": "Complex Portal", "tax_id": "9606", "scope": "manually_curated"}:
                raise ValueError("Complex Portal SQLite metadata is invalid.")
        for name in TABLE_COLUMNS:
            if not (snapshot / "tables" / f"{name}.csv").is_file():
                raise ValueError(f"Complex Portal normalized table is missing: {name}")
        return True

    def _query(self, sql, params=()):
        database = self.database_path()
        if database is None:
            raise RuntimeError("Complex Portal database is not installed.")
        with sqlite3.connect(database) as conn:
            conn.row_factory = sqlite3.Row
            return [dict(row) for row in conn.execute(sql, params)]

    def get_complex(self, complex_id):
        rows = self._query("SELECT * FROM complexes WHERE complex_id=? OR complex_accession_raw=?", (complex_id, complex_id))
        return rows[0] if rows else None

    def search_complexes(self, query, limit=100):
        pattern = f"%{str(query).strip().casefold()}%"
        return self._query("""SELECT DISTINCT c.* FROM complexes c
            LEFT JOIN complex_aliases a ON a.complex_id=c.complex_id
            LEFT JOIN protein_membership p ON p.complex_id=c.complex_id
            WHERE lower(c.complex_id) LIKE ? OR lower(c.recommended_name) LIKE ?
               OR lower(a.alias) LIKE ? OR lower(p.uniprot_accession) LIKE ?
            ORDER BY c.complex_id LIMIT ?""", (pattern, pattern, pattern, pattern, int(limit)))

    def get_complexes_for_uniprot(self, accession):
        normalization = canonicalize_uniprot_for_complex_portal(accession)
        rows = self._query("""SELECT c.*, p.stoichiometry_raw, p.stoichiometry, p.stoichiometry_known
            FROM protein_membership p JOIN complexes c ON c.complex_id=p.complex_id
            WHERE p.uniprot_accession=? ORDER BY c.complex_id""", (normalization["lookup_accession"],))
        return [{**row, **normalization} for row in rows]

    def get_proteins_for_complex(self, complex_id):
        return self._query("SELECT * FROM protein_membership WHERE complex_id=? ORDER BY uniprot_accession", (complex_id,))

    def get_direct_participants(self, complex_id):
        return self._query("SELECT * FROM direct_participants WHERE complex_id=? ORDER BY CAST(occurrence_index AS INTEGER)", (complex_id,))

    def get_nonprotein_participants(self, complex_id):
        return self._query("SELECT * FROM nonprotein_participants WHERE complex_id=? ORDER BY CAST(occurrence_index AS INTEGER)", (complex_id,))

    def get_nested_complexes(self, complex_id):
        return self._query("SELECT * FROM nested_complexes WHERE parent_complex_id=? ORDER BY CAST(occurrence_index AS INTEGER)", (complex_id,))
