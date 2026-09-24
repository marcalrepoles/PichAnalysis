"""Offline presentation bundles built only from explicitly selected frozen results."""
from __future__ import annotations

import csv
import hashlib
import html
import json
import os
import shutil
import uuid
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

from PySide6.QtCore import QMarginsF, QUrl
from PySide6.QtGui import QFont, QFontDatabase, QPageLayout, QPageSize, QPdfWriter, QTextDocument


SCHEMA_VERSION = 1
GENERATOR_VERSION = "1.0"
NOT_AVAILABLE = "Not available in this historical run"
MAPPING_LEGACY = "Historical mapping data are not available for safe reconstruction."
PREVIEW_LIMITS = {10, 20, 50, 100}


class ReportError(RuntimeError):
    def __init__(self, code: str, message: str, details: dict | None = None):
        super().__init__(message)
        self.code, self.details = code, details or {}


@dataclass(frozen=True)
class SourceAdapter:
    module_id: str
    label: str
    list_runs: Callable
    load_run: Callable
    root_for: Callable
    artifact_dirs: tuple[str, ...]


@dataclass(frozen=True)
class SourceArtifact:
    artifact_id: str
    kind: str
    label: str
    path: Path


@dataclass(frozen=True)
class SourceRun:
    module: str
    label: str
    run_id: str
    run_date: str
    input_lineage: str
    snapshot: str
    status: str
    source_path: Path
    metadata: dict
    artifacts: tuple[SourceArtifact, ...] = ()
    reason: str = ""


@dataclass(frozen=True)
class ReportSection:
    module: str
    run_id: str
    artifact_ids: tuple[str, ...]
    preview_rows: int = 20
    user_notes: str = ""


@dataclass(frozen=True)
class ReportConfig:
    report_title: str = "PichAnalysis Consolidated Report"
    subtitle: str = ""
    author: str = ""
    description: str = ""
    sections: tuple[ReportSection, ...] = ()
    user_discussion: str = ""
    include_provenance_appendix: bool = True
    include_attachments: bool = True


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _safe_relative(root: Path, path: Path) -> str:
    try:
        relative = path.resolve(strict=True).relative_to(root.resolve(strict=True))
    except (OSError, ValueError) as error:
        raise ReportError("invalid_source_path", "Selected artifact is outside its frozen run.") from error
    return relative.as_posix()


def _metadata(outputs) -> dict:
    value = outputs.get("metadata", {}) if isinstance(outputs, dict) else getattr(outputs, "metadata", {})
    return value if isinstance(value, dict) else {}


def _run_root(project, module_id: str, run_id: str, outputs) -> Path:
    if module_id == "mapping":
        return Path(project.root) / "mapping"
    if isinstance(outputs, dict) and outputs.get("run_root"):
        return Path(outputs["run_root"])
    for name in ("run_root", "root"):
        path = getattr(outputs, name, None)
        if path:
            return Path(path)
    folders = {
        "presence_absence": "analyses/presence_absence/runs", "go": "analyses/GO/runs",
        "kegg": "analyses/KEGG/runs", "reactome": "analyses/Reactome/runs",
        "mitocarta": "analyses/MitoCarta/runs", "domains": "analyses/InterPro_Pfam/runs",
        "string": "analyses/STRING/runs", "complexes": "analyses/Complexes/runs",
        "mtdna_evidence": "analyses/mtDNA_Evidence/runs",
        "differential_preparation": "analyses/Differential/Preparation/runs",
    }
    return Path(project.root) / folders[module_id] / run_id


def default_sources() -> tuple[SourceAdapter, ...]:
    from .cross_module_integration import default_registry
    from . import differential_preparation

    labels = {
        "mapping": "Mapping", "presence_absence": "Presence / Absence",
        "go": "Gene Ontology", "kegg": "KEGG", "reactome": "Reactome",
        "mitocarta": "MitoCarta", "domains": "InterPro / Pfam",
        "string": "STRING", "complexes": "Complexes",
        "mtdna_evidence": "mtDNA Evidence", "proteomics_qc": "Proteomics QC",
        "differential": "Differential Statistics",
    }
    # These are the modules' declared output directories, never arbitrary project files.
    directories = {
        "mapping": ("tables",),
        "presence_absence": ("tables", "graphs"),
        "go": ("annotation", "frequency", "enrichment", "graphs"),
        "kegg": ("mapping", "frequency", "enrichment", "pathways", "graphs"),
        "reactome": ("mapping", "pathways", "frequency", "enrichment", "graphs", "diagrams"),
        "mitocarta": ("mapping", "membership", "subcompartments", "mitopathways", "enrichment", "graphs"),
        "domains": ("mapping", "frequency", "enrichment", "architecture", "plots"),
        "string": ("mapping", "tables", "plots"),
        "complexes": ("mapping", "tables", "plots", "coverage", "enrichment"),
        "mtdna_evidence": ("mapping", "tables", "plots", "evidence", "enrichment"),
        "proteomics_qc": ("tables", "plots"),
        "differential": ("tables", "results", "plots"),
    }
    adapters = []
    for adapter in default_registry().values():
        module_id = adapter.module_id
        adapters.append(SourceAdapter(module_id, labels[module_id], adapter.list_runs,
            adapter.load_run,
            lambda project, run_id, outputs, name=module_id: _run_root(project, name, run_id, outputs),
            directories[module_id]))
    adapters.append(SourceAdapter("differential_preparation", "Differential Preparation",
        differential_preparation.list_runs, differential_preparation.load_run,
        lambda project, run_id, outputs: _run_root(project, "differential_preparation", run_id, outputs),
        ("tables", "plots", "input")))
    return tuple(adapters)


def _discover_artifacts(root: Path, adapter: SourceAdapter) -> tuple[SourceArtifact, ...]:
    paths = []
    if adapter.module_id != "mapping" and (root / "summary.csv").is_file():
        paths.append(root / "summary.csv")
    for relative in adapter.artifact_dirs:
        directory = root / relative
        if not directory.is_dir():
            continue
        for path in sorted(directory.rglob("*")):
            if path.is_file() and path.suffix.lower() in {".csv", ".png"}:
                paths.append(path)
    artifacts = []
    for path in paths:
        artifact_id = _safe_relative(root, path)
        kind = "figure" if path.suffix.lower() == ".png" else "table"
        artifacts.append(SourceArtifact(artifact_id, kind,
            path.stem.replace("_", " ").title(), path))
    return tuple(artifacts)


def _lineage(project, module_id: str, outputs, metadata: dict) -> str:
    if module_id in {"proteomics_qc", "differential", "presence_absence"}:
        from .cross_module_integration import _project_id, _quantitative_lineage
        try:
            value = _quantitative_lineage(outputs, module_id, _project_id(project)).frozen_input_hash
            if value:
                return value
        except (OSError, ValueError, KeyError, TypeError):
            pass
    for key in ("input_lineage", "input_hash", "input_sha256", "matrix_hash"):
        value = metadata.get(key)
        if isinstance(value, str) and value:
            return value
    return NOT_AVAILABLE


def _snapshot(metadata: dict) -> str:
    for key in ("snapshot_id", "annotation_set_id", "database_snapshot_id",
                "reactome_release", "mitocarta_version", "complex_portal_release",
                "string_version", "go_db_version"):
        value = metadata.get(key)
        if value:
            return str(value)
    return NOT_AVAILABLE


def discover_runs(project, sources: tuple[SourceAdapter, ...] | None = None) -> list[SourceRun]:
    sources = sources or default_sources()
    runs = []
    for adapter in sources:
        for run_id in adapter.list_runs(project):
            run_id = str(run_id)
            try:
                outputs = adapter.load_run(project, run_id)
                metadata = _metadata(outputs)
                if str(metadata.get("run_id", run_id)) != run_id:
                    raise ReportError("run_mismatch", "Loaded artifact belongs to another run.")
                root = adapter.root_for(project, run_id, outputs)
                if not root.is_dir():
                    raise FileNotFoundError(root)
                artifacts = _discover_artifacts(root, adapter)
                run_date = str(metadata.get("created_at") or metadata.get("started_at") or
                    metadata.get("finished_at") or NOT_AVAILABLE)
                runs.append(SourceRun(adapter.module_id, adapter.label, run_id, run_date,
                    _lineage(project, adapter.module_id, outputs, metadata),
                    _snapshot(metadata), "Ready", root, metadata, artifacts))
            except (OSError, ValueError, RuntimeError, KeyError, TypeError) as error:
                reason = MAPPING_LEGACY if adapter.module_id == "mapping" else str(error)
                runs.append(SourceRun(adapter.module_id, adapter.label, run_id,
                    NOT_AVAILABLE, NOT_AVAILABLE, NOT_AVAILABLE, "Incomplete",
                    Path(project.root), {}, (), reason))
    return sorted(runs, key=lambda run: (run.label, run.run_date, run.run_id))


def _validate_config(config: ReportConfig) -> None:
    if not config.report_title.strip():
        raise ReportError("invalid_title", "Report title is required.")
    if not config.sections:
        raise ReportError("empty_selection", "Select at least one analysis run.")
    seen = set()
    for section in config.sections:
        key = (section.module, section.run_id)
        if key in seen:
            raise ReportError("duplicate_run", "Each selected run must appear in one section.")
        seen.add(key)
        if section.preview_rows not in PREVIEW_LIMITS:
            raise ReportError("invalid_preview", "Choose a preview limit of 10, 20, 50 or 100 rows.")


def _csv_preview(path: Path, limit: int) -> tuple[str, int, int]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.reader(stream)
        header = next(reader, [])
        rows = []
        count = 0
        for row in reader:
            count += 1
            if len(rows) < limit:
                rows.append(row)
    columns = min(len(header), 8)
    head = "".join(f"<th>{html.escape(value)}</th>" for value in header[:columns])
    body = "".join("<tr>" + "".join(f"<td>{html.escape(row[i] if i < len(row) else '')}</td>"
        for i in range(columns)) + "</tr>" for row in rows)
    table = f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"
    if len(header) > columns:
        table += f"<p class='small'>Showing the first {columns} of {len(header)} persisted columns. The complete source table is in the attachment.</p>"
    return table, count, len(rows)


def _h(value) -> str:
    return html.escape(str(value), quote=True)


def _notes(value: str) -> str:
    return f"<div class='notes'>{_h(value)}</div>" if value.strip() else ""


def _metadata_box(run: SourceRun) -> str:
    metadata = run.metadata
    excluded = {"run_id", "created_at", "started_at", "finished_at"}
    parameters = [(key, value) for key, value in metadata.items()
        if key not in excluded and isinstance(value, (str, int, float, bool))][:20]
    rows = [("Module", run.label), ("Run ID", run.run_id), ("Run date", run.run_date),
            ("Input lineage", run.input_lineage), ("Database snapshot/version", run.snapshot)]
    rows += [(key.replace("_", " ").title(), value) for key, value in parameters]
    return "<dl class='metadata'>" + "".join(f"<dt>{_h(key)}</dt><dd>{_h(value)}</dd>"
        for key, value in rows) + "</dl>"


def _write_csv(path: Path, header: tuple[str, ...], rows: list[tuple]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)


def _render_pdf(html_path: Path, pdf_path: Path) -> None:
    if not QFontDatabase.families():
        for font_path in (Path("C:/Windows/Fonts/arial.ttf"),
                          Path("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"),
                          Path("/Library/Fonts/Arial.ttf")):
            if font_path.is_file() and QFontDatabase.addApplicationFont(str(font_path)) >= 0:
                break
    if not QFontDatabase.families():
        raise ReportError("pdf_font_unavailable", "No local font is available for readable PDF output.")
    document = QTextDocument()
    document.setDefaultFont(QFont("Arial" if "Arial" in QFontDatabase.families() else QFontDatabase.families()[0]))
    document.setBaseUrl(QUrl.fromLocalFile(str(html_path.parent.resolve()) + os.sep))
    document.setHtml(html_path.read_text(encoding="utf-8"))
    writer = QPdfWriter(str(pdf_path))
    writer.setPageSize(QPageSize(QPageSize.PageSizeId.A4))
    writer.setPageMargins(QMarginsF(16, 18, 16, 18), QPageLayout.Unit.Millimeter)
    writer.setResolution(120)
    document.print_(writer)
    del document, writer
    if (not pdf_path.is_file() or pdf_path.stat().st_size < 100 or
            pdf_path.open("rb").read(5) != b"%PDF-"):
        raise ReportError("pdf_failed", "Local PDF generation failed.")


def _html_document(config: ReportConfig, project, report_id: str, generated: str,
                   sections_html: list[str], runs: list[SourceRun], artifacts: list[dict]) -> str:
    modules = list(dict.fromkeys(run.label for run in runs))
    lineages = sorted({run.input_lineage for run in runs if run.input_lineage != NOT_AVAILABLE})
    snapshots = sorted({run.snapshot for run in runs if run.snapshot != NOT_AVAILABLE})
    warning = "<p class='warning'>This report contains results from multiple input lineages.</p>" if len(lineages) > 1 else ""
    toc = "".join(f"<li><a href='#section-{i}'>{_h(run.label)} - {_h(run.run_id)}</a></li>"
        for i, run in enumerate(runs, 1))
    provenance = ""
    if config.include_provenance_appendix:
        provenance = "<section id='provenance'><h2>Provenance</h2>" + "".join(
            f"<h3>{_h(run.label)} - {_h(run.run_id)}</h3>{_metadata_box(run)}"
            for run in runs) + "</section>"
    manifest_rows = "".join(f"<tr><td>{_h(item['module'])}</td><td>{_h(item['run_id'])}</td>"
        f"<td>{_h(item['source_file'])}</td><td>{_h(item['source_hash'])}</td></tr>"
        for item in artifacts)
    discussion = f"<section><h2>User discussion</h2>{_notes(config.user_discussion)}</section>" if config.user_discussion.strip() else ""
    css = """@page{size:A4;margin:18mm 16mm}body{font-family:Arial,sans-serif;color:#182232;font-size:10pt;line-height:1.45}h1{font-size:25pt;color:#193956}h2{font-size:17pt;color:#193956;border-bottom:1px solid #9eb4c5;padding-bottom:5px;margin-top:24px}h3{font-size:12pt;color:#234f6e}a{color:#1b5d8b}section{page-break-inside:avoid;margin-bottom:20px}.cover{min-height:160px;padding-top:45px}.subtitle{font-size:15pt;color:#567}.small{font-size:8pt;color:#556}.warning{padding:9px;background:#fff0d8;border-left:4px solid #bb7a1d}.metadata{display:grid;grid-template-columns:150px 1fr;gap:4px 12px;background:#f3f6f8;padding:12px}.metadata dt{font-weight:bold}.metadata dd{margin:0;overflow-wrap:anywhere}table{border-collapse:collapse;width:100%;font-size:8pt;table-layout:fixed}th,td{border:1px solid #c8d4dd;padding:4px;overflow-wrap:anywhere;vertical-align:top}th{background:#e9f0f5}.artifact{page-break-inside:avoid;margin:15px 0 22px}.artifact img{max-width:100%;max-height:175mm;object-fit:contain}.caption{font-size:8pt;color:#526172}.notes{white-space:pre-wrap;background:#fff9e8;padding:12px;border-left:3px solid #d9a83e}footer{margin-top:30px;color:#667;font-size:8pt}"""
    return f"""<!doctype html><html lang='en'><head><meta charset='utf-8'><title>{_h(config.report_title)}</title><style>{css}</style></head><body>
<section class='cover'><h1>{_h(config.report_title)}</h1><p class='subtitle'>{_h(config.subtitle)}</p><p>{_h(project.name)}</p><p>Generated: {_h(generated)}<br>Report ID: {_h(report_id)}</p><p>{_h(config.author)}</p></section>
<section><h2>Report overview</h2><p>{_h(config.description)}</p><p>Selected analysis runs: {len(runs)}<br>Modules represented: {_h(', '.join(modules))}<br>Input lineages represented: {len(lineages)}<br>Database snapshots represented: {len(snapshots)}</p>{warning}<h3>Contents</h3><ol>{toc}<li><a href='#provenance'>Provenance</a></li></ol></section>
{''.join(sections_html)}{discussion}{provenance}<section><h2>Artifact manifest</h2><p>Complete file paths and SHA-256 hashes are also recorded in report_manifest.json.</p><table><thead><tr><th>Module</th><th>Run</th><th>Source artifact</th><th>SHA-256</th></tr></thead><tbody>{manifest_rows}</tbody></table></section>
<footer>Generated by PichAnalysis Consolidated Report {GENERATOR_VERSION}. This report presents existing results; it does not rerun or reinterpret analyses.</footer></body></html>"""


def generate_report(project, config: ReportConfig,
                    sources: tuple[SourceAdapter, ...] | None = None) -> Path:
    _validate_config(config)
    available = {(run.module, run.run_id): run for run in discover_runs(project, sources)}
    selected = []
    for section in config.sections:
        run = available.get((section.module, section.run_id))
        if run is None or run.status != "Ready":
            raise ReportError("incomplete_run", f"Selected run {section.module}/{section.run_id} is incomplete or unavailable.")
        ids = {artifact.artifact_id for artifact in run.artifacts}
        if not section.artifact_ids or any(item not in ids for item in section.artifact_ids):
            raise ReportError("invalid_artifact", "Selected report artifact is no longer available.")
        selected.append(run)
    report_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ") + "_" + uuid.uuid4().hex[:8]
    base = Path(project.root) / "reports/consolidated/runs"
    base.mkdir(parents=True, exist_ok=True)
    stage = base / (report_id + ".part")
    final = base / report_id
    stage.mkdir(exist_ok=False)
    (stage / "assets").mkdir()
    (stage / "attachments").mkdir()
    generated = datetime.now(timezone.utc).isoformat()
    try:
        (stage / "report_status.json").write_text(json.dumps({"report_id": report_id,
            "status": "Generating"}), encoding="utf-8")
        (stage / "report_config.json").write_text(json.dumps(asdict(config), ensure_ascii=False,
            indent=2), encoding="utf-8")
        source_rows, artifact_rows, artifact_manifest, sections_html = [], [], [], []
        for number, (section, run) in enumerate(zip(config.sections, selected), 1):
            source_rows.append((run.module, run.run_id, run.run_date, run.input_lineage,
                run.snapshot, run.status, str(run.source_path)))
            by_id = {artifact.artifact_id: artifact for artifact in run.artifacts}
            blocks = []
            for artifact_id in section.artifact_ids:
                artifact = by_id[artifact_id]
                source = artifact.path
                if not source.is_file() or _safe_relative(run.source_path, source) != artifact_id:
                    raise ReportError("missing_artifact", "Selected report artifact is no longer available.")
                source_hash = _sha256(source)
                attachment = None
                if config.include_attachments:
                    attachment = stage / "attachments" / f"section-{number}" / artifact_id
                    attachment.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copyfile(source, attachment)
                    if _sha256(attachment) != source_hash:
                        raise ReportError("hash_mismatch", "Copied artifact SHA-256 does not match source.")
                bundle_path = attachment
                if artifact.kind == "figure":
                    asset = stage / "assets" / f"section-{number}" / artifact_id
                    asset.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copyfile(source, asset)
                    if _sha256(asset) != source_hash:
                        raise ReportError("hash_mismatch", "Copied figure SHA-256 does not match source.")
                    bundle_path = asset
                    blocks.append(f"<div class='artifact'><h3>{_h(artifact.label)}</h3><img src='{_h(asset.relative_to(stage).as_posix())}' alt='{_h(artifact.label)}'><p class='caption'>Source: {_h(run.label)} · Run {_h(run.run_id)} · {_h(artifact_id)} · Snapshot: {_h(run.snapshot)}</p></div>")
                    companion = source.with_suffix(".pdf")
                    if config.include_attachments and companion.is_file():
                        _safe_relative(run.source_path, companion)
                        pdf_copy = attachment.with_suffix(".pdf")
                        shutil.copyfile(companion, pdf_copy)
                        if _sha256(pdf_copy) != _sha256(companion):
                            raise ReportError("hash_mismatch", "Plot PDF attachment hash mismatch.")
                        companion_item = {"section": number, "module": run.module, "run_id": run.run_id,
                            "artifact_type": "figure_pdf", "source_file": str(companion),
                            "source_hash": _sha256(companion),
                            "bundle_file": pdf_copy.relative_to(stage).as_posix(),
                            "bundle_hash": _sha256(pdf_copy), "display_mode": "attachment",
                            "snapshot": run.snapshot, "input_lineage": run.input_lineage}
                        artifact_manifest.append(companion_item)
                        artifact_rows.append(tuple(companion_item[key] for key in ("section", "module", "run_id",
                            "artifact_type", "source_file", "source_hash", "bundle_file", "bundle_hash", "display_mode")))
                else:
                    table_html, total, shown = _csv_preview(source, section.preview_rows)
                    description = (f"Showing the first {shown} rows in the persisted source order. "
                        + ("The complete source table is included in the report attachments." if attachment else
                           "The complete source table was not selected as an attachment."))
                    blocks.append(f"<div class='artifact'><h3>{_h(artifact.label)}</h3><p class='small'>{_h(description)} Total persisted rows: {total}.</p>{table_html}<p class='caption'>Source: {_h(run.label)} · Run {_h(run.run_id)} · {_h(artifact_id)} · Snapshot: {_h(run.snapshot)}</p></div>")
                copied = bundle_path.relative_to(stage).as_posix() if bundle_path else ""
                copied_hash = _sha256(bundle_path) if bundle_path else ""
                item = {"section": number, "module": run.module, "run_id": run.run_id,
                    "artifact_type": artifact.kind, "source_file": str(source),
                    "source_hash": source_hash, "bundle_file": copied,
                    "bundle_hash": copied_hash, "display_mode": "image" if artifact.kind == "figure" else "preview",
                    "snapshot": run.snapshot, "input_lineage": run.input_lineage}
                artifact_manifest.append(item)
                artifact_rows.append(tuple(item[key] for key in ("section", "module", "run_id",
                    "artifact_type", "source_file", "source_hash", "bundle_file", "bundle_hash", "display_mode")))
            sections_html.append(f"<section id='section-{number}'><h2>{_h(run.label)} - {_h(run.run_id)}</h2>{_metadata_box(run)}{''.join(blocks)}"
                + (f"<h3>User notes</h3>{_notes(section.user_notes)}" if section.user_notes.strip() else "") + "</section>")
        _write_csv(stage / "sources.csv", ("module", "run_id", "run_date", "input_lineage",
            "snapshot", "status", "source_path"), source_rows)
        _write_csv(stage / "artifacts.csv", ("section", "module", "run_id", "artifact_type",
            "source_file", "source_hash", "bundle_file", "bundle_hash", "display_mode"), artifact_rows)
        html_path = stage / "report.html"
        html_path.write_text(_html_document(config, project, report_id, generated,
            sections_html, selected, artifact_manifest), encoding="utf-8")
        _render_pdf(html_path, stage / "report.pdf")
        manifest = {"schema_version": SCHEMA_VERSION, "generator_version": GENERATOR_VERSION,
            "report_id": report_id, "status": "Ready", "generation_timestamp": generated,
            "project": str(project.root), "report_configuration": asdict(config),
            "selected_runs": [dict(module=run.module, run_id=run.run_id,
                run_date=run.run_date, input_lineage=run.input_lineage,
                database_snapshot=run.snapshot, metadata=run.metadata) for run in selected],
            "selected_artifacts": artifact_manifest,
            "bundle_hashes": {name: _sha256(stage / name) for name in
                ("report.html", "report.pdf", "report_config.json", "sources.csv", "artifacts.csv")}}
        (stage / "report_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False,
            indent=2, default=str), encoding="utf-8")
        (stage / "report_status.json").write_text(json.dumps({"report_id": report_id,
            "status": "Ready"}), encoding="utf-8")
        os.replace(stage, final)
        return final
    except Exception as error:
        (stage / "report_status.json").write_text(json.dumps({"report_id": report_id,
            "status": "Incomplete", "error": str(error)}), encoding="utf-8")
        raise


def list_reports(project) -> list[dict]:
    base = Path(project.root) / "reports/consolidated/runs"
    if not base.is_dir():
        return []
    reports = []
    for folder in sorted((path for path in base.iterdir() if path.is_dir()), reverse=True):
        manifest_path = folder / "report_manifest.json"
        status_path = folder / "report_status.json"
        try:
            status = json.loads(status_path.read_text(encoding="utf-8")) if status_path.is_file() else {}
            manifest = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.is_file() else {}
        except (OSError, ValueError):
            continue
        config = manifest.get("report_configuration", {})
        runs = manifest.get("selected_runs", [])
        reports.append({"report_id": status.get("report_id", folder.name),
            "status": status.get("status", "Incomplete"),
            "date": manifest.get("generation_timestamp", ""),
            "title": config.get("report_title", ""),
            "modules": ", ".join(dict.fromkeys(run.get("module", "") for run in runs)),
            "runs": ", ".join(run.get("run_id", "") for run in runs),
            "path": folder})
    return reports


def load_report(project, report_id: str) -> dict:
    base = (Path(project.root) / "reports/consolidated/runs").resolve()
    folder = (base / report_id).resolve()
    if folder.parent != base:
        raise ReportError("invalid_report_id", "Invalid report ID.")
    manifest = json.loads((folder / "report_manifest.json").read_text(encoding="utf-8"))
    status = json.loads((folder / "report_status.json").read_text(encoding="utf-8"))
    if (manifest.get("status") != "Ready" or manifest.get("report_id") != report_id
            or status.get("status") != "Ready") :
        raise ReportError("incomplete_report", "Historical report is not Ready.")
    for relative, digest in manifest["bundle_hashes"].items():
        if _sha256(folder / relative) != digest:
            raise ReportError("bundle_corrupt", f"Historical report artifact failed SHA-256: {relative}")
    for item in manifest["selected_artifacts"]:
        if item["bundle_file"] and _sha256(folder / item["bundle_file"]) != item["bundle_hash"]:
            raise ReportError("bundle_corrupt", "Historical attachment failed SHA-256.")
    return {"path": folder, "manifest": manifest,
        "html": folder / "report.html", "pdf": folder / "report.pdf"}


def export_report(project, report_id: str, destination: Path, kind: str) -> Path:
    report = load_report(project, report_id)
    destination = Path(destination)
    if destination.exists():
        raise ReportError("destination_exists", "Destination already exists; report export never overwrites.")
    if kind == "bundle":
        shutil.copytree(report["path"], destination)
    elif kind == "html":
        asset_dir = destination.parent / (destination.stem + "_assets")
        if asset_dir.exists():
            raise ReportError("destination_exists", "Destination assets already exist; export never overwrites.")
        shutil.copytree(report["path"] / "assets", asset_dir)
        html_text = report["html"].read_text(encoding="utf-8")
        html_text = html_text.replace("src='assets/", f"src='{asset_dir.name}/")
        destination.write_text(html_text, encoding="utf-8")
    elif kind == "pdf":
        shutil.copy2(report["pdf"], destination)
    else:
        raise ReportError("invalid_export", "Choose HTML, PDF or bundle export.")
    return destination
