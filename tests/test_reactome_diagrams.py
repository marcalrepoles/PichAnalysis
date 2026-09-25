from __future__ import annotations

import base64
import hashlib
import io
import tarfile
import urllib.error
from pathlib import Path

import pytest
from PySide6.QtWidgets import QApplication, QFileDialog

from pichanalysis.core.database_manager import DatabaseManager
from pichanalysis.core.database_registry import DatabaseState
from pichanalysis.core.databases.reactome import DIAGRAM_ARCHIVE, ReactomeDownloadCancelled, ReactomeProvider, safe_extract_tar
from pichanalysis.ui.database_manager_page import ReactomeDiagramDownloadWorker
from pichanalysis.ui.reactome_page import ReactomePage

PNG = base64.b64decode("iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mNk+A8AAQUBAScY42YAAAAASUVORK5CYII=")


def write_core(root: Path) -> Path:
    root.mkdir(parents=True)
    values = {
        "ReactomePathways.txt": "R-HSA-1\tPathway one\tHomo sapiens\n",
        "ReactomePathwaysRelation.txt": "R-HSA-1\tR-HSA-2\n",
        "UniProt2Reactome.txt": "P1\tR-HSA-1\turl\tPathway one\tIEA\tHomo sapiens\n",
        "NCBI2Reactome.txt": "1\tR-HSA-1\turl\tPathway one\tIEA\tHomo sapiens\n",
        "humanPathwaysWithDiagrams.txt": "R-HSA-1\tPathway one\n",
        "pathway2summation.txt": "R-HSA-1\tSummary\n",
    }
    for name, value in values.items(): (root / name).write_text(value, encoding="utf-8")
    return root


def test_staging_snapshot_ids_remain_unique_when_clock_collides(tmp_path, monkeypatch):
    from datetime import datetime, timezone
    from pichanalysis.core import reactome_database as module
    from pichanalysis.core.reactome_database import ReactomeDatabase

    class FrozenClock:
        @staticmethod
        def now(_timezone):
            return datetime(2026, 1, 1, tzinfo=timezone.utc)

    monkeypatch.setattr(module, "datetime", FrozenClock)
    database = ReactomeDatabase(tmp_path)
    first = database.create_staging_snapshot()
    second = database.create_staging_snapshot()
    assert first != second
    assert first.is_dir() and second.is_dir()
    assert database.manifest(second)["snapshot_id"] == second.name


def write_tar(path: Path, members: list[tuple[str, bytes]], *, kind=None, link="") -> Path:
    with tarfile.open(path, "w:gz") as bundle:
        for name, content in members:
            info = tarfile.TarInfo(name); info.size = len(content)
            if kind is not None: info.type = kind; info.linkname = link; info.size = 0
            bundle.addfile(info, None if kind is not None else io.BytesIO(content))
    return path


@pytest.mark.parametrize("name", ["../evil", "../../evil", "/absolute", "C:\\evil", "nested/../../../evil"])
def test_safe_tar_rejects_traversal_absolute_and_drive_paths(tmp_path, name):
    archive = write_tar(tmp_path / "unsafe.tgz", [(name, PNG)])
    with pytest.raises(ValueError, match="Unsafe"): safe_extract_tar(archive, tmp_path / "out")
    assert not (tmp_path / "escape.png").exists()


def test_diagram_archive_download_success_hash_and_atomic_part(tmp_path, monkeypatch):
    payload = b"official archive fixture"
    class Response:
        headers = {"Content-Length": str(len(payload))}
        def __enter__(self): return self
        def __exit__(self, *_args): return False
        def read(self, _size): value=getattr(self,"value",payload);self.value=b"";return value
    destination=tmp_path/DIAGRAM_ARCHIVE; observed=[]
    monkeypatch.setattr("urllib.request.urlopen",lambda *_a,**_k:Response())
    info=ReactomeProvider(retries=0).download(DIAGRAM_ARCHIVE,destination,progress=lambda *_:observed.append(destination.with_name(destination.name+".part").is_file()))
    assert destination.read_bytes()==payload and all(observed) and not destination.with_name(destination.name+".part").exists()
    assert info["sha256"]==hashlib.sha256(payload).hexdigest() and info["file_size"]==len(payload)


def test_diagram_archive_download_retries_bounded(tmp_path, monkeypatch):
    calls=[]
    class Response:
        headers={"Content-Length":"2"}
        def __enter__(self):return self
        def __exit__(self,*_args):return False
        def read(self,_size):value=getattr(self,"value",b"ok");self.value=b"";return value
    def open_url(*_a,**_k):
        calls.append(1)
        if len(calls)==1:raise urllib.error.URLError("planned")
        return Response()
    monkeypatch.setattr("urllib.request.urlopen",open_url);monkeypatch.setattr("time.sleep",lambda *_:None)
    destination=tmp_path/DIAGRAM_ARCHIVE;ReactomeProvider(retries=1).download(DIAGRAM_ARCHIVE,destination)
    assert len(calls)==2 and destination.read_bytes()==b"ok"


@pytest.mark.parametrize("kind", [tarfile.SYMTYPE, tarfile.LNKTYPE])
def test_safe_tar_rejects_symbolic_and_hard_links(tmp_path, kind):
    archive = write_tar(tmp_path / "link.tgz", [("R-HSA-1.png", b"")], kind=kind, link="outside")
    with pytest.raises(ValueError, match="member type"): safe_extract_tar(archive, tmp_path / "out")


def test_diagram_snapshot_index_manifest_lookup_and_core_preservation(tmp_path):
    manager = DatabaseManager(tmp_path / "db"); original = manager.reactome.install_from_directory(write_core(tmp_path / "core"))
    archive = write_tar(tmp_path / DIAGRAM_ARCHIVE, [("z/R-MMU-9.png", PNG), ("b/R-HSA-2.png", PNG), ("a/R-HSA-1.png", PNG)])
    staging = manager.reactome.create_diagram_staging(); completed = manager.reactome.install_diagrams(staging, archive)
    assert completed != original and manager.reactome.active_snapshot() == completed
    assert manager.reactome.diagram_index() == {"R-HSA-1": "a/R-HSA-1.png", "R-HSA-2": "b/R-HSA-2.png"}
    assert manager.reactome.diagram_path("R-HSA-1").read_bytes() == PNG
    assert manager.reactome.is_diagrams_ready() and manager.reactome.source_hashes()
    info = manager.reactome.diagram_manifest()
    assert info["license"] == "CC BY 4.0" and info["status"] == DatabaseState.READY and info["diagram_count"] == 2


def test_failed_or_canceled_diagram_update_keeps_previous_snapshot(tmp_path):
    manager = DatabaseManager(tmp_path / "db"); active = manager.reactome.install_from_directory(write_core(tmp_path / "core"))
    bad = write_tar(tmp_path / "bad.tgz", [("R-MMU-1.png", PNG)]); staging = manager.reactome.create_diagram_staging()
    with pytest.raises(ValueError): manager.reactome.install_diagrams(staging, bad)
    assert manager.reactome.active_snapshot() == active and manager.reactome.is_core_ready()
    canceled = manager.reactome.create_diagram_staging(); good = write_tar(tmp_path / "good.tgz", [("R-HSA-1.png", PNG)])
    with pytest.raises(ReactomeDownloadCancelled): manager.reactome.install_diagrams(canceled, good, cancel_requested=lambda: True)
    assert manager.reactome.active_snapshot() == active


def test_diagram_worker_downloads_stages_and_activates(tmp_path):
    manager = DatabaseManager(tmp_path / "db"); manager.reactome.install_from_directory(write_core(tmp_path / "core"))
    class Provider:
        def download(self, name, destination, progress=None, cancel_requested=None):
            write_tar(destination, [("R-HSA-1.png", PNG)]); progress(name, destination.stat().st_size, destination.stat().st_size)
            return {"source": "official-test", "sha256": "test"}
    worker = ReactomeDiagramDownloadWorker(manager, Provider()); events=[]; worker.progress.connect(events.append); worker.run()
    assert manager.reactome.is_diagrams_ready() and {event["stage"] for event in events} == {"download", "extract", "index"}


def test_local_viewer_zoom_limits_fit_export_and_original_immutable(tmp_path, monkeypatch):
    app = QApplication.instance() or QApplication([])
    manager = DatabaseManager(tmp_path / "db"); manager.reactome.install_from_directory(write_core(tmp_path / "core"))
    archive = write_tar(tmp_path / DIAGRAM_ARCHIVE, [("R-HSA-1.png", PNG)]); staging = manager.reactome.create_diagram_staging(); manager.reactome.install_diagrams(staging, archive)
    source = manager.reactome.diagram_path("R-HSA-1"); before = __import__("hashlib").sha256(source.read_bytes()).hexdigest()
    page = ReactomePage(manager); page.pathway.addItem("R-HSA-1", "R-HSA-1"); page._show_diagram("R-HSA-1")
    assert not page._diagram_pixmap.isNull() and page.diagram_export.isEnabled()
    page._set_diagram_scale(1.0); assert page._diagram_scale == 1.0
    page._set_diagram_scale(99); assert page._diagram_scale == 8.0
    page._set_diagram_scale(0.001); assert page._diagram_scale == 0.1
    page._fit_diagram(); assert 0.1 <= page._diagram_scale <= 8.0
    destination = tmp_path / "export.png"; monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: (str(destination), "PNG (*.png)")); page._export_diagram()
    assert destination.read_bytes() == source.read_bytes()
    assert __import__("hashlib").sha256(source.read_bytes()).hexdigest() == before
    page.close(); app.processEvents()


def test_page_handles_unavailable_diagram_without_blocking_analysis(tmp_path):
    app = QApplication.instance() or QApplication([]); manager = DatabaseManager(tmp_path / "db"); manager.reactome.install_from_directory(write_core(tmp_path / "core"))
    page = ReactomePage(manager); page._show_diagram("R-HSA-404")
    assert "No local diagram" in page.diagram_image.text() and not page.diagram_export.isEnabled()
    assert manager.reactome.is_core_ready() and not manager.reactome.is_diagrams_ready(); page.close(); app.processEvents()
