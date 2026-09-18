"""Manual official-data Reactome smoke. Not collected by pytest."""
from __future__ import annotations

import hashlib
import io
import json
import os
import shutil
import sys
import tarfile
import tempfile
import urllib.request
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pandas as pd
from PySide6.QtWidgets import QApplication

from pichanalysis.core.database_manager import DatabaseManager
from pichanalysis.core.databases.reactome import CORE_FILES, ReactomeProvider
from pichanalysis.core.organism import set_organism
from pichanalysis.core.project import create_project
from pichanalysis.core.r_runtime import RRuntime
from pichanalysis.core.reactome_analysis import ReactomeParameters, read_reactome_outputs, run_reactome_analysis
from pichanalysis.ui.database_manager_page import DatabaseManagerPage
from pichanalysis.ui.reactome_page import ReactomePage


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(1024 * 1024): digest.update(chunk)
    return digest.hexdigest()


def download_png(ids: list[str], destination: Path) -> tuple[str, str]:
    for identifier in ids[:30]:
        url = f"https://reactome.org/ContentService/exporter/diagram/{identifier}.png"
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "PichAnalysis/0.1", "Accept": "image/png"})
            with urllib.request.urlopen(request, timeout=90) as response:
                payload = response.read()
            if payload.startswith(b"\x89PNG\r\n\x1a\n"):
                destination.write_bytes(payload); return identifier, url
        except Exception:
            continue
    raise RuntimeError("Official Reactome Content Service did not return a PNG for the candidate pathways.")


def main() -> int:
    temporary = Path(tempfile.mkdtemp(prefix="pichanalysis-reactome-real-"))
    report: dict = {"temporary_path": str(temporary)}
    try:
        manager = DatabaseManager(temporary / "databases")
        snapshot = manager.reactome.create_staging_snapshot(); raw = snapshot / "raw"; provider = ReactomeProvider(retries=2, timeout=120)
        downloads = []
        for name in CORE_FILES: downloads.append(provider.download(name, raw / name))
        snapshot = manager.reactome.finalize_snapshot(snapshot)
        pathways = pd.read_csv(raw / "ReactomePathways.txt", sep="\t", header=None, names=["id", "name", "species"])
        diagrams = pd.read_csv(raw / "humanPathwaysWithDiagrams.txt", sep="\t", header=None)
        diagram_ids = sorted(set(diagrams.iloc[:, 0].astype(str)).intersection(pathways.loc[pathways.species.eq("Homo sapiens"), "id"].astype(str)))
        uniprot = pd.read_csv(raw / "UniProt2Reactome.txt", sep="\t", header=None, dtype=str)
        human = uniprot[uniprot.iloc[:, 1].str.startswith("R-HSA-", na=False)].sort_values([0, 1])
        chosen = list(dict.fromkeys(human.iloc[:, 0]))[:12]
        if not chosen: raise RuntimeError("No deterministic human UniProt identifiers were found.")
        project = create_project(temporary / "projects", "OfficialReactomeSmoke"); set_organism(project, "Homo sapiens", "9606")
        catalog = pd.DataFrame({"source_row": range(1, len(chosen)+1), "original_id": chosen, "mapping_status": "mapped_unique", "uniprot_accession": chosen, "ncbi_gene_id": "", "gene_symbol": chosen, "protein_name": chosen})
        catalog.to_csv(project.root / "mapping" / "tables" / "protein_catalog.csv", index=False)
        associated = sorted(set(human.loc[human.iloc[:, 0].isin(chosen), 1]).intersection(diagram_ids))
        png = temporary / "official.png"; diagram_id, png_url = download_png(associated, png)
        archive = temporary / "one-official-diagram.tgz"
        with tarfile.open(archive, "w:gz") as bundle:
            info = tarfile.TarInfo(f"official/{diagram_id}.png"); payload = png.read_bytes(); info.size = len(payload); bundle.addfile(info, io.BytesIO(payload))
        diagram_snapshot = manager.reactome.create_diagram_staging(); diagram_snapshot = manager.reactome.install_diagrams(diagram_snapshot, archive, {"source": png_url, "sha256": sha(archive)})
        before = sha(manager.reactome.diagram_path(diagram_id))

        original_urlopen = urllib.request.urlopen
        urllib.request.urlopen = lambda *_a, **_k: (_ for _ in ()).throw(RuntimeError("network disabled after acquisition"))
        try:
            outputs = run_reactome_analysis(project, manager, RRuntime(), run_id="official-real", parameters=ReactomeParameters(minimum_overlap=1, top_n=10), timeout=300)
            loaded = read_reactome_outputs(project, "official-real")
            app = QApplication.instance() or QApplication([])
            db_page = DatabaseManagerPage(manager); page = ReactomePage(manager); page.set_project(project); page.show_outputs(loaded)
            index = page.pathway.findData(diagram_id)
            if index < 0:
                for candidate in manager.reactome.diagram_index():
                    index = page.pathway.findData(candidate)
                    if index >= 0: diagram_id = candidate; break
            page.pathway.setCurrentIndex(max(index, 0)); page._show_pathway_members(); page._set_diagram_scale(1.25)
            exported = temporary / "exported-official.png"; shutil.copy2(manager.reactome.diagram_path(diagram_id), exported)
            app.processEvents(); page.close(); db_page.close()
        finally:
            urllib.request.urlopen = original_urlopen
        after = sha(manager.reactome.diagram_path(diagram_id))
        pixmap = page._diagram_pixmap
        statuses = outputs.mapping.drop_duplicates("source_row")["mapping_status"].value_counts().to_dict()
        report.update({
            "snapshot_path": str(diagram_snapshot), "snapshot_id": diagram_snapshot.name,
            "release": manager.reactome.manifest().get("release_version", "unknown"),
            "source_files": downloads, "selected_identifiers": chosen,
            "mapped_unique": int(statuses.get("mapped_unique", 0)), "ambiguous": int(statuses.get("ambiguous", 0)), "unmapped": int(statuses.get("unmapped", 0)),
            "memberships": len(outputs.membership), "hierarchy": len(outputs.hierarchy), "frequency": len(outputs.frequency), "enrichment": len(outputs.enrichment),
            "plots": len(outputs.graphs), "workbook": str(outputs.workbook), "diagram_id": diagram_id, "png_source": png_url,
            "png_dimensions": [pixmap.width(), pixmap.height()], "png_sha256_before": before, "png_sha256_after": after,
            "export_sha256": sha(exported), "offline_analysis_and_gui": True,
        })
        assert report["mapped_unique"] > 0 and report["memberships"] > 0 and report["hierarchy"] > 0 and report["frequency"] > 0
        assert outputs.workbook.is_file() and outputs.graphs and before == after == sha(exported)
        print("REACTOME_REAL_SMOKE " + json.dumps(report, indent=2)); return 0
    except Exception as error:
        report["error"] = repr(error); print("REACTOME_REAL_SMOKE_FAILED " + json.dumps(report, indent=2), file=sys.stderr); return 1


if __name__ == "__main__": raise SystemExit(main())
