from __future__ import annotations

import hashlib
import os
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pandas as pd
import pytest
from PySide6.QtWidgets import QApplication, QGroupBox

from pichanalysis.core.database_manager import DatabaseManager
from pichanalysis.core.database_registry import DATABASE_REGISTRY, DatabaseState
from pichanalysis.core.databases.mitocarta import FILES, MitoCartaDownloadCancelled, MitoCartaProvider
from pichanalysis.core.mitocarta_database import parse_gmx, parse_workbook
from pichanalysis.ui.database_manager_page import DatabaseManagerPage, MitoCartaDownloadWorker


def write_fixture(root: Path, invalid=False) -> Path:
    root.mkdir()
    workbook = root / FILES[0]
    temporary = root / "fixture.xlsx"
    if invalid:
        with pd.ExcelWriter(temporary, engine="openpyxl") as writer: pd.DataFrame({"wrong": [1]}).to_excel(writer, sheet_name="Wrong", index=False)
    else:
        columns = {
            "HumanGeneID": [1, 2, 3], "Symbol": ["GENE1", "GENE2", "GENE3"], "Description": ["one", "two", "three"],
            "MitoCarta3.0_List": ["MitoCarta3.0", "MitoCarta3.0", "0"], "UniProt": ["P1", "P2", "P3"],
            "Synonyms": ["A|B", None, "C"], "MitoCarta2.0_Score": [10.0, 9.0, 1.0], "MitoCarta3.0_Evidence": ["literature", "MS/MS", None],
            "MitoCarta3.0_SubMitoLocalization": ["matrix", "MIM|IMS", None],
            "MitoCarta3.0_MitoPathways": ["Metabolism > TCA cycle", "OXPHOS > Complex I; Signaling", None], "Tissues": ["all 14", "heart", None],
        }
        with pd.ExcelWriter(temporary, engine="openpyxl") as writer:
            pd.DataFrame({"notes": ["fixture"]}).to_excel(writer, sheet_name="Description", index=False)
            pd.DataFrame({key: value[:2] for key, value in columns.items()}).to_excel(writer, sheet_name="A Human MitoCarta3.0", index=False)
            pd.DataFrame(columns).to_excel(writer, sheet_name="B Human All Genes", index=False)
    temporary.replace(workbook)
    (root / FILES[1]).write_text("Path_A\tPath_B\nRoot > Child\tOther\nGENE1\tGENE2\nGENE1\tGENE2\nGENE2\t\n", encoding="utf-8")
    return root


def test_registry_and_workbook_detection_membership(tmp_path):
    assert DATABASE_REGISTRY["mitocarta"].taxonomy_id == "9606"
    genes, compartments, metadata = parse_workbook(write_fixture(tmp_path / "source") / FILES[0])
    assert metadata["workbook_sheet_used"] == "B Human All Genes" and metadata["total_workbook_rows"] == 3
    assert genes.is_mitocarta.tolist() == [True, True, False] and metadata["mitocarta_gene_count"] == 2
    assert set(map(tuple, compartments.values)) == {("GENE1", "matrix"), ("GENE2", "MIM"), ("GENE2", "IMS")}


def test_invalid_workbook_headers_fail_clearly(tmp_path):
    with pytest.raises(ValueError, match="required official headers"): parse_workbook(write_fixture(tmp_path / "bad", True) / FILES[0])


def test_gmx_membership_dedup_and_hierarchy(tmp_path):
    path=write_fixture(tmp_path/"source")/FILES[1];members,hierarchy,metadata=parse_gmx(path)
    assert len(members)==3 and metadata["mitopathway_count"]==2
    row=hierarchy[hierarchy.mitopathway.eq("Path_A")].iloc[0];assert row.pathway_full=="Root > Child" and row.pathway_level_1=="Root" and row.pathway_level_2=="Child"


def test_snapshot_manifest_api_and_atomic_updates(tmp_path):
    manager=DatabaseManager(tmp_path/"db");first=manager.mitocarta.install_from_directory(write_fixture(tmp_path/"a"));before=(first/"tables/mitocarta_genes.csv").read_bytes()
    assert manager.mitocarta.is_ready() and manager.mitocarta.version()=="3.0" and manager.mitocarta.gene_count()==2 and manager.mitocarta.pathway_count()==2
    assert all(path.is_file() for path in (manager.mitocarta.genes_path(),manager.mitocarta.pathway_membership_path(),manager.mitocarta.subcompartment_membership_path()))
    assert len(manager.mitocarta.source_hashes())==2 and manager.mitocarta.manifest()["workbook_sheet_used"]=="B Human All Genes"
    second=manager.mitocarta.install_from_directory(write_fixture(tmp_path/"b"));assert second!=first and manager.mitocarta.active_snapshot()==second and first.is_dir()
    with pytest.raises(ValueError):manager.mitocarta.install_from_directory(write_fixture(tmp_path/"bad",True))
    assert manager.mitocarta.active_snapshot()==second and (first/"tables/mitocarta_genes.csv").read_bytes()==before


def test_provider_download_hash_part_retry_and_cancel(tmp_path,monkeypatch):
    payload=b"fixture"
    class Response:
        headers={"Content-Length":str(len(payload))}
        def __enter__(self):return self
        def __exit__(self,*_):return False
        def read(self,_):value=getattr(self,"value",payload);self.value=b"";return value
    calls=[]
    def urlopen(*a,**k):
        calls.append(1)
        if len(calls)==1:raise __import__("urllib.error").error.URLError("planned")
        return Response()
    monkeypatch.setattr("urllib.request.urlopen",urlopen);monkeypatch.setattr("time.sleep",lambda *_:None);destination=tmp_path/FILES[1];seen=[]
    info=MitoCartaProvider(retries=1).download(FILES[1],destination,progress=lambda *_:seen.append(destination.with_name(destination.name+".part").exists()))
    assert len(calls)==2 and all(seen) and info["sha256"]==hashlib.sha256(payload).hexdigest() and not destination.with_name(destination.name+".part").exists()
    checks=iter((False,True));monkeypatch.setattr("urllib.request.urlopen",lambda *_a,**_k:Response())
    with pytest.raises(MitoCartaDownloadCancelled):MitoCartaProvider(retries=0).download(FILES[1],tmp_path/"cancel.gmx",cancel_requested=lambda:next(checks))


class FixtureProvider:
    def __init__(self,source,fail=False,cancel=None):self.source=source;self.fail=fail;self.cancel=cancel
    def download(self,name,destination,progress=None,cancel_requested=None):
        if self.fail and name==FILES[1]:raise RuntimeError("planned failure")
        data=(self.source/name).read_bytes();progress(name,len(data),len(data));destination.write_bytes(data)
        if self.cancel:self.cancel()
        if cancel_requested and cancel_requested():raise MitoCartaDownloadCancelled("MitoCarta3.0 download was canceled.")
        return {"filename":name,"source_url":"fixture","size":len(data),"sha256":hashlib.sha256(data).hexdigest()}


def test_worker_success_failure_cancel_and_previous_active(tmp_path):
    source=write_fixture(tmp_path/"source");manager=DatabaseManager(tmp_path/"db");worker=MitoCartaDownloadWorker(manager,FixtureProvider(source));success=[];worker.succeeded.connect(success.append);worker.run();active=manager.mitocarta.active_snapshot();assert success and active
    failed=[];worker=MitoCartaDownloadWorker(manager,FixtureProvider(source,fail=True));worker.failed.connect(failed.append);worker.run();assert failed and manager.mitocarta.active_snapshot()==active
    canceled=[];provider=FixtureProvider(source);worker=MitoCartaDownloadWorker(manager,provider);provider.cancel=worker.cancel;worker.canceled.connect(canceled.append);worker.run();assert canceled and manager.mitocarta.active_snapshot()==active


def test_database_manager_card_states_and_gui_smoke(tmp_path,monkeypatch):
    app=QApplication.instance() or QApplication([]);manager=DatabaseManager(tmp_path/"db");page=DatabaseManagerPage(manager);page.show();app.processEvents()
    assert "MitoCarta3.0 — Homo sapiens" in {group.title() for group in page.findChildren(QGroupBox)} and page.mitocarta_state_label.text()=="Status: Not installed" and page.mitocarta_primary_button.text()=="Download"
    manager.mitocarta.install_from_directory(write_fixture(tmp_path/"source"));page.refresh();assert page.mitocarta_state_label.text()=="Status: Ready" and page.mitocarta_primary_button.text()=="Update"
    opened=[];monkeypatch.setattr("pichanalysis.ui.database_manager_page.QDesktopServices.openUrl",lambda url:opened.append(url) or True);page._open_mitocarta_folder();assert opened;page.close()
    import pichanalysis.ui.database_manager_page as database_page_module
    from pichanalysis.ui.main_window import MainWindow
    monkeypatch.setattr(database_page_module,"DatabaseManager",lambda:manager);window=MainWindow();window.show();window.navigation.setCurrentRow(3);app.processEvents();assert window.database_manager_page.mitocarta_state_label.text()=="Status: Ready";window.close();app.processEvents()
