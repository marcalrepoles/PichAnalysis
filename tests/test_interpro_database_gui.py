from __future__ import annotations

import os
import pytest
from PySide6.QtWidgets import QApplication

from pichanalysis.core.database_manager import DatabaseManager
from pichanalysis.core.interpro_database import InterProCancelled
from pichanalysis.ui.database_manager_page import DatabaseManagerPage, InterProAnnotationWorker, InterProMetadataWorker
from test_interpro_database import FakeProvider

os.environ.setdefault("QT_QPA_PLATFORM","offscreen")

@pytest.fixture(scope="module")
def app():yield QApplication.instance() or QApplication([])

def test_database_manager_not_installed_ready_update_error_cancel_and_folder(app,tmp_path,monkeypatch):
    manager=DatabaseManager(tmp_path/"db");manager.interpro.provider=FakeProvider();page=DatabaseManagerPage(manager)
    assert "Not installed" in page.interpro_state_label.text() and page.interpro_primary_button.text()=="Download metadata" and "Cached proteins: 0" in page.interpro_details.text()
    manager.interpro.download_metadata();page.refresh();assert "Ready" in page.interpro_state_label.text() and page.interpro_primary_button.text()=="Update metadata"
    page._interpro_attempt_error="planned";page.refresh();assert "planned" in page.interpro_error_label.text()
    opened=[];monkeypatch.setattr("pichanalysis.ui.database_manager_page.QDesktopServices.openUrl",lambda url:opened.append(url) or True);page._open_interpro_folder();assert opened
    worker=InterProMetadataWorker(manager,FakeProvider());events=[];worker.succeeded.connect(lambda value:events.append("ready"));worker.run();assert events==["ready"]
    canceled=InterProMetadataWorker(manager,FakeProvider());canceled.cancel();canceled.canceled.connect(lambda value:events.append("cancelled"));canceled.run();assert events[-1]=="cancelled";page.close()
    broken=DatabaseManager(tmp_path/"broken");broken.interpro.provider=FakeProvider(fail="unused")
    broken.interpro.provider.metadata=lambda cancel_requested=None:(_ for _ in ()).throw(RuntimeError("planned metadata error"))
    with pytest.raises(RuntimeError):broken.interpro.download_metadata()
    broken_page=DatabaseManagerPage(broken);assert "Error" in broken_page.interpro_state_label.text();broken_page.close()

def test_annotation_worker_success_error_and_cancel(app,tmp_path):
    manager=DatabaseManager(tmp_path/"db");manager.interpro.provider=FakeProvider();manager.interpro.download_metadata();events=[]
    worker=InterProAnnotationWorker(manager.interpro,tmp_path/"project",["P1"],"worker");worker.progress.connect(lambda value:events.append("progress"));worker.succeeded.connect(lambda value:events.append("success"));worker.run();assert "progress" in events and events[-1]=="success"
    error=InterProAnnotationWorker(manager.interpro,tmp_path/"project",["P1"],"worker");error.failed.connect(lambda value:events.append("failed"));error.run();assert events[-1]=="failed"
    cancel=InterProAnnotationWorker(manager.interpro,tmp_path/"project2",["P1"],"cancel");cancel.cancel();cancel.canceled.connect(lambda value:events.append("cancelled"));cancel.run();assert events[-1]=="cancelled"
