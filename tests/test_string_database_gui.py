import os
os.environ.setdefault("QT_QPA_PLATFORM","offscreen")
import pytest
from PySide6.QtWidgets import QApplication
from pichanalysis.core.database_manager import DatabaseManager
from pichanalysis.ui.database_manager_page import DatabaseManagerPage
from pichanalysis.ui.string_database_card import StringDownloadWorker
from test_string_database import fixture

@pytest.fixture(scope="module")
def app():return QApplication.instance() or QApplication([])

class FixtureProvider:
    def __init__(self,source,mode="ok"):self.source=source;self.mode=mode
    def detect_version(self):return {"string_version":"fixture","string_stable_address":"fixture://official"}
    def filenames(self,version):return tuple(x.name for x in sorted(self.source.glob("*.txt.gz")))
    def download(self,address,name,target,progress=None,cancel_requested=None):
        if self.mode=="error":raise RuntimeError("fixture failure")
        if self.mode=="cancel" or (cancel_requested and cancel_requested()):
            from pichanalysis.core.databases.string import StringDownloadCancelled
            raise StringDownloadCancelled("STRING download was canceled.")
        import shutil,hashlib
        shutil.copyfile(self.source/name,target);data=target.read_bytes()
        if progress:progress(name,len(data),len(data))
        return {"filename":name,"source":str(self.source/name),"size":len(data),"sha256":hashlib.sha256(data).hexdigest()}

def test_database_manager_card_states(tmp_path,app):
    source=tmp_path/"source";fixture(source);manager=DatabaseManager(tmp_path/"db");page=DatabaseManagerPage(manager)
    card=page.string_card;assert "Not installed" in card.state_label.text() and card.primary_button.text()=="Download"
    worker=StringDownloadWorker(manager,FixtureProvider(source));success=[];worker.succeeded.connect(success.append);worker.run();page.refresh()
    assert success and "Ready" in card.state_label.text() and card.primary_button.text()=="Update" and "Functional network" in card.details.text()
    worker=StringDownloadWorker(manager,FixtureProvider(source,"error"));errors=[];worker.failed.connect(errors.append);worker.run();assert errors and manager.string.active_snapshot()
    worker=StringDownloadWorker(manager,FixtureProvider(source,"cancel"));canceled=[];worker.canceled.connect(canceled.append);worker.run();assert canceled and manager.string.active_snapshot()
    page.close()
