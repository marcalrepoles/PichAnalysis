"""Main-window Complexes wiring, safeguard choices, and close guard."""
import os

os.environ.setdefault("QT_QPA_PLATFORM","offscreen")

from PySide6.QtGui import QCloseEvent
from PySide6.QtWidgets import QApplication, QMessageBox

from pichanalysis.ui.main_window import MainWindow


def test_complexes_tab_and_safeguard_choices(monkeypatch):
    app=QApplication.instance() or QApplication([])
    window=MainWindow()
    assert window.analyses_page.complex_page is not None
    captured=[]
    monkeypatch.setattr(window,"_run_complex",lambda params, allowed=False:captured.append((params,allowed)))
    details={"entities_outside_background":["P12345"],"initial_target_size":2,"initial_background_size":1}
    original_exec=QMessageBox.exec
    def choose_continue(box):
        buttons=box.buttons()
        assert {b.text() for b in buttons}=={"Continue","Cancel"}
        assert "target size: 2" in box.text() and "background size: 1" in box.text()
        box._selected=next(b for b in buttons if b.text()=="Continue")
        return 0
    monkeypatch.setattr(QMessageBox,"exec",choose_continue)
    monkeypatch.setattr(QMessageBox,"clickedButton",lambda box:getattr(box,"_selected",None))
    window._complex_target_outside(details,{"minimum_overlap":3})
    assert captured==[({"minimum_overlap":3},True)]
    captured.clear()
    def choose_cancel(box):
        box._selected=next(b for b in box.buttons() if b.text()=="Cancel")
        return 0
    monkeypatch.setattr(QMessageBox,"exec",choose_cancel)
    window._complex_target_outside(details,{"minimum_overlap":3})
    assert captured==[]
    class Running:
        def isRunning(self): return True
    window.complex_worker=Running()
    notices=[]
    monkeypatch.setattr(QMessageBox,"information",lambda *args:notices.append(args))
    event=QCloseEvent()
    window.closeEvent(event)
    assert not event.isAccepted() and notices
    window.complex_worker=None
    monkeypatch.setattr(QMessageBox,"exec",original_exec)
    window.close(); app.processEvents()
