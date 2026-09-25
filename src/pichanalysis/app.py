from __future__ import annotations

import sys
import time

from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QIcon, QPixmap
from PySide6.QtWidgets import QApplication, QSplashScreen

from .core.resources import branding_asset
from .ui.main_window import MainWindow


def show_branded_window(application: QApplication, minimum_splash_ms: int = 3000) -> tuple[MainWindow, QSplashScreen]:
    """Keep the splash visible for at least the minimum or until startup ends."""
    icon = QIcon(str(branding_asset("icon.ico")))
    if icon.isNull():
        raise RuntimeError("PichAnalysis icon could not be loaded.")
    application.setWindowIcon(icon)
    pixmap = QPixmap(str(branding_asset("splash.png")))
    if pixmap.isNull():
        raise RuntimeError("PichAnalysis splash image could not be loaded.")
    pixmap = pixmap.scaled(1000, 600, Qt.KeepAspectRatio, Qt.SmoothTransformation)
    splash = QSplashScreen(pixmap)
    splash.setWindowIcon(icon)
    splash.show()
    application.processEvents()
    shown_at = time.monotonic()

    window = MainWindow()
    window.setWindowIcon(icon)

    def reveal() -> None:
        window.show()
        splash.finish(window)

    remaining_ms = max(0, minimum_splash_ms - int((time.monotonic() - shown_at) * 1000))
    QTimer.singleShot(remaining_ms, reveal)
    return window, splash

def main() -> int:
    if len(sys.argv) >= 3 and sys.argv[1] == "--release-smoke":
        from pathlib import Path
        from .release_smoke import run_packaged_smoke

        return run_packaged_smoke(Path(sys.argv[2]),
            simulate_missing_r="--simulate-missing-r" in sys.argv[3:])
    if sys.platform == "win32":
        import ctypes
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("PichAnalysis.Desktop")
    application = QApplication.instance() or QApplication(sys.argv)
    application.setApplicationName("PichAnalysis")
    window, splash = show_branded_window(application)
    return application.exec()

