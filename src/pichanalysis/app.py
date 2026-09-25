from __future__ import annotations

import sys

from PySide6.QtWidgets import QApplication

from .ui.main_window import MainWindow


def main() -> int:
    if len(sys.argv) >= 3 and sys.argv[1] == "--release-smoke":
        from pathlib import Path
        from .release_smoke import run_packaged_smoke

        return run_packaged_smoke(Path(sys.argv[2]),
            simulate_missing_r="--simulate-missing-r" in sys.argv[3:])
    application = QApplication.instance() or QApplication(sys.argv)
    application.setApplicationName("PichAnalysis")
    window = MainWindow()
    window.show()
    return application.exec()

