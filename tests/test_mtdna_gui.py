"""Offscreen mtDNA Evidence navigation and historical-run smoke."""
from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from smoke_mtdna_analysis_gui import main


def test_mtdna_python_r_gui_history():
    main()
