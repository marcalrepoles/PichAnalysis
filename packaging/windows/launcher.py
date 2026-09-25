"""PyInstaller entry point; application code remains in the package."""

if __name__ == "__main__":
    try:
        from pichanalysis.app import main

        raise SystemExit(main())
    except Exception:
        import os
        import traceback
        from pathlib import Path

        diagnostic = os.environ.get("PICHANALYSIS_BOOT_ERROR_FILE")
        if diagnostic:
            Path(diagnostic).write_text(traceback.format_exc(), encoding="utf-8")
        raise
