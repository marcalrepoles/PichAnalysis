from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QThread, Signal

from ..core.r_runtime import RRuntime


class MappingWorker(QThread):
    succeeded = Signal(str, str)
    failed = Signal(str, str)

    def __init__(self, runtime: RRuntime, script: Path, arguments: list[str]) -> None:
        super().__init__()
        self.runtime, self.script, self.arguments = runtime, script, arguments

    def run(self) -> None:
        try:
            result = self.runtime.run(self.script, *self.arguments, timeout=1800)
            if result.returncode == 0:
                self.succeeded.emit(result.stdout, result.stderr)
            else:
                self.failed.emit(f"R terminou com código {result.returncode}.", result.stderr)
        except Exception as error:
            self.failed.emit(str(error), "")

