from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class RResult:
    command: tuple[str, ...]
    returncode: int
    stdout: str
    stderr: str


class RRuntime:
    def __init__(self, executable: str | None = None) -> None:
        self.executable = executable or shutil.which("Rscript")

    @property
    def available(self) -> bool:
        return self.executable is not None

    def run(self, script: Path, *arguments: str, timeout: int = 60) -> RResult:
        if not self.executable:
            raise RuntimeError("Rscript não foi encontrado no PATH.")
        command = (self.executable, str(Path(script)), *(str(value) for value in arguments))
        try:
            completed = subprocess.run(
                command,
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=timeout,
                check=False,
                shell=False,
            )
        except (OSError, subprocess.TimeoutExpired) as error:
            raise RuntimeError(f"Falha ao executar Rscript: {error}") from error
        return RResult(command, completed.returncode, completed.stdout, completed.stderr)

    def version(self) -> str | None:
        if not self.executable:
            return None
        try:
            completed = subprocess.run(
                (self.executable, "--version"),
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=10,
                check=False,
                shell=False,
            )
        except (OSError, subprocess.TimeoutExpired):
            return None
        output = (completed.stdout or completed.stderr).strip()
        return output.splitlines()[0] if output else None

