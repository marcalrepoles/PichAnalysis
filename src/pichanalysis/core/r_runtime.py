from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable


@dataclass(frozen=True)
class RResult:
    command: tuple[str, ...]
    returncode: int
    stdout: str
    stderr: str


def _version_key(path: Path) -> tuple[int, ...]:
    match = re.search(r"R-(\d+(?:\.\d+)*)", str(path), re.IGNORECASE)
    return tuple(int(part) for part in match.group(1).split(".")) if match else ()


def discover_rscript(configured: str | None = None, *, platform: str | None = None,
                     program_files: str | None = None,
                     which: Callable[[str], str | None] = shutil.which) -> str | None:
    in_path = which("Rscript")
    if in_path and Path(in_path).is_file():
        return str(Path(in_path))
    if configured and Path(configured).is_file():
        return str(Path(configured))
    current_platform = platform or sys.platform
    if current_platform.startswith("win"):
        root = Path(program_files or os.environ.get("ProgramFiles", r"C:\Program Files")) / "R"
        candidates = [path for path in root.glob("R-*/bin/Rscript.exe") if path.is_file()]
        if candidates:
            return str(max(candidates, key=_version_key))
    if current_platform == "darwin":
        for candidate in (
            Path("/Library/Frameworks/R.framework/Resources/bin/Rscript"),
            Path("/opt/homebrew/bin/Rscript"), Path("/usr/local/bin/Rscript"),
        ):
            if candidate.is_file():
                return str(candidate)
    return None


class RRuntime:
    def __init__(self, configured_executable: str | None = None) -> None:
        self.executable = discover_rscript(configured_executable)

    @property
    def available(self) -> bool:
        return self.executable is not None

    def command(self, script: Path, *arguments: str) -> tuple[str, ...]:
        if not self.executable:
            raise RuntimeError("Rscript não foi encontrado no PATH nem nas instalações padrão.")
        return (self.executable, str(Path(script)), *(str(value) for value in arguments))

    def run(self, script: Path, *arguments: str, timeout: int = 60) -> RResult:
        command = self.command(script, *arguments)
        try:
            completed = subprocess.run(command, capture_output=True, text=True, encoding="utf-8",
                errors="replace", timeout=timeout, check=False, shell=False)
        except (OSError, subprocess.TimeoutExpired) as error:
            raise RuntimeError(f"Falha ao executar Rscript: {error}") from error
        return RResult(command, completed.returncode, completed.stdout, completed.stderr)

    def version(self) -> str | None:
        if not self.executable:
            return None
        try:
            completed = subprocess.run((self.executable, "--version"), capture_output=True,
                text=True, encoding="utf-8", errors="replace", timeout=10,
                check=False, shell=False)
        except (OSError, subprocess.TimeoutExpired):
            return None
        output = (completed.stdout or completed.stderr).strip()
        return output.splitlines()[0] if output else None
