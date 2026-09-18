from __future__ import annotations

import hashlib
import os
import time
import urllib.error
import urllib.request
from pathlib import Path

BASE_URL = "https://personal.broadinstitute.org/scalvo/MitoCarta_Download"
VERSION = "3.0"
FILES = ("Human.MitoCarta3.0.xls", "Human.MitoPathways3.0.gmx")
URLS = {name: f"{BASE_URL}/{name}" for name in FILES}


class MitoCartaDownloadCancelled(RuntimeError):
    pass


class MitoCartaProvider:
    def __init__(self, retries: int = 3, timeout: int = 90):
        self.retries, self.timeout = retries, timeout

    def download(self, name: str, destination: Path, progress=None, cancel_requested=None) -> dict:
        if name not in URLS:
            raise ValueError(f"Unsupported MitoCarta file: {name}")
        destination.parent.mkdir(parents=True, exist_ok=True)
        error = None
        for attempt in range(self.retries + 1):
            part = destination.with_name(destination.name + ".part")
            try:
                if cancel_requested and cancel_requested():
                    raise MitoCartaDownloadCancelled("MitoCarta3.0 download was canceled.")
                request = urllib.request.Request(URLS[name], headers={"User-Agent": "PichAnalysis/0.1"})
                digest, size = hashlib.sha256(), 0
                with urllib.request.urlopen(request, timeout=self.timeout) as response, part.open("wb") as output:
                    total = int(response.headers.get("Content-Length", 0))
                    while True:
                        if cancel_requested and cancel_requested():
                            raise MitoCartaDownloadCancelled("MitoCarta3.0 download was canceled.")
                        chunk = response.read(1024 * 1024)
                        if not chunk:
                            break
                        output.write(chunk); digest.update(chunk); size += len(chunk)
                        if progress: progress(name, size, total)
                os.replace(part, destination)
                return {"filename": name, "source_url": URLS[name], "size": size, "sha256": digest.hexdigest()}
            except MitoCartaDownloadCancelled:
                part.unlink(missing_ok=True); raise
            except (OSError, urllib.error.URLError) as exc:
                error = exc; part.unlink(missing_ok=True)
                if attempt < self.retries: time.sleep(min(2 ** attempt, 4))
        raise RuntimeError(f"MitoCarta download failed: {error}")


def validate_downloads(raw: Path) -> None:
    workbook, gmx = raw / FILES[0], raw / FILES[1]
    if not workbook.is_file() or workbook.stat().st_size < 8:
        raise ValueError("Invalid MitoCarta workbook: file is missing or empty.")
    signature = workbook.read_bytes()[:8]
    if signature != bytes.fromhex("D0CF11E0A1B11AE1") and not signature.startswith(b"PK"):
        raise ValueError("Invalid MitoCarta workbook: expected an XLS or OOXML workbook.")
    if not gmx.is_file() or gmx.stat().st_size == 0:
        raise ValueError("Invalid MitoPathways GMX: file is missing or empty.")
    lines = gmx.read_text(encoding="utf-8-sig").splitlines()
    if len(lines) < 3 or len(lines[0].split("\t")) < 1 or len(lines[0].split("\t")) != len(lines[1].split("\t")):
        raise ValueError("Invalid MitoPathways GMX structure.")
