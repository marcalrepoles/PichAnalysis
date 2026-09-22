from __future__ import annotations

import hashlib
import json
import os
import time
import urllib.request
from pathlib import Path

VERSION_URL = "https://string-db.org/api/json/version"
TAX_ID = "9606"


class StringDownloadCancelled(RuntimeError):
    pass


class StringProvider:
    def __init__(self, timeout: float = 60, retries: int = 3):
        self.timeout, self.retries = timeout, retries

    def detect_version(self) -> dict:
        request = urllib.request.Request(VERSION_URL, headers={"User-Agent": "PichAnalysis/0.1"})
        data = json.loads(self._read(request).decode("utf-8"))
        if isinstance(data, list): data = data[0]
        version = str(data.get("string_version", "")).strip().removeprefix("v")
        address = str(data.get("string_stable_address") or data.get("stable_address") or "").strip().rstrip("/")
        if not version or not address: raise ValueError("The official STRING version endpoint returned incomplete metadata.")
        return {"string_version": version, "string_stable_address": address}

    @staticmethod
    def filenames(version: str) -> tuple[str, ...]:
        return tuple(f"{TAX_ID}.{stem}.v{version}.txt.gz" for stem in ("protein.links.full", "protein.physical.links.full", "protein.info", "protein.aliases"))

    def download(self, stable_address: str, filename: str, target: Path, progress=None, cancel_requested=None) -> dict:
        collection = filename.removeprefix(f"{TAX_ID}.").removesuffix(".txt.gz")
        url = f"https://stringdb-downloads.org/download/{collection}/{filename}"
        target.parent.mkdir(parents=True, exist_ok=True); part = target.with_name(target.name + ".part"); error = None
        for attempt in range(self.retries):
            try:
                request = urllib.request.Request(url, headers={"User-Agent": "PichAnalysis/0.1"})
                with urllib.request.urlopen(request, timeout=self.timeout) as response, part.open("wb") as output:
                    total = int(response.headers.get("Content-Length") or 0); done = 0; digest = hashlib.sha256()
                    while True:
                        if cancel_requested and cancel_requested(): raise StringDownloadCancelled("STRING download was canceled.")
                        chunk = response.read(1024 * 1024)
                        if not chunk: break
                        output.write(chunk); digest.update(chunk); done += len(chunk)
                        if progress: progress(filename, done, total)
                os.replace(part, target)
                return {"filename": filename, "source": url, "size": target.stat().st_size, "sha256": digest.hexdigest()}
            except StringDownloadCancelled:
                part.unlink(missing_ok=True); raise
            except Exception as exc:
                error = exc; part.unlink(missing_ok=True)
                if attempt + 1 < self.retries: time.sleep(min(2 ** attempt, 4))
        raise RuntimeError(f"Could not download {filename}: {error}") from error

    def _read(self, request) -> bytes:
        error = None
        for attempt in range(self.retries):
            try:
                with urllib.request.urlopen(request, timeout=self.timeout) as response: return response.read()
            except Exception as exc:
                error = exc
                if attempt + 1 < self.retries: time.sleep(min(2 ** attempt, 4))
        raise RuntimeError(f"Could not query the official STRING version endpoint: {error}") from error


