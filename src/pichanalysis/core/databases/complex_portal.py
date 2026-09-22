"""Official, manually curated human Complex Portal ComplexTab acquisition."""
from __future__ import annotations

import hashlib
import os
import time
import urllib.request
from pathlib import Path
from urllib.parse import urlparse

CURATED_HUMAN_URL = "https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/9606.tsv"
TERMS_URL = "https://www.ebi.ac.uk/complexportal/about#license_privacy"
CITATION_URL = "https://doi.org/10.1093/nar/gkae1085"


class ComplexPortalCancelled(RuntimeError):
    pass


class ComplexPortalProvider:
    def __init__(self, timeout: float = 60, retries: int = 3):
        self.timeout = timeout
        self.retries = retries

    def download(self, target: Path, progress=None, cancel_requested=None) -> dict:
        """Stream exactly the official curated 9606 file into a staged .part file."""
        target = Path(target)
        target.parent.mkdir(parents=True, exist_ok=True)
        part = target.with_name(target.name + ".part")
        error = None
        for attempt in range(self.retries):
            try:
                request = urllib.request.Request(CURATED_HUMAN_URL, headers={"User-Agent": "PichAnalysis/0.1"})
                with urllib.request.urlopen(request, timeout=self.timeout) as response:
                    final = response.geturl()
                    if urlparse(final).hostname not in {"ftp.ebi.ac.uk", "www.ebi.ac.uk"}:
                        raise ValueError("Complex Portal download redirected outside EMBL-EBI.")
                    total = int(response.headers.get("Content-Length") or 0)
                    metadata = {
                        "source_url": CURATED_HUMAN_URL,
                        "resolved_url": final,
                        "last_modified": response.headers.get("Last-Modified"),
                        "etag": response.headers.get("ETag"),
                        "content_length": total or None,
                    }
                    digest = hashlib.sha256()
                    done = 0
                    with part.open("wb") as output:
                        while True:
                            if cancel_requested and cancel_requested():
                                raise ComplexPortalCancelled("Complex Portal download was canceled.")
                            chunk = response.read(1024 * 1024)
                            if not chunk:
                                break
                            output.write(chunk)
                            digest.update(chunk)
                            done += len(chunk)
                            if progress:
                                progress(done, total)
                    if not done or (total and done != total):
                        raise ValueError("Complex Portal download is empty or incomplete.")
                os.replace(part, target)
                return {**metadata, "filename": target.name, "size": done, "sha256": digest.hexdigest()}
            except ComplexPortalCancelled:
                part.unlink(missing_ok=True)
                raise
            except Exception as exc:
                error = exc
                part.unlink(missing_ok=True)
                if attempt + 1 < self.retries:
                    time.sleep(min(2 ** attempt, 4))
        raise RuntimeError(f"Could not download curated human ComplexTab: {error}") from error
