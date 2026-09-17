from __future__ import annotations

import threading
import time
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Callable


KEGG_REST_BASE = "https://rest.kegg.jp"
KEGG_USAGE_URL = "https://www.kegg.jp/kegg/rest/"
PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


class RateLimiter:
    """Thread-safe start-rate limiter. The default permits at most 3 calls/s."""

    def __init__(self, calls_per_second: float = 3.0, *, clock=time.monotonic, sleeper=time.sleep):
        self.interval = 1.0 / calls_per_second
        self.clock = clock
        self.sleeper = sleeper
        self._last: float | None = None
        self._lock = threading.Lock()

    def wait(self) -> None:
        with self._lock:
            now = self.clock()
            if self._last is not None:
                delay = self.interval - (now - self._last)
                if delay > 0:
                    self.sleeper(delay)
                    now = self.clock()
            self._last = now


class KEGGProvider:
    base_url = KEGG_REST_BASE

    def __init__(self, *, limiter: RateLimiter | None = None, retries: int = 3, timeout: int = 45):
        self.limiter = limiter or RateLimiter()
        self.retries = retries
        self.timeout = timeout

    def fetch(self, route: str) -> bytes:
        last_error: Exception | None = None
        for attempt in range(self.retries + 1):
            self.limiter.wait()
            try:
                request = urllib.request.Request(
                    f"{self.base_url}/{route.lstrip('/')}",
                    headers={"User-Agent": "PichAnalysis/0.1 (academic database downloader)"},
                )
                with urllib.request.urlopen(request, timeout=self.timeout) as response:
                    return response.read()
            except urllib.error.HTTPError as error:
                last_error = error
                if error.code != 429 and not 500 <= error.code < 600:
                    raise
            except (urllib.error.URLError, TimeoutError, OSError) as error:
                last_error = error
            if attempt < self.retries:
                time.sleep(min(2 ** attempt, 4))
        raise RuntimeError(f"KEGG request failed after {self.retries + 1} attempts: {last_error}")

    @staticmethod
    def pathway_routes(pathway_id: str) -> dict[str, str]:
        return {
            "entries": f"get/{pathway_id}",
            "kgml": f"get/{pathway_id}/kgml",
            "images": f"get/{pathway_id}/image",
        }


def decode_text(data: bytes) -> str:
    return data.decode("utf-8")


def parse_pathways(data: bytes) -> list[dict[str, str]]:
    rows = []
    for line in decode_text(data).splitlines():
        if not line.strip():
            continue
        identifier, name = line.split("\t", 1)
        pathway_id = identifier.removeprefix("path:")
        rows.append({"pathway_id": pathway_id, "name": name, "organism_code": "hsa"})
    if not rows:
        raise ValueError("The KEGG pathway list is empty.")
    return rows


def parse_genes(data: bytes) -> list[dict[str, str]]:
    rows = []
    for line in decode_text(data).splitlines():
        if not line.strip():
            continue
        identifier, description = line.split("\t", 1)
        gene_id = identifier.removeprefix("hsa:")
        first = description.split(";", 1)[0]
        symbol = first.split(",", 1)[0].strip() if first else ""
        rows.append({"gene_id": gene_id, "kegg_gene_id": f"hsa:{gene_id}", "symbol": symbol, "raw_description": description})
    if not rows:
        raise ValueError("The KEGG gene list is empty.")
    return rows


def parse_gene_pathway_links(data: bytes) -> list[tuple[str, str]]:
    rows = []
    for line in decode_text(data).splitlines():
        if not line.strip():
            continue
        gene, pathway = line.split("\t", 1)
        rows.append((gene.removeprefix("hsa:"), pathway.removeprefix("path:")))
    if not rows:
        raise ValueError("The KEGG gene-to-pathway table is empty.")
    return rows

def parse_conversion(data: bytes, external_prefix: str) -> list[tuple[str, str]]:
    rows=[]
    for line in decode_text(data).splitlines():
        if not line.strip(): continue
        left,right=line.split("\t",1)
        if left.startswith("hsa:"): kegg,external=left,right
        else: external,kegg=left,right
        rows.append((external.removeprefix(external_prefix),kegg.removeprefix("hsa:")))
    if not rows: raise ValueError("The KEGG identifier conversion table is empty.")
    return rows


def validate_text(data: bytes) -> None:
    if not data.strip() or b"\t" not in data:
        raise ValueError("The downloaded KEGG text table is empty or malformed.")


def validate_png(data: bytes) -> None:
    if not data.startswith(PNG_SIGNATURE):
        raise ValueError("The downloaded pathway image is not a valid PNG file.")


def validate_kgml(data: bytes) -> None:
    try:
        root = ET.fromstring(data)
    except ET.ParseError as error:
        raise ValueError(f"The downloaded KGML is invalid XML: {error}") from error
    if root.tag.rsplit("}", 1)[-1] != "pathway":
        raise ValueError("The downloaded KGML does not have a pathway root element.")


VALIDATORS: dict[str, Callable[[bytes], None]] = {
    "entries": lambda data: None if b"ENTRY" in data else (_ for _ in ()).throw(ValueError("The KEGG entry is malformed.")),
    "kgml": validate_kgml,
    "images": validate_png,
}
