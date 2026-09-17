import pytest
from pichanalysis.core.database_manager import DatabaseManager
from pichanalysis.core.database_registry import DatabaseState
from pichanalysis.core.databases.kegg import RateLimiter,parse_gene_pathway_links,parse_genes,parse_pathways,validate_kgml,validate_png
PNG=b"\x89PNG\r\n\x1a\nmock";KGML=b'<pathway name="path:hsa00010" org="hsa" number="00010" />'
class FakeProvider:
    base_url="https://rest.kegg.jp"
    def __init__(self,fail=None):self.calls=[];self.fail=fail
    @staticmethod
    def pathway_routes(pid):return {"entries":f"get/{pid}","kgml":f"get/{pid}/kgml","images":f"get/{pid}/image"}
    def fetch(self,route):
        self.calls.append(route)
        if route==self.fail:raise RuntimeError("planned failure")
        return {"list/pathway/hsa":b"path:hsa00010\tGlycolysis\npath:hsa00020\tTCA cycle\n","list/hsa":b"hsa:1\tA1BG; description\n","link/pathway/hsa":b"hsa:1\tpath:hsa00010\n"}.get(route,KGML if route.endswith("/kgml") else PNG if route.endswith("/image") else b"ENTRY       hsa00010\n///\n")
def test_parsers_and_validators():
    assert parse_pathways(b"path:hsa00010\tName\n")[0]["pathway_id"]=="hsa00010"
    assert parse_genes(b"hsa:1\tA1BG; description\n")[0]["symbol"]=="A1BG"
    assert parse_gene_pathway_links(b"hsa:1\tpath:hsa00010\n")==[("1","hsa00010")]
    validate_png(PNG);validate_kgml(KGML)
    with pytest.raises(ValueError):validate_png(b"bad")
    with pytest.raises(ValueError):validate_kgml(b"<html />")
def test_rate_limiter_never_exceeds_three_calls_per_second():
    now=[0.0];sleeps=[]
    def sleep(delay):sleeps.append(delay);now[0]+=delay
    limiter=RateLimiter(clock=lambda:now[0],sleeper=sleep)
    for _ in range(4):limiter.wait()
    assert sum(sleeps)==pytest.approx(1.0)
def test_download_creates_complete_atomic_snapshot(tmp_path):
    manager=DatabaseManager(tmp_path,FakeProvider());snapshot=manager.download(["entries","kgml","images"],pathway_subset=["hsa00010"])
    assert manager.manifest(snapshot)["status"]==DatabaseState.READY and manager.active_snapshot()==snapshot
    assert (snapshot/"tables/gene_to_pathway.tsv").is_file() and not list(snapshot.rglob("*.part"))
def test_resume_skips_valid_completed_files(tmp_path):
    provider=FakeProvider("get/hsa00010/kgml");manager=DatabaseManager(tmp_path,provider);first=manager.download(["entries","kgml"],pathway_subset=["hsa00010"])
    assert manager.manifest(first)["status"]==DatabaseState.INCOMPLETE
    provider.fail=None;provider.calls.clear();second=manager.download(["entries","kgml"],resume=True,pathway_subset=["hsa00010"])
    assert second==first and "get/hsa00010" not in provider.calls
    assert "get/hsa00010/kgml" in provider.calls and manager.state()==DatabaseState.READY
def test_failed_update_preserves_previous_active_snapshot(tmp_path):
    manager=DatabaseManager(tmp_path,FakeProvider());active=manager.download([],pathway_subset=[]);manager.provider=FakeProvider("list/hsa")
    with pytest.raises(RuntimeError):manager.download([],pathway_subset=[])
    assert manager.active_snapshot()==active
def test_cancel_marks_snapshot_and_does_not_activate(tmp_path):
    manager=DatabaseManager(tmp_path,FakeProvider())
    def progress(payload):
        if payload["done"]>=3:manager.cancel()
    snapshot=manager.download(["entries"],progress,pathway_subset=["hsa00010","hsa00020"])
    assert manager.manifest(snapshot)["status"]==DatabaseState.CANCELLED and manager.active_snapshot() is None
