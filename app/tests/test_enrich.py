"""enrich_py.py 的离线单元测试（不联网、不需要 gseapy）。"""
from __future__ import annotations
import json
from pathlib import Path

import pytest

from lib import enrich_py as ep


# ---------------------------------------------------------------- 方向标签（RED-01 单一来源）
def test_direction_label_two_parts():
    assert ep.direction_label("LPS-C", "up") == "LPS高于C"
    assert ep.direction_label("LPS-C", "down") == "LPS低于C"
    assert ep.direction_label("LPS-C", "all") == "全部"


def test_direction_label_fallback():
    # 非两段式比较名回退中文方向名
    assert ep.direction_label("A-B-C", "up") == "上调"
    assert ep.direction_label("weird", "down") == "下调"


def test_direction_weight_order():
    assert ep.direction_weight("X高于Y") == 0
    assert ep.direction_weight("up") == 0
    assert ep.direction_weight("X低于Y") == 1
    assert ep.direction_weight("其他") == 2


def test_safe_name():
    assert ep._safe_name("A/B:C") == "A_B_C"
    assert ep._safe_name("...") == "x"


def test_upper_lib():
    assert ep._upper_lib({"p": ["tp53", "TP53", "brca1"]}) == {"p": ["BRCA1", "TP53"]}


# ---------------------------------------------------------------- GTF 解析
def test_parse_gtf_symbols(tmp_path: Path):
    gtf = tmp_path / "anno.gtf"
    gtf.write_text(
        "#!genome-build GRCh38\n"
        'chr1\tHAVANA\tgene\t11869\t14409\t.\t+\t.\t'
        'gene_id "ENSG00000223972"; gene_name "DDX11L1";\n'
        'chr1\tHAVANA\tgene\t14404\t29570\t.\t-\t.\t'
        'gene_id "ENSG00000227232"; gene_name "WASH7P";\n'
        # 无 gene_name 的行应被跳过
        'chr1\tHAVANA\tgene\t30000\t31000\t.\t+\t.\tgene_id "ENSG999";\n',
        encoding="utf-8")
    m = ep.parse_gtf_symbols(gtf)
    assert m == {"ENSG00000223972": "DDX11L1", "ENSG00000227232": "WASH7P"}


def test_parse_gtf_missing_file(tmp_path: Path):
    assert ep.parse_gtf_symbols(tmp_path / "nope.gtf") == {}


# ---------------------------------------------------------------- diff_genes 归组
def test_collect_deg_sets_split(tmp_path: Path):
    d = tmp_path / "diff_genes"
    d.mkdir()
    (d / "LPS-C.txt").write_text("G1\nG2\nG3\n")      # 全部=并集，有拆分时忽略
    (d / "LPS-C_up.txt").write_text("G1\nG2\n")
    (d / "LPS-C_down.txt").write_text("G3\n")
    out = ep.collect_deg_sets(d)
    assert out == {"LPS-C": {"up": ["G1", "G2"], "down": ["G3"]}}


def test_collect_deg_sets_all_only(tmp_path: Path):
    d = tmp_path / "diff_genes"
    d.mkdir()
    (d / "TTP-C.txt").write_text("G1\n# comment\nG2\n")
    out = ep.collect_deg_sets(d)
    assert out == {"TTP-C": {"all": ["G1", "G2"]}}


def test_collect_deg_sets_empty(tmp_path: Path):
    assert ep.collect_deg_sets(tmp_path / "missing") == {}


# ---------------------------------------------------------------- 基因集库缓存 + sha256（SEC-05）
def test_ensure_library_uses_cache_offline(tmp_path: Path, monkeypatch):
    cache = tmp_path / "cache"
    lib = {"PATH_A": ["tp53", "brca1"]}
    raw = json.dumps(lib, ensure_ascii=False)
    cache.mkdir()
    cache_file = cache / "GO_Test__Human.json"
    cache_file.write_text(raw, encoding="utf-8")
    import hashlib
    (cache / "GO_Test__Human.json.sha256").write_text(
        hashlib.sha256(raw.encode()).hexdigest() + "\n", encoding="utf-8")

    def boom(*a, **k):
        raise AssertionError("有缓存时不该联网")
    monkeypatch.setattr(ep.gp, "get_library", boom) if ep.gp else None
    out = ep._ensure_library("GO_Test", "hsapiens", cache)
    assert out == {"PATH_A": ["BRCA1", "TP53"]}


def test_ensure_library_tampered_cache_redownloads(tmp_path: Path, monkeypatch):
    """SEC-05 回归：sha256 不匹配视为被篡改，删缓存重新下载。"""
    cache = tmp_path / "cache"
    cache.mkdir()
    cache_file = cache / "GO_Test__Human.json"
    cache_file.write_text('{"EVIL": ["gene"]}', encoding="utf-8")
    (cache / "GO_Test__Human.json.sha256").write_text("0" * 64 + "\n", encoding="utf-8")

    calls = []

    class FakeGp:
        @staticmethod
        def get_library(name, organism):
            calls.append((name, organism))
            return {"PATH_A": ["good"]}

    monkeypatch.setattr(ep, "gp", FakeGp)
    out = ep._ensure_library("GO_Test", "hsapiens", cache)
    assert calls == [("GO_Test", "Human")]
    assert out == {"PATH_A": ["GOOD"]}
    # 新缓存与校验和都已写回
    assert json.loads(cache_file.read_text(encoding="utf-8")) == {"PATH_A": ["good"]}


def test_require_gp_without_gseapy(monkeypatch):
    monkeypatch.setattr(ep, "gp", None)
    with pytest.raises(RuntimeError, match="gseapy"):
        ep._require_gp()
