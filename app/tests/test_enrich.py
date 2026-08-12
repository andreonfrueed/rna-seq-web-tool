"""enrich_py.py 的离线单元测试（不联网、不需要 gseapy）。"""
from __future__ import annotations
import json
import math
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


# ---------------------------------------------------------------- 原子提交路径改写
def test_rewrite_produced_paths(tmp_path: Path):
    tmp = tmp_path / "enrich_py.partial"
    final = tmp_path / "enrich_py"
    (tmp / "LPS-C" / "up").mkdir(parents=True)
    p = tmp / "LPS-C" / "up" / "GO_result.csv"
    p.write_text("x")
    produced = {"LPS-C/GO_result.csv": p}
    out = ep._rewrite_produced_paths(produced, tmp, final)
    assert out["LPS-C/GO_result.csv"] == final / "LPS-C" / "up" / "GO_result.csv"


# ---------------------------------------------------------------- GSEA（RED-06）
def test_collect_deseq_tables(tmp_path: Path):
    d = tmp_path / "diff"
    d.mkdir()
    (d / "DESeq2_LPS_vs_C.csv").write_text("x")
    (d / "DESeq2_T_vs_C.csv").write_text("x")
    (d / "other.csv").write_text("x")
    out = ep.collect_deseq_tables(d)
    assert set(out) == {"LPS-C", "T-C"}
    assert ep.collect_deseq_tables(tmp_path / "nope") == {}


def test_build_ranking_signed_pvalue(tmp_path: Path):
    csv = tmp_path / "DESeq2_LPS_vs_C.csv"
    csv.write_text(
        "Gene,Symbol,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj\n"
        "ENSG1,TP53,100,2.0,0.1,20,1e-10,1e-9\n"      # 上调显著 → 正分
        "ENSG2,BRCA1,100,-1.5,0.1,-15,1e-8,1e-7\n"    # 下调显著 → 负分
        "ENSG3,no_padj,10,3.0,0.1,30,0.01,\n"          # padj 缺失 → 过滤
        "ENSG4,bad_lfc,10,abc,0.1,1,0.5,0.6\n"         # lfc 非数字 → 过滤
        "ENSG5,padj_zero,10,1.0,0.1,10,0.0,0.0\n",     # padj=0 → 过滤
        encoding="utf-8")
    rows = ep.build_ranking(csv)
    assert rows == [("TP53", -math.log10(1e-9)), ("BRCA1", -math.log10(1e-7) * -1)]


def test_run_enrichment_method_dispatch(monkeypatch, tmp_path: Path):
    """分发器：method=gsea 走 run_gsea，ora 走 run_ora，非法值抛错。"""
    calls = []

    def fake_gsea(*a, **k):
        calls.append("gsea")
        return {}, {}, []

    def fake_ora(*a, **k):
        calls.append("ora")
        return {}, {}, []

    monkeypatch.setattr(ep, "run_gsea", fake_gsea)
    monkeypatch.setattr(ep, "run_ora", fake_ora)
    ep.run_enrichment(tmp_path, tmp_path, "hsapiens", tmp_path, tmp_path,
                      method="gsea")
    ep.run_enrichment(tmp_path, tmp_path, "hsapiens", tmp_path, tmp_path,
                      method="ora")
    assert calls == ["gsea", "ora"]
    with pytest.raises(ValueError, match="不支持的富集方法"):
        ep.run_enrichment(tmp_path, tmp_path, "hsapiens", tmp_path, tmp_path,
                          method="foo")


def test_run_gsea_requires_deseq_tables(monkeypatch, tmp_path: Path):
    """没有 DESeq2 差异表时 GSEA 直接报错（UI 会提示改用 ORA）。"""
    monkeypatch.setattr(ep, "_require_gp", lambda: None)  # 跳过 gseapy 检查
    with pytest.raises(ValueError, match="DESeq2 差异表"):
        ep.run_gsea(tmp_path, tmp_path / "x.gtf", "hsapiens",
                    tmp_path / "cache", tmp_path / "out")
