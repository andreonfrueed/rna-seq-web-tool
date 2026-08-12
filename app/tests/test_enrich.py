"""enrich_py.py 的离线单元测试（不联网、不需要 gseapy）。"""
from __future__ import annotations
import json
import math
import types
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


def _make_deseq_csv(tmp_path: Path, name: str = "DESeq2_LPS_vs_C.csv") -> Path:
    """迷你 DESeq2 差异表（真实列名/格式）。"""
    diff = tmp_path / "output" / "4.Differential_Expression"
    diff.mkdir(parents=True, exist_ok=True)
    p = diff / name
    p.write_text(
        "Gene,Symbol,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj\n"
        "ENSG1,TNF,100,2.0,0.1,20,1e-10,1e-9\n"
        "ENSG2,IL6,100,1.8,0.1,18,1e-9,1e-8\n"
        "ENSG3,TP53,100,-2.0,0.1,-20,1e-10,1e-9\n"
        "ENSG4,PTEN,100,-1.5,0.1,-15,1e-8,1e-7\n"
        "ENSG5,BRCA1,100,-1.2,0.1,-12,1e-7,1e-6\n"
        "ENSG6,FAKE1,10,0.5,0.1,5,0.01,0.02\n",
        encoding="utf-8")
    return p


def _fake_prerank_main(res2d_df):
    """主跑 mock：返回带真实 gseapy 列名（FDR q-val）的 res2d。"""
    import types

    def _prerank(rnk=None, gene_sets=None, outdir=None, **kwargs):
        od = str(outdir or "")
        if "GSEA_plots" in od:
            # 曲线图模式：必须显式 format="png"（gseapy>=1.2 默认 pdf）
            assert kwargs.get("format") == "png", \
                f"曲线图必须 format='png'，实际 {kwargs.get('format')}"
            # 模拟 gseapy 真实行为：图输出到 outdir/prerank/ 子目录
            import pathlib
            p = pathlib.Path(od) / "prerank" / "gsea_IL17.png"
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_bytes(b"\x89PNG\r\n\x1a\nfakepng")
            return types.SimpleNamespace(res2d=None)
        return types.SimpleNamespace(res2d=res2d_df)

    return _prerank


def test_run_gsea_full_flow_real_column_names(monkeypatch, tmp_path: Path):
    """回归（RED-06）：真实 gseapy 的 res2d 列名是 'FDR q-val' 而非 'fdr'。

    此前代码只认 'fdr' 列 → 显著通路统计恒为 0、富集曲线图永不生成。
    本测试用真实列名跑通全流程：GSEA_result.csv / NES 条形图 / sig_terms /
    top_term / GSEA_plots 曲线图产物。
    """
    import pandas as pd

    _make_deseq_csv(tmp_path)
    monkeypatch.setattr(ep, "_require_gp", lambda: None)

    # 迷你基因集库：TNF/IL6 上调相关 + TP53/PTEN/BRCA1 下调相关
    monkeypatch.setattr(
        ep, "_ensure_library",
        lambda name, species, cache_dir: {
            "IL17_SIGNALING": ["TNF", "IL6"],
            "P53_PATHWAY": ["TP53", "PTEN", "BRCA1"],
        })

    # 真实 gseapy 列名构造的显著结果（fdr < 0.25）
    res2d = pd.DataFrame({
        "Name": ["prerank", "prerank"],
        "Term": ["IL17_SIGNALING", "P53_PATHWAY"],
        "ES": [0.9, -0.85],
        "NES": [1.6, -1.5],
        "NOM p-val": [0.0, 0.01],
        "FDR q-val": [0.12, 0.20],
        "FWER p-val": [0.05, 0.10],
        "Tag %": [50, 60],
        "Gene %": [10, 12],
        "Lead_genes": ["TNF;IL6", "TP53;PTEN"],
    })
    monkeypatch.setattr(ep, "gp", types.SimpleNamespace(
        prerank=_fake_prerank_main(res2d)))

    out = tmp_path / "out"
    produced, stats, skipped = ep.run_gsea(
        tmp_path, tmp_path / "x.gtf", "hsapiens",
        tmp_path / "cache", out)

    assert not skipped, f"不应有 skipped: {skipped}"
    # 主产物
    assert "LPS-C/GSEA_result.csv" in produced
    assert "LPS-C/GSEA_NES_barplot.png" in produced
    # 显著统计（此前 bug：sig_terms 恒为 0）
    entry = stats["LPS-C"]
    assert entry["sig_terms"] == 2, f"sig_terms 应为 2，实际 {entry}"
    assert entry["top_term"] == "IL17_SIGNALING"
    assert entry["matched"] >= 5
    # 曲线图产物（此前 bug：永不生成）
    curve_keys = [k for k in produced if "GSEA_plots" in k]
    assert len(curve_keys) >= 1, f"应生成富集曲线图，produced={list(produced)}"
    # 原子提交：无 .partial 残留
    assert not (tmp_path / "out.partial").exists()


def test_find_fdr_column_variants():
    """列名解析 helper：兼容 gseapy 不同版本列名（fdr / FDR q-val）。"""
    import pandas as pd
    df1 = pd.DataFrame({"Term": ["a"], "FDR q-val": [0.1]})
    assert ep._find_fdr_col(df1) == "FDR q-val"
    df2 = pd.DataFrame({"Term": ["a"], "fdr": [0.1]})
    assert ep._find_fdr_col(df2) == "fdr"
    df3 = pd.DataFrame({"Term": ["a"]})
    assert ep._find_fdr_col(df3) is None


def test_finalize_retries_drvfs_rename(monkeypatch, tmp_path: Path):
    """回归（RED-06）：WSL 挂载盘（drvfs）下移动目录偶发 PermissionError。

    原子提交应短重试而非直接失败——否则真实 WSL 工作流里富集结果明明
    已完整产出，却因延迟写时序报错、网页显示失败。
    """
    import shutil as sh

    tmp_out = tmp_path / "out.partial"
    tmp_out.mkdir()
    (tmp_out / "a.csv").write_text("x", encoding="utf-8")
    final = tmp_path / "out"
    produced = {"k": tmp_out / "a.csv"}

    real_move = sh.move
    calls = {"n": 0}

    def flaky_move(src, dst):
        calls["n"] += 1
        if calls["n"] <= 2:  # 前两次模拟 drvfs 延迟写失败
            raise PermissionError(13, "Permission denied (simulated drvfs)")
        return real_move(src, dst)

    monkeypatch.setattr(sh, "move", flaky_move)
    out = ep._finalize(tmp_out, final, produced)
    assert calls["n"] == 3, f"应重试 3 次，实际 {calls['n']}"
    assert final.exists() and (final / "a.csv").exists()
    assert not tmp_out.exists()
    assert out["k"] == final / "a.csv"


def test_finalize_falls_back_to_copy_when_rename_always_fails(
        monkeypatch, tmp_path: Path):
    """回归（RED-06）：drvfs 目录级 rename 持续被锁时，fallback 逐文件复制。

    Windows 侧实时扫描可能锁住整个目录的 rename 数秒到数十秒，短重试
    不够；此时必须能退化为复制 + 删源，保证结果可用。
    """
    import os

    tmp_out = tmp_path / "out.partial"
    tmp_out.mkdir()
    (tmp_out / "a.csv").write_text("x", encoding="utf-8")
    final = tmp_path / "out"
    produced = {"k": tmp_out / "a.csv"}

    real_rename = os.rename

    def always_fail(src, dst):
        raise PermissionError(13, "Permission denied (simulated drvfs)")

    monkeypatch.setattr(os, "rename", always_fail)
    try:
        out = ep._finalize(tmp_out, final, produced)
    finally:
        monkeypatch.setattr(os, "rename", real_rename)
    assert final.exists() and (final / "a.csv").exists()
    assert not tmp_out.exists()
    assert out["k"] == final / "a.csv"
