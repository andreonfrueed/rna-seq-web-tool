"""plots.py 的离线单元测试（纯集合运算与数据解析，不依赖 matplotlib/sklearn）。"""
from __future__ import annotations
from pathlib import Path

import pytest

from lib import plots


# ---------------------------------------------------------------- read_vst_matrix
def test_read_vst_matrix(tmp_path: Path):
    csv = tmp_path / "vst.csv"
    csv.write_text(
        "Gene,Symbol,S1,S2,S3\n"
        "ENSG1,TP53,1.0,2.0,3.0\n"
        "ENSG2,BRCA1,4.0,5.0,6.0\n"
        "ENSG3,bad,7.0,x,9.0\n",  # 非数字样本值 → 该基因跳过
        encoding="utf-8")
    samples, genes = plots.read_vst_matrix(csv)
    assert samples == ["S1", "S2", "S3"]
    assert genes == {"ENSG1": [1.0, 2.0, 3.0], "ENSG2": [4.0, 5.0, 6.0]}


def test_read_vst_matrix_empty(tmp_path: Path):
    assert plots.read_vst_matrix(tmp_path / "nope.csv") == ([], {})


# ---------------------------------------------------------------- read_sample_conditions
def test_read_sample_conditions(tmp_path: Path):
    ts = tmp_path / "samples.tsv"
    ts.write_text("SampleName\tReplication\tIdentifier\tFile1\tFile2\n"
                  "S1\tS1\tC\ta_R1.fq.gz\ta_R2.fq.gz\n"
                  "S2\tS2\tLPS\tb_R1.fq.gz\tb_R2.fq.gz\n", encoding="utf-8")
    assert plots.read_sample_conditions(ts) == {"S1": "C", "S2": "LPS"}
    assert plots.read_sample_conditions(tmp_path / "nope.tsv") == {}


# ---------------------------------------------------------------- read_deg_sets
def test_read_deg_sets_thresholds(tmp_path: Path):
    d = tmp_path / "diff"
    d.mkdir()
    (d / "DESeq2_LPS_vs_C.csv").write_text(
        "Gene,Symbol,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj\n"
        "G1,TP53,100,2.0,0.1,20,1e-10,1e-9\n"     # 显著上调 ✓
        "G2,BRCA1,100,-1.5,0.1,-15,1e-8,1e-7\n"   # 显著下调 ✓
        "G3,NoSig,100,0.5,0.1,5,0.01,0.4\n"        # padj 不显著 ✗
        "G4,LowFC,100,0.8,0.1,8,1e-5,1e-4\n",      # |lfc|<1 ✗
        encoding="utf-8")
    out = plots.read_deg_sets(d)
    assert out == {"LPS vs C": {"TP53", "BRCA1"}}


# ---------------------------------------------------------------- venn / upset 纯逻辑
def test_venn_regions_two_sets():
    sets = {"A": {"1", "2", "3"}, "B": {"2", "3", "4"}}
    assert plots.venn_regions(sets) == {"A": 1, "B": 1, "AB": 2}


def test_venn_regions_three_sets():
    sets = {"A": {"1", "2", "3"}, "B": {"2", "3", "4"}, "C": {"3", "4", "5"}}
    regions = plots.venn_regions(sets)
    # B-only 与 AC-only 为空 → 不收录；其余区域元素互斥
    assert regions == {"A": 1, "C": 1, "AB": 1, "BC": 1, "ABC": 1}
    assert sum(regions.values()) == len({"1", "2", "3", "4", "5"})


def test_venn_regions_too_few():
    assert plots.venn_regions({"A": {"1"}}) == {}


def test_upset_combinations_exclusive():
    sets = {"A": {"1", "2", "3"}, "B": {"2", "3", "4"}, "C": {"3", "4", "5"}}
    combos = plots.upset_combinations(sets)
    # A∩B（不含 C）= {2}；A∩C（不含 B）= ∅；B∩C（不含 A）= {4}；A∩B∩C = {3}
    assert ("A", "B") in [c for c, _ in combos]
    assert ("A", "B", "C") in [c for c, _ in combos]
    sizes = {c: n for c, n in combos}
    assert sizes[("A", "B")] == 1
    assert sizes[("A", "B", "C")] == 1
    assert ("A", "C") not in sizes  # 空交集不出现


def test_upset_combinations_top_n():
    sets = {f"G{i}": {f"g{i}", "shared"} for i in range(5)}
    combos = plots.upset_combinations(sets, top_n=3)
    assert len(combos) <= 3
    counts = [n for _, n in combos]
    assert counts == sorted(counts, reverse=True)  # 按大小降序


def test_venn_regions_uses_letter_keys_not_names():
    """回归（BUG-27）：集合名是长比较名时，区域键必须是字母 A/B/AB…，不是名字拼接。"""
    sets = {"X vs Y": {"1", "2", "3"}, "Y vs Z": {"2", "3", "4"},
            "X vs Z": {"3", "4", "5"}}
    regions = plots.venn_regions(sets)
    assert set(regions) <= {"A", "B", "C", "AB", "AC", "BC", "ABC"}
    # 与三集合"字母命名"用例结果一致：A=1, C=1, AB=1, BC=1, ABC=1
    assert regions == {"A": 1, "C": 1, "AB": 1, "BC": 1, "ABC": 1}


def test_render_venn_anchors_cover_all_regions():
    """回归（BUG-27）：render_venn 的锚点表必须覆盖 venn_regions 产生的所有键，
    否则未知键会落到原点 (0,0) 叠在一起。"""
    for n in (2, 3):
        sets = {f"S{i}": {str(x) for x in range(i, i + 3)} for i in range(n)}
        regions = plots.venn_regions(sets)
        anchors = plots._VENN3_ANCHORS if n == 3 else plots._VENN2_ANCHORS
        missing = set(regions) - set(anchors)
        assert not missing, f"{n} 集：区域键 {missing} 不在锚点表里"


def _capture_figs(monkeypatch):
    """让 _finalize_figure 只捕获 figure 不落盘，供渲染回归测试测量布局。"""
    captured = {}

    def fake_finalize(fig, path, issues_map, name):
        captured[name] = fig
        return path

    monkeypatch.setattr(plots, "_finalize_figure", fake_finalize)
    return captured


def test_render_venn_and_upset_no_overlap(tmp_path, monkeypatch):
    """回归（BUG-27/28）：真实渲染 Venn/UpSet，自检必须报不出任何标签重叠。"""
    pytest.importorskip("matplotlib")
    from lib import figure_qa as fqa
    import matplotlib.pyplot as plt

    captured = _capture_figs(monkeypatch)
    sets = {
        "C vs LPS": {f"g{i}" for i in range(200)},
        "C vs TTP": {f"g{i}" for i in range(100, 300)},
        "LPS vs TTP": {f"g{i}" for i in range(250, 400)},
    }
    assert plots.render_venn(tmp_path, sets) is not None
    assert plots.render_upset(tmp_path, sets) is not None

    assert set(captured) == {"Venn", "UpSet"}
    for name, fig in captured.items():
        issues = fqa.audit_layout(fig)
        overlaps = [msg for _sev, msg in issues if "重叠" in msg]
        assert not overlaps, f"{name} 图仍有标签重叠: {overlaps}"
        plt.close(fig)


# ---------------------------------------------------------------- render_tsne（sklearn 参数兼容）
def _make_tsne_inputs(tmp_path: Path) -> tuple[Path, Path]:
    """6 样本 × 60 基因的 VST 表 + samples.tsv。"""
    csv = tmp_path / "vst.csv"
    lines = ["Gene,Symbol,S1,S2,S3,S4,S5,S6"]
    for i in range(60):
        lines.append("ENSG%05d,ENSG%05d,%s"
                     % (i, i, ",".join(str((i + j) % 7 + 1) for j in range(6))))
    csv.write_text("\n".join(lines), encoding="utf-8")
    samples = tmp_path / "samples.tsv"
    rows = ["SampleName\tReplication\tIdentifier\tFile1\tFile2"]
    for i in range(1, 7):
        rows.append("S%d\t1\t%s\tf1\tf2" % (i, "CT" if i <= 3 else "LPS"))
    samples.write_text("\n".join(rows), encoding="utf-8")
    return csv, samples


class _FakeFig:
    def savefig(self, p, **kw):
        Path(p).parent.mkdir(parents=True, exist_ok=True)
        Path(p).write_bytes(b"png")

    def tight_layout(self):
        pass


class _FakeAx:
    def scatter(self, *a, **k):
        pass

    def annotate(self, *a, **k):
        pass

    def set_xlabel(self, *a, **k):
        pass

    def set_ylabel(self, *a, **k):
        pass

    def set_title(self, *a, **k):
        pass

    def legend(self, *a, **k):
        pass

    def axis(self, *a, **k):
        pass

    def grid(self, *a, **k):
        pass

    def set_axisbelow(self, *a, **k):
        pass


def _patch_tsne_env(monkeypatch, tsne_cls):
    """用假 sklearn.manifold.TSNE + 假 matplotlib 替换真实依赖。"""
    import sys, types
    import numpy as np

    fake_manifold = types.ModuleType("sklearn.manifold")
    fake_manifold.TSNE = tsne_cls
    fake_sklearn = types.ModuleType("sklearn")
    fake_sklearn.manifold = fake_manifold
    monkeypatch.setitem(sys.modules, "sklearn", fake_sklearn)
    monkeypatch.setitem(sys.modules, "sklearn.manifold", fake_manifold)

    class _FakePlt:
        def subplots(self, **kw):
            return _FakeFig(), _FakeAx()

        def close(self, fig):
            pass

    monkeypatch.setattr(plots, "_try_matplotlib", lambda: _FakePlt())
    return np


def test_render_tsne_new_sklearn_uses_max_iter(monkeypatch, tmp_path: Path):
    """回归（RED-06）：sklearn>=1.6 用 max_iter（n_iter 已移除），t-SNE 必须成功。"""
    calls = {"max_iter": 0, "n_iter": 0}

    class NewTSNE:
        def __init__(self, **kw):
            calls["max_iter"] += int("max_iter" in kw)
            calls["n_iter"] += int("n_iter" in kw)

        def fit_transform(self, mat):
            return np.zeros((mat.shape[0], 2))

    np = _patch_tsne_env(monkeypatch, NewTSNE)
    csv, samples = _make_tsne_inputs(tmp_path)
    out = tmp_path / "out"
    p = plots.render_tsne(out, csv, samples)
    assert p is not None and p.exists()
    assert calls["max_iter"] == 1 and calls["n_iter"] == 0


def test_render_tsne_old_sklearn_falls_back_to_n_iter(monkeypatch, tmp_path: Path):
    """回归：老版 sklearn 只有 n_iter（不认 max_iter）时自动降级，t-SNE 仍成功。"""
    calls = {"n_iter": 0}

    class OldTSNE:
        def __init__(self, **kw):
            if "max_iter" in kw:
                raise TypeError("unexpected keyword 'max_iter'")
            calls["n_iter"] += 1

        def fit_transform(self, mat):
            return np.zeros((mat.shape[0], 2))

    np = _patch_tsne_env(monkeypatch, OldTSNE)
    csv, samples = _make_tsne_inputs(tmp_path)
    out = tmp_path / "out"
    p = plots.render_tsne(out, csv, samples)
    assert p is not None and p.exists()
    assert calls["n_iter"] == 1
