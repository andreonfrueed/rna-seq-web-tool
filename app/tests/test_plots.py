"""plots.py 的离线单元测试（纯集合运算与数据解析，不依赖 matplotlib/sklearn）。"""
from __future__ import annotations
from pathlib import Path

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
