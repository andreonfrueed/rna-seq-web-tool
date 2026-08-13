"""样本表 / run.ini / 全局配置 的单元测试。"""
from __future__ import annotations
from pathlib import Path

import pytest

from lib import config, config_builder, sample_sheet


# ---------------------------------------------------------------- sample_sheet
def test_sample_sheet_basic(tmp_path: Path):
    samples = [
        {"id": "S1", "r1": "S1_R1.fq.gz", "r2": "S1_R2.fq.gz"},
        {"id": "S2", "r1": "S2_R1.fq.gz", "r2": None},
    ]
    group_of = {"S1": "Ctrl", "S2": "Treat"}
    out = sample_sheet.build_sample_sheet(samples, group_of, tmp_path / "s.tsv")
    lines = out.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "SampleName\tReplication\tIdentifier\tFile1\tFile2"
    assert lines[1] == "S1\tS1\tCtrl\tS1_R1.fq.gz\tS1_R2.fq.gz"
    # Replication 列填样本名（pyseqrna 用它当 sample_id），File2 可空
    assert lines[2] == "S2\tS2\tTreat\tS2_R1.fq.gz\t"


def test_sample_sheet_ungrouped_raises(tmp_path: Path):
    with pytest.raises(ValueError):
        sample_sheet.build_sample_sheet(
            [{"id": "S1", "r1": "a.fq"}], {}, tmp_path / "s.tsv")


# ---------------------------------------------------------------- config_builder
def _params(tmp_path: Path) -> dict:
    return {
        "sample_sheet": tmp_path / "samples.tsv",
        "fastq_dir": tmp_path / "uploads",
        "genome": tmp_path / "ref/genome.fa",
        "gtf": tmp_path / "ref/anno.gtf",
        "outdir": tmp_path / "output",
        "species": "hsapiens",
        "fold_threshold": 2.0,
        "pvalue_threshold": 0.05,
        "threads": 8,
        "memory": 12,
    }


def test_build_ini_contains_key_sections(tmp_path: Path):
    ini = config_builder.build_ini(_params(tmp_path))
    for section in ("[General]", "[Alignment]", "[DifferentialExpression]",
                    "[Computational]", "[FunctionalAnnotation]"):
        assert section in ini
    assert "species = hsapiens" in ini
    assert "fold_threshold = 2.0" in ini
    assert "threads = 8" in ini


def test_build_ini_rejects_unknown_species(tmp_path: Path):
    p = _params(tmp_path)
    p["species"] = "celegans"
    with pytest.raises(ValueError):
        config_builder.build_ini(p)


def test_build_ini_skip_trim_switches_quality_trim(tmp_path: Path):
    ini = config_builder.build_ini({**_params(tmp_path), "skip_trim": True})
    assert "skip_trim = True" in ini
    assert "quality_trim = False" in ini
    ini2 = config_builder.build_ini(_params(tmp_path))
    assert "skip_trim = False" in ini2
    assert "quality_trim = True" in ini2


def test_build_ini_diffexp_tool_defaults_to_deseq2(tmp_path: Path):
    ini = config_builder.build_ini(_params(tmp_path))
    assert "diffexp_tool = deseq2" in ini
    ini2 = config_builder.build_ini(
        {**_params(tmp_path), "diffexp_tool": "pydiffexpress"})
    assert "diffexp_tool = pydiffexpress" in ini2


def test_build_ini_alignment_stats_source_is_logs(tmp_path: Path):
    """BUG-21 回归：比对统计源固定为 logs，避免 pyseqrna 在 auto 模式下
    无视 HISAT2 现成日志、回退去逐个 pysam 扫 9 个约 2GB 的 BAM（实测 5.8 小时）。
    logs 模式直接解析 .err 摘要，秒级完成。"""
    ini = config_builder.build_ini(_params(tmp_path))
    assert "alignment_stats_source = logs" in ini


def test_build_ini_volcano_plot_depends_on_engine(tmp_path: Path):
    # deseq2 引擎关闭 pyseqrna 旧火山图（R 后处理自绘 SCI 版）；
    # pydiffexpress 引擎不经 R 后处理，保留 pyseqrna 自带火山图
    ini_de = config_builder.build_ini(_params(tmp_path))
    assert "volcano_plot = False" in ini_de
    ini_py = config_builder.build_ini(
        {**_params(tmp_path), "diffexp_tool": "pydiffexpress"})
    assert "volcano_plot = True" in ini_py


# ---------------------------------------------------------------- config
def test_load_config_defaults():
    cfg = config.load_config()
    # web_config.yaml 覆盖的键
    assert cfg["workspace_dir"] == "~/rna_web_workspace"
    assert cfg["threads"] == 8
    assert cfg["memory"] == 12
    # config.py 的兜底默认键
    assert cfg["fold_threshold"] == 2.0
    assert cfg["pvalue_threshold"] == 0.05


def test_build_ini_thresholds_have_defaults(tmp_path: Path):
    """RED-10 回归：fold/pvalue 阈值漏传时用安全默认值，而非 KeyError。"""
    p = _params(tmp_path)
    del p["fold_threshold"]
    del p["pvalue_threshold"]
    ini = config_builder.build_ini(p)
    assert "fold_threshold = 2.0" in ini
    assert "pvalue_threshold = 0.05" in ini
