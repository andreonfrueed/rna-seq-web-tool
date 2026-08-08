"""preflight.py 与 results.py 的单元测试。"""
from __future__ import annotations
import gzip
from pathlib import Path

import pytest

from lib import preflight, results


# ---------------------------------------------------------------- preflight
def test_check_fastq_ok(tmp_path: Path):
    f = tmp_path / "s.fq.gz"
    with gzip.open(f, "wb") as g:
        g.write(b"@read1\nACGT\n+\nIIII\n")
    assert preflight.check_fastq(f) is None


def test_check_fastq_truncated_gz(tmp_path: Path):
    f = tmp_path / "bad.fq.gz"
    with gzip.open(f, "wb") as g:
        g.write(b"@read1\nACGT\n")
    data = bytearray(f.read_bytes())
    data = data[: len(data) // 2]  # 截断 → gzip 解压必失败
    f.write_bytes(bytes(data))
    assert preflight.check_fastq(f) is not None


def test_check_fastq_plain_not_at(tmp_path: Path):
    f = tmp_path / "s.fastq"
    f.write_bytes(b">not a fastq\n")
    assert preflight.check_fastq(f) is not None


def test_check_fastq_plain_ok(tmp_path: Path):
    f = tmp_path / "s.fastq"
    f.write_bytes(b"@read1\nACGT\n+\nIIII\n")
    assert preflight.check_fastq(f) is None


def test_check_fastq_missing_and_empty(tmp_path: Path):
    assert preflight.check_fastq(tmp_path / "nope.fq") is not None
    empty = tmp_path / "empty.fq"
    empty.write_bytes(b"")
    assert preflight.check_fastq(empty) is not None


# ---------------------------------------------------------------- results
def _make_output(root: Path) -> Path:
    out = root / "output"
    (out / "3.Quantification").mkdir(parents=True)
    (out / "3.Quantification" / "counts.csv").write_text("g,c\n")
    (out / "hisat2_results").mkdir()
    (out / "hisat2_results" / "big.bam").write_bytes(b"x")
    (out / "4.Differential_Expression").mkdir()
    (out / "4.Differential_Expression" / "degs.xlsx").write_bytes(b"x")
    # 参考文件残留也要被排除
    (out / "3.Quantification" / "genome.fa").write_text(">1")
    return out


def test_find_outputs_excludes_intermediates(tmp_path: Path):
    out = _make_output(tmp_path)
    groups = results.find_outputs(out)
    all_files = [p.name for v in groups.values() for p in v]
    assert "counts.csv" in all_files
    assert "degs.xlsx" in all_files
    assert "big.bam" not in all_files      # 中间目录排除
    assert "genome.fa" not in all_files    # 参考后缀排除


def test_make_zip_and_signature_cache(tmp_path: Path):
    out = _make_output(tmp_path)
    zp = tmp_path / "r.zip"
    z1, sig1 = results.make_zip(out, zp)
    assert z1.exists()
    mtime1 = zp.stat().st_mtime_ns
    # 内容没变：复用旧包，不重写
    z2, sig2 = results.make_zip(out, zp)
    assert sig2 == sig1
    assert zp.stat().st_mtime_ns == mtime1


def test_make_zip_with_extra_dirs(tmp_path: Path):
    out = _make_output(tmp_path)
    enrich = tmp_path / "enrich_py"
    (enrich / "LPS-C" / "up").mkdir(parents=True)
    (enrich / "LPS-C" / "up" / "GO_result.csv").write_text("a")
    (enrich / "_stats.json").write_text("{}")  # _ 开头应被跳过
    zp = tmp_path / "r.zip"
    results.make_zip(out, zp, extra_dirs=[(enrich, "GO_KEGG_富集")])
    import zipfile
    names = zipfile.ZipFile(zp).namelist()
    assert "GO_KEGG_富集/LPS-C/up/GO_result.csv" in names
    assert not any("_stats.json" in n for n in names)


def test_cleanup_intermediates_frees(tmp_path: Path):
    out = _make_output(tmp_path)
    freed = results.cleanup_intermediates(out)
    assert freed > 0
    assert not (out / "hisat2_results").exists()
    assert not (out / "3.Quantification" / "genome.fa").exists()
    assert (out / "3.Quantification" / "counts.csv").exists()  # 结果保留
