"""reference.py 的多参考文件消歧 + CHECKSUMS 校验单元测试（不联网）。"""
from __future__ import annotations
from pathlib import Path

import pytest

from lib import reference as ref


# ---------------------------------------------------------------- _find_existing
def test_find_existing_single(tmp_path: Path):
    d = tmp_path / "hsapiens"
    d.mkdir()
    fa = d / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    gtf = d / "Homo_sapiens.GRCh38.113.gtf"
    fa.write_text(">1\nACGT\n")
    gtf.write_text("chr1\tHAVANA\tgene\t1\t10\t.\t+\t.\tgene_id \"A\";\n")
    got_fa, got_gtf = ref._find_existing(d, "hsapiens")
    assert got_fa == fa
    assert got_gtf == gtf


def test_find_existing_multi_fa_prefers_expected(tmp_path: Path):
    """BUG-02 回归：多个 .fa 时优先选官方名，绝不随机取。"""
    d = tmp_path / "hsapiens"
    d.mkdir()
    expected = d / ref._expected_fa_name("hsapiens")
    other = d / "old_copy.fa"
    expected.write_text("x")
    other.write_text("x")
    fa, _ = ref._find_existing(d, "hsapiens")
    assert fa == expected


def test_find_existing_multi_fa_no_expected_raises(tmp_path: Path):
    """BUG-02 回归：多个 .fa 且都不匹配官方名 → 报错要求清理。"""
    d = tmp_path / "hsapiens"
    d.mkdir()
    (d / "a.fa").write_text("x")
    (d / "b.fa").write_text("x")
    with pytest.raises(RuntimeError, match="无法确定"):
        ref._find_existing(d, "hsapiens")


def test_find_existing_multi_gtf_picks_newest_release(tmp_path: Path):
    d = tmp_path / "hsapiens"
    d.mkdir()
    (d / "Homo_sapiens.GRCh38.112.gtf").write_text("x")
    newer = d / "Homo_sapiens.GRCh38.113.gtf"
    newer.write_text("x")
    _, gtf = ref._find_existing(d, "hsapiens")
    assert gtf == newer


def test_find_existing_multi_gtf_no_match_raises(tmp_path: Path):
    d = tmp_path / "hsapiens"
    d.mkdir()
    (d / "weird.gtf").write_text("x")
    (d / "other.gtf").write_text("x")
    with pytest.raises(RuntimeError, match="无法确定|不匹配"):
        ref._find_existing(d, "hsapiens")


def test_find_existing_empty(tmp_path: Path):
    d = tmp_path / "hsapiens"
    d.mkdir()
    assert ref._find_existing(d, "hsapiens") == (None, None)


# ---------------------------------------------------------------- discover_gtf_url
def test_discover_gtf_url():
    listing = '<a href="Homo_sapiens.GRCh38.113.gtf.gz">x</a>'
    url = ref.discover_gtf_url("hsapiens", "https://ftp.ensembl.org/pub", listing)
    assert url == ("https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/"
                   "Homo_sapiens.GRCh38.113.gtf.gz")


def test_discover_gtf_url_no_match():
    assert ref.discover_gtf_url("hsapiens", "https://x", "<a href='none'>") is None


# ---------------------------------------------------------------- checksum 校验
def test_bsd_sum_deterministic(tmp_path: Path):
    f = tmp_path / "data.bin"
    f.write_bytes(b"hello world" * 1000)
    s1, n1 = ref._bsd_sum(f)
    s2, n2 = ref._bsd_sum(f)
    assert (s1, n1) == (s2, n2)
    assert n1 == f.stat().st_size


def test_verify_checksum_pass_and_fail(tmp_path: Path):
    f = tmp_path / "Homo_sapiens.fa.gz"
    f.write_bytes(b"payload" * 500)
    s, n = ref._bsd_sum(f)
    blocks = (n + 1023) // 1024

    # 匹配：不抛错、文件保留
    ref._verify_checksum(f, {f.name: (s, blocks)})
    assert f.exists()

    # 不匹配：删文件并抛错（SEC-02 回归）
    f.write_bytes(b"corrupted" * 500)
    with pytest.raises(RuntimeError, match="完整性校验失败"):
        ref._verify_checksum(f, {"Homo_sapiens.fa.gz": (s, blocks)})
    assert not f.exists()


def test_verify_checksum_skips_when_absent(tmp_path: Path):
    f = tmp_path / "x.gz"
    f.write_bytes(b"data")
    ref._verify_checksum(f, {})          # 空表跳过
    ref._verify_checksum(f, {"y.gz": (1, 1)})  # 无本文件条目跳过
    assert f.exists()
