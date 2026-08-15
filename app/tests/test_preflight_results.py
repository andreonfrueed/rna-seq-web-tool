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


def test_find_outputs_includes_quality_reports_excludes_big_files(tmp_path: Path):
    """质控目录：FastQC 报告（html/zip）进结果，修剪后的 fastq 大文件仍排除。"""
    out = _make_output(tmp_path)
    q = out / "1.Quality_and_trimming"
    q.mkdir()
    (q / "trimmed.fastq.gz").write_bytes(b"x")      # 大文件 → 排除
    (q / "S1_fastqc.html").write_text("<html></html>")  # 报告 → 纳入
    (q / "S1_fastqc.zip").write_bytes(b"PK")
    all_files = [p.name for v in results.find_outputs(out).values() for p in v]
    assert "trimmed.fastq.gz" not in all_files
    assert "S1_fastqc.html" in all_files
    assert "S1_fastqc.zip" in all_files


def test_cleanup_intermediates_keeps_quality_reports(tmp_path: Path):
    """清理中间文件：删质控目录的 fastq 大文件，保留 FastQC 报告。"""
    out = _make_output(tmp_path)
    q = out / "1.Quality_and_trimming"
    q.mkdir()
    (q / "trimmed.fastq.gz").write_bytes(b"x" * 1000)
    (q / "S1_fastqc.html").write_text("<html></html>")
    results.cleanup_intermediates(out)
    assert not (q / "trimmed.fastq.gz").exists()   # 大文件被删
    assert (q / "S1_fastqc.html").exists()         # 报告保留
    assert q.exists()                              # 目录不被整个删除


# ---------------------------------------------------------------- pyseqrna 旧版图归档
def _make_out_with_legacy(tmp_path: Path) -> Path:
    """构造含 pyseqrna 旧图 + R 新图混合的输出目录。"""
    out = tmp_path / "output"
    vis = out / "5.Visualization"
    # pyseqrna 旧图（目录级）
    (vis / "Volcano_Plots").mkdir(parents=True)
    (vis / "Volcano_Plots" / "LPS-C_volcano.png").write_bytes(b"x")
    (vis / "MA_Plots").mkdir()
    (vis / "MA_Plots" / "LPS-C_ma.png").write_bytes(b"x")
    (vis / "Venn_Plots").mkdir()
    (vis / "Venn_Plots" / "Venn_1.png").write_bytes(b"x")
    # R 新图（应保留）——Windows 文件系统大小写不敏感，MA_plots 与 MA_Plots
    # 实为同一物理目录，exist_ok 兼容（WSL/Linux 下是独立目录）
    (vis / "Volcano").mkdir()
    (vis / "Volcano" / "C_vs_LPS_volcano.png").write_bytes(b"x")
    (vis / "MA_plots").mkdir(exist_ok=True)
    (vis / "MA_plots" / "C_vs_LPS_MA.png").write_bytes(b"x")
    # Sample_Plots 混合：pyseqrna *_plot.* vs R *_vst.*
    sp = vis / "Sample_Plots"
    sp.mkdir()
    (sp / "All_Samples_PCA_plot.png").write_bytes(b"x")   # 旧
    (sp / "All_Samples_t-SNE_plot.png").write_bytes(b"x")  # 旧
    (sp / "All_Samples_PCA_vst.png").write_bytes(b"x")     # 新
    # Heatmaps 混合：All_* 旧 vs DEG_heatmap_vst 新
    hm = vis / "Heatmaps"
    hm.mkdir()
    (hm / "All_Top_50_genes_heatmap_clustered.png").write_bytes(b"x")  # 旧
    (hm / "DEG_heatmap_vst.png").write_bytes(b"x")                     # 新
    # 5.Clustering 混合
    cl = out / "5.Clustering"
    cl.mkdir()
    (cl / "sample_clustering_samples_dendrogram.png").write_bytes(b"x")  # 旧
    (cl / "sample_clustering_vst_heatmap.png").write_bytes(b"x")         # 新
    return out


def test_is_legacy_rules():
    assert results._is_legacy(Path("5.Visualization/Volcano_Plots/a.png"))
    assert results._is_legacy(Path("5.Visualization/MA_Plots/a.png"))
    assert results._is_legacy(Path("5.Visualization/Sample_Plots/All_Samples_PCA_plot.png"))
    assert results._is_legacy(Path("5.Visualization/Heatmaps/All_Top_50_genes_heatmap_clustered.png"))
    assert results._is_legacy(Path("5.Clustering/sample_clustering_samples_dendrogram.png"))
    # 新版图不算旧图
    assert not results._is_legacy(Path("5.Visualization/Volcano/C_vs_LPS_volcano.png"))
    assert not results._is_legacy(Path("5.Visualization/MA_plots/C_vs_LPS_MA.png"))
    assert not results._is_legacy(Path("5.Visualization/Sample_Plots/All_Samples_PCA_vst.png"))
    assert not results._is_legacy(Path("5.Visualization/Heatmaps/DEG_heatmap_vst.png"))
    assert not results._is_legacy(Path("5.Clustering/sample_clustering_vst_heatmap.png"))


def test_find_outputs_legacy_separate_group(tmp_path: Path):
    out = _make_out_with_legacy(tmp_path)
    groups = results.find_outputs(out)
    legacy_label = "🗄️ PySeqRNA 旧版图（已归档，可忽略）"
    legacy_names = [p.name for p in groups.get(legacy_label, [])]
    # 旧图全在归档组
    assert "LPS-C_volcano.png" in legacy_names
    assert "All_Samples_PCA_plot.png" in legacy_names
    assert "All_Top_50_genes_heatmap_clustered.png" in legacy_names
    assert "sample_clustering_samples_dendrogram.png" in legacy_names
    # 新图不在归档组，且在正常分组里
    normal = [p.name for k, v in groups.items() if k != legacy_label for p in v]
    assert "C_vs_LPS_volcano.png" in normal
    assert "All_Samples_PCA_vst.png" in normal
    assert "DEG_heatmap_vst.png" in normal
    assert "sample_clustering_vst_heatmap.png" in normal
    # 归档组在最后
    assert list(groups)[-1] == legacy_label


def test_make_zip_legacy_in_archive_folder(tmp_path: Path):
    out = _make_out_with_legacy(tmp_path)
    zp = tmp_path / "r.zip"
    results.make_zip(out, zp)
    import zipfile
    names = zipfile.ZipFile(zp).namelist()
    # 旧图收进归档前缀，新图保持原路径
    assert "5.Visualization/PySeqRNA旧版图/Volcano_Plots/LPS-C_volcano.png" in names
    assert "5.Visualization/PySeqRNA旧版图/Sample_Plots/All_Samples_PCA_plot.png" in names
    assert "5.Visualization/PySeqRNA旧版图/5.Clustering/sample_clustering_samples_dendrogram.png" in names
    assert "5.Visualization/Volcano/C_vs_LPS_volcano.png" in names
    assert "5.Visualization/Sample_Plots/All_Samples_PCA_vst.png" in names
    # 正常结果区（归档前缀之外）不再出现旧图目录
    normal_names = [n for n in names if "PySeqRNA旧版图" not in n]
    assert not any("Volcano_Plots/" in n for n in normal_names)
    assert not any("MA_Plots/" in n for n in normal_names)
    assert not any("Venn_Plots/" in n for n in normal_names)


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


def test_ensure_readme_idempotent(tmp_path: Path):
    """BUG-22 回归：《结果说明.txt》幂等写入，解释 4/5 双目录命名。"""
    out = tmp_path / "output"
    out.mkdir()
    p1 = results.ensure_readme(out)
    assert p1.exists()
    content = p1.read_text(encoding="utf-8")
    assert "4.Normalization" in content
    assert "5.Visualization" in content
    assert "不是出错" in content
    mtime1 = p1.stat().st_mtime_ns
    # 幂等：内容已存在时绝不重写
    p2 = results.ensure_readme(out)
    assert p2 == p1
    assert p2.stat().st_mtime_ns == mtime1


def test_ensure_readme_does_not_create_output_dir(tmp_path: Path):
    """BUG-26 回归：output 目录不存在时绝不抢先创建（引擎要求自己建，抢先会拒启）。"""
    out = tmp_path / "output"
    assert not out.exists()
    results.ensure_readme(out)
    assert not out.exists()  # 引擎还没建，绝不能 mkdir


def test_ensure_readme_writes_after_output_exists(tmp_path: Path):
    """BUG-26 回归：引擎建好 output 后，说明文件正常补写。"""
    out = tmp_path / "output"
    out.mkdir()
    p = results.ensure_readme(out)
    assert p.exists()
    assert "不是出错" in p.read_text(encoding="utf-8")


# ---------------------------------------------------------------- collect_vector_images
def test_collect_vector_images_mirrors_structure(tmp_path: Path):
    """矢量图收集：把 output 下各目录的 PDF 按相对路径收到 6.图片源码，PNG 不收。"""
    out = tmp_path / "output"
    (out / "5.Visualization" / "Volcano").mkdir(parents=True)
    (out / "4.Differential_Expression").mkdir(parents=True)
    (out / "5.Visualization" / "Volcano" / "a_volcano.pdf").write_bytes(b"%PDF-1.4")
    (out / "5.Visualization" / "Volcano" / "a_volcano.png").write_bytes(b"png")
    (out / "4.Differential_Expression" / "Filtered_DEG.pdf").write_bytes(b"%PDF-1.4")

    collected = results.collect_vector_images(out)

    vdir = out / "6.图片源码"
    assert (vdir / "5.Visualization" / "Volcano" / "a_volcano.pdf").exists()
    assert (vdir / "4.Differential_Expression" / "Filtered_DEG.pdf").exists()
    assert not (vdir / "5.Visualization" / "Volcano" / "a_volcano.png").exists()
    assert len(collected) == 2


def test_collect_vector_images_no_output_dir(tmp_path: Path):
    """输出目录不存在时不报错，返回空列表。"""
    assert results.collect_vector_images(tmp_path / "nope") == []


def test_collect_vector_images_idempotent_and_no_nesting(tmp_path: Path):
    """幂等：重复收集结果一致，且 6.图片源码 不会把自己再收进去。"""
    out = tmp_path / "output"
    (out / "5.Visualization").mkdir(parents=True)
    (out / "5.Visualization" / "x.pdf").write_bytes(b"%PDF-1.4")

    first = results.collect_vector_images(out)
    second = results.collect_vector_images(out)
    assert len(first) == len(second) == 1
    assert not (out / "6.图片源码" / "6.图片源码").exists()
