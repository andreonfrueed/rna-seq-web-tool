"""pyseqrna 输出产物发现与打包：返回全部结果文件（排除中间大文件）。

ZIP 带签名缓存：内容（文件数/总大小/最新修改时间）没变就直接复用旧包，
避免每次打开结果页都把几百 MB 重新压一遍。
"""
from __future__ import annotations
import shutil
import zipfile
from pathlib import Path

_EXCLUDE_DIR_PARTS = {
    "bam_preparation", "hisat2_results", "hisat2_index",
    "star_index", "star_results", "_STARtmp", "logs",
}
_EXCLUDE_SUFFIX = {
    ".bam", ".sam", ".fastq", ".fastq.gz", ".fq", ".fq.gz",
    ".fa", ".fasta", ".gtf", ".gff3", ".idx", ".bai",
}

# 输出子目录 → 展示名（按顺序展示）
_GROUP_LABELS = [
    ("8.Enhanced_Results", "增强结果 Enhanced（可读表/富集/投稿图）"),
    ("3.Quantification", "表达矩阵 Counts"),
    ("4.Normalization", "标准化表达量 Normalized"),
    ("4.Differential_Expression", "差异表达 DEGs"),
    ("5.Visualization/Volcano", "火山图 Volcano"),
    ("5.Visualization", "可视化图 Plots"),
    ("5.Clustering", "聚类分析 Clustering"),
    ("6.Functional_Annotation", "功能富集 Annotation"),
    ("7.Report", "报告 Report"),
    ("2.Alignment/alignment_stats", "比对统计 Alignment Stats"),
]


def _is_result_file(f: Path) -> bool:
    for part in f.parts:
        if part in _EXCLUDE_DIR_PARTS:
            return False
    if "1.Quality_and_trimming" in f.parts:
        return False  # 质控/修剪中间产物，不放进结果
    if f.name.lower().endswith(tuple(_EXCLUDE_SUFFIX)):
        return False
    return True


def find_outputs(outdir: Path) -> dict[str, list[Path]]:
    """返回 {展示名: [结果文件...]}，排除中间大文件。"""
    outdir = Path(outdir)
    if not outdir.exists():
        return {}
    groups: dict[str, list[Path]] = {}
    covered: set[Path] = set()
    for rel, label in _GROUP_LABELS:
        d = outdir / rel
        if not d.exists():
            continue
        files = sorted(
            f for f in d.rglob("*")
            if f.is_file() and _is_result_file(f) and f.resolve() not in covered
        )
        if files:
            groups[label] = files
            covered.update(f.resolve() for f in files)
    # 兜底：没被分组覆盖的结果文件
    extra = sorted(
        f for f in outdir.rglob("*")
        if f.is_file() and _is_result_file(f) and f.resolve() not in covered
    )
    if extra:
        groups["其他结果 Other"] = extra
    return groups


def _signature(paths: list[Path]) -> str:
    """文件清单签名：数量 + 总大小 + 最新修改时间。变了才需要重新打包。"""
    total = 0
    newest = 0
    for p in paths:
        st = p.stat()
        total += st.st_size
        newest = max(newest, st.st_mtime_ns)
    return f"{len(paths)}|{total}|{newest}"


def _zip_files(files: list[tuple[Path, Path]], zip_path: Path) -> tuple[Path, str]:
    """对已收集的 (源文件, ZIP 内名) 列表统一签名并打包。"""
    sig = _signature([src for src, _ in files])
    sig_path = zip_path.with_suffix(zip_path.suffix + ".sig")
    if (zip_path.exists() and sig_path.exists()
            and sig_path.read_text(encoding="utf-8") == sig):
        return zip_path, sig  # 内容没变，直接复用旧包

    zip_path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for src, arc in files:
            zf.write(src, arc)
    sig_path.write_text(sig, encoding="utf-8")
    return zip_path, sig


def make_zip(outdir: Path, zip_path: Path,
             extra_dirs: list[tuple[Path, str]] | None = None) -> tuple[Path, str]:
    """打包全部结果文件（排除中间大文件），带签名缓存。

    extra_dirs: 额外塞进 ZIP 的目录列表，每项 (目录路径, ZIP 内前缀)，
    递归收集（跳过 "_" 开头的临时目录，如 gseapy 草稿、统计缓存），
    用于把网页自己做的 GO/KEGG 富集结果（按 比较/方向 分目录）一起打包。
    返回 (zip 路径, 签名)，签名可用于调用方缓存 ZIP 字节。
    """
    outdir, zip_path = Path(outdir), Path(zip_path)
    files: list[tuple[Path, Path]] = [
        (f, f.relative_to(outdir))
        for f in outdir.rglob("*")
        if f.is_file() and _is_result_file(f)
    ]
    for d, prefix in (extra_dirs or []):
        d = Path(d)
        if d.exists():
            for f in sorted(d.rglob("*")):
                rel = f.relative_to(d)
                if f.is_file() and not any(part.startswith("_") for part in rel.parts):
                    files.append((f, Path(prefix) / rel))

    return _zip_files(files, zip_path)


def zip_folder(src: Path, zip_path: Path, prefix: str = "") -> tuple[Path, str]:
    """把单个目录整个打成 ZIP（跳过 "_" 开头的临时条目），带签名缓存。

    用于富集结果的独立一键打包：ZIP 内按 比较/方向 分文件夹。
    """
    src, zip_path = Path(src), Path(zip_path)
    files: list[tuple[Path, Path]] = []
    if src.exists():
        for f in sorted(src.rglob("*")):
            rel = f.relative_to(src)
            if f.is_file() and not any(part.startswith("_") for part in rel.parts):
                files.append((f, Path(prefix) / rel if prefix else rel))
    return _zip_files(files, zip_path)


def cleanup_intermediates(outdir: Path) -> int:
    """删除中间大文件（BAM/索引/修剪后的 fastq 等），返回释放的字节数。

    只删 _EXCLUDE_DIR_PARTS 目录和 _EXCLUDE_SUFFIX 后缀的文件，
    结果表格和图片不受影响。
    """
    outdir = Path(outdir)
    freed = 0
    if not outdir.exists():
        return 0
    for d in sorted(outdir.rglob("*"), reverse=True):
        if d.is_dir() and d.name in _EXCLUDE_DIR_PARTS:
            for f in d.rglob("*"):
                if f.is_file():
                    try:
                        freed += f.stat().st_size
                    except OSError:
                        pass
            shutil.rmtree(d, ignore_errors=True)
    for f in outdir.rglob("*"):
        if f.is_file() and f.name.lower().endswith(tuple(_EXCLUDE_SUFFIX)):
            try:
                freed += f.stat().st_size
                f.unlink()
            except OSError:
                pass
    return freed
