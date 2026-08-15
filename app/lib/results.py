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
    ("6.图片源码", "图片源码（矢量图 Vector）"),
    ("6.Functional_Annotation", "功能富集 Annotation"),
    ("7.Report", "报告 Report"),
    ("2.Alignment/alignment_stats", "比对统计 Alignment Stats"),
]


def _is_result_file(f: Path) -> bool:
    # 质控目录（1.Quality_and_trimming）保留 FastQC/修剪报告（html/zip/png），
    # 修剪后的 fastq 大文件由 _EXCLUDE_SUFFIX 排除
    for part in f.parts:
        if part in _EXCLUDE_DIR_PARTS:
            return False
    if f.name.lower().endswith(tuple(_EXCLUDE_SUFFIX)):
        return False
    return True


# ---------------------------------------------------------------- pyseqrna 旧版图归档
# pyseqrna 自带的旧样式图与 R 论文级版本并存，按用户要求归档到单独文件夹，
# 不混入正常结果（预览/下载分组/ZIP 内统一收进 PySeqRNA旧版图/）。
_LEGACY_DIR_NAMES = {"Volcano_Plots", "MA_Plots", "Venn_Plots"}
_LEGACY_PREFIX = Path("5.Visualization") / "PySeqRNA旧版图"


def _is_legacy(f: Path) -> bool:
    """是否 pyseqrna 旧样式图/表（R 版已替代）。

    目录级：Volcano_Plots/MA_Plots/Venn_Plots 整目录是 pyseqrna 的；
    文件级（混合目录）：Sample_Plots 的 *_plot.*（pyseqrna PCA/t-SNE）、
    Heatmaps 的 All_*（pyseqrna 热图）、5.Clustering 非 vst_heatmap 的
    聚类图/表（R 版是 sample_clustering_vst_heatmap.*）。
    """
    parts = f.parts
    if any(p in _LEGACY_DIR_NAMES for p in parts):
        return True
    name = f.name
    if "_plot." in name:
        return True
    if "Heatmaps" in parts and name.startswith("All_"):
        return True
    if "5.Clustering" in parts and "vst_heatmap" not in name:
        return True
    return False


def find_outputs(outdir: Path) -> dict[str, list[Path]]:
    """返回 {展示名: [结果文件...]}，排除中间大文件；pyseqrna 旧版图单独归档分组。"""
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
            if f.is_file() and _is_result_file(f) and not _is_legacy(f)
            and f.resolve() not in covered
        )
        if files:
            groups[label] = files
            covered.update(f.resolve() for f in files)
    # 兜底：没被分组覆盖的结果文件
    extra = sorted(
        f for f in outdir.rglob("*")
        if f.is_file() and _is_result_file(f) and not _is_legacy(f)
        and f.resolve() not in covered
    )
    if extra:
        groups["其他结果 Other"] = extra
    # pyseqrna 旧版图单独归档（放最后，不污染正常结果）
    legacy = sorted(f for f in outdir.rglob("*") if f.is_file() and _is_legacy(f))
    if legacy:
        groups["🗄️ PySeqRNA 旧版图（已归档，可忽略）"] = legacy
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
        if f.is_file() and _is_result_file(f) and not _is_legacy(f)
    ]
    # pyseqrna 旧版图：ZIP 内统一收进 5.Visualization/PySeqRNA旧版图/，不污染正常结果。
    # 开头的 5.Visualization/ 用归档前缀替换（5.Clustering 等则直接挂到前缀下）
    def _legacy_arc(f: Path) -> Path:
        rel = f.relative_to(outdir)
        rest = rel.parts[1:] if rel.parts[0] == "5.Visualization" else rel.parts
        return _LEGACY_PREFIX / Path(*rest)

    files += [
        (f, _legacy_arc(f))
        for f in outdir.rglob("*")
        if f.is_file() and _is_legacy(f)
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


_README_NAME = "结果说明.txt"
_README_CONTENT = (
    "RNA-seq 分析结果说明\n"
    "====================\n\n"
    "1. 为什么有「4.」和「5.」两个编号的文件夹？\n"
    "   这是分析引擎（PySeqRNA）的固定命名规则，不是出错：\n"
    "   - 编号 4 下面包含「标准化」和「差异表达」两个独立模块，\n"
    "     所以出现 4.Normalization 和 4.Differential_Expression 两个文件夹；\n"
    "   - 编号 5 下面包含「聚类」和「可视化」两个独立模块，\n"
    "     所以出现 5.Clustering 和 5.Visualization 两个文件夹。\n"
    "   这个编号是引擎写死的目录结构，改名会导致网页无法正常读取结果，故保持不变。\n\n"
    "2. 各文件夹含义：\n"
    "   1.Quality_and_trimming       数据质检（FastQC 报告）\n"
    "   2.Alignment                  比对与比对统计\n"
    "   3.Quantification             表达定量（原始计数）\n"
    "   4.Normalization              标准化表达量（RPKM / VST）\n"
    "   4.Differential_Expression    差异表达分析（DESeq2 差异表、差异基因列表）\n"
    "   5.Clustering                 样本聚类（VST 聚类热图）\n"
    "   5.Visualization              可视化图（火山图、MA 图、热图、PCA）\n"
)


def ensure_readme(outdir: Path) -> Path:
    """在结果目录写《结果说明.txt》（幂等），解释 4/5 双目录命名等结构。

    说明文件落在 output 根目录，会自然进入 results.zip 和结果页的「其他结果」分组，
    用户下载后打开文件夹第一眼就能看到。幂等：内容已存在就不重写。

    竞态守卫（BUG-26）：引擎（pyseqrna）非交互模式要求 output 目录必须由它自己
    先创建——若网页侧抢先 mkdir(output)，引擎启动时会因「Output directory already
    exists」抛 FileExistsError 拒启。触发机制：Streamlit 每次 rerun 会执行全部 tab
    的函数体（不只当前激活的 tab），用户点「开始分析」→ st.rerun() → tab_results()
    无条件调用本函数建 output，而引擎冷启动要几秒，轮到它 create_main_output_directory
    时目录已存在。所以 output 不存在时绝不 mkdir，等引擎建好目录后再补写说明文件。
    """
    outdir = Path(outdir)
    p = outdir / _README_NAME
    if p.exists():
        return p
    if not outdir.exists():
        return p  # 引擎还没建 output，绝不抢先创建（见上方竞态说明）
    try:
        p.write_text(_README_CONTENT, encoding="utf-8")
    except OSError:
        pass  # 说明文件写失败不影响结果页
    return p


_VECTOR_DIR = "6.图片源码"


def collect_vector_images(outdir: Path) -> list[Path]:
    """把结果里的矢量图（PDF）收集到「6.图片源码」，保持相对目录结构。

    顾客可无限放大、并据此验证图确由代码对真实数据绘制（非 AI 图像生成）。
    幂等：每次全量重建该目录（避免残留上一轮的旧图）。返回收集到的目标路径。
    """
    outdir = Path(outdir)
    dest_root = outdir / _VECTOR_DIR
    collected: list[Path] = []
    if not outdir.exists():
        return collected
    if dest_root.exists():
        shutil.rmtree(dest_root, ignore_errors=True)
    for pdf in sorted(outdir.rglob("*.pdf")):
        if _VECTOR_DIR in pdf.parts:  # 跳过 6.图片源码 自身，避免递归收集
            continue
        rel = pdf.relative_to(outdir)
        dest = dest_root / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(pdf, dest)
        collected.append(dest)
    return collected
