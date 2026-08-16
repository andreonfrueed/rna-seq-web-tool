"""论文级补充图（SCI 风格 300dpi）：VST 版 t-SNE、差异基因 Venn、UpSet。

与 pyseqrna 自带图**并存**（不覆盖）：pyseqrna 的版本基于 RPKM/默认样式，
本模块基于 VST 标准化矩阵与 DESeq2 差异表重新绘制，用于论文展示。

- t-SNE：VST 数据 + scikit-learn（可选依赖，未安装时跳过并提示）
- Venn：≤3 个比较的显著基因重叠（matplotlib 手绘，不用 matplotlib-venn）
- UpSet：多个比较的交集组合柱状图（手绘）

集合运算部分抽成纯函数，便于单元测试；绘图失败只跳过，绝不影响结果页。
"""
from __future__ import annotations

import csv
import itertools
import json
from pathlib import Path

# 与 R 侧一致的低饱和学术配色（DEG 热图/火山图同款）
_PALETTE = ["#C1666B", "#6B8EAE", "#7FA886", "#C9A227", "#8E7CA8", "#B07D4F"]

# 期刊规范预设（借鉴 scipilot-figure-skill 的"按最终尺寸出图/字号可读"原则）：
# 默认 SCI 通用——300 dpi、单栏 3.5 in / 双栏 7.2 in、正文最小字号 6 pt。
# figsize 直接用英寸定最终尺寸，导出后不在 Word/LaTeX 二次缩放。
_JOURNAL_SPECS: dict[str, dict] = {
    "sci_general": {
        "dpi": 300,
        "single_col_in": 3.5,
        "double_col_in": 7.2,
        "min_font_pt": 6,
    },
}


def get_journal_spec(name: str = "sci_general") -> dict:
    """取期刊规范预设（缺失时回退 sci_general）。"""
    return _JOURNAL_SPECS.get(name, _JOURNAL_SPECS["sci_general"])


def _finalize_figure(fig, path: Path, issues_map: dict | None, name: str) -> Path:
    """出图收尾闭环：布局自检 → 保存 PNG → PNG 有效性审计。

    借鉴 scipilot-figure-skill 的视觉自检：程序先抓缺字/裁切/刻度重叠，
    再验落盘文件有效性；结果记入 issues_map（None 则丢弃），供结果页
    汇总成 _figure_qa.json。任何检查失败都不影响出图本身（先出图，
    问题记录在案，绝不让坏图无声进结果页）。
    """
    from lib import figure_qa as fqa

    spec = get_journal_spec()
    issues = fqa.audit_layout(fig)
    try:
        fig.savefig(path, dpi=spec["dpi"], bbox_inches="tight")
    except Exception as e:
        issues.append(("FAIL", f"保存 PNG 失败: {e}"))
        if issues_map is not None:
            issues_map[name] = issues
        return path
    try:
        # 矢量版（PDF）：与 PNG 同名同目录，供「6.图片源码」收集。矢量图无
        # dpi 概念、可无损放大；导出失败不影响 PNG 主结果（仍返回 PNG 路径）。
        fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    except Exception:
        pass
    try:
        # SVG 源码版：svg.fonttype='none' 让文字保留为可编辑的 <text> 元素
        # （matplotlib 默认 'path' 会把文字转成曲线，改不了字）。顾客可直接
        # 改颜色/字号/位置等属性值，无需重算数据。导出失败不影响 PNG/PDF。
        import matplotlib as mpl
        with mpl.rc_context({"svg.fonttype": "none"}):
            fig.savefig(path.with_suffix(".svg"), bbox_inches="tight")
    except Exception:
        pass
    try:
        fig.close()
    except Exception:
        pass
    issues.extend(fqa.audit_png(path))
    if issues_map is not None:
        issues_map[name] = issues
    return path

# 输出位置（相对 outdir = run_dir/output）
_TSNE_REL = "5.Visualization/Sample_Plots/All_Samples_tSNE_vst.png"
_VENN_REL = "5.Visualization/Venn/deg_venn.png"
_UPSET_REL = "5.Visualization/Upset/deg_upset.png"


# ---------------------------------------------------------------- 纯逻辑（可测）

def read_vst_matrix(vst_csv: Path) -> tuple[list[str], dict[str, list[float]]]:
    """读 VST 表，返回 (样本列名, {基因: 每样本值})。

    VST_normalized_counts.csv 首列 Gene、次列 Symbol，其余为样本。
    文件缺失/损坏时返回空，由调用方跳过对应图。
    """
    samples: list[str] = []
    genes: dict[str, list[float]] = {}
    try:
        with open(vst_csv, encoding="utf-8", errors="replace", newline="") as f:
            reader = csv.reader(f)
            header = next(reader, None)
            if not header:
                return samples, genes
            samples = header[2:]  # 跳过 Gene、Symbol
            for row in reader:
                if len(row) < 3:
                    continue
                gid = row[0].strip()
                if not gid:
                    continue
                try:
                    values = [float(x) for x in row[2:]]
                except ValueError:
                    continue
                if len(values) == len(samples):
                    genes[gid] = values
    except OSError:
        pass
    return samples, genes


def read_sample_conditions(sample_sheet: Path) -> dict[str, str]:
    """读 samples.tsv，返回 {SampleName: Identifier(分组)}。"""
    conds: dict[str, str] = {}
    try:
        with open(sample_sheet, encoding="utf-8", errors="replace", newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for r in reader:
                name = (r.get("SampleName") or "").strip()
                ident = (r.get("Identifier") or "").strip()
                if name:
                    conds[name] = ident
    except OSError:
        pass
    return conds


def read_deg_sets(diff_dir: Path, padj_cut: float = 0.05,
                  lfc_cut: float = 1.0) -> dict[str, set[str]]:
    """从 DESeq2 差异表收集每个比较的显著基因集：{比较名: {基因符号大写}}。

    阈值与 R 侧火山图一致（padj<0.05 且 |log2FC|>=1），保证图与表口径相同。
    """
    out: dict[str, set[str]] = {}
    d = Path(diff_dir)
    if not d.exists():
        return out
    for f in sorted(d.glob("DESeq2_*_vs_*.csv")):
        import re
        m = re.match(r"^DESeq2_(.*)_vs_(.*)\.csv$", f.name)
        if not m:
            continue
        genes: set[str] = set()
        with open(f, encoding="utf-8", errors="replace", newline="") as fh:
            for r in csv.DictReader(fh):
                sym = (r.get("Symbol") or "").strip()
                if not sym:
                    continue
                try:
                    lfc = float(r["log2FoldChange"])
                    padj = float(r["padj"])
                except (KeyError, TypeError, ValueError):
                    continue
                if padj < padj_cut and abs(lfc) >= lfc_cut:
                    genes.add(sym.upper())
        if genes:
            out[f"{m.group(1)} vs {m.group(2)}"] = genes
    return out


def venn_regions(sets: dict[str, set[str]]) -> dict[str, int]:
    """把 N 个集合拆成 2^N 个互斥区域（纯集合运算，供 Venn 图标注）。

    返回 {"A", "AB", "ABC", ...} → 元素个数。区域键用字母 A/B/C…，
    按集合在 dict 中的顺序编号（与 render_venn 的圆圈/锚点一一对应）。
    BUG-27：旧实现用集合名拼接当键（如 "C vs LPSC vs TTP"），而
    render_venn 的 anchor 用字母键，两者对不上 → 所有数字叠在原点。
    """
    names = list(sets)
    n = len(names)
    if n < 2:
        return {}
    letters = [chr(ord("A") + i) for i in range(n)]
    regions: dict[str, int] = {}
    for mask in range(1, 1 << n):
        in_idx = [i for i in range(n) if mask & (1 << i)]
        in_names = [names[i] for i in in_idx]
        elem = set.intersection(*(sets[x] for x in in_names))
        for i in range(n):
            if i not in in_idx:
                elem -= sets[names[i]]
        if elem:
            regions["".join(letters[i] for i in in_idx)] = len(elem)
    return regions


def upset_combinations(sets: dict[str, set[str]], top_n: int = 15) -> list[tuple[tuple[str, ...], int]]:
    """UpSet 交集组合：[(成员元组, 交集大小)]，按大小降序取 top_n。"""
    names = list(sets)
    combos: list[tuple[tuple[str, ...], int]] = []
    for k in range(2, len(names) + 1):
        for combo in itertools.combinations(names, k):
            inter = set.intersection(*(sets[c] for c in combo))
            # 只保留「不属于任何其他集合」的独有交集，避免重复计数
            others = set(names) - set(combo)
            for o in others:
                inter -= sets[o]
            if inter:
                combos.append((combo, len(inter)))
    combos.sort(key=lambda x: x[1], reverse=True)
    return combos[:top_n]


# ---------------------------------------------------------------- 绘图

def _try_matplotlib():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def render_tsne(outdir: Path, vst_csv: Path, sample_sheet: Path,
                qa_map: dict | None = None) -> Path | None:
    """VST 版 t-SNE（300dpi，按分组着色）；sklearn 缺失时返回 None。"""
    samples, genes = read_vst_matrix(vst_csv)
    if len(samples) < 4 or len(genes) < 50:
        return None  # 样本/基因太少，t-SNE 不可靠
    try:
        from sklearn.manifold import TSNE
    except ImportError:
        return None
    import numpy as np

    mat = np.array([genes[g] for g in sorted(genes)], dtype=float).T  # 样本×基因
    perplexity = min(30, max(2, len(samples) - 1))
    # sklearn>=1.6 把 n_iter 改名为 max_iter（1.9 已移除旧名），兼容两者
    tsne_kwargs = dict(n_components=2, perplexity=perplexity, random_state=42,
                       init="pca", learning_rate="auto")
    try:
        tsne = TSNE(**tsne_kwargs, max_iter=1000)
    except TypeError:
        tsne = TSNE(**tsne_kwargs, n_iter=1000)  # 老版本 sklearn
    xy = tsne.fit_transform(mat)

    conds = read_sample_conditions(sample_sheet)
    labels = [conds.get(s, "未分组") for s in samples]
    uniq = list(dict.fromkeys(labels))
    color_of = {c: _PALETTE[i % len(_PALETTE)] for i, c in enumerate(uniq)}

    plt = _try_matplotlib()
    fig, ax = plt.subplots(figsize=(7, 5.5))
    for i, s in enumerate(samples):
        ax.scatter(xy[i, 0], xy[i, 1], s=70, color=color_of[labels[i]],
                   edgecolors="#333333", linewidths=0.4, alpha=0.92, zorder=3)
        ax.annotate(s, (xy[i, 0], xy[i, 1]), fontsize=7, color="#444444",
                    xytext=(5, 5), textcoords="offset points")
    for c in uniq:
        ax.scatter([], [], s=60, color=color_of[c], label=c)
    ax.legend(frameon=False, fontsize=9)
    ax.set_xlabel("t-SNE 1", fontsize=10)
    ax.set_ylabel("t-SNE 2", fontsize=10)
    ax.grid(True, color="#E3E3E3", linewidth=0.8, alpha=0.6)
    ax.set_axisbelow(True)
    fig.tight_layout()

    p = outdir / _TSNE_REL
    p.parent.mkdir(parents=True, exist_ok=True)
    return _finalize_figure(fig, p, qa_map, "t-SNE")


# ---------------------------------------------------------------- Venn 几何常量
# 三圆 Venn（A/B 在上，C 在下）区域锚点 = 各互斥区域的质心（离线采样计算），
# 保证数字落在各自区域中央、互不重叠（BUG-27 配套：键用字母，锚点也按字母对齐）。
_VENN3_CENTERS = [(-0.5, 0.30), (0.5, 0.30), (0.0, -0.52)]
_VENN3_R = 0.95
_VENN3_ANCHORS = {
    "A": (-0.87, 0.52), "B": (0.87, 0.52), "C": (0.0, -0.97),
    "AB": (0.0, 0.67), "AC": (-0.55, -0.29), "BC": (0.55, -0.29),
    "ABC": (0.0, 0.04),
}
# BUG-29：上方两圆（A/B）的组别名原位置 (-1.15,1.02)/(1.15,1.02) 到各自圆心
# 距离约 0.97 ≈ 圆半径 0.95，文字一半压在圈里。改为沿左上/右上射线外移到
# 圈外（到圆心约 1.33），与下方 C 的「全在圈外」效果一致；ylim 上界同步抬高。
_VENN3_NAME_POS = {"A": (-1.5, 1.30), "B": (1.5, 1.30), "C": (0.0, -1.62)}
_VENN3_LIMITS = ((-2.0, 2.0), (-1.8, 1.55))

_VENN2_CENTERS = [(-0.55, 0.0), (0.55, 0.0)]
_VENN2_R = 0.85
_VENN2_ANCHORS = {"A": (-0.95, 0.0), "B": (0.95, 0.0), "AB": (0.0, 0.0)}
_VENN2_NAME_POS = {"A": (-1.35, 0.75), "B": (1.35, 0.75)}
_VENN2_LIMITS = ((-1.9, 1.9), (-1.2, 1.2))


def _venn_letters(n: int) -> list[str]:
    return [chr(ord("A") + i) for i in range(n)]


def render_venn(outdir: Path, sets: dict[str, set[str]],
                qa_map: dict | None = None) -> Path | None:
    """≤3 个比较的显著基因 Venn（手绘，300dpi）。"""
    n = len(sets)
    if n < 2 or n > 3:
        return None
    regions = venn_regions(sets)
    names = list(sets)
    letters = _venn_letters(n)

    plt = _try_matplotlib()
    from matplotlib.patches import Circle
    fig, ax = plt.subplots(figsize=(7, 6))

    if n == 2:
        centers, r = _VENN2_CENTERS, _VENN2_R
        anchor, name_pos = _VENN2_ANCHORS, _VENN2_NAME_POS
        xlim, ylim = _VENN2_LIMITS
    else:
        centers, r = _VENN3_CENTERS, _VENN3_R
        anchor, name_pos = _VENN3_ANCHORS, _VENN3_NAME_POS
        xlim, ylim = _VENN3_LIMITS

    for i, (cx, cy) in enumerate(centers):
        ax.add_patch(Circle((cx, cy), r, alpha=0.32, color=_PALETTE[i % len(_PALETTE)]))

    for key, cnt in regions.items():
        x, y = anchor.get(key, (0, 0))
        ax.text(x, y, str(cnt), ha="center", va="center", fontsize=11,
                fontweight="bold", color="#333333")
    for i, name in enumerate(names):
        x, y = name_pos[letters[i]]
        ax.text(x, y, name, ha="center", va="center", fontsize=11, color="#222222")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.tight_layout()

    p = outdir / _VENN_REL
    p.parent.mkdir(parents=True, exist_ok=True)
    return _finalize_figure(fig, p, qa_map, "Venn")


def render_upset(outdir: Path, sets: dict[str, set[str]], top_n: int = 15,
                 qa_map: dict | None = None) -> Path | None:
    """UpSet 图：顶部交集大小柱状 + 底部集合成员矩阵（手绘，300dpi）。"""
    if len(sets) < 2:
        return None
    combos = upset_combinations(sets, top_n)
    if not combos:
        return None
    names = list(sets)
    n_combo = len(combos)

    plt = _try_matplotlib()
    fig = plt.figure(figsize=(max(8, 0.55 * n_combo), 6.5))
    # 布局（figure 分数坐标）：左侧「集合大小」窄面板，右侧上下两块留足缝隙，
    # 让交集柱状图的纵轴刻度数字（可能到 4~5 位）落在缝隙里、不压左侧条形图
    # （BUG-28：旧位置 [0.28, ...] 与左面板右缘 0.24 只差 0.04，刻度数字溢出压图）。
    ax_bar = fig.add_axes([0.32, 0.52, 0.60, 0.36])
    ax_mat = fig.add_axes([0.32, 0.12, 0.60, 0.32])
    ax_size = fig.add_axes([0.04, 0.12, 0.14, 0.72])

    counts = [n for _, n in combos]
    ax_bar.bar(range(n_combo), counts, color="#6B8EAE", edgecolor="none")
    ax_bar.set_ylabel("Intersection size", fontsize=9)
    ax_bar.tick_params(axis="x", labelbottom=False)
    for i, v in enumerate(counts):
        ax_bar.text(i, v + max(counts) * 0.01, str(v), ha="center",
                    fontsize=7, color="#333333")
    ax_bar.spines[["top", "right"]].set_visible(False)

    for j, name in enumerate(names):
        for i, (combo, _) in enumerate(combos):
            in_set = name in combo
            ax_mat.plot(i, j, "o" if in_set else "o", color=("#2F618C" if in_set else "#D9D9D9"),
                        markersize=7 if in_set else 4)
        # 同一组合的点用线段连接（UpSet 惯例）
        for i in range(n_combo - 1):
            combo = combos[i][0]
            nxt = combos[i + 1][0]
            if combo == nxt:
                continue
            for j, name in enumerate(names):
                if name in combo and name in nxt:
                    ax_mat.plot([i, i + 1], [j, j], color="#B0C4DE", linewidth=2,
                                zorder=0)
    ax_mat.set_yticks(range(len(names)))
    ax_mat.set_yticklabels(names, fontsize=8)
    ax_mat.set_xlim(-0.5, n_combo - 0.5)
    ax_mat.set_ylim(-0.5, len(names) - 0.5)
    ax_mat.tick_params(axis="x", labelbottom=False)
    ax_mat.spines[["top", "right", "left", "bottom"]].set_visible(False)

    sizes = [len(sets[nm]) for nm in names]
    ax_size.barh(range(len(names)), sizes, color="#7FA886", edgecolor="none")
    ax_size.set_yticks(range(len(names)))
    ax_size.set_yticklabels(names, fontsize=8)
    ax_size.invert_yaxis()
    ax_size.set_xlabel("Set size", fontsize=9)
    ax_size.spines[["top", "right"]].set_visible(False)

    p = outdir / _UPSET_REL
    p.parent.mkdir(parents=True, exist_ok=True)
    return _finalize_figure(fig, p, qa_map, "UpSet")


def ensure_aux_plots(run_dir: Path, progress_cb=None) -> dict[str, Path]:
    """确保补充图存在（缺啥补啥），返回 {类别: 路径}。

    数据源：run_dir/output 的 VST 表 + DESeq2 差异表 + run_dir/samples.tsv。
    任何一步失败只跳过对应图，不影响结果页。
    出图同时做自检闭环，结果写 outdir/_figure_qa.json（结果页可展示）。
    """
    run_dir = Path(run_dir)
    outdir = run_dir / "output"
    made: dict[str, Path] = {}
    qa_map: dict[str, list[tuple[str, str]]] = {}

    def cb(p, msg):
        if progress_cb:
            progress_cb(min(max(p, 0.0), 1.0), msg)

    vst_csv = outdir / "4.Normalization" / "VST_normalized_counts.csv"
    if vst_csv.exists() and not (outdir / _TSNE_REL).exists():
        cb(0.1, "绘制 t-SNE（VST 数据）…")
        try:
            p = render_tsne(outdir, vst_csv, run_dir / "samples.tsv", qa_map)
            if p:
                made["t-SNE"] = p
        except Exception:
            pass

    diff_dir = outdir / "4.Differential_Expression"
    try:
        sets = read_deg_sets(diff_dir)
    except Exception:
        sets = {}
    if len(sets) >= 2:
        if not (outdir / _VENN_REL).exists():
            cb(0.55, "绘制差异基因 Venn…")
            try:
                p = render_venn(outdir, sets, qa_map)
                if p:
                    made["Venn"] = p
            except Exception:
                pass
        if not (outdir / _UPSET_REL).exists():
            cb(0.75, "绘制差异基因 UpSet…")
            try:
                p = render_upset(outdir, sets, qa_map=qa_map)
                if p:
                    made["UpSet"] = p
            except Exception:
                pass

    if qa_map:
        try:
            from lib import figure_qa as fqa
            fqa.save_report(qa_map, outdir / "_figure_qa.json")
        except Exception:
            pass

    cb(1.0, "补充图就绪")
    return made
