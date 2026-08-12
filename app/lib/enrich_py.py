"""GO/KEGG 富集（纯 Python：gseapy + Enrichr 基因集库）。

按"比较 × X高于Y/X低于Y"分别富集——上调和下调基因涉及的通路往往相反，
混在一起做会互相稀释，什么也看不出来。

流程：
1. 读取 pyseqrna 的 diff_genes 目录：每个比较有 比较_up.txt（上调）、
   比较_down.txt（下调）、比较.txt（全部 = 两者并集，有拆分时忽略它）。
2. 从 GTF 离线解析 ENSEMBL → 基因符号（逐行读，避免大文件吃内存）。
3. 首次从 Enrichr 拉取 GO/KEGG 基因集库并缓存到本地（之后离线）。
   缓存文件名带物种，避免人/小鼠共用一份缓存导致张冠李戴；
   写入时附 sha256 旁路校验，防止缓存被篡改后悄悄影响富集结论。
4. 匹配前把基因名和基因库统一转成大写——人基因全大写、小鼠首字母大写，
   大小写不一致会让小鼠 GO 富集静默匹配不到任何结果。
5. gseapy.enrich 做富集 → 表格 CSV + 气泡图 PNG。
6. 统计与跳过原因写进 outdir/_stats.json，网页刷新后仍能展示。

原子提交（v2.1 新增）：富集先写进 <outdir>.partial 临时目录，全部跑完
再整体 rename 为正式目录。中途被打断（浏览器刷新/断连杀掉脚本）只会
留下 .partial 目录，正式目录要么是上一次的完整结果、要么不存在——
网页永远不会把半成品当成正常结果展示。
"""
from __future__ import annotations
import hashlib
import json
import re
import shutil
from pathlib import Path

# gseapy 只在真正跑富集时才需要（顶层导入会让没装 gseapy 的环境
# 连本模块都 import 不了，网页其他功能也跟着挂）——惰性导入。
try:
    import gseapy as gp
except ImportError:  # pragma: no cover - 装了 gseapy 的正式环境不会走到这
    gp = None

GO_LIBS = [
    "GO_Biological_Process_2023",
    "GO_Cellular_Component_2023",
    "GO_Molecular_Function_2023",
]

# 方向 → 排序权重 / 中文名（网页展示用）。这是唯一定义处，
# app.py 直接从这里导入，不再各自维护一份。
DIRECTION_ORDER = {"up": 0, "down": 1, "all": 2}
DIRECTION_CN = {"up": "上调", "down": "下调", "all": "全部"}


def direction_weight(name: str) -> int:
    """方向目录排序权重：up/高于 最前，down/低于 次之（网页展示排序用）。"""
    if "高于" in name or name == "up":
        return 0
    if "低于" in name or name == "down":
        return 1
    return 2


def direction_label(cmp_name: str, direction: str) -> str:
    """方向目录显示名：比较名 X-Y 时用自解释格式（X高于Y / X低于Y），
    旧数据/多横杠名回退到中文方向名（上调/下调/全部）。"""
    parts = [p for p in cmp_name.split("-") if p]
    if len(parts) == 2:
        x, y = _safe_name(parts[0]), _safe_name(parts[1])
        if direction == "up":
            return f"{x}高于{y}"
        if direction == "down":
            return f"{x}低于{y}"
        if direction == "all":
            return "全部"
    return DIRECTION_CN.get(direction, direction)


def _require_gp():
    """确认 gseapy 可用；没装时给出明确提示（惰性导入的配套守卫）。"""
    if gp is None:
        raise RuntimeError(
            "富集分析需要 gseapy 库，但当前环境没有安装。\n"
            "修复：在 WSL 里运行  pip install gseapy  （或重新双击『一键安装.bat』）。")


def _kegg_lib(species: str) -> str:
    # 注意：Enrichr 无 KEGG_2021_Mouse，小鼠库实际名为 KEGG_2019_Mouse
    return "KEGG_2021_Human" if species == "hsapiens" else "KEGG_2019_Mouse"


def _organism(species: str) -> str:
    return "Human" if species == "hsapiens" else "Mouse"


def _parse_attrs(attr_str: str) -> dict[str, str]:
    d: dict[str, str] = {}
    for item in attr_str.split(";"):
        item = item.strip()
        if not item or " " not in item:
            continue
        k, _, v = item.partition(" ")
        d[k] = v.strip('"')
    return d


def parse_gtf_symbols(gtf: Path) -> dict[str, str]:
    """解析 GTF，返回 {ENSEMBL 基因ID: 基因符号}。逐行读，不全量进内存。"""
    mapping: dict[str, str] = {}
    gtf_path = Path(gtf)
    if not gtf_path.exists():
        return mapping
    with open(gtf_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            attrs = _parse_attrs(parts[8])
            gid, gname = attrs.get("gene_id"), attrs.get("gene_name")
            if gid and gname:
                mapping.setdefault(gid, gname)
    return mapping


def _read_gene_list(f: Path) -> set[str]:
    return {ln.strip() for ln in f.read_text(encoding="utf-8", errors="replace").splitlines()
            if ln.strip() and not ln.strip().startswith("#")}


def collect_deg_sets(diff_genes_dir: Path) -> dict[str, dict[str, list[str]]]:
    """把 diff_genes 目录归组成 {比较名: {"up"/"down"/"all": [基因ID]}}。

    pyseqrna 每个比较写三个 txt：比较.txt（全部）、比较_up.txt（上调）、
    比较_down.txt（下调）。有上/下调拆分时忽略同比较的"全部"文件（它就是
    并集，放进去会重复计算）；没有拆分的老结果按 "all" 处理。
    """
    d = Path(diff_genes_dir)
    groups: dict[str, dict[str, set[str]]] = {}
    if not d.exists():
        return {}
    for f in sorted(d.glob("*.txt")):
        stem = f.stem
        if stem.endswith("_up"):
            cmp_name, direction = stem[:-3], "up"
        elif stem.endswith("_down"):
            cmp_name, direction = stem[:-5], "down"
        else:
            cmp_name, direction = stem, "all"
        groups.setdefault(cmp_name, {}).setdefault(direction, set()).update(_read_gene_list(f))
    out: dict[str, dict[str, list[str]]] = {}
    for cmp_name, dirs in groups.items():
        if "up" in dirs or "down" in dirs:
            dirs.pop("all", None)  # "全部"就是 up∪down，有拆分时丢弃
        out[cmp_name] = {k: sorted(v) for k, v in dirs.items() if v}
    return {k: v for k, v in out.items() if v}


def _safe_name(name: str) -> str:
    """比较名会进目录名，去掉路径不友好字符。"""
    cleaned = re.sub(r'[\\/:*?"<>|\s]+', "_", name).strip("_")
    if not cleaned or set(cleaned) <= {".", "_"}:
        return "x"
    return cleaned


def _upper_lib(lib: dict[str, list[str]]) -> dict[str, list[str]]:
    """把基因库里所有基因名转成大写（查询侧同样转大写，做到大小写不敏感匹配）。"""
    return {k: sorted({g.upper() for g in v}) for k, v in lib.items()}


def _ensure_library(name: str, species: str, cache_dir: Path) -> dict[str, list[str]]:
    """返回 Enrichr 基因集库字典（已转大写）；已缓存则离线加载。

    缓存文件名带物种：GO 库人/小鼠同名，不带物种会出现"先跑哪个物种，
    另一个物种就一直用错库"的静默错误。

    完整性（v2.1 新增）：缓存旁路记录 sha256（<缓存名>.sha256）。
    读取时校验，不匹配视为被篡改/损坏，删掉重新联网下载——
    防止有人改动本地缓存后悄悄影响富集结论。
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    organism = _organism(species)
    cache_file = cache_dir / f"{name}__{organism}.json"
    sum_file = cache_file.with_suffix(cache_file.suffix + ".sha256")

    lib = None
    if cache_file.exists():
        try:
            raw = cache_file.read_text(encoding="utf-8")
            if sum_file.exists():
                expected = sum_file.read_text(encoding="utf-8").strip().split()[0]
                actual = hashlib.sha256(raw.encode("utf-8")).hexdigest()
                if expected != actual:
                    cache_file.unlink(missing_ok=True)
                    sum_file.unlink(missing_ok=True)
                    raw = None
            if raw is not None:
                lib = json.loads(raw)
        except Exception:
            lib = None  # 缓存损坏时降级为重新联网下载
    if lib is None:
        lib = gp.get_library(name=name, organism=organism)
        raw = json.dumps(lib, ensure_ascii=False)
        tmp = cache_file.with_suffix(cache_file.suffix + ".tmp")
        tmp.write_text(raw, encoding="utf-8")
        tmp.replace(cache_file)
        sum_file.write_text(
            hashlib.sha256(raw.encode("utf-8")).hexdigest() + "\n", encoding="utf-8")
    return _upper_lib(lib)


def _run_ora(symbols: list[str], gene_sets: dict[str, list[str]], outdir: Path):
    """gseapy 富集（ORA），返回结果 DataFrame。"""
    res = gp.enrich(
        gene_list=symbols,
        gene_sets=gene_sets,
        no_plot=True,
        outdir=str(outdir),
        verbose=False,
    )
    return res.results


def _save_dotplot(df, png_path: Path, title: str, column: str = "Adjusted P-value") -> None:
    n = min(15, len(df))
    if n == 0:
        return
    try:
        import math
        import textwrap

        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import pandas as pd
        from matplotlib.colors import LinearSegmentedColormap

        top = df.sort_values(column).head(n)
        adjusted = pd.to_numeric(top[column], errors="coerce")
        valid = adjusted.notna()
        if not valid.any():
            return
        top = top.loc[valid].copy()
        adjusted = adjusted.loc[valid]

        def _overlap_size(value) -> int:
            try:
                return int(str(value).split("/", 1)[0])
            except (TypeError, ValueError):
                return 40

        sizes = [max(40, min(260, _overlap_size(v))) for v in top["Overlap"]]
        x = [-math.log10(max(float(v), 1e-300)) for v in adjusted]
        y = list(range(len(top)))[::-1]
        terms = [textwrap.fill(str(v), width=45) for v in top["Term"]]

        fig, ax = plt.subplots(figsize=(8, max(4.5, 0.5 * n)))
        cmap = LinearSegmentedColormap.from_list(
            "soft_blue", ["#C7DBEF", "#2F618C"])
        scatter = ax.scatter(x, y, s=sizes, c=adjusted, cmap=cmap,
                             edgecolors="#6B6B6B", linewidths=0.4, alpha=0.9)
        ax.set_yticks(y)
        ax.set_yticklabels(terms, fontsize=8)
        ax.set_xlabel("-log10(Adjusted P-value)", fontsize=10)
        ax.set_ylabel("Term", fontsize=10)
        ax.grid(True, color="#E3E3E3", linewidth=0.8, alpha=0.7)
        ax.set_axisbelow(True)
        cbar = fig.colorbar(scatter, ax=ax, pad=0.02)
        cbar.set_label("Adjusted P-value", fontsize=9)
        # title 参数保留兼容；新版气泡图按要求不渲染标题。
        fig.tight_layout()
        fig.savefig(png_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
    except Exception:
        return


def _enrich_one(symbols: list[str], go_sets: dict[str, list[str]],
                kegg_sets: dict[str, list[str]], go_universe: set[str],
                kegg_universe: set[str], od: Path,
                produced: dict[str, Path], label: str,
                entry: dict[str, int], skipped: list[str]) -> None:
    """对一组基因（某比较·某方向）跑 GO + KEGG，产物写进 od。

    单组失败（网络、匹配不上等）只记入 skipped，不影响其他组。
    """
    entry["matched_go"] = len(set(symbols) & go_universe)
    if entry["matched_go"] > 0:
        try:
            go_df = _run_ora(symbols, go_sets, od / "_gseapy_go")
            if go_df is not None and len(go_df) > 0:
                go_df.to_csv(od / "GO_result.csv", index=False)
                produced[f"{label}/GO_result.csv"] = od / "GO_result.csv"
                _save_dotplot(go_df, od / "GO_dotplot.png", "GO Enrichment")
                if (od / "GO_dotplot.png").exists():
                    produced[f"{label}/GO_dotplot.png"] = od / "GO_dotplot.png"
        except Exception as e:
            skipped.append(f"{label} 的 GO 富集出错：{str(e)[:60]}")

    entry["matched_kegg"] = len(set(symbols) & kegg_universe)
    if entry["matched_kegg"] > 0:
        try:
            kegg_df = _run_ora(symbols, kegg_sets, od / "_gseapy_kegg")
            if kegg_df is not None and len(kegg_df) > 0:
                kegg_df.to_csv(od / "KEGG_result.csv", index=False)
                produced[f"{label}/KEGG_result.csv"] = od / "KEGG_result.csv"
                _save_dotplot(kegg_df, od / "KEGG_dotplot.png", "KEGG Pathway Enrichment")
                if (od / "KEGG_dotplot.png").exists():
                    produced[f"{label}/KEGG_dotplot.png"] = od / "KEGG_dotplot.png"
        except Exception as e:
            skipped.append(f"{label} 的 KEGG 富集出错：{str(e)[:60]}")


def _rewrite_produced_paths(produced: dict[str, Path], tmp_root: Path, final_root: Path) -> dict[str, Path]:
    """原子提交后把 produced 里的路径从 .partial 临时目录改指正式目录。"""
    out: dict[str, Path] = {}
    for key, p in produced.items():
        try:
            rel = p.relative_to(tmp_root)
            out[key] = final_root / rel
        except ValueError:
            out[key] = p
    return out


def run_enrichment(run_dir: Path, gtf: Path, species: str, cache_dir: Path,
                   outdir: Path, progress_cb=None) -> tuple[dict[str, Path], dict[str, dict[str, int]], list[str]]:
    """对某次分析按 比较×上调/下调 分别跑 GO/KEGG 富集。

    返回 (产物文件 dict, 统计 dict, 跳过原因列表)。统计与跳过原因同时写进
    outdir/_stats.json，网页刷新后仍能展示匹配情况。

    原子提交：全部产物先写进 <outdir>.partial，成功后 rename 为 outdir；
    中途被打断不会留下被网页当成正常结果的半成品。
    progress_cb(进度0-1, 提示语) 可选，网页用它展示富集进度。
    """
    run_dir, gtf, cache_dir, outdir = Path(run_dir), Path(gtf), Path(cache_dir), Path(outdir)
    if species not in ("hsapiens", "mmusculus"):
        raise ValueError(f"不支持的物种: {species}")
    _require_gp()

    def cb(p, msg):
        if progress_cb:
            progress_cb(min(max(p, 0.0), 1.0), msg)

    diff_genes_dir = run_dir / "output" / "4.Differential_Expression" / "diff_genes"
    deg_sets = collect_deg_sets(diff_genes_dir)
    if not deg_sets:
        raise ValueError("没有找到差异基因（diff_genes 目录为空或不存在）")

    cb(0.05, "解析基因注释（GTF）…")
    mapping = parse_gtf_symbols(gtf)
    cb(0.10, "准备 GO/KEGG 基因集库（首次需联网）…")
    go_sets: dict[str, list[str]] = {}
    for n in GO_LIBS:
        go_sets.update(_ensure_library(n, species, cache_dir))
    kegg_sets = _ensure_library(_kegg_lib(species), species, cache_dir)
    go_universe = {g for s in go_sets.values() for g in s}
    kegg_universe = {g for s in kegg_sets.values() for g in s}

    # 临时目录：跑完整体 rename 为正式目录（原子提交）
    tmp_outdir = outdir.with_name(outdir.name + ".partial")
    if tmp_outdir.exists():
        shutil.rmtree(tmp_outdir, ignore_errors=True)  # 上次被打断的残留
    tmp_outdir.mkdir(parents=True, exist_ok=True)

    produced: dict[str, Path] = {}
    stats: dict[str, dict[str, int]] = {}
    skipped: list[str] = []

    tasks: list[tuple[str, str, list[str]]] = []
    for cmp_name in sorted(deg_sets):
        dirs = deg_sets[cmp_name]
        for direction in sorted(dirs, key=lambda k: DIRECTION_ORDER.get(k, 9)):
            tasks.append((cmp_name, direction, dirs[direction]))

    for i, (cmp_name, direction, genes) in enumerate(tasks):
        safe_cmp = _safe_name(cmp_name)
        dir_label = direction_label(cmp_name, direction)
        label = f"{safe_cmp} · {dir_label}"
        cb(0.15 + 0.80 * i / max(len(tasks), 1), f"富集分析：{label}（{i + 1}/{len(tasks)}）")
        mapped = [g for g in genes if mapping.get(g)]
        # 统一大写后再去重（小鼠 Trp53 → TRP53，与人版基因库写法对齐）
        symbols = sorted({mapping[g].upper() for g in mapped})
        entry = {"total": len(genes), "mapped": len(mapped),
                 "unmapped": len(genes) - len(mapped),
                 "matched_go": 0, "matched_kegg": 0}
        stats[label] = entry
        if len(symbols) < 3:
            skipped.append(f"{label}：能映射上基因名的差异基因只有 {len(symbols)} 个，太少")
            continue
        od = tmp_outdir / safe_cmp / direction
        od.mkdir(parents=True, exist_ok=True)
        _enrich_one(symbols, go_sets, kegg_sets, go_universe, kegg_universe,
                    od, produced, label, entry, skipped)
        if entry["matched_go"] == 0 and entry["matched_kegg"] == 0:
            skipped.append(f"{label}：{len(symbols)} 个基因名在 GO/KEGG 库里一个都没匹配上")

    (tmp_outdir / "_stats.json").write_text(
        json.dumps({"stats": stats, "skipped": skipped}, ensure_ascii=False, indent=1),
        encoding="utf-8")

    if not produced:
        detail = "；".join(skipped) if skipped else "无差异基因可用"
        shutil.rmtree(tmp_outdir, ignore_errors=True)
        raise RuntimeError(
            f"所有比较的基因都没能在基因库里匹配上，没有产出任何富集结果（{detail}）。"
            "通常是基因名风格不一致，请把此提示截图反馈。")

    # 原子提交：删旧目录（若存在）→ rename 临时目录为正式目录
    if outdir.exists():
        shutil.rmtree(outdir)
    tmp_outdir.rename(outdir)
    produced = _rewrite_produced_paths(produced, tmp_outdir, outdir)
    cb(1.0, "富集完成")
    return produced, stats, skipped
