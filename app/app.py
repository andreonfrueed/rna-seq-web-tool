"""RNA-seq 分析网页外壳（Streamlit）。

启动（在 WSL 内）：
  cd ~/rna_web_app && bash run.sh
"""
from __future__ import annotations
import configparser
import hashlib
import io
import itertools
import json
import os
import re
import shutil
import time
import zipfile
from pathlib import Path

import streamlit as st

try:
    from streamlit_autorefresh import st_autorefresh
except Exception:  # 没装自动刷新包时退回手动刷新
    st_autorefresh = None

from lib import (config, config_builder, enrich_py, env_check, plots,
                 preflight, reference, results, runner, sample_sheet)

VERSION = "2.1"

st.set_page_config(page_title="RNA 分析小助手", page_icon="🧬", layout="wide")

WORKSPACE = config.workspace_dir()
UPLOAD_DIR = WORKSPACE / "uploads"
RUNS_DIR = WORKSPACE / "runs"
REF_DIR = WORKSPACE / "references"
for d in (UPLOAD_DIR, RUNS_DIR, REF_DIR):
    d.mkdir(parents=True, exist_ok=True)

SPECIES_LABEL = {"hsapiens": "人 (GRCh38)", "mmusculus": "小鼠 (GRCm39)"}

_RE_R1 = re.compile(r"(.+)_R1(_\d+)?\.(fastq\.gz|fq\.gz|fastq|fq)$", re.I)
_RE_R2 = re.compile(r"(.+)_R2(_\d+)?\.(fastq\.gz|fq\.gz|fastq|fq)$", re.I)
# BUG-18：_1/_2 命名此前不支持 lane 后缀（sample_1_001.fastq.gz 识别失败），
# 与 _R1/_R2 对齐补上 (_\d+)?
_RE_1 = re.compile(r"(.+)_1(_\d+)?\.(fastq\.gz|fq\.gz|fastq|fq)$", re.I)
_RE_2 = re.compile(r"(.+)_2(_\d+)?\.(fastq\.gz|fq\.gz|fastq|fq)$", re.I)
_RE_SINGLE = re.compile(r"(.+)\.(fastq\.gz|fq\.gz|fastq|fq)$", re.I)
# ZIP 下载的两级阈值（v2.1 修复 BUG-01）：
# - 缓存阈值：ZIP 字节进 session_state 内存的上限
# - 下载阈值：download_button 的数据经 WebSocket 下发，Streamlit 默认
#   单消息上限约 200MB（config.toml server.maxMessageSize=200），
#   超限浏览器端必然失败。超过下载阈值的包改走「复制到上传文件夹」，
#   用户用『打开上传文件夹.bat』自取，不再硬塞给浏览器。
_MAX_ZIP_CACHE_BYTES = 150 * 1024 * 1024
_MAX_ZIP_DOWNLOAD_BYTES = 180 * 1024 * 1024

# ---------------------------------------------------------------- UI 样式

_CSS = """
<style>
.block-container {padding-top: 1.5rem; max-width: 1180px;}
.rna-banner {
  background: linear-gradient(120deg, #1d4ed8 0%, #2563eb 48%, #0ea5e9 100%);
  border-radius: 16px; padding: 26px 30px; color: #fff; margin-bottom: 4px;
  box-shadow: 0 4px 18px rgba(37, 99, 235, .18);
}
.rna-banner h1 {margin: 0; font-size: 1.75rem; color: #fff; letter-spacing: .5px;}
.rna-banner .sub {color: #dbeafe; margin-top: 6px; font-size: .95rem;}
.rna-steps {display: flex; gap: 8px; margin: 12px 0 20px; flex-wrap: wrap;}
.rna-step {
  flex: 1; min-width: 110px; text-align: center; font-size: .85rem;
  background: #f1f5fb; color: #64748b; border: 1px solid #e2e8f0;
  border-radius: 999px; padding: 7px 10px;
}
.rna-step.done {background: #e0ecff; color: #1d4ed8; border-color: #bfdbfe;}
.rna-step.active {background: #2563eb; color: #fff; border-color: #2563eb; font-weight: 600;}
div[data-testid="stMetric"] {
  background: #f8fafc; border: 1px solid #e2e8f0;
  border-radius: 12px; padding: 12px 14px;
}
section[data-testid="stSidebar"] {background: #f8fafc;}
section[data-testid="stSidebar"] .block-container {padding-top: 2rem;}
</style>
"""


def _inject_css() -> None:
    st.markdown(_CSS, unsafe_allow_html=True)


def _banner() -> None:
    st.markdown(
        '<div class="rna-banner"><h1>🧬 RNA 分析小助手</h1>'
        '<div class="sub">上传测序数据 → 自动比对 → 差异基因 → 富集分析，'
        '全程在你自己的电脑上完成，数据不外传</div></div>',
        unsafe_allow_html=True,
    )


def _stepper() -> None:
    """顶部步骤条：数据 → 分组 → 运行 → 结果，点亮已完成的步骤。"""
    species = st.session_state.get("species", "hsapiens")
    has_refs = bool(st.session_state.get("refs_map", {}).get(species))
    has_samples = bool(st.session_state.get("samples_sel"))
    has_groups = bool(st.session_state.get("group_of"))
    running = runner.find_active_run(RUNS_DIR) is not None
    done = bool(st.session_state.get("done_run"))
    states = [has_refs and has_samples, has_groups, running or done, done]
    labels = ["📤 数据与参考", "👥 分组参数", "▶️ 运行", "📥 结果"]
    current = next((i for i, ok in enumerate(states) if not ok), len(states))
    pills = []
    for i, lab in enumerate(labels):
        cls = "done" if states[i] else ("active" if i == current else "")
        pills.append(f'<div class="rna-step {cls}">{lab}</div>')
    st.markdown(f'<div class="rna-steps">{"".join(pills)}</div>', unsafe_allow_html=True)


def _sidebar() -> None:
    with st.sidebar:
        st.markdown("## 🧬 RNA 分析小助手")
        st.caption(f"v{VERSION} · 本地运行 · 数据不上传")
        st.divider()
        try:
            free, where = preflight.disk_free_gb(WORKSPACE)
            st.metric(f"磁盘剩余（{where}）", f"{free:.0f} GB")
        except Exception:
            pass
        st.caption(f"工作区：`{WORKSPACE}`")
        active = runner.find_active_run(RUNS_DIR)
        if active:
            st.info(f"▶️ 正在分析：{active[0]}")
        st.divider()
        with st.expander("💡 使用小贴士"):
            st.markdown(
                "- 分析要跑几小时，**先把电脑睡眠关掉**\n"
                "- 大文件拖进 uploads 文件夹后点『扫描已有文件』\n"
                "- 每组建议至少 **3 个样本**\n"
                "- 跑完几次后记得在结果页**清理中间文件**"
            )


# ---------------------------------------------------------------- 样本识别

def _classify_pair(files: list[Path], paired: bool = True) -> tuple[list[dict], list[str]]:
    """把文件名配对成样本；返回 (samples, messages)。

    样本名里的空格/中文等特殊字符自动改成下划线（只改样本名，不改文件名）。
    双端模式下缺 R1 或 R2 的样本会被排除并明确提示，不再静默丢弃。
    """
    raw_to_sid: dict[str, str] = {}
    renames: list[str] = []

    def sid_of(raw: str) -> str:
        if raw in raw_to_sid:
            return raw_to_sid[raw]
        sid = re.sub(r"[^A-Za-z0-9_.-]+", "_", raw).strip("_") or "sample"
        base, i = sid, 2
        while sid in raw_to_sid.values():
            sid = f"{base}_{i}"
            i += 1
        raw_to_sid[raw] = sid
        if sid != raw:
            renames.append(f"样本名「{raw}」含空格或特殊字符，已自动改为「{sid}」（文件名不变）")
        return sid

    samples: dict[str, dict] = {}
    errors: list[str] = []

    def slot(raw: str, which: str, fname: str) -> None:
        sid = sid_of(raw)
        s = samples.setdefault(sid, {"id": sid, "r1": "", "r2": ""})
        if s[which]:
            errors.append(f"样本 {sid} 已有 {which.upper()} 文件（{s[which]}），{fname} 被忽略")
        else:
            s[which] = fname

    for f in files:
        name = f.name
        m1 = _RE_R1.match(name) or _RE_1.match(name)
        m2 = _RE_R2.match(name) or _RE_2.match(name)
        if paired:
            if m1:
                slot(m1.group(1), "r1", name)
            elif m2:
                slot(m2.group(1), "r2", name)
            elif _RE_SINGLE.match(name):
                errors.append(f"{name} 找不到配对的 R1/R2 命名，双端模式下被忽略")
            else:
                errors.append(f"{name} 无法识别，请用 sample_R1/sample_R2 或 sample_1/sample_2 命名")
        else:
            if m2:
                errors.append(f"{name} 是 R2 文件，单端模式下被忽略")
                continue
            m = m1 or _RE_SINGLE.match(name)
            if m:
                slot(m.group(1), "r1", name)
            else:
                errors.append(f"{name} 无法识别，请检查文件名")

    if paired:
        ok, bad = [], []
        for s in samples.values():
            (ok if s["r1"] and s["r2"] else bad).append(s)
        for s in bad:
            miss = "R2" if s["r1"] else "R1"
            errors.append(f"样本 {s['id']} 缺少 {miss} 文件，已被排除（双端模式需要 R1+R2 成对）")
        final = ok
    else:
        final = [s for s in samples.values() if s["r1"]]

    lane_groups: dict[str, list[str]] = {}
    lane_tags: dict[str, list[str]] = {}
    for s in final:
        m = re.search(r"_L\d+$", s["id"], flags=re.I)
        if m:
            base = s["id"][:m.start()]
            lane_groups.setdefault(base, []).append(s["id"])
            lane_tags.setdefault(base, []).append(m.group(0)[1:])
    for base, lanes in lane_groups.items():
        if len(lanes) > 1:
            tags = lane_tags[base]
            desc = "/".join(tags[:2]) + ("…" if len(tags) > 2 else "")
            errors.append(
                f"检测到样本 {base} 疑似多个 lane（{desc}），当前会被当作独立样本；"
                "测序仪分 lane 的数据请先用 cat 合并成一个 R1/R2 文件再上传")
    return final, renames + errors


def _upload_digest(f) -> str:
    """文件内容指纹：首尾各 1MB + 大小做 md5（BUG-03 修复）。

    旧签名只用 (文件名, 大小)：重传同名同大小但内容不同的文件时签名不变，
    不会覆盖，静默沿用旧数据。加入内容指纹后，内容变了就一定会重新保存。
    只读首尾 2MB，几 GB 的文件也几乎不花时间。
    """
    try:
        f.seek(0)
        head = f.read(1024 * 1024)
        f.seek(0, 2)
        size = f.tell()
        f.seek(max(0, size - 1024 * 1024))
        tail = f.read(1024 * 1024)
        f.seek(0)
        return hashlib.md5(head + tail + str(size).encode()).hexdigest()
    except Exception:
        return ""


def _check_upload_head(f) -> bool:
    """上传内容快速校验（SEC-03）：.gz 验 gzip 魔数 1f 8b，
    未压缩 fastq 验首字符 '@'。挡住拖错的/恶意的非测序文件。"""
    try:
        f.seek(0)
        head = f.read(4)
        f.seek(0)
    except Exception:
        return True  # 读不了就放行，开跑前还有完整检查兜底
    if not head:
        return False  # 空文件
    if f.name.lower().endswith(".gz"):
        return head[:2] == b"\x1f\x8b"
    return head[:1] == b"@"


def _disk_room_for(total_bytes: int) -> str | None:
    """上传前磁盘余量检查（SEC-03）：总量超过剩余空间 90% 就拒绝，
    防止误拖超大文件把磁盘写满。"""
    if total_bytes <= 0:
        return None
    free, where = preflight.disk_free_gb(UPLOAD_DIR)
    if total_bytes > free * 1e9 * 0.9:
        return (f"上传文件总共约 {total_bytes / 1e9:.1f} GB，但磁盘（{where}）"
                f"只剩 {free:.0f} GB，写不下这么多数据。"
                "请先清理磁盘，或改用小批量上传。")
    return None


def _save_uploads(uploaded) -> tuple[list[Path], list[str], list[str], list[str]]:
    """把网页上传的文件写入 uploads 目录；同名文件直接覆盖（避免改了数据重传却被静默忽略）。

    先写 .part 临时文件再 os.replace 原子替换，失败时保留旧文件。
    返回 (已保存路径, 被覆盖名, 运行中未更新名, 保存失败名)。

    BUG-04 修复：分析运行中同名文件跳过覆盖时，不再把旧文件冒充为
    「本次上传的有效样本」返回——旧逻辑会让用户误以为更正后的数据已生效。
    现在单独放进 skipped 列表并明确提示。
    """
    paths, replaced, skipped, failed = [], [], [], []
    active = runner.find_active_run(RUNS_DIR)
    for f in uploaded:
        dest = UPLOAD_DIR / Path(f.name).name  # 只取文件名，防路径穿越
        tmp = dest.parent / (dest.name + ".part")
        if dest.exists() and active:
            skipped.append(f.name)
            continue
        if not _check_upload_head(f):
            st.error(f"文件 {f.name} 内容不对（不是有效的 fastq/fastq.gz），已拒收")
            failed.append(f.name)
            continue
        try:
            if dest.exists():
                replaced.append(f.name)
            f.seek(0)
            with open(tmp, "wb") as out:
                shutil.copyfileobj(f, out)
            os.replace(tmp, dest)
            paths.append(dest)
        except Exception as e:
            tmp.unlink(missing_ok=True)
            st.error(f"保存文件 {f.name} 失败：{e}")
            failed.append(f.name)
    return paths, replaced, skipped, failed


def _show_classification(files: list[Path], paired: bool) -> None:
    """识别样本 → 存 session → 展示并让用户勾选本次要分析的样本。"""
    samples, msgs = _classify_pair(files, paired)
    st.session_state["classified_paired"] = paired
    st.session_state["files"] = files
    for m in msgs:
        st.warning(m)
    if not samples:
        st.error("没有识别到可用样本，请根据上面的提示检查文件名。")
        st.session_state["samples_sel"] = []
        return

    ids = [s["id"] for s in samples]
    prev = st.session_state.get("sel_samples")
    default = [i for i in (prev or ids) if i in ids] or ids
    sel = st.multiselect(
        f"识别到 **{len(samples)}** 个样本，勾选本次要分析的（不想用的旧样本可以去掉勾，或到下方删除文件）",
        ids, default=default, key="sel_samples",
    )
    st.session_state["samples_sel"] = [s for s in samples if s["id"] in sel]
    for s in samples:
        mark = "✅" if s["id"] in sel else "⬜"
        r2 = f"  ｜  R2: `{s['r2']}`" if s["r2"] else "（单端）"
        st.write(f"{mark} **{s['id']}**　R1: `{s['r1']}`{r2}")


# ---------------------------------------------------------------- 五个页面

def tab_env() -> None:
    st.header("🔍 环境体检")
    st.caption("检查分析要用的程序是否就绪。缺什么，页面会告诉你怎么办。")
    if st.button("🔁 重新检查") or "env_results" not in st.session_state:
        with st.spinner("检查中…"):
            st.session_state["env_results"] = env_check.check_all()
    results = st.session_state["env_results"]
    cols = st.columns(2)
    for i, r in enumerate(results):
        with cols[i % 2]:
            icon = "✅" if r["ok"] else "❌"
            st.write(f"{icon} **{r['name']}** — {r['version']}")
            if not r["ok"] and r["hint"]:
                st.code(r["hint"], language="bash")
    if all(r["ok"] for r in results):
        st.success("全部就绪，可以去传数据了。")


def tab_upload() -> None:
    st.header("📤 数据与参考文件")

    c1, c2 = st.columns(2)
    species = c1.selectbox(
        "物种",
        ["hsapiens", "mmusculus"],
        format_func=lambda s: SPECIES_LABEL[s],
        key="species_sel",
    )
    st.session_state["species"] = species
    paired_label = c2.radio(
        "测序类型",
        ["双端（R1 + R2 成对，最常见）", "单端"],
        horizontal=True,
        key="paired_radio",
    )
    paired = paired_label.startswith("双端")
    st.session_state["paired"] = paired

    # ---- 参考文件：按物种各自记录，切换物种不会串 ----
    refs_map = st.session_state.setdefault("refs_map", {})
    refs = refs_map.get(species)
    with st.container(border=True):
        st.subheader("① 参考文件（基因组 + 基因注释）")
        if refs:
            st.success(f"✅ {SPECIES_LABEL[species]} 参考文件已就绪："
                       f"`{Path(refs['genome']).name}` ｜ `{Path(refs['gtf']).name}`")
            if st.button("重新检查", key="recheck_ref"):
                refs_map.pop(species, None)
                st.rerun()
        else:
            st.caption("第一次使用需要联网下载（几个 GB，约几分钟）；下过一次以后永久离线可用。")
            if st.button("⬇️ 检查 / 下载参考文件", type="primary"):
                with st.status("处理参考文件…", expanded=True) as status:
                    bar = st.progress(0.0)

                    def cb(p, msg):
                        bar.progress(min(p, 1.0))
                        status.update(label=msg)

                    try:
                        refs_map[species] = reference.ensure_reference(species, REF_DIR, cb)
                        status.update(label="参考文件就绪", state="complete")
                        st.rerun()
                    except Exception as e:
                        status.update(label="参考文件准备失败", state="error")
                        st.error(f"参考文件准备失败：{e}\n\n"
                                 f"也可以手动准备：在 WSL 里把 .fa 和 .gtf 放到 "
                                 f"~/rna_web_workspace/references/{species}/ 下，再点上方按钮。")

    st.divider()

    # ---- 测序文件 ----
    with st.container(border=True):
        st.subheader("② 测序文件（fastq）")
        uploaded = st.file_uploader(
            "上传 fastq 文件（可多选；双端需命名 sample_R1/sample_R2 或 sample_1/sample_2）",
            type=["gz", "fastq", "fq"], accept_multiple_files=True,
        )
        files = st.session_state.get("files")
        blocked = st.session_state.setdefault("blocked_uploads", set())
        in_uploader = {f.name for f in uploaded} if uploaded else set()
        blocked &= in_uploader
        if uploaded:
            active_uploaded = [f for f in uploaded if f.name not in blocked]
            if blocked:
                st.warning("已删除的文件仍在上传框中，已拦截回写；"
                           "要重新上传请先在组件里移除再拖入：" + "、".join(sorted(blocked)))
            # BUG-03 修复：签名加入内容指纹，同名同大小但内容不同的文件也会触发重新保存
            sig = [(f.name, f.size, _upload_digest(f)) for f in active_uploaded]
            if sig != st.session_state.get("upload_sig"):
                room_err = _disk_room_for(sum(f.size for f in active_uploaded))
                if room_err:
                    st.error(room_err)
                else:
                    with st.spinner("保存文件…"):
                        paths, replaced, skipped, failed = _save_uploads(active_uploaded)
                    if replaced:
                        st.info(f"已覆盖同名旧文件：{', '.join(replaced)}")
                    if skipped:
                        st.warning("分析运行中，以下文件未更新（仍用磁盘上的旧文件）；"
                                   "如需使用新数据，请等分析结束后重新上传："
                                   + "、".join(skipped))
                        paths += [UPLOAD_DIR / Path(n).name for n in skipped
                                  if (UPLOAD_DIR / Path(n).name).exists()]
                    if not failed:
                        st.session_state["upload_sig"] = sig
                    if paths:
                        st.session_state["files"] = paths
                        files = paths
                    else:
                        st.session_state.pop("files", None)
                        files = None
            else:
                missing_uploaded = [f.name for f in active_uploaded
                                    if not (UPLOAD_DIR / Path(f.name).name).exists()]
                if missing_uploaded:
                    st.warning("这些文件还在上传框中，请在组件里点 x 移除，"
                               "否则下次操作会把它写回来：" + "、".join(missing_uploaded))
                paths = [UPLOAD_DIR / Path(f.name).name for f in active_uploaded
                         if (UPLOAD_DIR / Path(f.name).name).exists()]
                if paths:
                    st.session_state["files"] = paths
                    files = paths
                else:
                    st.session_state.pop("files", None)
                    files = None
        else:
            st.info("大文件（超过 1GB 或一次很多个）建议：双击『打开上传文件夹.bat』，"
                    "把文件拖进弹出的文件夹，然后回来点下面的按钮。")
            if st.button("🔍 扫描已有文件"):
                existing = sorted(set(
                    list(UPLOAD_DIR.glob("*.gz"))
                    + list(UPLOAD_DIR.glob("*.fastq"))
                    + list(UPLOAD_DIR.glob("*.fq"))
                ))
                if existing:
                    st.session_state["files"] = existing
                    files = existing
                else:
                    st.warning("uploads 目录里还没有文件。")
        old_paired = st.session_state.get("classified_paired")
        if files and old_paired is not None and old_paired != paired:
            st.info("测序类型已切换，已重新识别")
        if files:
            _show_classification(files, paired)

        # ---- 文件管理：旧文件可以在这里删掉，避免混进新分析 ----
        all_files = sorted(p.name for p in UPLOAD_DIR.iterdir() if p.is_file())
        if all_files:
            with st.expander(f"🗂️ 管理上传文件夹（共 {len(all_files)} 个文件，可删除旧文件）"):
                doomed = st.multiselect("选择要删除的文件", all_files, key="del_files")
                active = runner.find_active_run(RUNS_DIR)
                if st.button("🗑️ 删除选中文件", disabled=not doomed or bool(active)):
                    for n in doomed:
                        (UPLOAD_DIR / n).unlink(missing_ok=True)
                    still = [n for n in doomed
                             if any(u.name == n for u in (uploaded or []))]
                    st.success(f"已删除 {len(doomed)} 个文件")
                    for key in ("samples", "samples_sel", "files",
                                "sel_samples", "group_of"):
                        st.session_state.pop(key, None)
                    if still:
                        st.session_state.setdefault("blocked_uploads", set()).update(still)
                        st.warning("已删除的文件仍在上传框中，已拦截回写；"
                                   "要重新上传请先在组件里移除再拖入：" + "、".join(still))
                    time.sleep(0.5)
                    st.rerun()
                if active:
                    st.caption("分析运行中，结束后才能删除文件")


def _render_direction_preview(order: list[str]) -> None:
    """预告比较方向：样本表逆序后，vs 对照 = 处理÷对照。"""
    comps = list(itertools.combinations(reversed(order), 2))
    if not comps:
        return
    st.caption(f"本次将计算 {len(comps)} 个比较：")
    for x, y in comps:
        st.caption(f"{x} vs {y}（差异倍数为正 = {x} 组更高）")


def tab_groups() -> None:
    st.header("👥 分组与参数")
    samples = st.session_state.get("samples_sel") or []
    if not samples:
        st.session_state.pop("group_order", None)
        st.warning("请先在『数据与参考文件』页上传或扫描样本。")
        return

    group_names = st.text_input(
        "分组名（用逗号分隔，建议对照组写第一个）", "C, LPS",
        help="差异比较的方向以结果报告中的比较名为准（如 A_vs_B 表示 A 相对于 B）；"
             "通常把对照组写在第一个",
    )
    invalid: list[str] = []
    names: list[str] = []
    for g in (x.strip() for x in group_names.split(",")):
        if not g:
            continue
        if any(ch in g for ch in "\t\r\n/\\-") or "_vs_" in g or not g.strip("."):
            invalid.append(g)
            continue
        if g not in names:
            names.append(g)
    if invalid:
        st.error("非法分组名：" + "、".join(invalid)
                 + "（不能包含制表/换行/斜杠/横杠，也不能是 . 或 ..）")
        st.session_state["group_of"] = {}
        st.session_state.pop("group_order", None)
        return
    if len(names) < 2:
        st.error("至少需要两个不同的分组（例如 C, LPS），否则没法做差异比较。")
        st.session_state["group_of"] = {}
        st.session_state.pop("group_order", None)
        return

    st.write("为每个样本选择分组（**必须逐个选**，不选不会开跑——防止忘了选导致全被当成对照组）：")
    group_of: dict[str, str] = {}
    unassigned: list[str] = []
    cols_per_row = 2
    for row_start in range(0, len(samples), cols_per_row):
        cols = st.columns(cols_per_row)
        for col, s in zip(cols, samples[row_start:row_start + cols_per_row]):
            with col:
                g = st.radio(f"样本 **{s['id']}**", names, index=None,
                             key=f"grp_{s['id']}", horizontal=True)
            if g is None:
                unassigned.append(s["id"])
            else:
                group_of[s["id"]] = g

    ready = True
    if unassigned:
        st.warning(f"还有 {len(unassigned)} 个样本没分组：{', '.join(unassigned)}")
        ready = False
    else:
        counts = {n: sum(1 for v in group_of.values() if v == n) for n in names}
        used = [n for n in names if counts[n] > 0]
        st.write("　".join(f"**{n}**：{counts[n]} 个样本" for n in used))
        if len(used) < 2:
            st.error("所有样本都在同一组里，没法做差异比较。请把样本分到至少两个组。")
            ready = False
        for n in used:
            if counts[n] < 2:
                st.error(f"组「{n}」只有 {counts[n]} 个样本，每组至少要 2 个重复才能开跑。")
                ready = False
            elif counts[n] < 3:
                st.warning(f"组「{n}」只有 {counts[n]} 个样本，建议每组至少 3 个，"
                           "太少的话差异结果可信度低。")
    st.session_state["group_of"] = group_of if ready else {}
    if ready:
        st.session_state["group_order"] = used  # 只记实际用到的组，空组不进比较
        st.success("分组完成 ✅")
        _render_direction_preview(used)
    else:
        st.session_state.pop("group_order", None)

    st.divider()
    st.subheader("分析参数（默认值即可，不懂不用动）")
    cfg = config.load_config()  # RED-04：参数默认值统一来自 config.py
    c1, c2, c3, c4 = st.columns(4)
    fold = c1.number_input("差异倍数阈值", 1.0, 100.0,
                           float(cfg["fold_threshold"]), 0.5,
                           help="表达量相差多少倍才算差异基因，默认 2 倍；"
                                "DESeq2 引擎下对应 log2FC 阈值")
    pval = c2.number_input("P 值阈值", 0.0001, 0.5,
                           float(cfg["pvalue_threshold"]), 0.01,
                           help="越小越严格，默认 0.05；"
                                "DESeq2 引擎下对应 padj 阈值")
    threads = c3.number_input("线程数", 1, 32, cfg["threads"],
                              help="电脑 CPU 核心越多可以开越大，跑得快些")
    memory = c4.number_input("内存 (GB)", 4, 128, cfg["memory"],
                             help="别超过电脑实际内存，留 2-4GB 给系统")
    st.session_state["fold"] = fold
    st.session_state["pval"] = pval
    st.session_state["threads"] = threads
    st.session_state["memory"] = memory


def tab_run() -> None:
    st.header("▶️ 运行分析")

    # ---- 断线重连：后台有正在跑的分析就直接接上进度 ----
    active = runner.find_active_run(RUNS_DIR)
    if active:
        active_name, active_pid = active
        st.info(f"检测到正在运行的分析「{active_name}」，已自动连接进度。")
        _show_run_progress(RUNS_DIR / active_name, runner.RunningPid(active_pid),
                           active_name, active_pid)
        return

    samples = st.session_state.get("samples_sel") or []
    group_of = st.session_state.get("group_of") or {}
    species = st.session_state.get("species", "hsapiens")
    refs = st.session_state.get("refs_map", {}).get(species)

    missing = []
    if not samples:
        missing.append("样本（去『数据与参考文件』页）")
    if not group_of:
        missing.append("分组（去『分组与参数』页，每个样本都要选组）")
    if not refs:
        missing.append("参考文件（去『数据与参考文件』页下载）")
    if missing:
        st.warning("还缺：" + "、".join(missing))
        return

    # ---- 本次分析概览 ----
    paired = st.session_state.get("paired", True)
    c = st.columns(4)
    c[0].metric("样本数", len(samples))
    c[1].metric("分组", " vs ".join(sorted(set(group_of.values()))))
    c[2].metric("物种", SPECIES_LABEL[species])
    c[3].metric("测序类型", "双端" if paired else "单端")
    _render_direction_preview(st.session_state.get("group_order") or [])
    with st.expander("查看样本与分组明细"):
        st.table({"样本": [s["id"] for s in samples],
                  "分组": [group_of[s["id"]] for s in samples],
                  "R1 文件": [s["r1"] for s in samples]})

    run_name = st.text_input("本次分析名称（只能英文/数字/下划线）", "my_run")
    if not re.match(r"^[A-Za-z0-9_]+$", run_name):
        st.error("名称只能包含字母、数字、下划线")
        return

    c1, c2 = st.columns(2)
    skip_trim = c1.checkbox("数据已经清洗过（跳过修剪，更快更省空间）", value=False)
    alignment_tool = c2.selectbox(
        "比对引擎",
        ["hisat2", "star"],
        index=0,
        help="HISAT2 省内存（推荐本机）；STAR 更主流，但建人/鼠索引需要 30GB+ 内存，本机可能失败",
    )
    c3, c4 = st.columns(2)
    diffexp_label = c3.selectbox(
        "差异分析引擎",
        ["DESeq2（推荐，论文标准）", "pydiffexpress（旧引擎）"],
        index=0,
        help="DESeq2 是当前论文标准差异分析方法（输出 log2FC/padj，"
             "并生成 VST 标准化矩阵用于热图）；pydiffexpress 为旧引擎，"
             "仅作为装不上 R/DESeq2 时的回退",
    )
    volcano_gene_labels = c4.checkbox(
        "火山图标注 Top 差异基因名",
        value=True,
        help="基因很多时可关闭以避免文字堆叠",
    )
    diffexp_tool = "deseq2" if diffexp_label.startswith("DESeq2") else "pydiffexpress"
    st.session_state["diffexp_tool"] = diffexp_tool

    if st.button("🚀 开始分析", type="primary"):
        _start_analysis(run_name, samples, group_of, species, refs, paired,
                        skip_trim, alignment_tool, diffexp_tool,
                        volcano_gene_labels)


def _start_analysis(run_name: str, samples: list[dict], group_of: dict,
                    species: str, refs: dict, paired: bool,
                    skip_trim: bool, alignment_tool: str,
                    diffexp_tool: str, volcano_gene_labels: bool) -> None:
    """开跑前自检（磁盘 + 文件完整性）→ 生成配置 → 后台启动。"""
    if paired:
        missing_r2 = [s["id"] for s in samples if not s.get("r2")]
        if missing_r2:
            st.error("双端模式下以下样本缺少 R2 文件：" + "、".join(missing_r2))
            return

    if not Path(refs["genome"]).exists() or not Path(refs["gtf"]).exists():
        st.error("参考文件已被删除，请回『数据与参考文件』页重新检查")
        return

    errors: list[str] = []
    with st.status("开跑前自检…", expanded=True) as status:
        st.write("🖴 检查磁盘空间…")
        disk_err = preflight.check_disk(WORKSPACE, species)
        if disk_err:
            errors.append(disk_err)
        else:
            _free, _where = preflight.disk_free_gb(WORKSPACE)
            st.write(f"🖴 磁盘空间充足（{_where}剩余 {_free:.0f} GB）")

        st.write("🧰 检查分析环境…")
        env_ok = True
        for r in env_check.check_all():
            # pydiffexpress 回退模式不要求 R/DESeq2，避免装了旧引擎却开不了跑
            if diffexp_tool != "deseq2" and r["name"] == "R/DESeq2":
                continue
            if not r["ok"]:
                env_ok = False
                hint = r.get("hint") or "请回『环境体检』页按提示安装"
                errors.append(f"分析环境 {r['name']} 未就绪：{hint}")
        if env_ok:
            st.write("🧰 分析环境正常")

        files = [UPLOAD_DIR / s["r1"] for s in samples]
        files += [UPLOAD_DIR / s["r2"] for s in samples if s.get("r2")]
        bar = st.progress(0.0)

        def cb(p, msg):
            bar.progress(min(max(p, 0.0), 1.0), text=msg)

        fq_errors = preflight.check_fastqs(files, cb)
        errors.extend(fq_errors)
        if not fq_errors:
            st.write(f"🧬 {len(files)} 个测序文件完整无损")

        if errors:
            status.update(label="自检未通过", state="error")
        else:
            status.update(label="自检通过，启动分析", state="complete")

    if errors:
        for e in errors:
            st.error(e)
        return

    # 分组顺序与实际分组必须一致，否则比较方向预告会与实际不符
    # （放在抢占名称之前校验，避免校验失败留下空占位目录）
    order = st.session_state.get("group_order") or []
    if not order or set(order) != set(group_of.values()):
        st.error("分组信息有变化，请回『分组与参数』页重新确认分组。")
        return

    # 名称避让（BUG-10 + BUG-15 修复）：分配不冲突的目录名并原子占位。
    # 只要目录已存在（正在运行/已完成/残留）就换序号，绝不复用旧目录，
    # 否则同名重跑会把旧的 output 文件混进新结果。逻辑在 runner.reserve_run_name。
    base_name = run_name
    run_name = runner.reserve_run_name(base_name, RUNS_DIR)
    if run_name != base_name:
        st.info(f"名称「{base_name}」已使用过，为避免覆盖旧结果，本次自动改用「{run_name}」。")

    try:
        run_dir = RUNS_DIR / run_name
        # pyseqrna 按 Identifier 首次出现顺序生成比较 "X-Y"=X÷Y；
        # 逆序后 vs 对照的比较变成"处理÷对照"（logFC>0=处理组高，符合惯例）。
        rank = {g: i for i, g in enumerate(reversed(order))}
        samples = sorted(samples, key=lambda s: rank[group_of[s["id"]]])

        sample_tsv = run_dir / "samples.tsv"
        sample_sheet.build_sample_sheet(samples, group_of, sample_tsv)
        ini = config_builder.build_ini({
            "sample_sheet": sample_tsv,
            "fastq_dir": UPLOAD_DIR,
            "genome": refs["genome"], "gtf": refs["gtf"],
            "outdir": run_dir / "output",
            "species": species,
            "alignment_tool": alignment_tool,
            "diffexp_tool": diffexp_tool,
            "threads": int(st.session_state.get("threads", config.load_config()["threads"])),
            "memory": int(st.session_state.get("memory", config.load_config()["memory"])),
            "fold_threshold": st.session_state.get("fold", 2.0),
            "pvalue_threshold": st.session_state.get("pval", 0.05),
            "skip_trim": skip_trim,
            "paired": paired,
        })
        ini_path = run_dir / "run.ini"
        ini_path.write_text(ini, encoding="utf-8")
        runner.start_run(ini_path, run_dir, run_dir / "pyseqrna.log",
                         extra_env={"VOLCANO_GENE_LABELS":
                                    "1" if volcano_gene_labels else "0"})
        st.session_state.pop("done_run", None)
    except Exception as e:
        # 启动失败：清掉占位标记和半成品目录，下次可以同名重开
        (RUNS_DIR / run_name / ".active.json").unlink(missing_ok=True)
        shutil.rmtree(RUNS_DIR / run_name, ignore_errors=True)
        st.error(f"启动失败：{e}")
        return
    st.rerun()


def _show_run_progress(run_dir: Path, proc, run_name: str, pid: int | None = None) -> None:
    log_path = run_dir / "pyseqrna.log"
    checkpoint = run_dir / "output" / "pyseqrna_checkpoint.json"
    progress = runner.read_progress(log_path, proc, checkpoint)

    # 阶段只进不退（日志关键词偶尔会误判，导致进度条"倒退"）
    key = f"maxstage_{run_name}"
    prev = st.session_state.get(key, -1)
    stage_idx = max(progress["stage_index"], prev)
    st.session_state[key] = stage_idx

    pct = (stage_idx + 1) / progress["total_stages"] if stage_idx >= 0 else 0.0
    st.progress(min(pct, 1.0))
    head = st.columns([5, 1])
    head[0].write(f"当前阶段：**{progress['stage'] or '初始化'}**"
                  f"（第 {stage_idx + 1 if stage_idx >= 0 else 0}/{progress['total_stages']} 步）")
    if not progress["done"] and pid:
        stopping_key = f"stopping_{run_name}"
        stopping = st.session_state.get(stopping_key, False)
        label = "正在停止…" if stopping else "🛑 停止"
        if head[1].button(label, key=f"stop_{run_name}",
                          disabled=stopping,
                          help="中途停止本次分析（已算完的部分保留）"):
            st.session_state[stopping_key] = True
            runner.stop_run(pid)
            st.warning("已发送停止指令，正在收尾…")
            time.sleep(1)
            st.rerun()
    if progress.get("detail"):
        st.caption(f"最近动作：{progress['detail']}")
    st.text_area("最近日志", progress["tail"], height=200)

    if progress["done"]:
        st.session_state.pop(f"stopping_{run_name}", None)
        st.session_state.pop(key, None)
        runner.clear_active(run_dir)
        if progress["success"]:
            st.success("🎉 分析完成！去『结果下载』页取结果。")
            st.session_state["done_run"] = run_name
        else:
            outdir = RUNS_DIR / run_name / "output"
            if results.find_outputs(outdir):
                note = progress.get("partial_note") or "常见原因是 GO/KEGG 富集需要联网"
                st.warning(
                    f"分析主体已完成，但**部分阶段失败**（{note}）。"
                    "核心结果（差异基因表、火山图、热图、报告）仍可在『结果下载』页获取。")
                st.session_state["done_run"] = run_name
            else:
                rc = progress["returncode"]
                # BUG-05 修复：断线重连场景（RunningPid）拿不到真实退出码，
                # 旧逻辑显示"退出码 None"误导用户；明确区分"未知"与真实码
                rc_text = "未知（进程已不在，常见于网页刷新后重连）" if rc is None else str(rc)
                st.error(f"分析失败（退出码 {rc_text}）。"
                         "可以把下面的日志下载下来，发给懂技术的人排查。")
        if log_path.exists():
            st.download_button("📄 下载运行日志", data=log_path.read_bytes(),
                               file_name=f"{run_name}.log", key=f"log_{run_name}")
            # 诊断包：日志 + 配置 + 样本表 + checkpoint，报错排查一次拿全
            buf = io.BytesIO()
            with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
                zf.write(log_path, "pyseqrna.log")
                for name in ("run.ini", "samples.tsv",
                             "output/pyseqrna_checkpoint.json"):
                    p = run_dir / name
                    if p.exists():
                        zf.write(p, name)
            st.download_button("🧰 下载诊断包（日志+配置+样本表）",
                               data=buf.getvalue(),
                               file_name=f"{run_name}_diagnostics.zip",
                               key=f"diag_{run_name}")
    else:
        st.caption("每 5 秒自动刷新…")
        if st_autorefresh is not None:
            st_autorefresh(interval=5000, key=f"ar_{run_name}")
        else:
            # BUG-14 修复：streamlit-autorefresh 已在 requirements.txt 里，
            # 缺失说明安装不完整——旧逻辑退回 sleep(5)+rerun 轮询，每个看进度的
            # 浏览器会话都占一个工作线程，多人同时看会耗尽线程池。
            # 现在明确提示安装问题，改为用户手动点按钮刷新（不占线程）。
            st.warning("自动刷新组件（streamlit-autorefresh）未安装，进度不会自动更新。\n\n"
                       "修复：在 WSL 里运行 `pip install streamlit-autorefresh`，"
                       "或点下面按钮手动刷新。")
            if st.button("🔄 手动刷新进度", key=f"manualrefresh_{run_name}"):
                st.rerun()


def _read_run_ini(run_dir: Path) -> configparser.ConfigParser:
    """用 configparser 读 run.ini（RED-03：替代逐行手写解析）。"""
    ini_path = run_dir / "run.ini"
    cp = configparser.ConfigParser()
    if ini_path.exists():
        cp.read(ini_path, encoding="utf-8")
    return cp


def _run_species(run_dir: Path) -> str:
    """从 run.ini 读物种代码。"""
    cp = _read_run_ini(run_dir)
    val = cp.get("Species", "species", fallback="")
    return val if val in ("hsapiens", "mmusculus") else "hsapiens"


def _run_feature_file(run_dir: Path) -> Path | None:
    """从 run.ini 读注释文件（GTF）路径。"""
    cp = _read_run_ini(run_dir)
    val = cp.get("General", "feature_file", fallback="")
    return Path(val) if val else None


# RED-01：方向显示名与排序权重统一从 enrich_py 导入，不再各自维护一份
direction_label = enrich_py.direction_label
direction_weight = enrich_py.direction_weight


def _zip_download_section(zp: Path, sig: str, key: str, label: str, file_name: str,
                          fallback_note: str = "") -> None:
    """统一的 ZIP 下载区（v2.1 抽取，消除 app.py 里两处复制粘贴）。

    带签名缓存：内容没变就不重新读字节。按 ZIP 大小分三档处理（修复 BUG-01）：
    - ≤ 下载阈值（180MB）：正常 download_button（必要时缓存字节）
    - > 下载阈值：不再硬塞给浏览器（超 Streamlit 消息上限必然失败），
      复制到工作区 zip_export/ 供用户自取
    """
    zp = Path(zp)
    size = zp.stat().st_size
    bkey = f"{key}_zipcache"
    if size <= _MAX_ZIP_DOWNLOAD_BYTES:
        cached = st.session_state.get(bkey)
        if not cached or cached[0] != sig:
            if size <= _MAX_ZIP_CACHE_BYTES:
                cached = (sig, zp.read_bytes())
                st.session_state[bkey] = cached
            else:
                cached = None  # 能下载但不值得占内存，每次现读
        data = cached[1] if cached else zp.read_bytes()
        st.download_button(label, data=data, file_name=file_name,
                           key=f"{key}_zipbtn", type="primary")
        return

    # 超大包：复制到 zip_export，提示用户自取（BUG-01 修复路径）
    st.session_state.pop(bkey, None)
    export_dir = WORKSPACE / "zip_export"
    export_dir.mkdir(parents=True, exist_ok=True)
    export_path = export_dir / file_name
    try:
        if export_path.exists():
            export_path.unlink()
        shutil.copy2(zp, export_path)
        size_gb = size / 1e9
        st.info(f"📦 {label.replace('📦 ', '')}\n\n"
                f"结果包约 **{size_gb:.1f} GB**，超过网页单次下载上限，已为你复制到：\n\n"
                f"`{export_path}`\n\n"
                "双击『打开上传文件夹.bat』旁的文件夹，或到该路径直接取走；"
                + (fallback_note or ""))
    except Exception as e:
        st.warning(f"结果包过大且复制失败（{e}）。可以手动复制：`{zp}`")


def _enrich_zip_button(enrich_out: Path, key_prefix: str) -> None:
    """一键打包下载全部富集结果（ZIP 内按 比较×上调/下调 分好文件夹）。"""
    try:
        zp, sig = results.zip_folder(enrich_out, enrich_out.parent / "enrich_results.zip",
                                     prefix="GO_KEGG_富集")
    except Exception as e:
        st.warning(f"富集打包失败：{e}")
        return
    _zip_download_section(
        zp, sig, key_prefix,
        "📦 一键下载全部富集结果（ZIP，已按 比较×上调/下调 分好文件夹）",
        f"{enrich_out.parent.name}_富集结果.zip")


def _render_enrich_files(enrich_out: Path, key_prefix: str) -> bool:
    """按富集方法分组展示产物（从磁盘读，刷新后仍在）。

    GSEA：目录结构为 enrich_out/比较名/GSEA_result.csv 等；
    ORA：enrich_out/比较名/up|down/GO_result.csv 等。
    "_" 开头的条目（gseapy 草稿、_stats.json）不展示。有内容返回 True。
    """
    if not enrich_out.exists():
        return False
    files = [p for p in enrich_out.rglob("*") if p.is_file()
             and not any(part.startswith("_") for part in p.relative_to(enrich_out).parts)]
    if not files:
        return False
    _enrich_zip_button(enrich_out, key_prefix)
    stats_map: dict[str, dict] = {}
    skipped: list[str] = []
    sj = enrich_out / "_stats.json"
    if sj.exists():
        try:
            d = json.loads(sj.read_text(encoding="utf-8"))
            stats_map, skipped = d.get("stats", {}), d.get("skipped", [])
        except Exception:
            pass
    if _enrich_out_is_gsea(enrich_out):
        return _render_gsea_files(enrich_out, files, stats_map, skipped, key_prefix)
    return _render_ora_files(enrich_out, files, stats_map, skipped, key_prefix)


def _render_gsea_files(enrich_out: Path, files: list[Path],
                       stats_map: dict[str, dict], skipped: list[str],
                       key_prefix: str) -> bool:
    """GSEA 产物展示：每比较一目录（结果表 + NES 条形图 + 富集曲线图）。"""
    groups: dict[str, dict[str, Path]] = {}
    for p in files:
        rel = p.relative_to(enrich_out)
        parts = rel.parts
        if len(parts) < 2:
            continue
        gkey, fname = parts[0], "/".join(parts[1:])
        groups.setdefault(gkey, {})[fname] = p
    for gkey in sorted(groups):
        name_map = groups[gkey]
        st.markdown(f"**{gkey}**")
        s = stats_map.get(gkey)
        if s:
            top = f"；最显著：{s['top_term']}" if s.get("top_term") else ""
            st.caption(f"{s.get('genes', '?')} 个基因参与排序，匹配基因集 {s.get('matched', '?')} 个；"
                       f"显著通路 {s.get('sig_terms', '?')} 个（fdr<0.25）{top}")
        csv = name_map.get("GSEA_result.csv")
        if csv:
            st.download_button("⬇️ 下载 GSEA 结果表 (CSV)", data=csv.read_bytes(),
                               file_name=f"{gkey}_GSEA_result.csv",
                               key=f"{key_prefix}_{gkey}_gsea_csv")
        bar = name_map.get("GSEA_NES_barplot.png")
        if bar:
            st.image(str(bar), caption="Top NES 条形图", use_container_width=True)
            st.download_button("⬇️ 下载条形图", data=bar.read_bytes(),
                               file_name=f"{gkey}_GSEA_NES_barplot.png",
                               key=f"{key_prefix}_{gkey}_bar")
        curves = [fname for fname in name_map
                  if fname.startswith("GSEA_plots/")
                  and fname.lower().endswith(".png")]  # 只展示 PNG（gseapy 还有 .gmt/.rnk/.log 辅助文件）
        if curves:
            with st.expander(f"📈 显著通路富集曲线（{len(curves)} 张）"):
                for fname in sorted(curves):
                    p = name_map[fname]
                    st.image(str(p), caption=fname.rsplit("/", 1)[-1],
                             use_container_width=True)
                    st.download_button(f"⬇️ 下载 {fname.rsplit('/', 1)[-1]}",
                                       data=p.read_bytes(),
                                       file_name=f"{gkey}_{fname.rsplit('/', 1)[-1]}",
                                       key=f"{key_prefix}_{gkey}_{fname}")
    if skipped:
        st.caption("已跳过：" + "；".join(skipped))
    return True


def _render_ora_files(enrich_out: Path, files: list[Path],
                      stats_map: dict[str, dict], skipped: list[str],
                      key_prefix: str) -> bool:
    """ORA 产物展示：按 比较×方向 分组（GO/KEGG 表 + 气泡图）。"""
    groups: dict[str, dict[str, Path]] = {}
    for p in files:
        rel = p.relative_to(enrich_out)
        gkey = "/".join(rel.parts[:-1])
        groups.setdefault(gkey, {})[p.name] = p

    def _direction_label(parts: list[str]) -> str:
        raw = parts[1] if len(parts) > 1 else ""
        if raw in ("up", "down", "all") and parts:
            return direction_label(parts[0], raw)
        return enrich_py.DIRECTION_CN.get(raw, raw)

    def _sort_key(gk: str):
        parts = gk.split("/")
        direction = parts[-1] if len(parts) > 1 else ""
        return (parts[0] if parts else "", direction_weight(direction))

    for gkey in sorted(groups, key=_sort_key):
        name_map = groups[gkey]
        parts = gkey.split("/")
        if len(parts) == 2:
            label = f"{parts[0]} · {_direction_label(parts)}"
        else:
            label = parts[0] if gkey else ""
        if label:
            st.markdown(f"**{label}**")
        s = stats_map.get(label)
        if s is None and len(parts) == 2 and parts[1] in enrich_py.DIRECTION_CN:
            # 旧版 _stats.json 的键是"比较 · 上调/下调"，显示名已改为自解释格式，回退旧键再查
            s = stats_map.get(f"{parts[0]} · {enrich_py.DIRECTION_CN[parts[1]]}")
        if s:
            st.caption(f"差异基因 {s.get('total', '?')} 个，映射基因名 {s.get('mapped', '?')} 个；"
                       f"GO 匹配 {s.get('matched_go', '?')}、KEGG 匹配 {s.get('matched_kegg', '?')}")
        cols = st.columns(4)
        for col, fn in zip(cols, ("GO_result.csv", "GO_dotplot.png",
                                  "KEGG_result.csv", "KEGG_dotplot.png")):
            if fn in name_map:
                col.download_button(f"下载 {fn}", data=name_map[fn].read_bytes(),
                                    file_name=f"{gkey.replace('/', '_')}_{fn}" if gkey else fn,
                                    key=f"{key_prefix}_{gkey}_{fn}")
    if skipped:
        st.caption("已跳过：" + "；".join(skipped))
    return True


def _enrich_out_is_gsea(enrich_out: Path) -> bool:
    """富集目录里有没有 GSEA 产物（用于展示与重跑逻辑）。"""
    return any(p.name == "GSEA_result.csv" for p in enrich_out.rglob("*"))


def _do_enrich(run_dir: Path, gtf: Path, species: str, method: str) -> None:
    """跑富集并展示结果（供"做富集"和"重新运行富集"两个入口复用）。

    进度条通过 enrich_py 的 progress_cb 实时更新（解析 GTF / 拉基因库 / 逐组富集）。
    """
    method_cn = "GSEA（全基因排序富集）" if method == "gsea" else "ORA（上调/下调分别富集）"
    with st.status(f"正在按比较跑{method_cn}（首次需联网拉取基因集库）…",
                   expanded=True) as status:
        bar = st.progress(0.0)

        def cb(p, msg):
            bar.progress(min(p, 1.0))
            status.update(label=msg)

        produced, stats, skipped = enrich_py.run_enrichment(
            run_dir, gtf, species,
            cache_dir=WORKSPACE / "enrich_cache",
            outdir=run_dir / "enrich_py",
            progress_cb=cb,
            method=method,
        )
        status.update(label="富集完成", state="complete")
    if method == "gsea":
        st.success("GSEA 完成！每个比较按全部基因排序做富集：NES 为正 = 该通路在高表达组富集，"
                   "NES 为负 = 在低表达组富集。")
    else:
        st.success("富集完成！已按每个比较的「上调」「下调」分别出结果：")
    for label, s in stats.items():
        if method == "gsea":
            top = f"；最显著：{s['top_term']}" if s.get("top_term") else ""
            st.caption(f"{label}：{s['genes']} 个基因参与排序，匹配基因集 {s['matched']} 个；"
                       f"显著通路 {s['sig_terms']} 个（fdr<0.25）{top}")
        else:
            st.caption(
                f"{label}：差异基因 {s['total']} 个，映射基因名 {s['mapped']} 个"
                f"（未映射 {s['unmapped']} 个）；GO 匹配 {s['matched_go']}、KEGG 匹配 {s['matched_kegg']}")
    if skipped:
        st.warning("以下分组被跳过：" + "；".join(skipped))
    if species == "mmusculus":
        st.caption("小鼠的 GO 富集用人版基因库做直系同源匹配，解读时注意。")
    for key in sorted(produced):
        st.download_button(f"下载 {key}", data=produced[key].read_bytes(),
                           file_name=key.replace("/", "_").replace(" · ", "_"),
                           key=f"rich_now_{key}_{run_dir.name}")
    _enrich_zip_button(run_dir / "enrich_py", f"rich_nowzip_{run_dir.name}")


def tab_results() -> None:
    st.header("📥 结果下载")
    run_names = sorted([p.name for p in RUNS_DIR.iterdir() if p.is_dir()]) \
        if RUNS_DIR.exists() else []
    if not run_names:
        st.info("还没有任何分析。")
        return

    done_run = st.session_state.get("done_run")
    # BUG-13 修复：用显式 key 管理选择状态。Streamlit 带 key 的 selectbox
    # 只在首次渲染采用 index 默认值，之后用户手动切换不会被重置回 done_run。
    default_idx = run_names.index(done_run) if done_run in run_names else 0
    run_name = st.selectbox("选择一次分析", run_names, index=default_idx,
                            key="selected_run")
    run_dir = RUNS_DIR / run_name
    outdir = run_dir / "output"

    active = runner.find_active_run(RUNS_DIR)
    if active and active[0] == run_name:
        st.info("这次分析还在运行中，目前只能看到部分结果。")

    groups = results.find_outputs(outdir)
    if not groups:
        st.warning("该次分析输出目录为空或尚未完成。")
        return
    total = sum(len(v) for v in groups.values())
    st.write(f"共找到 **{total}** 个结果文件，按类别分组，可单个下载或整体打包：")

    # ---- 打包下载（带缓存，内容没变不重新压缩；富集结果也一起打进 ZIP）----
    # RED-02 修复：ZIP 下载逻辑统一走 _zip_download_section，不再复制粘贴
    enrich_out = run_dir / "enrich_py"
    extra = [(enrich_out, "GO_KEGG_富集")] if enrich_out.exists() else []
    zip_path = run_dir / "results.zip"
    try:
        zp, sig = results.make_zip(outdir, zip_path, extra_dirs=extra)
        _zip_download_section(
            zp, sig, f"zipb_{run_name}",
            "📦 打包下载全部结果 (ZIP，含 GO/KEGG 富集)",
            f"{run_name}_results.zip")
    except Exception as e:
        st.warning(f"打包失败：{e}")

    if (outdir / "4.Normalization" / "VST_normalized_counts.csv").exists():
        st.caption("RPKM 表保留供浏览；本次分析另输出 VST 标准化矩阵"
                   "（VST_normalized_counts.csv），论文级热图/下游分析请用 VST。")

    # ---- 论文级补充图（VST t-SNE / 差异基因 Venn / UpSet），缺啥补啥 ----
    if (outdir / "4.Normalization" / "VST_normalized_counts.csv").exists():
        with st.spinner("生成论文级补充图（t-SNE / Venn / UpSet）…"):
            try:
                plots.ensure_aux_plots(run_dir)
            except Exception:
                pass  # 补充图失败不影响结果页

    # ---- 图片预览：火山图/热图等直接看，不用先下载 ----
    # 排除 pyseqrna 旧版图（已归档，单独分组下载，不进预览）
    pngs = [p for files in groups.values() for p in files
            if p.suffix.lower() == ".png" and not results._is_legacy(p)]
    volcano_dir = outdir / "5.Visualization" / "Volcano"
    if volcano_dir.exists():
        seen = {p.resolve() for p in pngs}
        pngs.extend(p for p in sorted(volcano_dir.glob("*"))
                    if p.is_file() and p.suffix.lower() == ".png"
                    and p.resolve() not in seen)
    pngs.sort(key=lambda p: (0 if "vst" in p.name.lower() else 1, p.name.lower()))
    if pngs:
        with st.expander(f"👀 图片预览（{len(pngs)} 张）", expanded=True):
            cols = st.columns(3)
            for i, p in enumerate(pngs[:12]):
                try:
                    cols[i % 3].image(str(p), caption=p.name, use_container_width=True)
                except Exception as e:
                    st.caption(f"图片 {p.name} 读取失败：{e}")
            if len(pngs) > 12:
                st.caption(f"……共 {len(pngs)} 张，其余请在下方分类下载。")

    for label, files in groups.items():
        with st.expander(f"{label}（{len(files)} 个文件）"):
            for p in files:
                try:
                    data = p.read_bytes()
                except Exception as e:
                    st.caption(f"文件 {p.name} 读取失败：{e}")
                    continue
                st.download_button(
                    f"⬇️ {p.name}",
                    data=data,
                    file_name=p.name,
                    key=f"dl_{run_name}_{p}",
                    use_container_width=True,
                )

    # ---- 磁盘清理 ----
    with st.expander("🧹 清理这次分析的中间文件（释放磁盘）"):
        st.caption("删除 BAM、比对索引、修剪后的 fastq 等中间大文件。"
                   "结果表格和图片都会保留，放心删。")
        running_this = bool(active and active[0] == run_name)
        if st.button("🗑️ 删除中间大文件", key=f"clean_{run_name}", disabled=running_this):
            freed = results.cleanup_intermediates(outdir)
            st.success(f"已释放 {freed / 1e9:.1f} GB 磁盘空间。")
        if running_this:
            st.caption("分析运行中，结束后才能操作")

    st.divider()
    st.subheader("🧬 GO/KEGG 富集")
    species = _run_species(run_dir)
    gtf = _run_feature_file(run_dir)

    # 富集方法：GSEA 需要 DESeq2 引擎的全基因差异表；否则自动回退 ORA
    diff_dir = outdir / "4.Differential_Expression"
    has_deseq = diff_dir.exists() and any(diff_dir.glob("DESeq2_*_vs_*.csv"))
    method_label = st.selectbox(
        "富集方法",
        ["GSEA（推荐：全基因排序富集，论文级）", "ORA（经典：上调/下调基因分别富集）"],
        index=0 if has_deseq else 1,
        key=f"enrich_method_{run_name}",
        help="GSEA 用全部基因按表达差异排序做富集（输出 NES/fdr），是当前论文主流方法；"
             "ORA 只对显著差异基因做经典过表达分析。",
    )
    enrich_method = "gsea" if method_label.startswith("GSEA") else "ora"
    if enrich_method == "gsea" and not has_deseq:
        st.warning("该次分析没有 DESeq2 差异表（旧结果或用了 pydiffexpress 引擎），"
                   "GSEA 不可用，已自动切到 ORA。")
        enrich_method = "ora"
    if enrich_method == "gsea":
        st.caption("GSEA 方向说明：NES 为正 = 通路在高表达组富集，为负 = 在低表达组富集；"
                   "显著阈值 fdr<0.25（论文通用标准）。")
    else:
        st.caption("ORA 方向说明：比较「X vs Y」= X 组相对于 Y 组；上调 = X 组更高。")

    has_results = _render_enrich_files(enrich_out, key_prefix=f"rich_{run_name}")
    if has_results:
        if gtf and gtf.exists():
            if st.button("🔁 重新运行富集", key=f"rerun_enrich_{run_name}",
                         disabled=running_this):
                shutil.rmtree(enrich_out, ignore_errors=True)
                try:
                    _do_enrich(run_dir, gtf, species, enrich_method)
                except Exception as e:
                    st.error(f"富集失败：{e}")
            if running_this:
                st.caption("分析运行中，结束后才能操作")
    else:
        if enrich_method == "gsea":
            st.caption("读取该次分析的 DESeq2 全基因差异表，按每个比较做 GSEA"
                       "（首次运行从 Enrichr 下载基因集库，缓存后离线）。")
        else:
            st.caption("读取该次分析的差异基因名单，按每个比较的「上调」「下调」分别跑 GO/KEGG 富集。"
                       "首次运行从 Enrichr 下载基因集库（需联网，缓存后离线）。")
        if not gtf or not gtf.exists():
            st.warning("找不到该次分析的注释文件（GTF），无法做基因名映射。")
        elif st.button("🚀 做 GO/KEGG 富集", type="primary", key=f"run_enrich_{run_name}"):
            try:
                _do_enrich(run_dir, gtf, species, enrich_method)
            except Exception as e:
                st.error(f"富集失败：{e}")


def main() -> None:
    _inject_css()
    _sidebar()
    _banner()
    _stepper()
    t1, t2, t3, t4, t5 = st.tabs(
        ["🔍 环境体检", "📤 数据与参考文件", "👥 分组与参数", "▶️ 运行分析", "📥 结果下载"])
    with t1:
        tab_env()
    with t2:
        tab_upload()
    with t3:
        tab_groups()
    with t4:
        tab_run()
    with t5:
        tab_results()


if __name__ == "__main__":
    main()
