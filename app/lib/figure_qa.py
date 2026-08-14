"""出图自检模块（借鉴 scipilot-figure-skill 的视觉自检思想，机器可判部分）。

期刊投稿级图形不能"画完就算"，必须过了确定性检查才能进结果页/论文：

1. **PNG 有效性**（对 R 和 Python 出的图都适用）：
   文件头必须是 \x89PNG、非空、宽高 > 0 —— 截断/损坏的图直接拦截。
2. **缺字拦截**（matplotlib 图）：渲染时同时挂 warnings 与 logging 两条
   告警通道，任一报 "missing from font / Glyph" 即判定成图会出方框/乱码
   （中文期刊最常见的翻车点）。
3. **文字裁切**（matplotlib 图）：Text 的 window_extent 超出画布边界。
4. **刻度标签重叠**（matplotlib 图）：相邻 tick label 的包围盒相交。

severity 约定：INFO < WARN < FAIL。审计结果可汇总成 JSON 报告，
供结果页展示（run 目录下 _figure_qa.json）。

matplotlib 缺失时 2-4 自动降级为空（不崩溃），PNG 有效性检查无依赖。
"""
from __future__ import annotations

import json
import logging
import struct
import warnings
from pathlib import Path

_SEVERITY = {"INFO": 0, "WARN": 1, "FAIL": 2}
_GLYPH_MARKERS = ("missing from", "Glyph", "findfont", "does not contain")


class _GlyphLogHandler(logging.Handler):
    """拦截 matplotlib logger 里关于缺字 / 找不到字体的记录。"""

    def __init__(self) -> None:
        super().__init__()
        self.messages: list[str] = []

    def emit(self, record: logging.LogRecord) -> None:
        msg = record.getMessage()
        if any(m in msg for m in _GLYPH_MARKERS):
            self.messages.append(msg)


def audit_png(path) -> list[tuple[str, str]]:
    """PNG 有效性检查：文件头 / 非空 / 宽高。返回 [(severity, msg), ...]。"""
    p = Path(path)
    issues: list[tuple[str, str]] = []
    if not p.exists():
        issues.append(("FAIL", f"PNG 不存在: {p}"))
        return issues
    size = p.stat().st_size
    if size == 0:
        issues.append(("FAIL", f"PNG 为空文件: {p}"))
        return issues
    try:
        head = p.read_bytes()[:24]
        if not head.startswith(b"\x89PNG\r\n\x1a\n"):
            issues.append(("FAIL", f"文件头不是 PNG 格式: {p}"))
            return issues
        # IHDR: 宽高在字节 16-23（大端 uint32 ×2）
        width, height = struct.unpack(">II", head[16:24])
        if width <= 0 or height <= 0:
            issues.append(("FAIL", f"PNG 尺寸非法 {width}x{height}: {p}"))
        elif size < 1000:
            issues.append(("WARN", f"PNG 偏小 ({size} B)，可能内容空洞: {p.name}"))
    except (OSError, struct.error) as e:
        issues.append(("FAIL", f"PNG 读取失败: {p} ({e})"))
    return issues


def audit_layout(fig) -> list[tuple[str, str]]:
    """matplotlib 图布局检查：缺字 / 文字裁切 / 刻度重叠。

    缺字走渲染时的 warnings + logging 双通道拦截；裁切与重叠用
    renderer 测量 Text 包围盒。matplotlib 不可用时返回空列表。
    """
    issues: list[tuple[str, str]] = []
    try:
        import matplotlib
        import matplotlib.text as mtext
        import matplotlib.pyplot as plt
    except ImportError:
        return issues

    # 1) 缺字拦截（渲染一次触发字体解析）
    handler = _GlyphLogHandler()
    mpl_logger = logging.getLogger("matplotlib")
    prev_level = mpl_logger.level
    mpl_logger.setLevel(logging.WARNING)
    mpl_logger.addHandler(handler)
    collected: list[str] = []
    try:
        with warnings.catch_warnings(record=True) as wlist:
            warnings.simplefilter("always")
            try:
                fig.canvas.draw()
            except Exception:
                pass
            for w in wlist:
                if any(m in str(w.message) for m in _GLYPH_MARKERS):
                    collected.append(str(w.message))
    finally:
        mpl_logger.setLevel(prev_level)
        mpl_logger.removeHandler(handler)
    for msg in collected[:3]:
        issues.append(("FAIL", f"缺字/乱码风险: {msg[:120]}"))
    for msg in handler.messages[:3]:
        issues.append(("FAIL", f"缺字/乱码风险: {msg[:120]}"))

    # 2) 文字裁切 + 3) 刻度重叠 + 4) 跨轴重叠（需要 renderer 测量）
    try:
        renderer = fig.canvas.get_renderer()
        bbox = fig.bbox

        def _overlap(b1, b2) -> bool:
            return (b1.x0 < b2.x1 and b2.x0 < b1.x1
                    and b1.y0 < b2.y1 and b2.y0 < b1.y1)

        for ax in fig.axes:
            for txt in ax.texts + [ax.title, ax.xaxis.label, ax.yaxis.label]:
                if not txt.get_text():
                    continue
                try:
                    ext = txt.get_window_extent(renderer=renderer)
                    if ext.x0 < bbox.x0 - 1 or ext.x1 > bbox.x1 + 1 \
                            or ext.y0 < bbox.y0 - 1 or ext.y1 > bbox.y1 + 1:
                        issues.append(("WARN", f"文字裁切风险: '{txt.get_text()[:20]}'"))
                except Exception:
                    continue
            # 相邻刻度标签重叠（x 轴）
            try:
                labels = [t for t in ax.get_xticklabels() if t.get_text()]
                prev_ext = None
                for t in labels:
                    ext = t.get_window_extent(renderer=renderer)
                    if prev_ext is not None and ext.x0 < prev_ext.x1:
                        issues.append(("WARN", "x 轴刻度标签重叠"))
                        break
                    prev_ext = ext
            except Exception:
                pass
            # 相邻刻度标签重叠（y 轴）——按竖直位置排序后比较（倒置轴也正确）
            try:
                exts = []
                for t in ax.get_yticklabels():
                    if not t.get_text():
                        continue
                    try:
                        exts.append(t.get_window_extent(renderer=renderer))
                    except Exception:
                        pass
                exts.sort(key=lambda e: e.y0)
                prev = None
                for ext in exts:
                    if prev is not None and ext.y0 < prev.y1:
                        issues.append(("WARN", "y 轴刻度标签重叠"))
                        break
                    prev = ext
            except Exception:
                pass

        # 跨轴重叠：某轴的刻度/轴标签落到另一轴的绘图区里
        # （UpSet 纵轴刻度数字压到左侧集合大小面板即此类，BUG-28）
        try:
            axes = fig.axes
            for a in range(len(axes)):
                for b in range(len(axes)):
                    if a == b:
                        continue
                    bb_b = axes[b].get_window_extent(renderer=renderer)
                    for label in list(axes[a].get_xticklabels()) \
                            + list(axes[a].get_yticklabels()) \
                            + [axes[a].xaxis.label, axes[a].yaxis.label]:
                        if not label.get_text().strip():
                            continue
                        bb = label.get_window_extent(renderer=renderer)
                        if _overlap(bb, bb_b):
                            issues.append(("WARN", "跨轴标签重叠: "
                                           f"'{label.get_text()[:12]}' 压到相邻子图"))
                            break
        except Exception:
            pass
    except Exception:
        pass
    return issues


def audit_figure(fig, png_path) -> list[tuple[str, str]]:
    """组合审计：matplotlib 图布局检查 + 落盘后 PNG 有效性检查。"""
    issues = audit_layout(fig)
    try:
        fig.savefig(png_path, dpi=300, bbox_inches="tight")
    except Exception as e:
        issues.append(("FAIL", f"保存 PNG 失败: {e}"))
        return issues
    try:
        fig.close()  # 关闭失败不影响审计结论
    except Exception:
        pass
    issues.extend(audit_png(png_path))
    return issues


def worst_severity(issues: list[tuple[str, str]]) -> str:
    """取问题列表的最高严重级：FAIL > WARN > INFO > CLEAN。"""
    if not issues:
        return "CLEAN"
    return max((s for s, _ in issues), key=lambda s: _SEVERITY.get(s, 0))


def save_report(issues_map: dict[str, list[tuple[str, str]]], out_path) -> None:
    """审计结果汇总成 JSON（{图名: {severity, issues}}），供结果页展示。"""
    out = {
        name: {
            "severity": worst_severity(issues),
            "issues": [{"severity": s, "msg": m} for s, m in issues],
        }
        for name, issues in sorted(issues_map.items())
    }
    Path(out_path).write_text(
        json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
