"""figure_qa.py 的单元测试（出图自检：PNG 有效性 / 缺字拦截 / 报告汇总）。"""
from __future__ import annotations
import json
import struct
import sys
import zlib
from pathlib import Path

import pytest

from lib import figure_qa as fqa


def _chunk(typ: bytes, data: bytes) -> bytes:
    c = struct.pack(">I", len(data)) + typ + data
    return c + struct.pack(">I", zlib.crc32(typ + data) & 0xFFFFFFFF)


def _make_png(path: Path, width: int = 8, height: int = 8) -> Path:
    """构造最小合法 PNG（RGB，无滤波）。"""
    ihdr = struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0)
    row = b"\x00" + b"\x80\x80\x80" * width
    png = (b"\x89PNG\r\n\x1a\n" + _chunk(b"IHDR", ihdr)
           + _chunk(b"IDAT", zlib.compress(row * height)) + _chunk(b"IEND", b""))
    path.write_bytes(png)
    return path


# ---------------------------------------------------------------- audit_png
def test_audit_png_valid(tmp_path: Path):
    p = _make_png(tmp_path / "ok.png", 50, 30)
    issues = fqa.audit_png(p)
    assert all(s != "FAIL" for s, _ in issues), issues


def test_audit_png_missing(tmp_path: Path):
    issues = fqa.audit_png(tmp_path / "nope.png")
    assert issues and issues[0][0] == "FAIL"


def test_audit_png_empty(tmp_path: Path):
    p = tmp_path / "empty.png"
    p.write_bytes(b"")
    issues = fqa.audit_png(p)
    assert issues[0][0] == "FAIL"


def test_audit_png_wrong_header(tmp_path: Path):
    p = tmp_path / "fake.png"
    p.write_bytes(b"not a png at all" * 10)
    issues = fqa.audit_png(p)
    assert issues[0][0] == "FAIL"


def test_audit_png_zero_size_dimensions(tmp_path: Path):
    """PNG 头合法但宽高为 0 → FAIL。"""
    ihdr = struct.pack(">IIBBBBB", 0, 0, 8, 2, 0, 0, 0)
    png = (b"\x89PNG\r\n\x1a\n" + _chunk(b"IHDR", ihdr)
           + _chunk(b"IEND", b""))
    p = tmp_path / "zero.png"
    p.write_bytes(png)
    issues = fqa.audit_png(p)
    assert issues[0][0] == "FAIL"


# ---------------------------------------------------------------- audit_layout（无 matplotlib 降级）
def test_audit_layout_degrades_without_matplotlib(monkeypatch):
    """matplotlib 不可用时 audit_layout 返回空（不崩溃）。"""
    monkeypatch.setitem(sys.modules, "matplotlib", None)

    class FakeFig:
        pass

    assert fqa.audit_layout(FakeFig()) == []


def test_audit_layout_detects_cross_axis_overlap():
    """回归（BUG-28）：左面板右缘贴近右面板左缘时，纵轴刻度数字压到左面板，
    audit_layout 必须报「跨轴标签重叠」。"""
    mpl = pytest.importorskip("matplotlib")
    mpl.use("Agg")
    import matplotlib.pyplot as plt

    # 复刻旧 UpSet 布局：左面板右缘 0.24，右面板左缘 0.28，缝隙仅 0.04
    fig = plt.figure(figsize=(8, 6.5))
    ax_bar = fig.add_axes([0.28, 0.52, 0.64, 0.36])
    ax_size = fig.add_axes([0.06, 0.12, 0.18, 0.72])
    ax_bar.bar([0, 1, 2], [2000, 1500, 1000])
    ax_bar.set_ylabel("Intersection size")
    ax_bar.set_yticks([0, 500, 1000, 1500, 2000, 2500])  # 强制 4 位刻度
    ax_size.barh([0, 1, 2], [3500, 3000, 1000])
    try:
        issues = fqa.audit_layout(fig)
    finally:
        plt.close(fig)
    assert any("跨轴标签重叠" in msg for _sev, msg in issues), issues


def test_audit_layout_clean_figure_no_cross_axis_overlap():
    """回归（BUG-28 配套）：修好布局后，自检不应误报跨轴重叠。"""
    mpl = pytest.importorskip("matplotlib")
    mpl.use("Agg")
    import matplotlib.pyplot as plt

    # 新 UpSet 布局：左面板右缘 0.18，右面板左缘 0.32，缝隙 0.14
    fig = plt.figure(figsize=(8, 6.5))
    ax_bar = fig.add_axes([0.32, 0.52, 0.60, 0.36])
    ax_size = fig.add_axes([0.04, 0.12, 0.14, 0.72])
    ax_bar.bar([0, 1, 2], [2000, 1500, 1000])
    ax_bar.set_ylabel("Intersection size")
    ax_bar.set_yticks([0, 500, 1000, 1500, 2000, 2500])
    ax_size.barh([0, 1, 2], [3500, 3000, 1000])
    try:
        issues = fqa.audit_layout(fig)
    finally:
        plt.close(fig)
    assert not any("跨轴标签重叠" in msg for _sev, msg in issues), issues


# ---------------------------------------------------------------- worst_severity / save_report
def test_worst_severity():
    assert fqa.worst_severity([]) == "CLEAN"
    assert fqa.worst_severity([("WARN", "a")]) == "WARN"
    assert fqa.worst_severity([("WARN", "a"), ("FAIL", "b")]) == "FAIL"


def test_save_report(tmp_path: Path):
    out = tmp_path / "_figure_qa.json"
    fqa.save_report({
        "fig1.png": [("WARN", "偏小"), ("FAIL", "缺字")],
        "fig2.png": [],
    }, out)
    data = json.loads(out.read_text(encoding="utf-8"))
    assert data["fig1.png"]["severity"] == "FAIL"
    assert data["fig1.png"]["issues"][0]["msg"] == "偏小"
    assert data["fig2.png"]["severity"] == "CLEAN"


def test_audit_figure_saves_and_audits_png(tmp_path: Path):
    """audit_figure：布局检查 + savefig + PNG 有效性一条龙。"""
    import types

    calls = {"saved": False}
    fake_renderer = types.SimpleNamespace()

    class FakeCanvas:
        def draw(self):
            pass

        def get_renderer(self):
            return fake_renderer

    class FakeAx:
        texts = []
        title = types.SimpleNamespace(get_text=lambda: "", get_window_extent=lambda **k: None)
        xaxis = types.SimpleNamespace(label=types.SimpleNamespace(get_text=lambda: "", get_window_extent=lambda **k: None))
        yaxis = types.SimpleNamespace(label=types.SimpleNamespace(get_text=lambda: "", get_window_extent=lambda **k: None))

        def get_xticklabels(self):
            return []

    class FakeFig:
        bbox = types.SimpleNamespace(x0=0, x1=100, y0=0, y1=100)
        axes = [FakeAx()]
        canvas = FakeCanvas()

        def savefig(self, path, **kw):
            calls["saved"] = True
            _make_png(Path(path), 40, 30)

        def close(self):
            pass

    png_path = tmp_path / "out.png"
    issues = fqa.audit_figure(FakeFig(), png_path)
    assert calls["saved"] is True
    assert png_path.exists()
    assert all(s != "FAIL" for s, _ in issues), issues
