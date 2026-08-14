"""runner.py 的日志解析与进度判定单元测试。"""
from __future__ import annotations
from pathlib import Path

import pytest

from lib import runner


# ---------------------------------------------------------------- stage_of
def test_stage_of_quality():
    assert runner.stage_of("Running FastQC on sample_1") == "quality"
    assert runner.stage_of("Trimming adapters") == "quality"


def test_stage_of_alignment():
    assert runner.stage_of("Building HISAT2 index") == "alignment"
    assert runner.stage_of("STAR alignment in progress") == "alignment"


def test_stage_of_diffexp():
    assert runner.stage_of("Differential expression analysis") == "diffexp"


def test_stage_of_none():
    assert runner.stage_of("some unrelated line") is None


# ---------------------------------------------------------------- env_with_bindir
def test_env_with_bindir_prepends_and_does_not_touch_global(monkeypatch):
    import os, sys
    original_path = os.environ.get("PATH", "")
    env = runner.env_with_bindir()
    bindir = str(Path(sys.executable).resolve().parent)
    assert env["PATH"].startswith(bindir)
    # 关键：全局 PATH 未被修改（BUG-07 回归测试）
    assert os.environ.get("PATH", "") == original_path


# ---------------------------------------------------------------- _diffexp_tool_from_ini
def test_diffexp_tool_from_ini(tmp_path: Path):
    ini = tmp_path / "run.ini"
    ini.write_text("[DifferentialExpression]\ndiffexp_tool = deseq2\n",
                   encoding="utf-8")
    assert runner._diffexp_tool_from_ini(ini) == "deseq2"
    ini.write_text("[DifferentialExpression]\ndiffexp_tool = pydiffexpress\n",
                   encoding="utf-8")
    assert runner._diffexp_tool_from_ini(ini) == "pydiffexpress"


def test_diffexp_tool_from_ini_missing_or_invalid(tmp_path: Path):
    assert runner._diffexp_tool_from_ini(tmp_path / "nope.ini") == ""
    bad = tmp_path / "bad.ini"
    bad.write_text("diffexp_tool = deseq2\n", encoding="utf-8")  # 无段名
    assert runner._diffexp_tool_from_ini(bad) == ""


# ---------------------------------------------------------------- _read_log_lines
def test_read_log_lines_missing_file(tmp_path: Path):
    assert runner._read_log_lines(tmp_path / "nope.log") == []


def test_read_log_lines_strips_ansi(tmp_path: Path):
    log = tmp_path / "run.log"
    log.write_bytes("\x1b[32mhello\x1b[0m\nworld\n".encode())
    lines = runner._read_log_lines(log)
    assert lines == ["hello", "world"]


# ---------------------------------------------------------------- read_progress
def _write(tmp_path: Path, content: str) -> Path:
    log = tmp_path / "run.log"
    log.write_text(content, encoding="utf-8")
    return log


class FakeProc:
    def __init__(self, rc):
        self._rc = rc

    def poll(self):
        return self._rc

    @property
    def returncode(self):
        return self._rc


def test_read_progress_running(tmp_path: Path):
    log = _write(tmp_path, "Starting alignment\n")
    ck = tmp_path / "output" / "checkpoint.json"
    p = runner.read_progress(log, FakeProc(None), ck)
    assert p["done"] is False
    assert p["stage"] == "alignment"


def test_read_progress_success_needs_checkpoint_and_report(tmp_path: Path):
    # 日志说结束、rc=0，但缺 checkpoint → 不算成功
    log = _write(tmp_path, "End of PySeqRNA 1.0.0 Session\n")
    out = tmp_path / "output"
    out.mkdir()
    ck = out / "pyseqrna_checkpoint.json"
    p = runner.read_progress(log, FakeProc(0), ck)
    assert p["done"] is True
    assert p["success"] is False

    # 有 checkpoint 但无 7.Report → 仍不算成功（部分阶段失败，BUG-05 回归）
    ck.write_text("{}", encoding="utf-8")
    p = runner.read_progress(log, FakeProc(0), ck)
    assert p["success"] is False
    assert "报告" in p["partial_note"]

    # checkpoint + 报告齐全 → 成功
    (out / "7.Report").mkdir()
    (out / "7.Report" / "index.html").write_text("<html></html>")
    p = runner.read_progress(log, FakeProc(0), ck)
    assert p["success"] is True
    assert p["returncode"] == 0


def test_skip_report_from_ini(tmp_path: Path):
    """读 run.ini 的 skip_report：读不到/False 时为 False，True 时为 True。"""
    # 没有 run.ini → 默认 False（保持原报告校验）
    assert runner._skip_report_from_ini(tmp_path / "output" / "ck.json") is False
    # 明确的 skip_report
    (tmp_path / "run.ini").write_text("[Report]\nskip_report = True\n",
                                      encoding="utf-8")
    assert runner._skip_report_from_ini(tmp_path / "output" / "ck.json") is True
    (tmp_path / "run.ini").write_text("[Report]\nskip_report = False\n",
                                      encoding="utf-8")
    assert runner._skip_report_from_ini(tmp_path / "output" / "ck.json") is False


def test_read_progress_skip_report_does_not_require_report(tmp_path: Path):
    """方案 B：skip_report=True 时 7.Report 不产出是预期，不再判为失败。"""
    out = tmp_path / "output"
    out.mkdir()
    (tmp_path / "run.ini").write_text("[Report]\nskip_report = True\n",
                                      encoding="utf-8")
    ck = out / "pyseqrna_checkpoint.json"
    ck.write_text("{}", encoding="utf-8")
    log = _write(tmp_path, "End of PySeqRNA 1.0.0 Session\n")
    p = runner.read_progress(log, FakeProc(0), ck)
    assert p["done"] is True
    assert p["success"] is True
    assert p["returncode"] == 0


def test_read_progress_failed_marker(tmp_path: Path):
    log = _write(tmp_path, "Pipeline execution failed\n")
    ck = tmp_path / "checkpoint.json"
    ck.write_text("{}", encoding="utf-8")
    p = runner.read_progress(log, FakeProc(1), ck)
    assert p["done"] is True
    assert p["success"] is False


def test_read_progress_none_proc_uses_log_end(tmp_path: Path):
    """无进程句柄（重连）时靠日志判结束，退出码为 None（显示'未知'）。"""
    log = _write(tmp_path, "End of PySeqRNA Session\n")
    ck = tmp_path / "checkpoint.json"
    p = runner.read_progress(log, None, ck)
    assert p["done"] is True
    assert p["returncode"] is None


# ---------------------------------------------------------------- find_active_run
def test_find_active_run_cleans_dead_marker(tmp_path: Path):
    runs = tmp_path / "runs"
    d = runs / "run1"
    d.mkdir(parents=True)
    # 用一个几乎肯定不存在的 pid
    marker = d / ".active.json"
    marker.write_text('{"pid": 999999999}', encoding="utf-8")
    assert runner.find_active_run(runs) is None
    assert not marker.exists()  # 过期标记被清理


# ---------------------------------------------------------------- reserve_run_name
def test_reserve_run_name_free_name(tmp_path: Path):
    runs = tmp_path / "runs"
    name = runner.reserve_run_name("my_run", runs)
    assert name == "my_run"
    # 占位：目录 + .active.json 已创建
    assert (runs / "my_run").is_dir()
    assert (runs / "my_run" / ".active.json").exists()


def test_reserve_run_name_completed_dir_gets_suffix(tmp_path: Path):
    """BUG-15 回归：已完成的同名分析（目录在、.active.json 已被 clear_active 删）
    也要换序号，不能复用旧目录——否则旧 output 文件会混进新结果。"""
    runs = tmp_path / "runs"
    (runs / "my_run" / "output").mkdir(parents=True)  # 已完成的旧结果，无 .active.json
    name = runner.reserve_run_name("my_run", runs)
    assert name == "my_run_1"
    assert (runs / "my_run_1" / ".active.json").exists()


def test_reserve_run_name_active_dir_gets_suffix(tmp_path: Path):
    """正在运行的同名（有 .active.json）也要换序号（BUG-10 原有语义保留）。"""
    runs = tmp_path / "runs"
    (runs / "my_run").mkdir(parents=True)
    (runs / "my_run" / ".active.json").write_text('{"pid": 1}', encoding="utf-8")
    name = runner.reserve_run_name("my_run", runs)
    assert name == "my_run_1"


def test_reserve_run_name_increments_until_free(tmp_path: Path):
    runs = tmp_path / "runs"
    (runs / "my_run").mkdir(parents=True)
    (runs / "my_run_1").mkdir(parents=True)
    (runs / "my_run_2").mkdir(parents=True)
    name = runner.reserve_run_name("my_run", runs)
    assert name == "my_run_3"
