"""app.py 的静态结构检查（依赖 streamlit 无法 import，用读源码兜底）。

覆盖 app.py 里被抽离/新增、却因无法 import 而没有行为级测试的纯逻辑。
与 test_r_scripts.py 同款「读源码 + 字符串断言」模式。
"""
from __future__ import annotations
from pathlib import Path

APP = Path(__file__).resolve().parent.parent / "app.py"


def _src() -> str:
    assert APP.exists(), f"app.py 不存在: {APP}"
    return APP.read_text(encoding="utf-8")


def test_app_re1_re2_support_lane_suffix():
    r"""BUG-18 回归：_1/_2 命名必须支持 lane 后缀（sample_1_001.fastq.gz）。

    旧正则 `(.+)_1\.(...)` 不认 lane 后缀，`sample_1_001.fastq.gz` 会被当成
    单个样本；修复后应与 _R1/_R2 对齐，含 `(_\d+)?`。
    """
    text = _src()
    assert re_compile_with_lane(text, "_RE_1"), "_RE_1 应含 (_\\d+)? 支持 lane 后缀"
    assert re_compile_with_lane(text, "_RE_2"), "_RE_2 应含 (_\\d+)? 支持 lane 后缀"


def re_compile_with_lane(text: str, var: str) -> bool:
    """变量 var 的正则定义行里是否含 `(_\\d+)?`。"""
    for line in text.splitlines():
        if line.strip().startswith(f"{var} = re.compile("):
            return r"(_\d+)?" in line
    return False


def test_app_uses_reserve_run_name():
    """BUG-15 回归：_start_analysis 用 runner.reserve_run_name，而非手写 while 循环。"""
    text = _src()
    assert "runner.reserve_run_name(" in text
    assert ".active.json" in text  # 占位标记仍在用（reserve_run_name 内部也写，但 app.py 不再手写 os.open 抢占）
