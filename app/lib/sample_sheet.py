"""分组 → pyseqrna 5 列样本表（tab 分隔）。

列：SampleName, Replication, Identifier, File1, File2

关于 Replication 列（已对照 pyseqrna 源码 utils/input_processor.py 核实）：
pyseqrna 把 Replication 列直接用作样本唯一 ID（sample_id），而不是统计学
意义上的"重复编号"——填样本名正是正确写法。若误填数字重复号，
反而会让不同样本共享 ID、互相覆盖。列名沿用上游习惯，特此说明。
"""
from __future__ import annotations
from pathlib import Path


def build_sample_sheet(samples, group_of, out_path) -> Path:
    """samples: list[dict]，每项含 id/r1/r2；group_of: {sample_id: group}。

    Replication 列填样本名（= pyseqrna 的样本唯一 ID），
    Identifier 列填分组名（= 差异比较的依据）。
    """
    out_path = Path(out_path)
    lines = ["SampleName\tReplication\tIdentifier\tFile1\tFile2"]
    for s in samples:
        sid = str(s["id"])
        if sid not in group_of:
            raise ValueError(f"样本 {sid} 未分组，无法生成样本表")
        r2 = str(s.get("r2") or "")
        lines.append(f"{sid}\t{sid}\t{group_of[sid]}\t{s['r1']}\t{r2}")
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out_path
