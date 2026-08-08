"""环境检查：确认分析所需工具可用。

注意：STAR 的可执行文件名是大写 `STAR`（小写 `star` 无软链），
因此 bin 字段单独标注实际可执行名。
直接按绝对路径调用 python 时 PATH 缺 env/bin：把解释器 bin 目录前插后的
环境副本通过 subprocess 的 env 参数传入（runner.env_with_bindir 提供），
不修改全局 os.environ——多个浏览器会话并发体检时不会互相踩 PATH。
"""
from __future__ import annotations
import shutil
import subprocess

from .runner import env_with_bindir

TOOLS = [
    {"name": "pyseqrna", "bin": "pyseqrna", "cmd": "pyseqrna --version",
     "hint": "在 WSL 里 conda activate pyseqrna && pip install -e ~/pyseqrna"},
    {"name": "hisat2", "bin": "hisat2", "cmd": "hisat2 --version",
     "hint": "conda install -c conda-forge -c bioconda hisat2"},
    {"name": "star", "bin": "STAR", "cmd": "STAR --version",
     "hint": "conda install -c conda-forge -c bioconda star"},
    {"name": "samtools", "bin": "samtools", "cmd": "samtools --version",
     "hint": "conda install -c conda-forge -c bioconda samtools"},
    {"name": "fastqc", "bin": "fastqc", "cmd": "fastqc --version",
     "hint": "conda install -c conda-forge -c bioconda fastqc"},
    {"name": "trim_galore", "bin": "trim_galore", "cmd": "trim_galore --version",
     "hint": "conda install -c conda-forge -c bioconda trim-galore"},
]


def check_tool(tool: dict) -> dict:
    name = tool["name"]
    # which 也走扩展后的 PATH，否则解释器 bin 目录里的工具查不到
    if shutil.which(tool["bin"], path=env_with_bindir()["PATH"]) is None:
        return {"name": name, "ok": False, "version": "", "hint": tool["hint"]}
    try:
        out = subprocess.run(
            tool["cmd"].split(), capture_output=True, text=True, timeout=15,
            env=env_with_bindir(),  # PATH 前插 bin 目录，但不碰全局环境
        )
        version = (out.stdout or out.stderr).strip().splitlines()
        if out.returncode != 0:
            return {"name": name, "ok": False,
                    "version": version[0][:80] if version else "执行失败",
                    "hint": tool["hint"]}
        return {"name": name, "ok": True, "version": version[0][:80] if version else "ok",
                "hint": ""}
    except Exception as e:
        return {"name": name, "ok": False,
                "version": f"(无法读取版本: {e})", "hint": tool["hint"]}


def check_all() -> list[dict]:
    return [check_tool(t) for t in TOOLS]
