"""后台运行 pyseqrna + 日志进度解析。"""
from __future__ import annotations
import configparser
import os
import json
import re
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

_ALIGNMENT_KEYWORDS = (r"\balign", r"\bstar\b", "hisat")
_ALIGNMENT_RE = re.compile("|".join(_ALIGNMENT_KEYWORDS), re.IGNORECASE)

STAGE_ORDER = [
    "quality", "alignment", "quantification", "normalization",
    "diffexp", "clustering", "annotation", "report",
]

_STAGE_KEYWORDS = [
    ("quality", ("quality", "trim", "fastqc")),
    ("alignment", _ALIGNMENT_KEYWORDS),
    ("quantification", ("quantif", "count", "featurecount")),
    ("normalization", ("normaliz",)),
    ("diffexp", ("differential", "diffexp", "deseq")),
    ("clustering", ("cluster",)),
    ("annotation", ("annotation", "gene ontology", "kegg", "pathway")),
    ("report", ("report",)),
]

# 每次只读日志末尾这么多字节——跑十几小时后日志可能几百 MB，
# 全量读取会让进度页越跑越卡
_LOG_TAIL_BYTES = 512 * 1024


def stage_of(line: str) -> str | None:
    low = line.lower()
    for stage, kws in _STAGE_KEYWORDS:
        if stage == "alignment":
            if _ALIGNMENT_RE.search(low):
                return stage
        else:
            if any(k in low for k in kws):
                return stage
    return None


def env_with_bindir() -> dict:
    """返回环境变量副本，PATH 前插当前解释器的 bin 目录。

    直接按绝对路径调用 python（不 conda activate）时，PATH 缺 env/bin，
    pyseqrna 内部调用的 fastqc/STAR/trim_galore 等会找不到。
    env_check 与 runner 共用本函数（只构造副本传给 subprocess 的 env 参数，
    不修改全局 os.environ，避免多会话并发时的竞态）。
    """
    env = os.environ.copy()
    bindir = str(Path(sys.executable).resolve().parent)
    env["PATH"] = bindir + os.pathsep + env.get("PATH", "")
    # rpy2 需要显式 R_HOME（conda 的 R 不在默认搜索路径）
    r_home = Path(bindir).parent / "lib" / "R"
    if r_home.is_dir():
        env["R_HOME"] = str(r_home)
    return env


def pyseqrna_executable() -> str:
    """解析 pyseqrna 可执行文件路径，不依赖 PATH（直接调用环境 python 时 PATH 常缺 env/bin）。"""
    candidate = Path(sys.executable).resolve().parent / "pyseqrna"
    if candidate.exists():
        return str(candidate)
    found = shutil.which("pyseqrna")
    if found:
        return found
    return "pyseqrna"


def _diffexp_tool_from_ini(ini_path: Path) -> str:
    """从 run.ini 读差异分析引擎（configparser，与 app.py 的 _read_run_ini 一致）。"""
    cp = configparser.ConfigParser()
    try:
        cp.read(ini_path, encoding="utf-8")
    except (OSError, configparser.Error):
        return ""
    return cp.get("DifferentialExpression", "diffexp_tool", fallback="")


def start_run(ini_path: Path, cwd: Path, log_path: Path,
              extra_env: dict | None = None) -> subprocess.Popen:
    """启动 pyseqrna 后台进程，并写入活动标记便于断线重连。"""
    cwd = Path(cwd)
    marker = cwd / ".active.json"
    if marker.exists():
        try:
            old_pid = json.loads(marker.read_text(encoding="utf-8")).get("pid")
        except Exception:
            old_pid = None
        if old_pid and _pid_alive(int(old_pid)):
            raise RuntimeError(
                f"已有正在运行的分析「{cwd.name}」，请先停止或等它结束")
        marker.unlink(missing_ok=True)  # 旧标记已失效，清掉再继续
    log_path.parent.mkdir(parents=True, exist_ok=True)
    f = open(log_path, "wb", buffering=0)
    try:
        # deseq2 引擎：pyseqrna 不输出 VST 矩阵/热图，用 bash 包装层串 R 后处理
        if _diffexp_tool_from_ini(ini_path) == "deseq2":
            wrapper = Path(__file__).resolve().parent.parent / "r_scripts" / "run_pipeline.sh"
            cmd = ["bash", str(wrapper), "-c", str(ini_path)]
        else:
            cmd = [pyseqrna_executable(), "-c", str(ini_path)]
        env = env_with_bindir()
        if extra_env:
            env.update(extra_env)
        proc = subprocess.Popen(
            cmd,
            cwd=str(cwd),
            env=env,
            stdout=f,
            stderr=subprocess.STDOUT,
            start_new_session=True,  # 独立进程组：停止时可以整组杀掉（含 STAR/hisat2 子进程）
        )
    finally:
        # 子进程已继承该句柄，父进程侧必须显式关闭，
        # 否则多次启动会累积泄漏文件描述符
        f.close()
    mark_active(Path(cwd), proc.pid)
    return proc


def stop_run(pid: int) -> None:
    """停止一次分析：先礼后兵（SIGTERM → 5 秒不走 → SIGKILL），整组杀。

    动手前复核 pid 身份（cmdline 含 pyseqrna 或 run_pipeline.sh）：调用方查找活动运行与
    用户点停止之间存在时间差，极端情况下原进程已死、pid 被操作系统
    复用给无关程序，不复核可能误杀别人的进程组。
    """
    if not _pid_alive(pid):
        return

    def _kill(sig) -> None:
        try:
            os.killpg(os.getpgid(pid), sig)
        except Exception:
            try:
                os.kill(pid, sig)
            except Exception:
                pass

    _kill(signal.SIGTERM)
    deadline = time.time() + 5
    while time.time() < deadline and _pid_alive(pid):
        time.sleep(0.5)
    if _pid_alive(pid):
        _kill(signal.SIGKILL)
        deadline = time.time() + 2
        while time.time() < deadline and _pid_alive(pid):
            time.sleep(0.2)


def mark_active(run_dir: Path, pid: int) -> None:
    """写活动标记（pid），供服务器重启后重连。"""
    Path(run_dir).joinpath(".active.json").write_text(
        f'{{"pid": {pid}}}\n', encoding="utf-8")


def clear_active(run_dir: Path) -> None:
    Path(run_dir).joinpath(".active.json").unlink(missing_ok=True)


def reserve_run_name(base_name: str, runs_dir: Path) -> str:
    """分配一个不与现有 run 冲突的目录名，并原子占位（.active.json）。

    只要目标目录已存在（正在运行 / 已完成 / 残留）就换序号，**绝不复用旧目录**——
    否则同名重跑时旧的 output 文件会残留、混进新结果（BUG-15 修复：名字避让
    旧逻辑只防「正在运行」的同名，靠 .active.json 存在与否判断；但分析完成后
    .active.json 已被 clear_active 删掉，同名重跑就会复用旧目录、污染结果）。

    占位用 os.open 的 O_CREAT|O_EXCL 抢占：两个标签页同时点「开始分析」也不会撞名。
    返回最终 run_name（可能带 _1/_2 后缀）。
    """
    run_name = base_name
    counter = 0
    while True:
        target = Path(runs_dir) / run_name
        if target.exists():
            # 目录已存在（无论是否还在运行）→ 换序号，绝不覆盖旧结果
            counter += 1
            run_name = f"{base_name}_{counter}"
            continue
        target.mkdir(parents=True, exist_ok=True)
        try:
            fd = os.open(str(target / ".active.json"),
                         os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            os.close(fd)  # 占位只需文件存在；真正的 pid 由 start_run 写入
            return run_name
        except FileExistsError:
            # 并发抢占失败（另一个会话刚建了同名目录）→ 换序号重试
            counter += 1
            run_name = f"{base_name}_{counter}"


def _pid_alive(pid: int) -> bool:
    """pid 是否存活，且确实是 pyseqrna 进程。

    只看 pid 数字有一个坑：分析异常死掉后，操作系统可能把这个号码分给别的
    程序，网页就会永远"连接"到一个不相干的进程上（僵尸运行）。所以进一步
    看 /proc/<pid>/cmdline 里有没有 pyseqrna；deseq2 引擎外层是
    run_pipeline.sh，也需要认它，否则活动标记会被误判成死进程。
    """
    try:
        os.kill(pid, 0)
    except (OSError, ProcessLookupError):
        return False
    try:
        content = Path(f"/proc/{pid}/cmdline").read_bytes().replace(b"\x00", b" ").decode(
            errors="replace")
        return "pyseqrna" in content or "run_pipeline.sh" in content
    except OSError:
        return True  # 读不到 /proc（非 Linux）时退化为只查活


class RunningPid:
    """轻量"进程句柄"：只记录 pid，用于断线重连场景。

    与真正的 Popen 不同：进程消失时拿不到退出码（我们没有 wait 过它），
    所以 returncode 如实返回 None（未知），而不是伪装成 0——
    否则失败提示会显示"退出码 0"误导用户。
    """

    def __init__(self, pid: int):
        self.pid = int(pid)
        self.exited = False
        self._rc: int | None = None

    def poll(self):
        if not self.exited and not _pid_alive(self.pid):
            self.exited = True
        if self.exited:
            # 返回哨兵值表示"已退出"（语义同 Popen.poll 的非 None），
            # 真实退出码未知，由 returncode 属性如实给出 None
            return self._rc if self._rc is not None else -1
        return None

    @property
    def returncode(self):
        return self._rc


def find_active_run(runs_dir: Path) -> tuple[str, int] | None:
    """扫描 runs 目录，返回第一个 pid 仍存活的活动运行 (run_name, pid)。"""
    runs_dir = Path(runs_dir)
    if not runs_dir.exists():
        return None
    first_active = None
    for d in sorted(runs_dir.iterdir()):
        marker = d / ".active.json"
        if marker.exists():
            try:
                pid = json.loads(marker.read_text(encoding="utf-8")).get("pid")
            except Exception:
                continue
            if pid and _pid_alive(pid):
                if first_active is None:
                    first_active = (d.name, int(pid))
            else:
                marker.unlink(missing_ok=True)  # pid 已死，清掉过期标记
    return first_active


_ANSI_RE = re.compile(r"\x1b\[[0-9;]*m")


def _clean(line: str) -> str:
    return _ANSI_RE.sub("", line)


def _read_log_lines(log_path: Path) -> list[str]:
    """只读日志末尾一段（进度判断只需要尾部信息）。"""
    try:
        with open(log_path, "rb") as f:
            f.seek(0, os.SEEK_END)
            size = f.tell()
            f.seek(max(0, size - _LOG_TAIL_BYTES))
            data = f.read()
        text = data.decode("utf-8", errors="replace")
        if size > _LOG_TAIL_BYTES:
            text = text.split("\n", 1)[-1]  # 丢掉被截断的半行
        return [_clean(l) for l in text.splitlines()]
    except OSError:
        return []


def _skip_report_from_ini(checkpoint: Path) -> bool:
    """读 run.ini 的 [Report] skip_report 配置。

    checkpoint 在 <run>/output/ 下，run.ini 在 <run>/ 下（父目录的父目录）。
    跳过 report 阶段（skip_report=True）时本就不产出 7.Report 目录，
    此时不应再把「无 7.Report」判为失败。读不到 run.ini 时按 False 处理，
    保持原有的报告校验行为（未跳过却缺报告 = 阶段中途失败）。
    """
    ini_path = Path(checkpoint).parent.parent / "run.ini"
    cp = configparser.ConfigParser()
    try:
        cp.read(ini_path, encoding="utf-8")
    except (OSError, configparser.Error):
        return False
    return cp.getboolean("Report", "skip_report", fallback=False)


def read_progress(log_path: Path, proc, checkpoint: Path) -> dict:
    """从日志与进程状态汇总进度。proc 可为 None（重连场景/单测）。

    returncode 语义：真 Popen 给出真实退出码；重连用的 RunningPid
    退出码未知时为 None（界面应显示"未知"而非 0）。
    """
    lines = _read_log_lines(Path(log_path))
    tail = "\n".join(lines[-20:])
    detail = ""
    stage = None
    for line in reversed(lines):
        if line.strip():
            detail = line.strip()
            break
    for line in reversed(lines):
        s = stage_of(line)
        if s:
            stage = s
            break

    # 横幅实际带版本号（"End of PySeqRNA 1.0.0 Session"），整串匹配会漏
    log_ended = any("End of PySeqRNA" in l and "Session" in l for l in lines)
    failed_log = any("Pipeline execution failed" in l for l in lines)
    # DESeq2 模式下 wrapper 在 pyseqrna 结束后还要跑 R/VST 后处理，
    # 仅靠 pyseqrna 的结束横幅判完成会把后处理阶段漏掉
    wrapper_mode = any("DESeq2 后处理开始" in l for l in lines)
    wrapper_finished = any("Wrapper pipeline finished" in l for l in lines)

    done = False
    rc = None
    if proc is not None and proc.poll() is not None:
        rc = proc.returncode  # RunningPid 重连场景下为 None（未知）
        done = True
    if not done and log_ended and (not wrapper_mode or wrapper_finished):
        done = True  # 无进程句柄时靠日志判结束；成败用失败标记 + checkpoint

    success = False
    partial_note = ""
    if done:
        ck_exists = Path(checkpoint).exists()
        if rc is not None:
            success = (rc == 0) and (not failed_log) and ck_exists
        else:
            success = (not failed_log) and ck_exists

        # report 是收尾阶段，但 pyseqrna 从不把它写进 checkpoint（已核对源码），
        # 所以直接验证报告产物：没有 7.Report 说明有阶段中途失败（如 KEGG 断网）
        # 却仍正常退出，旧逻辑会把这种运行误报为 🎉成功。
        # 方案 B：默认跳过功能注释与报告（上游 GO/KEGG 服务已失效），
        # 此时 7.Report 不产出是预期行为，不再判为失败。
        if success and not _skip_report_from_ini(checkpoint):
            report_dir = Path(checkpoint).parent / "7.Report"
            try:
                if not report_dir.is_dir() or not any(report_dir.iterdir()):
                    success = False
                    partial_note = "报告未生成（常见原因是 GO/KEGG 富集阶段需要联网）"
            except OSError:
                pass  # 读不了报告目录时跳过，不误伤

    idx = STAGE_ORDER.index(stage) if stage in STAGE_ORDER else -1
    total = len(STAGE_ORDER)
    return {
        "done": done,
        "success": success,
        "partial_note": partial_note,
        "stage": stage,
        "stage_index": idx,
        "total_stages": total,
        "tail": tail,
        "detail": detail,
        "returncode": rc,
    }
