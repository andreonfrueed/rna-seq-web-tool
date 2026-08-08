"""开跑前自检：磁盘空间 + fastq 文件完整性。

目的：把"跑了几小时才因为磁盘满了/文件坏了而失败"这种事，拦在开跑之前。
"""
from __future__ import annotations
import gzip
import os
import re
import shutil
import subprocess
from pathlib import Path

# 各物种建议的最少剩余磁盘空间（中间文件很占地方：索引 + BAM + 修剪文件）
_MIN_FREE_GB = {"hsapiens": 45.0, "mmusculus": 30.0}

# WSL 虚拟盘所在的 Windows 盘挂载点（如 /mnt/c），探测一次后缓存
_HOST_MOUNT: Path | None = None
_HOST_MOUNT_DONE = False


def _detect_host_mount() -> Path | None:
    """找到 WSL 虚拟盘文件实际所在的 Windows 盘的挂载点（如 /mnt/c）。

    为什么需要：WSL 里直接量工作区路径，得到的是虚拟盘的数字——
    微软默认给虚拟盘画 1TB 的"饼"，所以永远显示还剩约 1TB，永远"够用"，
    和 Windows 真实盘的剩余量完全脱钩。真实的约束是虚拟盘文件（ext4.vhdx）
    放在 Windows 的哪块盘上，那块盘的盘符从 Windows 注册表查。
    """
    me = os.environ.get("WSL_DISTRO_NAME", "").strip().lower()
    if not me:
        return None  # 不在 WSL 里（比如单元测试），无法探测，交给调用方兜底
    try:
        raw = subprocess.run(
            ["reg.exe", "query",
             r"HKCU\Software\Microsoft\Windows\CurrentVersion\Lxss", "/s"],
            capture_output=True, timeout=10,
        ).stdout
    except Exception:
        return None
    # reg.exe 重定向输出常见为 UTF-16LE（带 BOM），多试几种解码兜底
    text = ""
    for enc in ("utf-16", "utf-8", "gbk"):
        try:
            text = raw.decode(enc)
            break
        except (UnicodeDecodeError, ValueError):
            continue
    if not text:
        return None
    # 输出按 HKEY...\{guid} 分块，每块含 DistributionName 和 BasePath
    drive: str | None = None
    for block in re.split(r"(?=^HKEY)", text, flags=re.M):
        name = re.search(r"DistributionName\s+REG_SZ\s+(\S+)", block)
        base = re.search(r"BasePath\s+REG_SZ\s+([A-Za-z]):", block)
        if name and base and name.group(1).strip().lower() == me:
            drive = base.group(1).lower()
            break
    if drive is None:
        # 没匹配到当前发行版：抓第一个 BasePath 的盘符兜底（多半就是那块盘）
        m = re.search(r"BasePath\s+REG_SZ\s+([A-Za-z]):", text)
        drive = m.group(1).lower() if m else None
    if drive:
        mnt = Path(f"/mnt/{drive}")
        if mnt.exists():
            return mnt
    return None


def _host_mount() -> Path | None:
    global _HOST_MOUNT, _HOST_MOUNT_DONE
    if not _HOST_MOUNT_DONE:
        _HOST_MOUNT = _detect_host_mount()
        _HOST_MOUNT_DONE = True
    return _HOST_MOUNT


def disk_free_gb(path: Path) -> tuple[float, str]:
    """返回 (剩余GB, 说明标签如"C 盘")。

    优先报 WSL 虚拟盘实际所在的 Windows 宿主盘；探测失败才退回量
    工作区本身（此时是虚拟盘数字，标签会注明仅供参考）。
    """
    path = Path(path)
    if (path.parts and path.parts[0] in ("/", "\\") and "mnt" in path.parts):
        # 工作区直接放在 /mnt/<盘> 下时，量的是真实 Windows 盘，不再绕注册表
        idx = path.parts.index("mnt")
        drive = path.parts[idx + 1] if idx + 1 < len(path.parts) else None
        if drive:
            mnt = Path("/") / "mnt" / drive
            try:
                return shutil.disk_usage(str(mnt)).free / 1e9, f"{mnt.name.upper()} 盘"
            except OSError:
                pass
    mnt = _host_mount()
    if mnt is not None:
        try:
            return shutil.disk_usage(str(mnt)).free / 1e9, f"{mnt.name.upper()} 盘"
        except OSError:
            pass
    return shutil.disk_usage(str(path)).free / 1e9, "虚拟盘，仅供参考"


def check_disk(path: Path, species: str) -> str | None:
    """磁盘不够返回错误文案，够用返回 None。"""
    free, where = disk_free_gb(path)
    need = _MIN_FREE_GB.get(species, 30.0)
    if free < need:
        return (f"磁盘（{where}）剩余空间只有 {free:.0f}GB，低于建议的 {need:.0f}GB。"
                "分析中途会产生大量中间文件，空间不够会半路失败。"
                "请先清理磁盘（结果页有『清理中间文件』按钮可以释放旧分析的空间）。")
    return None


def check_fastq(path: Path) -> str | None:
    """单个 fastq 文件完整性检查；有问题返回错误文案，没问题返回 None。

    .gz 文件做完整解压试读（能揪出"拷贝到一半"的坏文件）；
    未压缩的 fastq 检查首字符是否为 '@'。
    """
    p = Path(path)
    if not p.exists():
        return f"{p.name} 文件不存在（是不是在上传文件夹里被删掉或改名了？）"
    if p.stat().st_size == 0:
        return f"{p.name} 大小为 0，不是有效的测序文件"
    try:
        if p.suffix.lower() == ".gz":
            with gzip.open(p, "rb") as f:
                while f.read(8 * 1024 * 1024):
                    pass
        else:
            with open(p, "rb") as f:
                if f.read(1) != b"@":
                    return f"{p.name} 不是有效的 fastq 文件（首行应以 @ 开头）"
    except Exception:
        return (f"{p.name} 已损坏（可能拷贝/下载不完整），"
                "请删除后重新上传或重新拷贝到上传文件夹")
    return None


def check_fastqs(paths: list[Path], cb=None) -> list[str]:
    """批量检查，cb(进度0-1, 提示语) 汇报进度；返回错误文案列表。"""
    errors: list[str] = []
    n = max(len(paths), 1)
    for i, p in enumerate(paths):
        if cb:
            cb(i / n, f"检查 {Path(p).name}…")
        err = check_fastq(p)
        if err:
            errors.append(err)
    if cb:
        cb(1.0, "文件检查完成")
    return errors
