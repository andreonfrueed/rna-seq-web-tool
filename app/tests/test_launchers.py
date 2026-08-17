"""Windows 启动脚本（*.bat）结构检查。

本机/CI 无 WSL 真跑 bat 的环境，用静态结构检查兜底，锁死
「打开上传文件夹.bat」的修复行为，防止将来改回脆弱的 interop 直传模式。
真实端到端验证（Explorer 打开正确目录）在用户本机手动双击完成。
"""
from __future__ import annotations
from pathlib import Path

from lib import config

ROOT = Path(__file__).resolve().parents[2]  # tests -> app -> 仓库根
BAT = ROOT / "打开上传文件夹.bat"


def _bat_text() -> str:
    assert BAT.exists(), f"启动脚本不存在: {BAT}"
    return BAT.read_text(encoding="utf-8")


def test_upload_folder_bat_resolves_wsl_path():
    """用 wslpath 解析真实路径，而非 interop 直传（直传丢参会退回系统文档目录）。"""
    assert "wslpath -w" in _bat_text()


def test_upload_folder_bat_opens_resolved_path_directly():
    """单独用 explorer 打开解析出的变量，不把路径内联进 wsl 命令。"""
    assert 'explorer "%UPLOAD_WIN%"' in _bat_text()


def test_upload_folder_bat_has_failure_guard():
    """WSL 不可用时明确报错退出，不静默打开错误的文件夹。"""
    text = _bat_text()
    assert "if not defined UPLOAD_WIN" in text
    assert "exit /b 1" in text


def test_upload_folder_bat_no_fragile_interop_forwarding():
    """旧 bug 特征（explorer.exe $(wslpath ...) 内联直传）必须绝迹。"""
    assert "explorer.exe $(" not in _bat_text()


def test_upload_folder_bat_matches_workspace_config():
    """bat 写死的工作区路径必须与配置一致，否则又跳回系统「文档」目录。"""
    ws = config.load_config()["workspace_dir"]
    assert ws in _bat_text()
