"""全局配置：工作区路径解析。"""
from __future__ import annotations
import os
from pathlib import Path
from typing import Any

import yaml

_DEFAULT_CONFIG = {
    "workspace_dir": "~/rna_web_workspace",
    "ensembl_base": "https://ftp.ensembl.org/pub",
    "ensembl_release": 113,
    "fold_threshold": 2.0,
    "pvalue_threshold": 0.05,
    "fdr_threshold": 0.05,
    "threads": 8,
    "memory": 32,
}


def _app_dir() -> Path:
    return Path(__file__).resolve().parent.parent


def load_config() -> dict[str, Any]:
    cfg = dict(_DEFAULT_CONFIG)
    p = _app_dir() / "web_config.yaml"
    if p.exists():
        data = yaml.safe_load(p.read_text(encoding="utf-8"))
        if data is None:
            data = {}
        if not isinstance(data, dict):
            raise ValueError(
                "web_config.yaml 内容必须是键值对格式（形如 workspace_dir: ...），请检查文件")
        cfg.update(data)
    return cfg


def workspace_dir() -> Path:
    cfg = load_config()
    override = os.environ.get("RNA_WEB_WORKSPACE")
    raw = override or str(cfg["workspace_dir"])
    p = Path(raw).expanduser()
    if p.exists() and not p.is_dir():
        raise RuntimeError(
            f"工作区路径「{p}」已存在且不是文件夹，请修改 web_config.yaml 或删除该文件")
    p.mkdir(parents=True, exist_ok=True)
    return p
