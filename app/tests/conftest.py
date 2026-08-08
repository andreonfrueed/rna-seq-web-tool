"""pytest 配置：把 app/ 目录放进 sys.path，测试里直接 `from lib import ...`。

运行（在 app/ 目录或其上级）：
    python -m pytest app/tests/ -q
"""
from __future__ import annotations
import sys
from pathlib import Path

APP_DIR = Path(__file__).resolve().parent.parent
if str(APP_DIR) not in sys.path:
    sys.path.insert(0, str(APP_DIR))
