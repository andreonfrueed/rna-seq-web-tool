"""安装后环境自检：工具链 + 迷你样例完整性 + 核心逻辑测试。"""
from __future__ import annotations

import gzip
import subprocess
import sys
from pathlib import Path

from lib import env_check


def _check_tools() -> bool:
    print("== [1/3] 分析工具检查 ==")
    ok = True
    for r in env_check.check_all():
        mark = "OK" if r["ok"] else "MISSING/BROKEN"
        print(f"  [{mark}] {r['name']}: {r['version'] or '-'}")
        ok = ok and r["ok"]
    return ok


def _check_example() -> bool:
    print("== [2/3] 迷你样例完整性 ==")
    up = Path.home() / "rna_web_workspace" / "uploads"
    expect = {f"{s}_{m}.fq.gz" for s in ("CON_1", "CON_2", "TREAT_1", "TREAT_2")
              for m in ("R1", "R2")}
    files = [up / n for n in sorted(expect) if (up / n).exists()]
    have = {f.name for f in files}
    if have != expect:
        print(f"  [FAIL] 上传目录缺样例文件：{sorted(expect - have)}")
        return False
    ok = True
    for f in files:
        try:
            with gzip.open(f, "rb") as gz:
                first = gz.readline()
            valid = first.startswith(b"@")
            print(f"  [{'OK' if valid else 'FAIL'}] {f.name}")
            ok = ok and valid
        except Exception as e:
            print(f"  [FAIL] {f.name}: {e}")
            ok = False
    return ok


def _check_tests() -> bool:
    print("== [3/3] 核心逻辑测试 ==")
    r = subprocess.run([sys.executable, "-m", "pytest", "tests/", "-q"],
                       capture_output=True, text=True)
    tail = (r.stdout or r.stderr).strip().splitlines()[-1]
    print(f"  pytest: {tail}")
    return r.returncode == 0


def main() -> int:
    tools = _check_tools()
    example = _check_example()
    tests = _check_tests()
    if tools and example and tests:
        print("\n环境自检全部通过 ✅ 可以开始分析了。")
        return 0
    print("\n环境自检未通过，请把以上输出截图发给懂技术的人。")
    return 1


if __name__ == "__main__":
    sys.exit(main())
