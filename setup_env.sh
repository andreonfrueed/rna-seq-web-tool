#!/bin/bash
# ============================================================
#  RNA 分析网页 一键安装脚本（在 WSL 内运行）
#  由 Claude Code 为实验室师兄师姐制作
#  作用：安装 miniconda + 分析工具 + 网页应用，全自动
# ============================================================
set -e

echo ""
echo "==============================================="
echo "  RNA 分析网页 · 一键安装"
echo "  全程自动，约需 30~60 分钟，请勿关闭窗口"
echo "==============================================="

# 1. 安装 miniconda
if [ ! -x "$HOME/miniconda3/bin/conda" ]; then
  echo ""
  echo "[1/5] 安装 miniconda（基础软件管理）..."
  curl -fsSL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/mc.sh
  bash /tmp/mc.sh -b -p "$HOME/miniconda3"
  rm -f /tmp/mc.sh
else
  echo "[1/5] miniconda 已存在，跳过"
fi

source "$HOME/miniconda3/etc/profile.d/conda.sh"

# 2. 创建分析环境并安装工具
if conda env list 2>/dev/null | grep -q "^pyseqrna "; then
  echo "[2/5] 分析环境已存在，跳过"
else
  echo "[2/5] 创建分析环境并安装比对工具（内存友好的 HISAT2 + STAR + 辅助工具）..."
  if ! conda create -y -n pyseqrna -c conda-forge -c bioconda python=3.10 \
      hisat2 star samtools fastqc trim-galore; then
    echo "     尝试分开安装..."
    conda create -y -n pyseqrna -c conda-forge python=3.10
    conda install -y -n pyseqrna -c conda-forge -c bioconda \
      hisat2 star samtools fastqc trim-galore
  fi
fi

conda activate pyseqrna

# 3. 安装 pyseqrna 本体，并记录版本号（方便日后排查"两台机器结果不一样"）
# 安全修复（SEC-01）：锁定到经过审计的 release tag v1.0.0，并在安装前校验
# commit 哈希——上游仓库若被投毒或 tag 被挪动，安装会直接失败而不是静默中招。
# 将来升级引擎：改这两行 + 人工审计新代码即可。
PYSEQRNA_TAG="v1.0.0"
PYSEQRNA_COMMIT="0006e8af51c0a940bc8f37ade65c3e3a23c79029"
echo "[3/5] 安装分析引擎 pyseqrna（锁定 ${PYSEQRNA_TAG}）..."
if [ ! -d "$HOME/pyseqrna/.git" ]; then
  rm -rf "$HOME/pyseqrna"
  git clone --depth 1 --branch "$PYSEQRNA_TAG" \
    https://github.com/navduhan/pyseqrna.git "$HOME/pyseqrna"
fi
# 校验 commit：与锁定值不一致就重新克隆；仍不一致则拒绝安装
actual_commit="$(git -C "$HOME/pyseqrna" rev-parse HEAD 2>/dev/null || echo '')"
if [ "$actual_commit" != "$PYSEQRNA_COMMIT" ]; then
  echo "     已存在的 pyseqrna 副本 commit 不符（期望 ${PYSEQRNA_COMMIT:0:12}，实际 ${actual_commit:0:12}），重新克隆锁定版本..."
  rm -rf "$HOME/pyseqrna"
  git clone --depth 1 --branch "$PYSEQRNA_TAG" \
    https://github.com/navduhan/pyseqrna.git "$HOME/pyseqrna"
  actual_commit="$(git -C "$HOME/pyseqrna" rev-parse HEAD 2>/dev/null || echo '')"
fi
if [ "$actual_commit" != "$PYSEQRNA_COMMIT" ]; then
  echo "     错误：pyseqrna 版本校验失败（可能被篡改或网络问题），已中止安装"
  exit 1
fi
pip install -e "$HOME/pyseqrna"

# 4. 复制网页应用，并按 requirements.txt 安装网页依赖（首次安装后用 requirements-lock.txt 快照固定）
echo "[4/5] 复制网页应用并安装依赖..."
APP_SRC="$(cd "$(dirname "$0")" && pwd)/app"
if [ -d "$APP_SRC" ]; then
  rm -rf "$HOME/rna_web_app"
  cp -r "$APP_SRC" "$HOME/rna_web_app"
  echo "     已复制到 ~/rna_web_app"
else
  echo "     错误：没找到 app 文件夹（脚本应在『工具』文件夹里运行）"
  exit 1
fi
# 依赖版本锁定（RED-05 修复）：优先用随分发包携带的 requirements-lock.txt
# （每台机器装出完全相同的版本组合）；只有分发包里没有 lock 文件时
# （首次制作分发包的场景）才按 requirements.txt 安装并生成新 lock。
# 旧逻辑只生成、从不使用 lock，版本锁定形同虚设。
if [ -f "$HOME/rna_web_app/requirements-lock.txt" ]; then
  echo "     使用 requirements-lock.txt 安装锁定版本的依赖..."
  pip install -r "$HOME/rna_web_app/requirements-lock.txt"
else
  echo "     分发包未携带 lock 文件，按 requirements.txt 安装并生成新 lock..."
  pip install -r "$HOME/rna_web_app/requirements.txt"
  pip freeze > "$HOME/rna_web_app/requirements-lock.txt" 2>/dev/null || true
fi
git -C "$HOME/pyseqrna" rev-parse HEAD > "$HOME/rna_web_app/PYSEQRNA_COMMIT.txt" 2>/dev/null || true

# 5. 预下载富集基因集库（GO/KEGG，让富集功能离线可用）
echo "[5/5] 预下载富集基因集库（GO/KEGG，约 1-2 分钟）..."
cd "$HOME/rna_web_app" && python - << 'PY'
import sys
from pathlib import Path
sys.path.insert(0, str(Path(".").resolve()))
from lib import enrich_py
cache = Path.home() / "rna_web_workspace/enrich_cache"
cache.mkdir(parents=True, exist_ok=True)
for sp in ("hsapiens", "mmusculus"):
    for name in enrich_py.GO_LIBS + [enrich_py._kegg_lib(sp)]:
        try:
            enrich_py._ensure_library(name, sp, cache)
        except Exception as e:
            print("  跳过", name, "（需联网，可稍后在网页里下载）", str(e)[:40])
print("富集库就绪")
PY

# 收尾
echo ""
echo "==============================================="
echo "  安装完成！"
echo "==============================================="
echo ""
echo "下一步（3 步）："
echo "  1. 回到 Windows，双击『打开分析网页.bat』启动网页"
echo "  2. 浏览器打开 http://localhost:8501"
echo "  3. 进『数据与参考文件』页：选人/小鼠 → 下载参考文件 → 传数据 → 开跑"
echo ""
echo "详细说明请看『使用说明.md』"
echo "安装结束，可以关闭本窗口了。"
