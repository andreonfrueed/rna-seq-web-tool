#!/usr/bin/env bash
# 启动 RNA 分析网页（在 WSL 内运行）
# 用法（Windows 侧）：双击『打开分析网页.bat』，或手动：
#   wsl bash -lc "cd ~/rna_web_app && bash run.sh"
# 然后浏览器打开 http://localhost:8501
set -e
cd "$(dirname "$0")"
export PYTHONPATH="$PWD"

# RED-06 修复：maxUploadSize 只保留在 .streamlit/config.toml 一处，
# 不再在命令行重复写一遍（两处维护容易改一漏一）。

# 记录 pid（供『更新网页.bat』精确停止本进程，避免 pkill 误杀同机其他 streamlit）
echo $$ > .web.pid
trap 'rm -f .web.pid' EXIT

exec ~/miniconda3/envs/pyseqrna/bin/python -m streamlit run app.py \
  --server.headless true \
  --server.port 8501 \
  --server.address localhost
