#!/bin/bash
# ============================================================
#  pyseqrna + DESeq2 VST 后处理包装脚本（在 WSL 内运行）
#  用法：bash run_pipeline.sh -c <run.ini>
#  先跑 pyseqrna 主流水线；差异引擎为 deseq2 时再接 R 后处理。
# ============================================================
set -u

INI=""
while [ $# -gt 0 ]; do
  case "$1" in
    -c)
      if [ $# -lt 2 ]; then
        echo "用法: bash run_pipeline.sh -c <run.ini>" >&2
        exit 2
      fi
      INI="$2"
      shift 2
      ;;
    *)
      echo "用法: bash run_pipeline.sh -c <run.ini>" >&2
      exit 2
      ;;
  esac
done
if [ -z "$INI" ]; then
  echo "用法: bash run_pipeline.sh -c <run.ini>" >&2
  exit 2
fi
if [ ! -f "$INI" ]; then
  echo "错误：run.ini 不存在: $INI" >&2
  exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_SCRIPT="$SCRIPT_DIR/DESeq2_vst.R"

ini_value() {
  local key="$1"
  grep -E "^[[:space:]]*${key}[[:space:]]*=" "$INI" \
    | head -n 1 \
    | sed -E "s/^[[:space:]]*${key}[[:space:]]*=[[:space:]]*//"
}

INPUT="$(ini_value input_file)"
FEATURE="$(ini_value feature_file)"
OUTDIR="$(ini_value outdir)"
DIFFEXP_TOOL="$(ini_value diffexp_tool)"

# 1. pyseqrna 主流水线，退出码原样保留
pyseqrna -c "$INI"
PYR=$?
if [ "$PYR" -ne 0 ]; then
  exit "$PYR"
fi

# 2. 仅 deseq2 引擎补 VST 矩阵与热图；pydiffexpress 保持原行为
RRC=0
if [ "$DIFFEXP_TOOL" = "deseq2" ]; then
  echo "DESeq2 后处理开始：生成 VST 矩阵与热图..."
  Rscript "$R_SCRIPT" \
    --counts "$OUTDIR/3.Quantification/Raw_Counts.xlsx" \
    --samples "$INPUT" \
    --gtf "$FEATURE" \
    --outdir "$OUTDIR"
  RRC=$?
  if [ "$RRC" -ne 0 ]; then
    echo "Pipeline execution failed: DESeq2 VST post-processing" >&2
  fi
fi

echo "Wrapper pipeline finished"
exit "$RRC"
