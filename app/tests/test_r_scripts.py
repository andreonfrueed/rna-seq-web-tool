"""DESeq2_vst.R 结构检查（本机无 R 环境：验证脚本结构完整性与双格式导出）。

R 脚本改动无法在 Windows 端真跑，用结构检查兜底：
- 括号平衡（语法级粗检）
- 双格式导出（PNG 预览 + PDF 矢量投稿）覆盖全部图型
- 关键函数定义齐全
真实运行验证在 WSL（.wsl-verify/r_verify.R）。
"""
from __future__ import annotations
from pathlib import Path

R_SCRIPT = Path(__file__).resolve().parent.parent / "r_scripts" / "DESeq2_vst.R"


def _r_text() -> str:
    assert R_SCRIPT.exists(), f"R 脚本不存在: {R_SCRIPT}"
    return R_SCRIPT.read_text(encoding="utf-8")


def test_r_script_balanced_parens():
    text = _r_text()
    assert text.count("(") == text.count(")"), "圆括号不平衡"
    assert text.count("{") == text.count("}"), "花括号不平衡"


def test_r_script_plot_functions_defined():
    text = _r_text()
    for fn in ["save_figure", "render_pca", "render_volcanoes",
               "render_ma_plots", "render_deg_heatmap"]:
        assert f"{fn} <- function" in text, f"缺少函数 {fn}"


def test_r_script_double_format_export():
    """期刊矢量导出：全部 5 类图都应有 PDF（矢量投稿）+ PNG（网页预览）。"""
    text = _r_text()
    # 双格式 helper 已定义且被调用
    assert "save_figure <- function" in text
    # base graphics 三图（PCA/火山/MA）走 save_figure
    assert "save_figure(pca_png, pca_pdf" in text
    assert "save_figure(volcano_png, volcano_pdf" in text
    assert "save_figure(ma_png, ma_pdf" in text
    # pheatmap 两图走 filename 双格式
    assert "draw_deg_heat(deg_png)" in text
    assert "draw_deg_heat(deg_pdf)" in text
    assert "draw_cluster_heat(cluster_png)" in text
    assert "draw_cluster_heat(cluster_pdf)" in text
    # PDF 输出路径齐全（5 类图）
    assert text.count(".pdf") >= 5
