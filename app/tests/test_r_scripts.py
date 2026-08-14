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


def test_r_script_deg_heatmap_prefers_deseq_csv():
    """BUG-19 回归：DEG 热图数据源首选 DESeq2 CSV，Filtered_DEGs.xlsx 仅作回退。

    结构检查（本机无 R，只能验证源码顺序，无法真跑）：
    CSV 首选分支必须出现在 xlsx 回退分支之前。
    """
    text = _r_text()
    csv_pos = text.find("length(csv_files) > 0L")   # render_deg_heatmap 的 CSV 首选分支
    xlsx_pos = text.find("file.exists(deg_xlsx)")   # xlsx 回退分支（全文件唯一）
    assert csv_pos != -1, "DESeq2 CSV 首选分支缺失"
    assert xlsx_pos != -1, "xlsx 回退分支缺失"
    assert csv_pos < xlsx_pos, "CSV 分支应在 xlsx 回退之前（数据源首选 CSV）"


def test_r_script_deg_heatmap_same_threshold_as_diffexp():
    """BUG-19 回归：CSV 分支显著标准与差异表/火山图一致（padj<0.05 且 |log2FC|>=1）。"""
    text = _r_text()
    assert "deg_df$padj < 0.05" in text
    assert "abs(deg_df$log2FoldChange) >= 1" in text


def test_r_script_all_heatmaps_fixed_sample_order_bug20_bug23():
    """BUG-20/23 回归：DEG 热图与样本聚类热图都固定样本列顺序，不按相似度打乱。"""
    text = _r_text()
    deg_start = text.find("draw_deg_heat <- function")
    deg_end = text.find("run_deseq_results <- function")
    cluster_start = text.find("draw_cluster_heat <- function")
    cluster_end = text.find("status <- tryCatch")
    # 防御：所有锚点都要命中且切片范围有效，避免 find 返回 -1 走负索引得空串、断言假绿
    assert deg_start != -1 and deg_end != -1 and deg_start < deg_end, "DEG 热图锚点定位失败"
    assert cluster_start != -1 and cluster_end != -1 and cluster_start < cluster_end, "聚类热图锚点定位失败"
    deg_seg = text[deg_start:deg_end]
    assert "cluster_cols = FALSE" in deg_seg, "DEG 热图应固定样本列顺序"
    cluster_seg = text[cluster_start:cluster_end]
    assert "cluster_cols = FALSE" in cluster_seg, "样本聚类热图也应固定样本列顺序"
