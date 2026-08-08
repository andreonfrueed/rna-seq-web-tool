from pathlib import Path
import json
import sys

import pandas as pd

ROOT = Path("/mnt/e/数据/new_run_results")
COMP = {
    "C-LPS": "LPS-C",
    "C-TTP": "TTP-C",
    "LPS-TTP": "TTP-LPS",
}
NEW_ORDER = ["LPS-C", "TTP-C", "TTP-LPS"]


def rename_comparisons(name: str) -> str:
    out = name
    for old, new in COMP.items():
        out = out.replace(old, new)
    return out


def swap_direction_label(name: str) -> str:
    if "上调" in name:
        return name.replace("上调", "下调")
    if "下调" in name:
        return name.replace("下调", "上调")
    return name


def write_workbook(out_path: Path, sheets):
    tmp = out_path.with_name(out_path.name + ".tmp_fix.xlsx")
    if tmp.exists():
        tmp.unlink()
    with pd.ExcelWriter(tmp, engine="openpyxl") as writer:
        for sheet_name, df in sheets.items():
            df.to_excel(writer, sheet_name=sheet_name, index=False)
    tmp.replace(out_path)


def negate_named(df, columns):
    for col in columns:
        df[col] = -df[col]


def plain_sheet(df):
    df = df.copy()
    negate_named(df, ["logFC", "stat"])
    return df


def suffixed_sheet(df, old_comp):
    new_comp = COMP[old_comp]
    df = df.copy()
    df = df.rename(
        columns={
            f"baseMean({old_comp})": f"baseMean({new_comp})",
            f"logFC({old_comp})": f"logFC({new_comp})",
            f"lfcSE({old_comp})": f"lfcSE({new_comp})",
            f"stat({old_comp})": f"stat({new_comp})",
            f"pvalue({old_comp})": f"pvalue({new_comp})",
            f"FDR({old_comp})": f"FDR({new_comp})",
        }
    )
    negate_named(df, [f"logFC({new_comp})", f"stat({new_comp})"])
    return df


def check_stat_relation(path):
    xl = pd.ExcelFile(path)
    for old_comp in COMP:
        if old_comp in xl.sheet_names:
            df = pd.read_excel(path, sheet_name=old_comp)
            if f"stat({old_comp})" in df.columns:
                stat_col = f"stat({old_comp})"
                logfc_col = f"logFC({old_comp})"
                lfcse_col = f"lfcSE({old_comp})"
            else:
                stat_col = "stat"
                logfc_col = "logFC"
                lfcse_col = "lfcSE"
        else:
            df = pd.read_excel(path, sheet_name="Sheet1")
            stat_col = f"stat({old_comp})"
            logfc_col = f"logFC({old_comp})"
            lfcse_col = f"lfcSE({old_comp})"
        ratio = df[stat_col] / (df[logfc_col] / df[lfcse_col])
        diff = (ratio - 1.0).abs()
        print(
            f"  stat relation [{old_comp}]: max_abs_dev={diff.max():.12g}, "
            f"median_dev={diff.median():.12g}, n={len(df)}"
        )
        if diff.max() > 1e-8:
            raise SystemExit(
                f"stat is not equal to logFC/lfcSE in {path.name} [{old_comp}]"
            )


def step_diff_genes():
    d = ROOT / "4.Differential_Expression" / "diff_genes"
    for old_comp, new_comp in COMP.items():
        base_src = d / f"{old_comp}.txt"
        base_dst = d / f"{new_comp}.txt"
        up_src = d / f"{old_comp}_up.txt"
        down_src = d / f"{old_comp}_down.txt"
        new_up = d / f"{new_comp}_up.txt"
        new_down = d / f"{new_comp}_down.txt"
        # Copy first so old source files survive until all new names exist.
        new_up.write_bytes(down_src.read_bytes())
        new_down.write_bytes(up_src.read_bytes())
        base_src.replace(base_dst)
        up_src.unlink()
        down_src.unlink()
    print("diff_genes done")


def step_all_gene_expression():
    path = ROOT / "4.Differential_Expression" / "All_gene_expression.xlsx"
    print("precheck All_gene_expression.xlsx")
    check_stat_relation(path)
    df = pd.read_excel(path, sheet_name="Sheet1")
    print(f"  rows before transform: {len(df)}")
    for old_comp, new_comp in COMP.items():
        df = df.rename(
            columns={
                f"baseMean({old_comp})": f"baseMean({new_comp})",
                f"logFC({old_comp})": f"logFC({new_comp})",
                f"lfcSE({old_comp})": f"lfcSE({new_comp})",
                f"stat({old_comp})": f"stat({new_comp})",
                f"pvalue({old_comp})": f"pvalue({new_comp})",
                f"FDR({old_comp})": f"FDR({new_comp})",
            }
        )
        negate_named(df, [f"logFC({new_comp})", f"stat({new_comp})"])
    write_workbook(path, {"Sheet1": df})
    print("All_gene_expression.xlsx done")


def step_all_gene_expression_sheet():
    path = ROOT / "4.Differential_Expression" / "All_gene_expression_sheet.xlsx"
    print("precheck All_gene_expression_sheet.xlsx")
    check_stat_relation(path)
    xl = pd.ExcelFile(path)
    sheets = {}
    for old_comp in COMP:
        sheets[COMP[old_comp]] = plain_sheet(
            pd.read_excel(path, sheet_name=old_comp)
        )
    write_workbook(path, sheets)
    print("All_gene_expression_sheet.xlsx done")


def step_filtered_degs():
    path = ROOT / "4.Differential_Expression" / "Filtered_DEGs.xlsx"
    print("precheck Filtered_DEGs.xlsx")
    check_stat_relation(path)
    xl = pd.ExcelFile(path)
    sheets = {}
    for old_comp in COMP:
        sheets[COMP[old_comp]] = suffixed_sheet(
            pd.read_excel(path, sheet_name=old_comp), old_comp
        )
    write_workbook(path, sheets)
    print("Filtered_DEGs.xlsx done")


def step_filtered_up_down():
    de = ROOT / "4.Differential_Expression"
    up_path = de / "Filtered_upDEGs.xlsx"
    down_path = de / "Filtered_downDEGs.xlsx"
    print("precheck Filtered_upDEGs.xlsx / Filtered_downDEGs.xlsx")
    check_stat_relation(up_path)
    check_stat_relation(down_path)
    xl_up = pd.ExcelFile(up_path)
    xl_down = pd.ExcelFile(down_path)
    old_up = {
        old_comp: pd.read_excel(up_path, sheet_name=old_comp)
        for old_comp in COMP
    }
    old_down = {
        old_comp: pd.read_excel(down_path, sheet_name=old_comp)
        for old_comp in COMP
    }
    new_up = {
        COMP[old_comp]: suffixed_sheet(old_down[old_comp], old_comp)
        for old_comp in COMP
    }
    new_down = {
        COMP[old_comp]: suffixed_sheet(old_up[old_comp], old_comp)
        for old_comp in COMP
    }
    write_workbook(up_path, new_up)
    write_workbook(down_path, new_down)
    print("Filtered_upDEGs.xlsx / Filtered_downDEGs.xlsx done")


def step_summary():
    path = ROOT / "4.Differential_Expression" / "Filtered_DEGs_summary.xlsx"
    df = pd.read_excel(path, sheet_name="Sheet1")
    df["Comparisons"] = df["Comparisons"].replace(COMP)
    df[["Up_DEGs", "Down_DEGs"]] = df[["Down_DEGs", "Up_DEGs"]].values
    write_workbook(path, {"Sheet1": df})
    print("Filtered_DEGs_summary.xlsx done")


def step_gene_ontology():
    base = ROOT / "6.Functional_Annotation" / "Gene_Ontology"
    for old_comp, new_comp in COMP.items():
        (base / old_comp).rename(base / new_comp)
    for comp_dir in NEW_ORDER:
        d = base / comp_dir
        for f in list(d.iterdir()):
            if not f.is_file():
                continue
            new_name = rename_comparisons(f.name)
            if new_name != f.name:
                f.rename(d / new_name)
    print("Gene_Ontology done")


def remap_json_keys(obj):
    if isinstance(obj, dict):
        return {
            rename_comparisons(str(k)): remap_json_keys(v)
            for k, v in obj.items()
        }
    if isinstance(obj, list):
        return [remap_json_keys(v) for v in obj]
    return obj


def step_go_kegg():
    base = ROOT / "GO.KEGG" / "上.下调" / "GO_KEGG_富集"
    for old_comp, new_comp in COMP.items():
        (base / old_comp).rename(base / new_comp)
    for comp in NEW_ORDER:
        d = base / comp
        up = d / "up"
        down = d / "down"
        tmp = d / ".up_tmp"
        up.rename(tmp)
        down.rename(up)
        tmp.rename(down)
        for sub in ("up", "down"):
            stats = d / sub / "_stats.json"
            if stats.exists():
                data = json.loads(stats.read_text(encoding="utf-8"))
                stats.write_text(
                    json.dumps(remap_json_keys(data), ensure_ascii=False, indent=2),
                    encoding="utf-8",
                )
    print("GO_KEGG_富集 done")


def step_legacy_go_kegg():
    src = ROOT / "GO.KEGG"
    dst = src / "旧版"
    dst.mkdir(exist_ok=True)
    for name in [
        "GO_result.csv",
        "KEGG_result.csv",
        "GO_dotplot.png",
        "KEGG_dotplot.png",
    ]:
        f = src / name
        if f.exists():
            f.rename(dst / name)
    print("GO.KEGG 旧版 done")


def step_images():
    d = ROOT / "图片"
    for f in list(d.iterdir()):
        if not f.is_file():
            continue
        new_name = swap_direction_label(rename_comparisons(f.name))
        if new_name != f.name:
            f.rename(d / new_name)
    print("图片 rename done")


def step_visualization():
    base = ROOT / "5.Visualization"
    for sub in ("Volcano_Plots", "MA_Plots", "Venn_Plots"):
        d = base / sub
        for f in list(d.iterdir()):
            if not f.is_file():
                continue
            new_name = rename_comparisons(f.name)
            if new_name != f.name:
                f.rename(d / new_name)
    print("5.Visualization rename done")


def main():
    start = sys.argv[1] if len(sys.argv) > 1 else None
    steps = [
        ("diff_genes", step_diff_genes),
        ("all_gene_expression", step_all_gene_expression),
        ("all_gene_expression_sheet", step_all_gene_expression_sheet),
        ("filtered_degs", step_filtered_degs),
        ("filtered_up_down", step_filtered_up_down),
        ("summary", step_summary),
        ("gene_ontology", step_gene_ontology),
        ("go_kegg", step_go_kegg),
        ("legacy_go_kegg", step_legacy_go_kegg),
        ("images", step_images),
        ("visualization", step_visualization),
    ]
    started = start is None
    for name, fn in steps:
        if name == start:
            started = True
        if started:
            fn()
    print("ALL STEPS COMPLETE")


if __name__ == "__main__":
    main()
