from pathlib import Path
import hashlib
import sys

import pandas as pd

ROOT = Path("/mnt/e/数据/new_run_results")
BACKUP = Path("/mnt/e/数据/new_run_results_原始备份")

COMP = {
    "C-LPS": "LPS-C",
    "C-TTP": "TTP-C",
    "LPS-TTP": "TTP-LPS",
}
NEW_ORDER = ["LPS-C", "TTP-C", "TTP-LPS"]
EXPECT = {
    "LPS-C": {"total": 3720, "up": 2310, "down": 1410},
    "TTP-C": {"total": 3578, "up": 2648, "down": 930},
    "TTP-LPS": {"total": 1071, "up": 850, "down": 221},
}


def rename_comparisons(name):
    for old, new in COMP.items():
        name = name.replace(old, new)
    return name


def rename_old_columns(df):
    df = df.copy()
    for old_comp, new_comp in COMP.items():
        for base in ["baseMean", "logFC", "lfcSE", "stat", "pvalue", "FDR"]:
            old_col = f"{base}({old_comp})"
            if old_col in df.columns:
                df = df.rename(columns={old_col: f"{base}({new_comp})"})
    return df


def swap_direction_label(name):
    if "上调" in name:
        return name.replace("上调", "下调")
    if "下调" in name:
        return name.replace("下调", "上调")
    return name


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def check(cond, msg):
    if not cond:
        print("FAIL:", msg)
        sys.exit(1)
    print("OK:", msg)


def compare_frame_values(new_df, old_df, neg_cols, unchanged_cols, label):
    max_neg = 0.0
    max_unchanged = 0.0
    for col in neg_cols:
        max_neg = max(max_neg, (new_df[col] + old_df[col]).abs().max())
    for col in unchanged_cols:
        max_unchanged = max(
            max_unchanged, (new_df[col] - old_df[col]).abs().max()
        )
    print(f"  {label}: max|new+old| on negated cols={max_neg:.12g}; "
          f"max|new-old| on unchanged cols={max_unchanged:.12g}")
    check(max_neg < 1e-6, f"{label} negated values match")
    check(max_unchanged < 1e-6, f"{label} unchanged values match")


def verify_diff_genes():
    d = ROOT / "4.Differential_Expression" / "diff_genes"
    files = list(d.iterdir())
    names = sorted(f.name for f in files if f.is_file())
    check(len(files) == 9, "diff_genes has exactly 9 files")
    for comp, exp in EXPECT.items():
        total = len((d / f"{comp}.txt").read_text(encoding="utf-8").splitlines())
        up = set((d / f"{comp}_up.txt").read_text(encoding="utf-8").splitlines())
        down = set((d / f"{comp}_down.txt").read_text(encoding="utf-8").splitlines())
        all_genes = set((d / f"{comp}.txt").read_text(encoding="utf-8").splitlines())
        check(total == exp["total"], f"{comp} total rows = {exp['total']}")
        check(len(up) == exp["up"], f"{comp} up rows = {exp['up']}")
        check(len(down) == exp["down"], f"{comp} down rows = {exp['down']}")
        check(all_genes == up | down, f"{comp} all = up | down")
        check(not (up & down), f"{comp} up and down are disjoint")
    il6 = "ENSMUSG00000025746"
    up = set((d / "LPS-C_up.txt").read_text(encoding="utf-8").splitlines())
    down = set((d / "LPS-C_down.txt").read_text(encoding="utf-8").splitlines())
    check(il6 in up, "Il6 in LPS-C_up.txt")
    check(il6 not in down, "Il6 not in LPS-C_down.txt")


def verify_all_gene_expression():
    new = ROOT / "4.Differential_Expression" / "All_gene_expression.xlsx"
    old = BACKUP / "4.Differential_Expression" / "All_gene_expression.xlsx"
    df_new = pd.read_excel(new, sheet_name="Sheet1")
    df_old = pd.read_excel(old, sheet_name="Sheet1")
    df_old = rename_old_columns(df_old)
    check(len(df_new) == 78298, "All_gene_expression.xlsx rows = 78298")
    for comp in NEW_ORDER:
        for base in ["baseMean", "logFC", "lfcSE", "stat", "pvalue", "FDR"]:
            col = f"{base}({comp})"
            check(col in df_new.columns, f"All_gene_expression column {col}")
    il6 = df_new[df_new["Gene"] == "ENSMUSG00000025746"].iloc[0]
    il6_lfc = float(il6["logFC(LPS-C)"])
    print(f"  Il6 logFC(LPS-C) = {il6_lfc:.6f}")
    check(il6_lfc > 0 and abs(il6_lfc - 10.83) < 0.01,
          "Il6 LPS-C logFC positive and near +10.83")
    neg_cols = []
    unchanged_cols = []
    for old_comp, new_comp in COMP.items():
        neg_cols += [f"logFC({new_comp})", f"stat({new_comp})"]
        unchanged_cols += [
            f"baseMean({new_comp})",
            f"lfcSE({new_comp})",
            f"pvalue({new_comp})",
            f"FDR({new_comp})",
        ]
    compare_frame_values(
        df_new, df_old, neg_cols, unchanged_cols, "All_gene_expression"
    )
    for comp in NEW_ORDER:
        ratio = df_new[f"stat({comp})"] / (
            df_new[f"logFC({comp})"] / df_new[f"lfcSE({comp})"]
        )
        check((ratio - 1.0).abs().max() < 1e-8,
              f"All_gene_expression stat relation after negation [{comp}]")


def verify_sheet_file(name):
    new_path = ROOT / "4.Differential_Expression" / name
    old_path = BACKUP / "4.Differential_Expression" / name
    xl_new = pd.ExcelFile(new_path)
    check(xl_new.sheet_names == NEW_ORDER, f"{name} sheet names")
    for old_comp, new_comp in COMP.items():
        df_new = pd.read_excel(new_path, sheet_name=new_comp)
        df_old = pd.read_excel(old_path, sheet_name=old_comp)
        df_old = rename_old_columns(df_old)
        check(len(df_new) == 78298, f"{name} {new_comp} rows = 78298")
        compare_frame_values(
            df_new, df_old, ["logFC", "stat"],
            ["baseMean", "lfcSE", "pvalue", "FDR"],
            f"{name} [{new_comp}]",
        )


def verify_filtered_degs():
    new_path = ROOT / "4.Differential_Expression" / "Filtered_DEGs.xlsx"
    old_path = BACKUP / "4.Differential_Expression" / "Filtered_DEGs.xlsx"
    xl_new = pd.ExcelFile(new_path)
    check(xl_new.sheet_names == NEW_ORDER, "Filtered_DEGs.xlsx sheet names")
    for old_comp, new_comp in COMP.items():
        df_new = pd.read_excel(new_path, sheet_name=new_comp)
        df_old = pd.read_excel(old_path, sheet_name=old_comp)
        df_old = rename_old_columns(df_old)
        check(len(df_new) == EXPECT[new_comp]["total"],
              f"Filtered_DEGs {new_comp} rows = {EXPECT[new_comp]['total']}")
        neg_cols = [f"logFC({new_comp})", f"stat({new_comp})"]
        unchanged_cols = [
            f"baseMean({new_comp})",
            f"lfcSE({new_comp})",
            f"pvalue({new_comp})",
            f"FDR({new_comp})",
        ]
        compare_frame_values(
            df_new, df_old, neg_cols, unchanged_cols,
            f"Filtered_DEGs [{new_comp}]",
        )


def verify_up_down():
    de = ROOT / "4.Differential_Expression"
    new_up_path = de / "Filtered_upDEGs.xlsx"
    new_down_path = de / "Filtered_downDEGs.xlsx"
    old_up_path = BACKUP / "4.Differential_Expression" / "Filtered_upDEGs.xlsx"
    old_down_path = BACKUP / "4.Differential_Expression" / "Filtered_downDEGs.xlsx"
    new_up = pd.ExcelFile(new_up_path)
    new_down = pd.ExcelFile(new_down_path)
    check(new_up.sheet_names == NEW_ORDER, "Filtered_upDEGs.xlsx sheet names")
    check(new_down.sheet_names == NEW_ORDER, "Filtered_downDEGs.xlsx sheet names")
    for old_comp, new_comp in COMP.items():
        up_new = pd.read_excel(new_up_path, sheet_name=new_comp)
        down_new = pd.read_excel(new_down_path, sheet_name=new_comp)
        old_up = pd.read_excel(old_up_path, sheet_name=old_comp)
        old_down = pd.read_excel(old_down_path, sheet_name=old_comp)
        old_up = rename_old_columns(old_up)
        old_down = rename_old_columns(old_down)
        check(len(up_new) == EXPECT[new_comp]["up"],
              f"Filtered_upDEGs {new_comp} rows = {EXPECT[new_comp]['up']}")
        check(len(down_new) == EXPECT[new_comp]["down"],
              f"Filtered_downDEGs {new_comp} rows = {EXPECT[new_comp]['down']}")
        neg_cols = [f"logFC({new_comp})", f"stat({new_comp})"]
        unchanged_cols = [
            f"baseMean({new_comp})",
            f"lfcSE({new_comp})",
            f"pvalue({new_comp})",
            f"FDR({new_comp})",
        ]
        compare_frame_values(
            up_new, old_down, neg_cols, unchanged_cols,
            f"Filtered_upDEGs [{new_comp}] vs old down",
        )
        compare_frame_values(
            down_new, old_up, neg_cols, unchanged_cols,
            f"Filtered_downDEGs [{new_comp}] vs old up",
        )


def verify_summary():
    new_path = ROOT / "4.Differential_Expression" / "Filtered_DEGs_summary.xlsx"
    df = pd.read_excel(new_path, sheet_name="Sheet1")
    rows = df.to_dict("records")
    check(len(rows) == 3, "summary has 3 rows")
    expected = [
        ("LPS-C", 3720, 2310, 1410),
        ("TTP-C", 3578, 2648, 930),
        ("TTP-LPS", 1071, 850, 221),
    ]
    for row, exp in zip(rows, expected):
        check(
            (row["Comparisons"], row["Total_DEGs"], row["Up_DEGs"], row["Down_DEGs"])
            == exp,
            f"summary row {exp}",
        )


def verify_go_ontology():
    base_new = ROOT / "6.Functional_Annotation" / "Gene_Ontology"
    base_old = BACKUP / "6.Functional_Annotation" / "Gene_Ontology"
    for old_comp, new_comp in COMP.items():
        check(not (base_new / old_comp).exists(), f"Gene_Ontology no {old_comp} dir")
        check((base_new / new_comp).exists(), f"Gene_Ontology {new_comp} dir exists")
        for f_old in (base_old / old_comp).iterdir():
            expected_name = rename_comparisons(f_old.name)
            f_new = base_new / new_comp / expected_name
            check(f_new.exists(), f"Gene_Ontology file {expected_name} exists")
            check(sha256(f_old) == sha256(f_new),
                  f"Gene_Ontology content unchanged {expected_name}")


def verify_go_kegg():
    base_new = ROOT / "GO.KEGG" / "上.下调" / "GO_KEGG_富集"
    base_old = BACKUP / "GO.KEGG" / "上.下调" / "GO_KEGG_富集"
    for old_comp, new_comp in COMP.items():
        check(not (base_new / old_comp).exists(), f"GO_KEGG no {old_comp} dir")
        check((base_new / new_comp).exists(), f"GO_KEGG {new_comp} dir exists")
        for old_sub, new_sub in [("up", "down"), ("down", "up")]:
            old_dir = base_old / old_comp / old_sub
            new_dir = base_new / new_comp / new_sub
            for f_old in old_dir.iterdir():
                f_new = new_dir / f_old.name
                check(f_new.exists(), f"GO_KEGG {new_comp}/{new_sub}/{f_old.name} exists")
                check(sha256(f_old) == sha256(f_new),
                      f"GO_KEGG content unchanged {new_comp}/{new_sub}/{f_old.name}")


def verify_legacy():
    old_dir = BACKUP / "GO.KEGG"
    new_dir = ROOT / "GO.KEGG" / "旧版"
    for name in [
        "GO_result.csv",
        "KEGG_result.csv",
        "GO_dotplot.png",
        "KEGG_dotplot.png",
    ]:
        check((new_dir / name).exists(), f"旧版/{name} exists")
        check(not (ROOT / "GO.KEGG" / name).exists(),
              f"top-level GO.KEGG/{name} moved away")
        check(sha256(old_dir / name) == sha256(new_dir / name),
              f"旧版/{name} content unchanged")


def expected_image_name(name):
    return swap_direction_label(rename_comparisons(name))


def verify_images():
    old_dir = BACKUP / "图片"
    new_dir = ROOT / "图片"
    for f_old in old_dir.iterdir():
        expected = expected_image_name(f_old.name)
        f_new = new_dir / expected
        check(f_new.exists(), f"图片/{expected} exists")
        check(sha256(f_old) == sha256(f_new), f"图片/{expected} content unchanged")


def verify_visualization():
    for sub in ["Volcano_Plots", "MA_Plots", "Venn_Plots"]:
        old_dir = BACKUP / "5.Visualization" / sub
        new_dir = ROOT / "5.Visualization" / sub
        for f_old in old_dir.iterdir():
            expected = rename_comparisons(f_old.name)
            f_new = new_dir / expected
            check(f_new.exists(), f"5.Visualization/{sub}/{expected} exists")
            check(sha256(f_old) == sha256(f_new),
                  f"5.Visualization/{sub}/{expected} content unchanged")
    for name in ["Filtered_DEG.png", "Filtered_DEG.pdf"]:
        old_f = BACKUP / "4.Differential_Expression" / name
        new_f = ROOT / "4.Differential_Expression" / name
        check(sha256(old_f) == sha256(new_f), f"4.Differential_Expression/{name} unchanged")


def verify_no_old_names():
    old_names = list(COMP)
    hits = []
    for p in ROOT.rglob("*"):
        if any(old in p.name for old in old_names):
            hits.append(str(p))
    check(not hits, "no file/dir name still contains old comparison names")
    tmp_files = list(ROOT.rglob("*.tmp_fix.xlsx"))
    check(not tmp_files, "no temporary xlsx files remain")


def main():
    verify_diff_genes()
    verify_all_gene_expression()
    verify_sheet_file("All_gene_expression_sheet.xlsx")
    verify_filtered_degs()
    verify_up_down()
    verify_summary()
    verify_go_ontology()
    verify_go_kegg()
    verify_legacy()
    verify_images()
    verify_visualization()
    verify_no_old_names()
    print("ALL VERIFICATION PASSED")


if __name__ == "__main__":
    main()
