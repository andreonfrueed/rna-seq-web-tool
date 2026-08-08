from pathlib import Path

import pandas as pd

root = Path("/mnt/e/数据/new_run_results")
de = root / "4.Differential_Expression"

files = [
    de / "All_gene_expression.xlsx",
    de / "All_gene_expression_sheet.xlsx",
    de / "Filtered_DEGs.xlsx",
    de / "Filtered_upDEGs.xlsx",
    de / "Filtered_downDEGs.xlsx",
    de / "Filtered_DEGs_summary.xlsx",
]

for f in files:
    print("=" * 80)
    print(f.name)
    xl = pd.ExcelFile(f)
    print("sheets:", xl.sheet_names)
    for sh in xl.sheet_names:
        df = pd.read_excel(f, sheet_name=sh, nrows=3)
        print(f"  [{sh}] shape_head={df.shape}")
        print("    columns:", list(df.columns))
        if len(df):
            print("    first row:", df.iloc[0].to_dict())

print("=" * 80)
print("Summary sheet full:")
df = pd.read_excel(de / "Filtered_DEGs_summary.xlsx", sheet_name=0)
print(df.to_string())
