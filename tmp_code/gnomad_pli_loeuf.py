import pandas as pd
import numpy as np

input_tsv = r"D:\下载\gnomad.v4.1.constraint_metrics.tsv"
output_tsv = "gnomad_v4.1_mane_NM_pli_loeuf.tsv"

cols = [
    "gene",
    "gene_id",
    "transcript",
    "mane_select",
    "canonical",
    "lof.pLI",
    "lof.oe_ci.upper",
]
df = pd.read_csv(input_tsv, sep="\t", usecols=cols, low_memory=False)

# 1. 清洗字段
df["gene_symbol"] = df["gene"].astype(str).str.strip()
df["raw_gene_id"] = df["gene_id"].astype(str).str.split(".").str[0].str.strip()
df["transcript"] = df["transcript"].astype(str).str.strip()
df["pLI"] = pd.to_numeric(df["lof.pLI"], errors="coerce")
df["LOEUF"] = pd.to_numeric(df["lof.oe_ci.upper"], errors="coerce")

# 过滤掉 gene_symbol 为空或无效的行
df = df[df["gene_symbol"].notna() & ~df["gene_symbol"].isin(["", "NA", "nan", "None"])].copy()

is_mane = df["mane_select"].astype(str).str.strip().str.lower() == "true"
is_nm = df["transcript"].str.startswith("NM")
is_canon = df["canonical"].astype(str).str.strip().str.lower() == "true"
has_score = df["pLI"].notna() | df["LOEUF"].notna()

# 2. 计算四级来源优先级
# 1: MANE_NM (mane_select==True 且 NM 开头，且有有效评分)
# 2: Canonical_fallback (canonical==True 且有有效评分)
# 3: Min_LOEUF_fallback (有有效评分，按 LOEUF 最小取)
# 4: Missing (该基因所有转录本均无有效评分)
conditions = [
    is_mane & is_nm & has_score,
    is_canon & has_score,
    has_score,
]
df["_tier_rank"] = np.select(conditions, [1, 2, 3], default=4)

# 同级排序辅助键：优先 MANE_NM > Canonical > ENSG 主记录，再按 LOEUF 最小排序
df["_sub_rank"] = np.select([is_mane & is_nm, is_canon], [0, 1], default=2)
df["_is_non_ensg"] = (~df["raw_gene_id"].str.startswith("ENSG")).astype(int)
df["_loeuf_sort"] = df["LOEUF"].fillna(1e9)

# 3. 关键改动：按 gene_symbol 分组排序并去重，确保每个基因只保留 1 行
df_sorted = df.sort_values(
    by=["gene_symbol", "_tier_rank", "_sub_rank", "_loeuf_sort", "_is_non_ensg", "transcript"],
    ascending=[True, True, True, True, True, True],
)

gene_df = df_sorted.drop_duplicates(subset=["gene_symbol"], keep="first").copy()

# 4. 映射 source_tier 标签，并区分 ENSG ID 与 Entrez ID
tier_label_map = {
    1: "MANE_NM",
    2: "Canonical_fallback",
    3: "Min_LOEUF_fallback",
    4: "Missing",
}
gene_df["source_tier"] = gene_df["_tier_rank"].map(tier_label_map)
gene_df["ensembl_gene_id"] = np.where(
    gene_df["raw_gene_id"].str.startswith("ENSG"), gene_df["raw_gene_id"], np.nan
)

# 如果选中的是 NM 行（其 raw_gene_id 是 Entrez ID），顺便把同基因的 ENSG ID 也补全方便后续按 ENSG 关联
ensg_map = (
    df[df["raw_gene_id"].str.startswith("ENSG")]
    .sort_values(["_sub_rank", "raw_gene_id"])
    .drop_duplicates(subset=["gene_symbol"], keep="first")
    .set_index("gene_symbol")["raw_gene_id"]
)
gene_df["ensembl_gene_id"] = gene_df["ensembl_gene_id"].fillna(gene_df["gene_symbol"].map(ensg_map))

final_cols = [
    "gene_symbol",
    "ensembl_gene_id",
    "transcript",
    "source_tier",
    "pLI",
    "LOEUF",
]
result = gene_df[final_cols].reset_index(drop=True)

print("=== 清洗去重统计 ===")
print(f"唯一 gene_symbol 总数: {len(result)}")
print(f"检查重复 gene_symbol 数: {result['gene_symbol'].duplicated().sum()}")
print("\n各 source_tier 分布:")
print(result["source_tier"].value_counts())
print("\nMATR3 结果验证:")
print(result[result["gene_symbol"] == "MATR3"])

result.to_csv(output_tsv, sep="\t", index=False)