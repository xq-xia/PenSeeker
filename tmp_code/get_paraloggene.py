from pathlib import Path
import pandas as pd

# 1. 文件路径配置
gene_list_path = Path(r"E:\外显率预测课题\基因背景网络\KGPEN\annotation\data_source\gene\gene_list.xlsx")
para_path = Path(r"E:\外显率预测课题\基因背景网络\KGPEN\annotation\data_source\gene\gene_paralog\paralog_ensembl_gene.txt")
out_dir = para_path.parent

# 2. 相似度过滤参数配置
SIM_THRESHOLD = 50.0
# "both": 双向相似度均 > 50%（推荐，防止长短基因片段局部匹配假阳性）
# "target": 仅要求 target %id > 50%
# "either": 任一方向 > 50%
SIM_MODE = "both"

# 3. 读取目标基因列表
genes_df = pd.read_excel(gene_list_path)
gene_col = "gene_symbol" if "gene_symbol" in genes_df.columns else genes_df.columns[0]
target_genes = set(
    genes_df[gene_col]
    .dropna()
    .astype(str)
    .str.strip()
    .loc[lambda s: s != ""]
)
print(f"目标基因总数: {len(target_genes)}")

# 4. 读取最新的 paralog_ensembl_gene.txt
try:
    para_df = pd.read_csv(para_path, sep="\t", low_memory=False)
except Exception:
    para_df = pd.read_csv(para_path, sep=None, engine="python")

col_map = {
    "Gene stable ID": "gene_id",
    "Gene name": "gene_symbol",
    "Human paralogue gene stable ID": "paralog_id",
    "Human paralogue associated gene name": "paralog_symbol",
    "Paralogue %id. target Human gene identical to query gene": "target_pct_id",
    "Paralogue %id. query gene identical to target Human gene": "query_pct_id",
}
df = para_df.rename(columns=col_map).copy()

# 5. 基础清洗与筛选
for col in ["gene_id", "gene_symbol", "paralog_id", "paralog_symbol"]:
    df[col] = df[col].fillna("").astype(str).str.strip()

df["target_pct_id"] = pd.to_numeric(df["target_pct_id"], errors="coerce")
df["query_pct_id"] = pd.to_numeric(df["query_pct_id"], errors="coerce")

# 保留在 gene_list 中的基因，且存在有效旁系同源基因（排除空值与自身）
df = df[
    df["gene_symbol"].isin(target_genes)
    & (df["paralog_id"] != "")
    & (df["paralog_symbol"] != "")
    & (df["gene_id"] != df["paralog_id"])
    & (df["gene_symbol"] != df["paralog_symbol"])
].copy()

# 6. 按序列相似度 > 50% 过滤
if SIM_MODE == "both":
    cond_sim = (df["target_pct_id"] > SIM_THRESHOLD) & (df["query_pct_id"] > SIM_THRESHOLD)
elif SIM_MODE == "target":
    cond_sim = df["target_pct_id"] > SIM_THRESHOLD
elif SIM_MODE == "either":
    cond_sim = (df["target_pct_id"] > SIM_THRESHOLD) | (df["query_pct_id"] > SIM_THRESHOLD)
else:
    raise ValueError(f"不支持的 SIM_MODE: {SIM_MODE}")

filtered_df = (
    df[cond_sim]
    .sort_values(["gene_symbol", "target_pct_id"], ascending=[True, False])
    .drop_duplicates(subset=["gene_symbol", "paralog_symbol"], keep="first")
    .reset_index(drop=True)
)

# 7. 保存配对明细表（供下一步关联 GTEx 脑组织表达使用）
pair_out_path = out_dir / f"gene_list_paralogs_sim50_{SIM_MODE}.tsv"
filtered_df.to_csv(pair_out_path, sep="\t", index=False)

# 8. 按 gene_symbol 聚合生成基因级候选代偿列表（方便直观查看每个基因有几个高相似旁系同源基因）
gene_summary = (
    filtered_df.groupby("gene_symbol", as_index=False)
    .agg(
        n_high_sim_paralogs=("paralog_symbol", "nunique"),
        high_sim_paralogs=("paralog_symbol", lambda s: ";".join(s)),
        max_target_pct_id=("target_pct_id", "max"),
        max_query_pct_id=("query_pct_id", "max"),
    )
)
summary_out_path = out_dir / f"gene_level_paralogs_sim50_{SIM_MODE}.tsv"
gene_summary.to_csv(summary_out_path, sep="\t", index=False)

print(f"\n过滤模式: SIM_MODE = '{SIM_MODE}' (> {SIM_THRESHOLD}%)")
print(f"保留高相似同源配对行数: {len(filtered_df)}")
print(f"拥有至少 1 个高相似旁系同源基因的目标基因数: {len(gene_summary)}")
print(f"明细表已保存至: {pair_out_path}")
print(f"基因级汇总表已保存至: {summary_out_path}")
print("\n前 15 行明细预览:")
print(filtered_df.head(15))