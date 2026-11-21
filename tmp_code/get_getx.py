from pathlib import Path
import pandas as pd

# 1. 路径配置（请把 gtex_path 改成你实际下载的文件名）
base_dir = Path(r"E:\外显率预测课题\基因背景网络\KGPEN\annotation\data_source\gene\gene_paralog")
para_pairs_path = base_dir / "gene_list_paralogs_sim50_both.tsv"
gtex_path = base_dir / "GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct"

# 2. 旁系同源基因脑组织“高表达”阈值（TPM）
PARALOG_BRAIN_TPM_THRESHOLD = 10.0

# 3. 读取 GTEx 中位数 TPM 表（跳过前 2 行 GCT 头部说明）
gtex = pd.read_csv(gtex_path, sep="\t", skiprows=2, low_memory=False)
gtex["ensembl_gene_id"] = gtex["Name"].astype(str).str.split(".").str[0].str.strip()
gtex["gene_symbol"] = gtex["Description"].astype(str).str.strip()

brain_cols = [c for c in gtex.columns if str(c).startswith("Brain_")]
if not brain_cols:
    raise ValueError("未在 GTEx 文件中找到以 'Brain_' 开头的脑组织列，请检查列名。")

for c in brain_cols:
    gtex[c] = pd.to_numeric(gtex[c], errors="coerce")

gtex["brain_max_tpm"] = gtex[brain_cols].max(axis=1)

gtex_by_id = (
    gtex.sort_values("brain_max_tpm", ascending=False)
    .drop_duplicates(subset=["ensembl_gene_id"], keep="first")
    .set_index("ensembl_gene_id")["brain_max_tpm"]
)
gtex_by_sym = (
    gtex.sort_values("brain_max_tpm", ascending=False)
    .drop_duplicates(subset=["gene_symbol"], keep="first")
    .set_index("gene_symbol")["brain_max_tpm"]
)

# 4. 读取双向 >50% 相似度的旁系同源配对表并匹配旁系基因脑表达
pairs_df = pd.read_csv(para_pairs_path, sep="\t", low_memory=False)
para_id_clean = pairs_df["paralog_id"].astype(str).str.split(".").str[0].str.strip()
para_sym_clean = pairs_df["paralog_symbol"].astype(str).str.strip()

pairs_df["paralog_brain_max_tpm"] = para_id_clean.map(gtex_by_id)
miss_mask = pairs_df["paralog_brain_max_tpm"].isna()
pairs_df.loc[miss_mask, "paralog_brain_max_tpm"] = para_sym_clean.loc[miss_mask].map(gtex_by_sym)

# 5. 筛选具有高表达旁系同源基因的配对与目标基因列表
valid_comp_pairs = pairs_df[
    pairs_df["paralog_brain_max_tpm"].fillna(0.0) >= PARALOG_BRAIN_TPM_THRESHOLD
].copy()

valid_comp_pairs = valid_comp_pairs.sort_values(
    ["gene_symbol", "paralog_brain_max_tpm", "target_pct_id"],
    ascending=[True, False, False],
)

comp_genes_df = (
    valid_comp_pairs.groupby("gene_symbol", as_index=False)
    .agg(
        compensating_paralogs=("paralog_symbol", lambda s: ";".join(dict.fromkeys(s.astype(str).str.strip()))),
        n_compensating_paralogs=("paralog_symbol", "nunique"),
        max_paralog_brain_tpm=("paralog_brain_max_tpm", "max"),
    )
    .sort_values("gene_symbol")
    .reset_index(drop=True)
)

# 纯基因列表文件（每行一个 gene_symbol）
gene_list_out = base_dir / "genes_with_brain_paralog_compensation.txt"
comp_genes_df["gene_symbol"].to_csv(gene_list_out, index=False, header=False)

# 带代偿基因明细表（方便核查对应哪些旁系基因）
detail_out = base_dir / "genes_with_brain_paralog_compensation_detail.tsv"
comp_genes_df.to_csv(detail_out, sep="\t", index=False)

print(f"具有脑高表达旁系同源基因的基因总数: {len(comp_genes_df)}")
print(f"纯基因列表已保存至: {gene_list_out}")
print(f"代偿明细表已保存至: {detail_out}")
print("\n前 20 个基因预览:")
print(comp_genes_df["gene_symbol"].head(20).tolist())