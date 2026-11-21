from pathlib import Path
import pandas as pd
import networkx as nx
import os

# ================= 1. 路径与参数配置 =================
DATA_DIR = Path(r"E:\外显率预测课题\基因背景网络")
ASD_GENE_FILE = DATA_DIR / "KGPEN" / "data" / "SFARI-Gene_genes_07-12-2026release_11-19-2026export.csv"         # ASD 基因集文件（Tab 分隔）
KG_FILE = DATA_DIR / "KGRD" / "kg" / "kg.csv"                                                                   # KGRD 宽表路径（Tab 分隔）
OUT_DIR = DATA_DIR / "KGPEN" / "data"                                                                           # 输出目录
HOPS = 2                                                                                                        # 第 1 层纯基因背景网跳数：1 或 2
ALLOWED_SCORES = None                                                                                           # 例如设为 {1, 2} 仅用高分 SFARI 基因作种子；None 为全量
GENE_RELATIONS = None                                                                                           # 例如设为 {"interact"}；None 为保留基因间所有关系
INCLUDE_TYPES_L2 = {"effect/phenotype", "disease", "biological_process", "pathway"}                             # 第 2 层纳入的异质实体类型（按需调整）


def normalize_str(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    l1_dir = OUT_DIR / "layer1_gene_bg"
    l2_dir = OUT_DIR / "layer2_pheno_disease"
    l1_dir.mkdir(parents=True, exist_ok=True)
    l2_dir.mkdir(parents=True, exist_ok=True)

    # ================= 2. 读取 ASD 种子基因 =================
    asd_df = pd.read_csv(ASD_GENE_FILE, sep=",", dtype=str)
    asd_df["gene_norm"] = normalize_str(asd_df["gene-symbol"]).str.upper()
    asd_df = asd_df[asd_df["gene_norm"] != ""]

    if ALLOWED_SCORES is not None and "gene-score" in asd_df.columns:
        valid_scores = {str(s).strip() for s in ALLOWED_SCORES}
        seed_df = asd_df[normalize_str(asd_df["gene-score"]).isin(valid_scores)].copy()
    else:
        seed_df = asd_df.copy()

    seed_genes = set(seed_df["gene_norm"].unique())
    print(f"[Info] ASD 种子基因数: {len(seed_genes)}")

    # ================= 3. 读取 KG 并规范化 =================
    use_cols = ["x_id", "y_id", "relation", "x_type", "x_name", "y_type", "y_name"]
    kg_df = pd.read_csv(KG_FILE, sep=",", usecols=lambda c: c in use_cols, dtype=str)

    kg_df["x_type_n"] = normalize_str(kg_df["x_type"])
    kg_df["y_type_n"] = normalize_str(kg_df["y_type"])
    kg_df["rel_n"] = normalize_str(kg_df["relation"])
    kg_df["x_name_n"] = normalize_str(kg_df["x_name"])
    kg_df["y_name_n"] = normalize_str(kg_df["y_name"])

    # ================= 4. 构建 Layer 1：纯基因背景网 =================
    g2g = kg_df[
        (kg_df["x_type_n"] == "gene/protein") &
        (kg_df["y_type_n"] == "gene/protein") &
        (kg_df["x_name_n"] != "") &
        (kg_df["y_name_n"] != "")
    ].copy()

    if GENE_RELATIONS is not None:
        g2g = g2g[g2g["rel_n"].isin(GENE_RELATIONS)]

    g2g["u_raw"] = g2g["x_name_n"].str.upper()
    g2g["v_raw"] = g2g["y_name_n"].str.upper()
    # 去自环
    g2g = g2g[g2g["u_raw"] != g2g["v_raw"]].copy()

    # 无向规范化并统计边权重
    g2g["u"] = g2g[["u_raw", "v_raw"]].min(axis=1)
    g2g["v"] = g2g[["u_raw", "v_raw"]].max(axis=1)
    g2g_dedup = (
        g2g.groupby(["u", "v", "rel_n"], as_index=False)
        .size()
        .rename(columns={"rel_n": "relation", "size": "weight"})
    )

    frontier = set(seed_genes)
    visited_genes = set(seed_genes)
    hop_chunks = []

    for hop in range(1, HOPS + 1):
        mask = g2g_dedup["u"].isin(frontier) | g2g_dedup["v"].isin(frontier)
        sub_e = g2g_dedup[mask].copy()
        sub_e["hop"] = hop
        hop_chunks.append(sub_e)

        next_g = set(sub_e["u"]).union(set(sub_e["v"])) - visited_genes
        visited_genes.update(next_g)
        frontier = next_g
        print(f"[Layer 1] Hop {hop}: 新增基因 {len(next_g)}，累计基因 {len(visited_genes)}")

    l1_edges = (
        pd.concat(hop_chunks, ignore_index=True)
        .sort_values("hop")
        .drop_duplicates(subset=["u", "v", "relation"], keep="first")
    )

    G1 = nx.Graph()
    for _, r in l1_edges.iterrows():
        if G1.has_edge(r["u"], r["v"]):
            G1[r["u"]][r["v"]]["weight"] += int(r["weight"])
        else:
            G1.add_edge(r["u"], r["v"], weight=int(r["weight"]), relation=r["relation"], hop=int(r["hop"]))

    asd_meta = seed_df.drop_duplicates("gene_norm").set_index("gene_norm")
    l1_nodes_rows = []
    for g in sorted(G1.nodes()):
        meta = asd_meta.loc[g] if g in asd_meta.index else None
        l1_nodes_rows.append({
            "gene_symbol": g,
            "is_asd_seed": int(g in seed_genes),
            "sfari_score": meta["gene-score"] if meta is not None and "gene-score" in meta else "",
            "syndromic": meta["syndromic"] if meta is not None and "syndromic" in meta else "",
            "degree_l1": G1.degree(g),
        })
    l1_nodes = pd.DataFrame(l1_nodes_rows)

    l1_edges.to_csv(l1_dir / f"l1_gene_edges_{HOPS}hop.tsv", sep="\t", index=False)
    l1_nodes.to_csv(l1_dir / f"l1_gene_nodes_{HOPS}hop.tsv", sep="\t", index=False)

    # 合并 Layer 1 的 node 与 edge 为完整网络表
    l1_node_lookup = l1_nodes.set_index("gene_symbol")
    l1_merged = l1_edges.copy()
    l1_merged["x_type"] = "gene/protein"
    l1_merged["x_name"] = l1_merged["u"]
    l1_merged["y_type"] = "gene/protein"
    l1_merged["y_name"] = l1_merged["v"]
    for col in ["is_asd_seed", "sfari_score", "syndromic", "degree_l1"]:
        l1_merged[f"x_{col}"] = l1_merged["u"].map(l1_node_lookup[col])
        l1_merged[f"y_{col}"] = l1_merged["v"].map(l1_node_lookup[col])
    l1_merged.to_csv(l1_dir / f"l1_merged_network_{HOPS}hop.tsv", sep="\t", index=False)

    # ================= 5. 构建 Layer 2：基因–表型/疾病关联层 =================
    l1_gene_set = set(G1.nodes())

    # 筛选一端为 L1 基因、另一端为目标异质类型（如 phenotype / disease）的边
    cond_x_gene_y_other = (
        (kg_df["x_type_n"] == "gene/protein") &
        (kg_df["x_name_n"].str.upper().isin(l1_gene_set)) &
        (kg_df["y_type_n"].isin(INCLUDE_TYPES_L2)) &
        (kg_df["y_name_n"] != "")
    )
    cond_y_gene_x_other = (
        (kg_df["y_type_n"] == "gene/protein") &
        (kg_df["y_name_n"].str.upper().isin(l1_gene_set)) &
        (kg_df["x_type_n"].isin(INCLUDE_TYPES_L2)) &
        (kg_df["x_name_n"] != "")
    )

    part_a = kg_df[cond_x_gene_y_other].assign(
        gene_symbol=lambda d: d["x_name_n"].str.upper(),
        other_name=lambda d: d["y_name_n"],
        other_type=lambda d: d["y_type_n"],
    )
    part_b = kg_df[cond_y_gene_x_other].assign(
        gene_symbol=lambda d: d["y_name_n"].str.upper(),
        other_name=lambda d: d["x_name_n"],
        other_type=lambda d: d["x_type_n"],
    )

    l2_raw = pd.concat([
        part_a[["gene_symbol", "other_name", "other_type", "rel_n"]],
        part_b[["gene_symbol", "other_name", "other_type", "rel_n"]],
    ], ignore_index=True).rename(columns={"rel_n": "relation"})

    l2_edges = (
        l2_raw.groupby(["gene_symbol", "other_name", "other_type", "relation"], as_index=False)
        .size()
        .rename(columns={"size": "weight"})
    )

    # 统计每个非基因实体（如 HPO 表型）连接的基因数（用于后续参考节点度衰减高频表型）
    other_deg = (
        l2_edges.groupby(["other_name", "other_type"], as_index=False)["gene_symbol"]
        .nunique()
        .rename(columns={"gene_symbol": "linked_gene_count_in_l2"})
    )
    l2_edges = l2_edges.merge(other_deg, on=["other_name", "other_type"], how="left")

    l2_edges.to_csv(l2_dir / "l2_gene_to_pheno_disease_edges.tsv", sep="\t", index=False)
    other_deg.to_csv(l2_dir / "l2_pheno_disease_nodes.tsv", sep="\t", index=False)
    l1_node_lookup = l1_nodes.set_index("gene_symbol")
    l2_merged = l2_edges.copy()
    for col in ["is_asd_seed", "sfari_score", "syndromic", "degree_l1"]:
        l2_merged[f"gene_{col}"] = l2_merged["gene_symbol"].map(l1_node_lookup[col])
    l2_merged.to_csv(l2_dir / "l2_merged_network.tsv", sep="\t", index=False)

    print(f"[Layer 2] 导出关联边 {len(l2_edges)} 条，关联表型/疾病节点 {len(other_deg)} 个")
    print(f"[Done] 双层网络已保存至: {OUT_DIR}")


if __name__ == "__main__":
    main()