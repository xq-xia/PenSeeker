from pathlib import Path
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

# ================= 1. 路径与参数 =================
BASE_DIR = Path(r"E:\外显率预测课题\基因背景网络\KGPEN\data")
L1_MERGED_FILE = BASE_DIR / "layer1_gene_bg" / "l1_merged_network_2hop.tsv"
L2_EDGES_FILE = BASE_DIR / "layer2_pheno_disease" / "l2_gene_to_pheno_disease_edges.tsv"
OUT_PNG = BASE_DIR / "vis_asd_local_subnetwork.png"

TARGET_GENES = ["CHD8", "SHANK3"]   # 聚焦查看的核心基因
TOP_L2_PER_GENE = 8                 # 每个目标基因最多展示多少条 Layer 2 关联边（防标签重叠）

plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial"]
plt.rcParams["axes.unicode_minus"] = False


def main():
    l1_df = pd.read_csv(L1_MERGED_FILE, sep="\t", dtype=str)
    targets = {g.strip().upper() for g in TARGET_GENES}

    # 筛选与目标基因直接相连的 Layer 1 边
    sub_l1 = l1_df[
        l1_df["x_name"].str.upper().isin(targets) | l1_df["y_name"].str.upper().isin(targets)
    ].copy()
    sub_l1["weight_num"] = pd.to_numeric(sub_l1["weight"], errors="coerce").fillna(1)
    # 边太多时可只看权重较高的前 60 条，保证图面清晰
    sub_l1 = sub_l1.sort_values("weight_num", ascending=False).head(60)

    G = nx.Graph()
    for _, r in sub_l1.iterrows():
        u = r["x_name"].strip().upper()
        v = r["y_name"].strip().upper()
        G.add_node(u, node_kind="seed" if u in targets else "gene")
        G.add_node(v, node_kind="seed" if v in targets else "gene")
        G.add_edge(u, v, layer="L1", weight=r["weight_num"])

    # 可选：叠加少量 Layer 2 节点（表型/疾病/通路等）
    if L2_EDGES_FILE.exists():
        l2_df = pd.read_csv(L2_EDGES_FILE, sep="\t", dtype=str)
        l2_df["weight_num"] = pd.to_numeric(l2_df["weight"], errors="coerce").fillna(1)
        sub_l2 = (
            l2_df[l2_df["gene_symbol"].str.upper().isin(targets)]
            .sort_values("weight_num", ascending=False)
            .groupby("gene_symbol", as_index=False)
            .head(TOP_L2_PER_GENE)
        )
        for _, r in sub_l2.iterrows():
            g = r["gene_symbol"].strip().upper()
            other = f"{r['other_name']} [{r['other_type']}]"
            if g in G:
                G.add_node(other, node_kind="l2_other")
                G.add_edge(g, other, layer="L2", weight=r["weight_num"])

    pos = nx.spring_layout(G, k=0.8, seed=42)

    seed_nodes = [n for n, d in G.nodes(data=True) if d.get("node_kind") == "seed"]
    gene_nodes = [n for n, d in G.nodes(data=True) if d.get("node_kind") == "gene"]
    l2_nodes = [n for n, d in G.nodes(data=True) if d.get("node_kind") == "l2_other"]

    fig, ax = plt.subplots(figsize=(12, 9), dpi=300)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    l1_edges_list = [(u, v) for u, v, d in G.edges(data=True) if d.get("layer") == "L1"]
    l2_edges_list = [(u, v) for u, v, d in G.edges(data=True) if d.get("layer") == "L2"]

    nx.draw_networkx_edges(G, pos, edgelist=l1_edges_list, edge_color="#4A90E2", width=1.4, alpha=0.65, ax=ax)
    nx.draw_networkx_edges(G, pos, edgelist=l2_edges_list, edge_color="#D9534F", style="dashed", width=1.1, alpha=0.55, ax=ax)

    nx.draw_networkx_nodes(G, pos, nodelist=gene_nodes, node_color="#5BC0DE", node_size=260, edgecolors="white", linewidths=1.0, ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=l2_nodes, node_color="#F0AD4E", node_size=220, edgecolors="white", linewidths=1.0, ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=seed_nodes, node_color="#C9302C", node_size=520, edgecolors="black", linewidths=1.5, ax=ax)

    nx.draw_networkx_labels(G, pos, font_size=8, font_color="#111111", ax=ax)

    ax.set_title("ASD Background Subnetwork (Layer 1 Genes + Top Layer 2 Associations)", fontsize=13, pad=12)
    ax.axis("off")
    plt.tight_layout()
    plt.savefig(OUT_PNG, dpi=300)
    print(f"[Done] 已保存可视化图片至: {OUT_PNG}")


if __name__ == "__main__":
    main()