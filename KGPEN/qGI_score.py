from pathlib import Path
import numpy as np
import pandas as pd

# ================= 1. 路径与参数配置 =================
QGI_FILE = Path(r"E:\外显率预测课题\基因背景网络\KGPEN\data\qGI_score.xlsx")
L1_EDGES_FILE = Path(r"E:\外显率预测课题\基因背景网络\KGPEN\data\layer1_gene_bg\l1_gene_edges_2hop.tsv")
L1_NODES_FILE = Path(r"E:\外显率预测课题\基因背景网络\KGPEN\data\layer1_gene_bg\l1_gene_nodes_2hop.tsv")
OUT_FILE = Path(r"E:\外显率预测课题\基因背景网络\KGPEN\data\layer1_gene_bg\l1_gene_edges_2hop_with_qgi.tsv")
CANDIDATE_EXTRA_FILE = Path(r"E:\外显率预测课题\基因背景网络\KGPEN\data\layer1_gene_bg\qgi_extra_edges_in_l1_nodes.tsv")

ONLY_STANDARD_GI = True    # True: 仅保留 GI_standard 非空 (|qGI|>0.3, FDR<0.1) 的高置信记录
ALPHA = 0.5                # qGI 强度对基础权重的调节系数（0 表示不调制）
ADD_EXTRA_TO_MAIN = False  # False: 保持原 KG 边结构不变；True: 仅将 L1 内部高置信 qGI 边并入主表


def norm_str(s: pd.Series) -> pd.Series:
    return s.fillna("").astype(str).str.strip()


def read_auto(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xls"}:
        df = pd.read_excel(path, dtype=str)
    else:
        with open(path, "rb") as fb:
            header_bytes = fb.read(4)
        if header_bytes.startswith(b"PK\x03\x04"):
            df = pd.read_excel(path, dtype=str)
        else:
            with open(path, "r", encoding="utf-8-sig", errors="replace") as f:
                first_line = f.readline()
            sep = "\t" if first_line.count("\t") >= first_line.count(",") else ","
            df = pd.read_csv(path, sep=sep, dtype=str, encoding="utf-8-sig")
    df.columns = [str(c).strip() for c in df.columns]
    return df


def main():
    l1_edges = read_auto(L1_EDGES_FILE)
    l1_nodes = read_auto(L1_NODES_FILE)
    qgi_raw = read_auto(QGI_FILE)

    l1_edges["u"] = norm_str(l1_edges["u"]).str.upper()
    l1_edges["v"] = norm_str(l1_edges["v"]).str.upper()
    l1_edges["weight_num"] = pd.to_numeric(l1_edges["weight"], errors="coerce").fillna(1.0)

    qgi_raw["g1"] = norm_str(qgi_raw["library_gene"]).str.upper()
    qgi_raw["g2"] = norm_str(qgi_raw["query_gene"]).str.upper()
    qgi_raw["qgi"] = pd.to_numeric(qgi_raw["qGI_score"], errors="coerce")
    qgi_raw = qgi_raw.dropna(subset=["qgi"])
    qgi_raw = qgi_raw[(qgi_raw["g1"] != "") & (qgi_raw["g2"] != "") & (qgi_raw["g1"] != qgi_raw["g2"])].copy()

    std_col = next((c for c in qgi_raw.columns if c.startswith("GI_standard")), None)
    if ONLY_STANDARD_GI and std_col is not None:
        qgi_raw = qgi_raw[norm_str(qgi_raw[std_col]) != ""].copy()

    # 无向基因对规范化
    qgi_raw["u"] = qgi_raw[["g1", "g2"]].min(axis=1)
    qgi_raw["v"] = qgi_raw[["g1", "g2"]].max(axis=1)
    qgi_raw["abs_qgi"] = qgi_raw["qgi"].abs()
    qgi_raw["neg_qgi"] = np.where(qgi_raw["qgi"] < 0, -qgi_raw["qgi"], 0.0)
    qgi_raw["pos_qgi"] = np.where(qgi_raw["qgi"] > 0, qgi_raw["qgi"], 0.0)

    # 聚合同一基因对在多组 screen 中的 qGI 信号
    qgi_pair = (
        qgi_raw.groupby(["u", "v"], as_index=False)
        .agg(
            qgi_n_screens=("qgi", "count"),
            qgi_mean=("qgi", "mean"),
            qgi_max_abs=("abs_qgi", "max"),
            qgi_mean_abs=("abs_qgi", "mean"),
            qgi_neg_max=("neg_qgi", "max"),
            qgi_pos_max=("pos_qgi", "max"),
        )
    )

    # 1. 合并到现有 L1 边表
    merged = l1_edges.merge(qgi_pair, on=["u", "v"], how="left")
    merged["has_qgi"] = merged["qgi_mean"].notna().astype(int)

    for col in ["qgi_n_screens", "qgi_mean", "qgi_max_abs", "qgi_mean_abs", "qgi_neg_max", "qgi_pos_max"]:
        merged[col] = merged[col].fillna(0.0)

    max_abs_ref = merged["qgi_max_abs"].max()
    norm_abs = (merged["qgi_max_abs"] / max_abs_ref) if max_abs_ref > 0 else 0.0

    base_w = np.log1p(merged["weight_num"])
    merged["weight_qgi_mod"] = base_w * (1.0 + ALPHA * norm_abs)

    # 2. 找出“两端都在 L1 基因集中、原 KG 无边、但有 qGI”的基因对（供消融/补边参考）
    l1_gene_set = set(norm_str(l1_nodes["gene_symbol"]).str.upper())
    existing_pairs = set(zip(l1_edges["u"], l1_edges["v"]))

    extra_qgi = qgi_pair[
        qgi_pair["u"].isin(l1_gene_set) & qgi_pair["v"].isin(l1_gene_set)
    ].copy()
    keep_mask = [not (pair in existing_pairs) for pair in zip(extra_qgi["u"], extra_qgi["v"])]
    extra_qgi = extra_qgi[keep_mask].copy()
    extra_qgi.to_csv(CANDIDATE_EXTRA_FILE, sep="\t", index=False)

    if ADD_EXTRA_TO_MAIN and not extra_qgi.empty:
        extra_rows = pd.DataFrame({
            "u": extra_qgi["u"],
            "v": extra_qgi["v"],
            "relation": "qgi_only",
            "weight": 1,
            "hop": -1,
            "weight_num": 1.0,
            "qgi_n_screens": extra_qgi["qgi_n_screens"],
            "qgi_mean": extra_qgi["qgi_mean"],
            "qgi_max_abs": extra_qgi["qgi_max_abs"],
            "qgi_mean_abs": extra_qgi["qgi_mean_abs"],
            "qgi_neg_max": extra_qgi["qgi_neg_max"],
            "qgi_pos_max": extra_qgi["qgi_pos_max"],
            "has_qgi": 1,
        })
        extra_norm_abs = (extra_rows["qgi_max_abs"] / max_abs_ref) if max_abs_ref > 0 else 0.0
        extra_rows["weight_qgi_mod"] = np.log1p(extra_rows["weight_num"]) * (1.0 + ALPHA * extra_norm_abs)
        merged = pd.concat([merged, extra_rows], ignore_index=True)

    merged.to_csv(OUT_FILE, sep="\t", index=False)
    print(f"[Done] L1 边数: {len(merged)}，匹配到 qGI 的边数: {int(merged['has_qgi'].sum())}")
    print(f"[Info] L1 内 KG 无边但有 qGI 的候选对数: {len(extra_qgi)} (已存至 {CANDIDATE_EXTRA_FILE})")


if __name__ == "__main__":
    main()