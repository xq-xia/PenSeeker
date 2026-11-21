import re
import requests
from dataclasses import dataclass
from typing import List, Optional

ENSEMBL_SERVER = "https://rest.ensembl.org"
HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}


@dataclass
class KozakHit:
    cds_pos_1based: int
    dist_downstream_bp: int
    kozak_window: str
    kozak_score: float
    is_strong_kozak: bool


def parse_snv(var_str: str) -> dict:
    m = re.match(r"^chr([0-9XYMT]+)-(\d+)-(\d+)-([ACGTN])-([ACGTN])$", var_str.strip(), re.IGNORECASE)
    if not m:
        raise ValueError(f"变异格式不符合预期: {var_str}")
    chrom, start, end, ref, alt = m.groups()
    return {
        "chrom": chrom.upper(),
        "start": int(start),
        "end": int(end),
        "ref": ref.upper(),
        "alt": alt.upper(),
    }


def fetch_vep_transcript_cds_pos(chrom: str, start: int, end: int, ref: str, alt: str) -> List[dict]:
    """调用 Ensembl VEP region 接口获取变异映射的转录本及 CDS 位置"""
    region = f"{chrom}:{start}-{end}"
    url = f"{ENSEMBL_SERVER}/vep/human/region/{region}/{alt}?content-type=application/json"
    resp = requests.get(url, headers=HEADERS, timeout=20)
    resp.raise_for_status()
    data = resp.json()
    if not data:
        return []

    hits = []
    for tc in data[0].get("transcript_consequences", []):
        if "cds_start" in tc:
            hits.append({
                "transcript_id": tc.get("transcript_id"),
                "gene_symbol": tc.get("gene_symbol"),
                "consequence_terms": tc.get("consequence_terms", []),
                "cds_start": int(tc["cds_start"]),
                "cds_end": int(tc.get("cds_end", tc["cds_start"])),
                "is_canonical": tc.get("canonical", 0) == 1,
            })
    return hits


def fetch_transcript_cds_seq(transcript_id: str) -> Optional[str]:
    """拉取指定转录本的 5'->3' CDS 序列"""
    url = f"{ENSEMBL_SERVER}/sequence/id/{transcript_id}?type=cds;content-type=application/json"
    resp = requests.get(url, headers=HEADERS, timeout=20)
    if resp.status_code != 200:
        return None
    return resp.json().get("seq")


def score_kozak(flank_seq: str) -> tuple[float, bool]:
    seq = flank_seq.upper()
    if len(seq) != 13 or seq[6:9] != "ATG":
        return 0.0, False

    base_minus3 = seq[3]
    base_plus4 = seq[9]

    score = 0.0
    if base_minus3 in ("A", "G"):
        score += 0.45
    if base_plus4 == "G":
        score += 0.35

    minor_weights = {
        0: ("G", 0.03),
        1: ("C", 0.03),
        2: ("C", 0.03),
        4: ("C", 0.03),
        5: ("C", 0.04),
        10: ("C", 0.02),
        11: ("T", 0.02),
    }
    for idx, (expected_base, w) in minor_weights.items():
        if seq[idx] == expected_base:
            score += w

    is_strong = (base_minus3 in ("A", "G")) and (base_plus4 == "G") and (score >= 0.80)
    return round(min(score, 1.0), 4), is_strong


def scan_downstream_kozak_atg(cds_seq: str, var_cds_pos_1based: int, window_bp: int = 100) -> dict:
    cds = cds_seq.upper().replace("U", "T")
    var_idx = var_cds_pos_1based - 1
    if var_idx < 0 or var_idx >= len(cds):
        return {
            "has_inframe_atg": False,
            "has_strong_kozak_rescue": False,
            "best_hit": None,
        }

    scan_start = var_idx + 1
    scan_end = min(len(cds), var_idx + 1 + window_bp)

    hits: List[KozakHit] = []
    for idx in range(scan_start, scan_end - 2):
        if (idx - var_idx) % 3 != 0:
            continue
        if cds[idx : idx + 3] == "ATG":
            left_start = idx - 6
            right_end = idx + 7
            if left_start >= 0 and right_end <= len(cds):
                flank = cds[left_start:right_end]
                k_score, is_strong = score_kozak(flank)
            else:
                flank = cds[max(0, left_start) : min(len(cds), right_end)]
                k_score, is_strong = 0.0, False

            hits.append(
                KozakHit(
                    cds_pos_1based=idx + 1,
                    dist_downstream_bp=idx - var_idx,
                    kozak_window=flank,
                    kozak_score=k_score,
                    is_strong_kozak=is_strong,
                )
            )

    strong_hits = [h for h in hits if h.is_strong_kozak]
    best_hit = strong_hits[0] if strong_hits else None
    return {
        "has_inframe_atg": len(hits) > 0,
        "has_strong_kozak_rescue": best_hit is not None,
        "best_hit": best_hit,
    }


def evaluate_single_variant(var_str: str, window_bp: int = 100):
    parsed = parse_snv(var_str)
    print(f"变异: {var_str} -> {parsed}")

    tx_hits = fetch_vep_transcript_cds_pos(
        parsed["chrom"], parsed["start"], parsed["end"], parsed["ref"], parsed["alt"]
    )
    if not tx_hits:
        print("该变异未落在任何转录本的 CDS 区域内（可能位于内含子、UTR 或基因间区），不触发 CDS 翻译重启动减分。")
        return

    for tx in tx_hits:
        tid = tx["transcript_id"]
        cds_seq = fetch_transcript_cds_seq(tid)
        if not cds_seq:
            continue
        res = scan_downstream_kozak_atg(cds_seq, tx["cds_start"], window_bp=window_bp)
        print(
            f"- 转录本 {tid} ({tx['gene_symbol']}, canonical={tx['is_canonical']}) | "
            f"CDS位置={tx['cds_start']} | 后果={','.join(tx['consequence_terms'])} | "
            f"下游同帧ATG={res['has_inframe_atg']} | 强Kozak重启动={res['has_strong_kozak_rescue']} | "
            f"最佳位点={res['best_hit']}"
        )


if __name__ == "__main__":
    evaluate_single_variant("chr20-63413494-63413494-G-A", window_bp=100)