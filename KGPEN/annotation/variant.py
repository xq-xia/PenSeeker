from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any, Callable

import pandas as pd


BASE_DIR = Path(r"E:\外显率预测课题\基因背景网络\KGPEN")
DEFAULT_INPUT_PATH = BASE_DIR / "data_source" / "variant" / "variant_list.xlsx"
DEFAULT_OUTPUT_PATH = BASE_DIR / "annotation" / "variant_annotation_results.xlsx"


ANNOTATION_COLUMNS: list[tuple[str, str]] = [
    ("protein_expression_impact", "变异对蛋白表达量的影响"),
    ("cold_hotspot_region", "变异冷热点区域"),
    ("active_site_localization", "活性中心定位"),
    ("functional_element_domain_overlap", "关键功能元件与结构域重叠"),
    ("protein_folding_ddg", "蛋白质折叠自由能变化"),
    ("grantham_distance", "氨基酸理化性质改变距离（Grantham 距离）"),
    ("splice_site_disruption", "剪切位点破坏预测"),
    ("nmd_trigger_probability", "NMD 触发概率"),
    ("transcript_impact_ratio", "转录本影响比率与主要转录本覆盖"),
]


EXTRA_RESULT_COLUMNS: dict[str, list[str]] = {}


@dataclass(slots=True)
class Variant:
    variant_id: str | None = None
    chrom: str | None = None
    pos: int | None = None
    ref: str | None = None
    alt: str | None = None
    gene_symbol: str | None = None
    transcript_id: str | None = None
    protein_change: str | None = None
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class AnnotationResult:
    name: str
    value: Any = None
    score: float | None = None
    evidence: list[str] = field(default_factory=list)
    details: dict[str, Any] = field(default_factory=dict)

    def to_cell(self) -> Any:
        if self.value is not None:
            return self.value
        if self.score is not None or self.evidence or self.details:
            return asdict(self)
        return None


class VariantAnnotation:
    def __init__(self, variant: Variant, **context: Any) -> None:
        self.variant = variant
        self.context = context
        self.results: dict[str, Any] = {}

    def annotate_protein_expression_impact(self) -> AnnotationResult | None:
        '''
        通过 UKB 的 pQTL 数据计算变异对蛋白表达量的影响大小.
        '''
        return None

    def annotate_cold_hotspot_region(self) -> AnnotationResult | None:
        '''
        计算该变异物理位置是否落于 HCSeeker 计算的变异热点或冷点区域.
        '''
        return None

    def annotate_active_site_localization(self) -> AnnotationResult | None:
        '''
        基于 InterPro、Pfam 数据库功能域注释判断变异是否位于蛋白质活性中心.
        '''
        return None

    def annotate_functional_element_domain_overlap(self) -> AnnotationResult | None:
        '''
        检索 Pfam/InterPro 结构域坐标与已知功能元件，判断变异是否位于核心结构域或关键功能元件内.
        '''
        return None

    def annotate_protein_folding_ddg(self) -> AnnotationResult | None:
        '''
        基于 AlphaFold2 预测结构运行 FoldX，计算变异前后折叠自由能变化量 ΔΔG.
        '''
        return None

    def annotate_grantham_distance(self) -> AnnotationResult | None:
        '''
        查表获取变异前后氨基酸的 Grantham Distance，用于衡量理化性质改变幅度.
        '''
        return None

    def annotate_splice_site_disruption(self) -> AnnotationResult | None:
        '''
        提取 SpliceAI 预测评分，判断变异是否破坏经典或隐蔽剪接位点.
        '''
        return None

    def annotate_nmd_trigger_probability(self) -> AnnotationResult | None:
        '''
        计算变异位置是否满足 NMD 触发规则，以判断是否高概率触发无义介导的 mRNA 降解.
        '''
        return None

    def annotate_transcript_impact_ratio(self) -> AnnotationResult | None:
        '''
        比对 GENCODE/MANE 数据库，计算变异对主要转录本及受影响转录本表达覆盖比例的影响.
        '''
        return None

    def annotate_all(self) -> dict[str, Any]:
        methods: dict[str, Callable[[], AnnotationResult | None]] = {
            key: getattr(self, f"annotate_{key}") for key, _ in ANNOTATION_COLUMNS
        }
        self.results = {name: method() for name, method in methods.items()}
        return self.results

    def set_context(self, **context: Any) -> None:
        self.context.update(context)

    def get_result(self, name: str) -> Any:
        return self.results.get(name)


@lru_cache(maxsize=1)
def load_variant_rows(input_path_str: str) -> list[Variant]:
    input_path = Path(input_path_str)
    if not input_path.exists():
        raise FileNotFoundError(f"Variant list file not found: {input_path}")

    suffix = input_path.suffix.lower()
    if suffix in {".xlsx", ".xls"}:
        df = pd.read_excel(input_path)
    elif suffix == ".csv":
        df = pd.read_csv(input_path)
    elif suffix == ".tsv":
        df = pd.read_csv(input_path, sep="\t")
    else:
        raise ValueError(f"Unsupported input format: {input_path.suffix}")

    if "variant_id" not in df.columns:
        df["variant_id"] = None
    if "chrom" not in df.columns:
        df["chrom"] = None
    if "pos" not in df.columns:
        df["pos"] = None
    if "ref" not in df.columns:
        df["ref"] = None
    if "alt" not in df.columns:
        df["alt"] = None
    if "gene_symbol" not in df.columns:
        df["gene_symbol"] = None
    if "transcript_id" not in df.columns:
        df["transcript_id"] = None
    if "protein_change" not in df.columns:
        df["protein_change"] = None

    variants: list[Variant] = []
    for _, row in df.iterrows():
        variant_id = None if pd.isna(row["variant_id"]) else str(row["variant_id"]).strip()
        chrom = None if pd.isna(row["chrom"]) else str(row["chrom"]).strip()
        pos = None if pd.isna(row["pos"]) else int(row["pos"])
        ref = None if pd.isna(row["ref"]) else str(row["ref"]).strip()
        alt = None if pd.isna(row["alt"]) else str(row["alt"]).strip()
        gene_symbol = None if pd.isna(row["gene_symbol"]) else str(row["gene_symbol"]).strip()
        transcript_id = None if pd.isna(row["transcript_id"]) else str(row["transcript_id"]).strip()
        protein_change = None if pd.isna(row["protein_change"]) else str(row["protein_change"]).strip()

        variants.append(
            Variant(
                variant_id=variant_id or None,
                chrom=chrom or None,
                pos=pos,
                ref=ref or None,
                alt=alt or None,
                gene_symbol=gene_symbol or None,
                transcript_id=transcript_id or None,
                protein_change=protein_change or None,
                metadata={},
            )
        )
    return variants


@lru_cache(maxsize=128)
def parse_variant_string(variant_string: str) -> Variant:
    parts = variant_string.strip().split("-")
    if len(parts) != 5:
        raise ValueError(
            "Variant string must be in the format chr-pos-pos-ref-alt, for example chr20-63413494-63413494-G-A"
        )

    chrom, start, end, ref, alt = parts
    if start != end:
        raise ValueError("This test parser expects start and end to be identical for SNV input")

    return Variant(
        variant_id=variant_string,
        chrom=chrom,
        pos=int(start),
        ref=ref,
        alt=alt,
        metadata={},
    )





def serialize_result(result: Any) -> Any:
    if result is None:
        return None
    if isinstance(result, AnnotationResult):
        return result.to_cell()
    if isinstance(result, dict):
        return result
    if isinstance(result, (list, tuple, set)):
        return list(result)
    return result


def annotate_variants(variants: list[Variant], **context: Any) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []

    for variant in variants:
        annotator = VariantAnnotation(variant, **context)
        annotation_map = annotator.annotate_all()

        row: dict[str, Any] = {
            "variant_id": variant.variant_id,
            "chrom": variant.chrom,
            "pos": variant.pos,
            "ref": variant.ref,
            "alt": variant.alt,
            "gene_symbol": variant.gene_symbol,
            "transcript_id": variant.transcript_id,
            "protein_change": variant.protein_change,
        }
        for column_name, _ in ANNOTATION_COLUMNS:
            result = serialize_result(annotation_map.get(column_name))
            extra_columns = EXTRA_RESULT_COLUMNS.get(column_name)
            if extra_columns and isinstance(result, dict):
                for extra_column in extra_columns:
                    row[extra_column] = result.get(extra_column)
            else:
                row[column_name] = result
        rows.append(row)

    columns = [
        "variant_id",
        "chrom",
        "pos",
        "ref",
        "alt",
        "gene_symbol",
        "transcript_id",
        "protein_change",
    ]
    
    for column_name, _ in ANNOTATION_COLUMNS:
        extra_columns = EXTRA_RESULT_COLUMNS.get(column_name)
        if extra_columns:
            columns.extend(extra_columns)
        else:
            columns.append(column_name)

    return pd.DataFrame(rows, columns=columns)


def save_annotations(df: pd.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if output_path.suffix.lower() == ".csv":
        df.to_csv(output_path, index=False, encoding="utf-8-sig")
    else:
        df.to_excel(output_path, index=False)


def run_pipeline(input_path: Path = DEFAULT_INPUT_PATH, output_path: Path = DEFAULT_OUTPUT_PATH) -> Path:
    variants = load_variant_rows(str(input_path))
    result_df = annotate_variants(variants)
    save_annotations(result_df, output_path)
    return output_path


def run_single_variant(variant_string: str) -> pd.DataFrame:
    variant = parse_variant_string(variant_string)
    result_df = annotate_variants([variant])
    return result_df


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run variant annotation on a variant list.")
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT_PATH,
        help="Path to the variant list file (xlsx/csv/tsv).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT_PATH,
        help="Path to the annotation output file.",
    )
    parser.add_argument(
        "--variant",
        type=str,
        default=None,
        help="Annotate a single variant string like chr20-63413494-63413494-G-A.",
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if args.variant:
        result_df = run_single_variant(args.variant)
        print(result_df.to_string(index=False))
        return 0

    if not args.input.exists():
        raise SystemExit(
            f"Input file not found: {args.input}. Use --variant chr20-63413494-63413494-G-A to test a single variant, or pass --input to a valid file."
        )

    output_path = run_pipeline(args.input, args.output)
    print(f"Annotation file written to: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
