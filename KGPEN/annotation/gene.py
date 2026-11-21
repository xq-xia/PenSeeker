from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any, Callable

import pandas as pd


BASE_DIR = Path(r"E:\外显率预测课题\基因背景网络\KGPEN")
DEFAULT_INPUT_PATH = BASE_DIR / "annotation" / "data_source" / "gene" / "gene_list.xlsx"
DEFAULT_OUTPUT_PATH = BASE_DIR / "annotation" / "gene_annotation_results.xlsx"
DEFAULT_INTRAGENIC_COMPLEMENTATION_PATH = (
    BASE_DIR / "annotation" / "data_source" / "gene" / "intragenic_complementation.xlsx"
)
DEFAULT_PARALOG_COMPENSATION_PATH = (
    BASE_DIR
    / "annotation"
    / "data_source"
    / "gene"
    / "gene_paralog"
    / "genes_with_brain_paralog_compensation.txt"
)
DEFAULT_PLI_LOEUF_PATH = (
    BASE_DIR / "annotation" / "data_source" / "gene" / "gnomad_v4.1_mane_NM_pli_loeuf.tsv"
)

ANNOTATION_COLUMNS: list[tuple[str, str]] = [
    ("cis_trans_regulation", "顺式/反式调控修饰"),
    ("paralog_compensation", "旁系同源代偿机制"),
    ("intragenic_complementation", "基因内互补"),
    ("translation_reinitiation", "翻译重启动"),
    ("splicing_rescue", "剪接救援评分"),
    ("pli_loeuf", "概率性单倍体不足评分（pLI / LOEUF）"),
    ("gene_epistasis", "基因间上位性修饰"),
    ("polygenic_risk_background", "多基因风险评分与背景负荷"),
    ("promoter_methylation_imprinting", "启动子甲基化与基因印迹状态"),
    ("histone_and_allele_specific_expression", "组蛋白修饰与单等位基因随机表达"),
    ("ncrna_and_3d_genome", "非编码 RNA 调控与 3D 基因组结构"),
]

EXTRA_RESULT_COLUMNS: dict[str, list[str]] = {
    "pli_loeuf": ["pli_loeuf_recorded", "pli_score", "loeuf_score"],
}


@dataclass(slots=True)
class Gene:
    gene_symbol: str
    gene_id: str | None = None
    chromosome: str | None = None
    start: int | None = None
    end: int | None = None
    strand: str | None = None
    transcript_id: str | None = None
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


@lru_cache(maxsize=1)
def load_paralog_compensation_gene_set(lookup_path_str: str) -> set[str]:
    lookup_path = Path(lookup_path_str)
    if not lookup_path.exists():
        raise FileNotFoundError(f"Paralog compensation file not found: {lookup_path}")

    genes: set[str] = set()
    with lookup_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            gene_symbol = line.strip()
            if gene_symbol:
                genes.add(gene_symbol)
    return genes


@lru_cache(maxsize=1)
def load_intragenic_complementation_gene_set(lookup_path_str: str) -> set[str]:
    lookup_path = Path(lookup_path_str)
    if not lookup_path.exists():
        raise FileNotFoundError(f"Intragenic complementation file not found: {lookup_path}")

    df = pd.read_excel(lookup_path, usecols=[0])
    first_column = df.columns[0]
    return {
        str(value).strip()
        for value in df[first_column].dropna().tolist()
        if str(value).strip()
    }


@lru_cache(maxsize=1)
def load_pli_loeuf_map(lookup_path_str: str) -> dict[str, tuple[float, float]]:
    lookup_path = Path(lookup_path_str)
    if not lookup_path.exists():
        raise FileNotFoundError(f"pLI/LOEUF file not found: {lookup_path}")

    df = pd.read_csv(lookup_path, sep="\t")
    if df.shape[1] < 3:
        raise ValueError(f"pLI/LOEUF file must have at least 3 columns: {lookup_path}")

    gene_column = df.columns[0]
    pli_column = df.columns[-2]
    loeuf_column = df.columns[-1]

    result: dict[str, tuple[float, float]] = {}
    for _, row in df[[gene_column, pli_column, loeuf_column]].dropna(subset=[gene_column]).iterrows():
        gene_symbol = str(row[gene_column]).strip()
        if not gene_symbol:
            continue

        pli_value = row[pli_column]
        loeuf_value = row[loeuf_column]
        pli = float(pli_value) if pd.notna(pli_value) else 0.0
        loeuf = float(loeuf_value) if pd.notna(loeuf_value) else 0.0
        result[gene_symbol] = (pli, loeuf)

    return result


class GeneAnnotation:
    def __init__(self, gene: Gene, **context: Any) -> None:
        self.gene = gene
        self.context = context
        self.results: dict[str, Any] = {}

    def annotate_cis_trans_regulation(self) -> AnnotationResult | None:
        return None

    def annotate_paralog_compensation(self) -> AnnotationResult | None:
        '''
        查询 Ensembl BioMart 旁系同源树，判断是否存在序列相似度 > 50%且 GTEx 脑组织高表达(>=10)的同源基因.
        '''
        lookup_path = Path(
            self.context.get(
                "paralog_compensation_lookup_path",
                DEFAULT_PARALOG_COMPENSATION_PATH,
            )
        )
        gene_set = load_paralog_compensation_gene_set(str(lookup_path))
        value = 1 if self.gene.gene_symbol in gene_set else 0
        return AnnotationResult(name="paralog_compensation", value=value)

    def annotate_intragenic_complementation(self) -> AnnotationResult | None:
        '''
        查询变异所处基因是否位于PNAS 2026这篇文章提供的Dataset S7中具有基因内互补机制的基因列表中.
        '''
        lookup_path = Path(
            self.context.get(
                "intragenic_complementation_lookup_path",
                DEFAULT_INTRAGENIC_COMPLEMENTATION_PATH,
            )
        )
        gene_set = load_intragenic_complementation_gene_set(str(lookup_path))
        value = 1 if self.gene.gene_symbol in gene_set else 0
        return AnnotationResult(name="intragenic_complementation", value=value)

    def annotate_translation_reinitiation(self) -> AnnotationResult | None:
        return None

    def annotate_splicing_rescue(self) -> AnnotationResult | None:
        return None

    def annotate_pli_loeuf(self) -> AnnotationResult | None:
        '''
        提取 gnomAD 数据库中的 pLI 及 LOEUF 官方评价指标.
        '''
        lookup_path = Path(
            self.context.get(
                "pli_loeuf_lookup_path",
                DEFAULT_PLI_LOEUF_PATH,
            )
        )
        pli_loeuf_map = load_pli_loeuf_map(str(lookup_path))

        if self.gene.gene_symbol in pli_loeuf_map:
            pli, loeuf = pli_loeuf_map[self.gene.gene_symbol]
            return AnnotationResult(
                name="pli_loeuf",
                value={
                    "pli_loeuf_recorded": 1,
                    "pli_score": pli,
                    "loeuf_score": loeuf,
                },
            )

        return AnnotationResult(
            name="pli_loeuf",
            value={
                "pli_loeuf_recorded": 0,
                "pli_score": 0,
                "loeuf_score": 0,
            },
        )

    def annotate_gene_epistasis(self) -> AnnotationResult | None:
        return None

    def annotate_polygenic_risk_background(self) -> AnnotationResult | None:
        return None

    def annotate_promoter_methylation_imprinting(self) -> AnnotationResult | None:
        return None

    def annotate_histone_and_allele_specific_expression(self) -> AnnotationResult | None:
        return None

    def annotate_ncrna_and_3d_genome(self) -> AnnotationResult | None:
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


def load_gene_list(input_path: Path) -> list[Gene]:
    df = pd.read_excel(input_path, usecols=["gene_symbol"])
    if "gene_symbol" not in df.columns:
        raise ValueError(f"Input file must contain a 'gene_symbol' column: {input_path}")

    genes: list[Gene] = []
    for raw_symbol in df["gene_symbol"].dropna().tolist():
        gene_symbol = str(raw_symbol).strip()
        if not gene_symbol:
            continue
        genes.append(Gene(gene_symbol=gene_symbol))
    return genes


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


def annotate_genes(genes: list[Gene], **context: Any) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []

    for gene in genes:
        annotator = GeneAnnotation(gene, **context)
        annotation_map = annotator.annotate_all()

        row: dict[str, Any] = {"gene_symbol": gene.gene_symbol}
        for column_name, _ in ANNOTATION_COLUMNS:
            result = serialize_result(annotation_map.get(column_name))
            extra_columns = EXTRA_RESULT_COLUMNS.get(column_name)
            if extra_columns and isinstance(result, dict):
                for extra_column in extra_columns:
                    row[extra_column] = result.get(extra_column)
            else:
                row[column_name] = result
        rows.append(row)

    columns = ["gene_symbol"]
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
    genes = load_gene_list(input_path)
    result_df = annotate_genes(genes)
    save_annotations(result_df, output_path)
    return output_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run gene annotation on a gene list.")
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT_PATH,
        help="Path to the gene list Excel file.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT_PATH,
        help="Path to the annotation output file.",
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    output_path = run_pipeline(args.input, args.output)
    print(f"Annotation file written to: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
