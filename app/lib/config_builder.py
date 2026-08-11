"""UI 参数 → pyseqrna run.ini。"""
from __future__ import annotations
from pathlib import Path

TEMPLATE = """[General]
input_file = {sample_sheet}
samples_path = {fastq_dir}
reference_genome = {genome}
feature_file = {gtf}
outdir = {outdir}
dryrun = False
force = False
resume_policy = skip
paired = {paired}

[Species]
source = ENSEMBL
species = {species}
organism_type = animals

[Quality]
skip_quality = False
quality_tool = fastqc
quality_trim = {quality_trim}
skip_trim = {skip_trim}
trimming_tool = trim_galore

[Alignment]
skip_alignment = False
alignment_tool = {alignment_tool}
alignment_stats = True
alignment_stats_source = auto

[Quantification]
skip_quantification = False
quant_method = genomic_overlaps
run_multimapped_groups = False

[Normalization]
skip_normalization = False
normalization_method = rpkm
skip_normalization_plots = False

[Clustering]
run_clustering = True
cluster_target = samples
cluster_method = hierarchical
cluster_metric = euclidean
cluster_linkage = average
cluster_top_variable = 1000
cluster_scale = row
cluster_no_log = False
cluster_no_heatmap = False
cluster_cmap = vlag

[Coexpression]
run_coexpression = False

[DifferentialExpression]
skip_diffexp = False
diffexp_tool = {diffexp_tool}
diffexp_normalization = median_ratio
diffexp_abundance = base_mean
diffexp_dispersion = map
diffexp_test = wald
fdr_threshold = {fdr_threshold}
fold_threshold = {fold_threshold}
pvalue_threshold = {pvalue_threshold}
add_gene_names = True
subset = False

[Visualization]
pca_plot = True
tsne_plot = True
volcano_plot = True
ma_plot = True
deg_heatmap = True
heatmap_top_genes = 50
venn = True
upset = True
venn_comparisons =
venn_label = updown

[FunctionalAnnotation]
skip_functional_annotation = False
gene_ontology = {gene_ontology}
kegg_pathway = {kegg_pathway}
go_pvalue_threshold = {pvalue_threshold}
kegg_pvalue_threshold = {pvalue_threshold}

[Report]
skip_report = False
report_formats = html,md,json
report_title = RNA-seq Analysis Report

[Computational]
threads = {threads}
memory = {memory}
local_jobs = 1
resume = all

[SLURM]
slurm = False
"""


def build_ini(params: dict) -> str:
    """params 见测试 BASE；species 仅限 hsapiens/mmusculus。"""
    species = params["species"]
    if species not in ("hsapiens", "mmusculus"):
        raise ValueError(f"不支持的物种: {species}")
    return TEMPLATE.format(
        sample_sheet=str(Path(params["sample_sheet"]).resolve()),
        fastq_dir=str(Path(params["fastq_dir"]).resolve()),
        genome=str(Path(params["genome"]).resolve()),
        gtf=str(Path(params["gtf"]).resolve()),
        outdir=str(Path(params["outdir"]).resolve()),
        diffexp_tool=params.get("diffexp_tool", "deseq2"),
        species=species,
        alignment_tool=params.get("alignment_tool", "hisat2"),
        paired="True" if params.get("paired", True) else "False",
        skip_trim="True" if params.get("skip_trim", False) else "False",
        quality_trim="False" if params.get("skip_trim", False) else "True",
        fdr_threshold=params.get("fdr_threshold", 0.05),
        fold_threshold=params["fold_threshold"],
        pvalue_threshold=params["pvalue_threshold"],
        gene_ontology="True" if params.get("gene_ontology", True) else "False",
        kegg_pathway="True" if params.get("kegg_pathway", True) else "False",
        threads=params["threads"],
        memory=params["memory"],
    )
