# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.6.0
#   kernelspec:
#     display_name: Python [conda env:mouse_hema_meth_py37] *
#     language: python
#     name: conda-env-mouse_hema_meth_py37-py
# ---

# # Setup

n_cores = 16

# +
# isort: off
import os

num_threads = str(n_cores)

# these need to be set prior to numpy import
os.environ["OMP_NUM_THREADS"] = num_threads
os.environ["OPENBLAS_NUM_THREADS"] = num_threads
os.environ["MKL_NUM_THREADS"] = num_threads
os.environ["VECLIB_MAXIMUM_THREADS"] = num_threads
os.environ["NUMEXPR_NUM_THREADS"] = num_threads

import numpy as np

# isort: on

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from dendritic_cells.config import paper_context
from IPython.display import Markdown, display
from pandas.api.types import CategoricalDtype

import mouse_hema_meth.utils as ut
# -

import dendritic_cells.analyses_for_manuscript_20201015_lib as lib

results_dir = "/icgc/dkfzlsdf/analysis/hs_ontogeny/notebook-data/gNs4xcMJscaLLwlt"
os.makedirs(results_dir, exist_ok=True)


def link_fn(s, markdown=True):
    s = str(s)
    curry_path = s.replace(
        "/icgc/dkfzlsdf/analysis/hs_ontogeny/",
        "https://currywurst.dkfz.de/hs-ontogeny/",
    )
    if markdown:
        link = f'[{curry_path}]({curry_path} "will only work within the DKFZ network")'
        # noinspection PyTypeChecker
        display(Markdown(link))
    else:
        print(curry_path)


# %matplotlib inline

# # Retrieve and format input data

# Get cluster ids and dmr boundaries to TSV format, load into python, adjust column names and convert Chromosome to string categorical

cluster_ids_rds = results_dir + "/nW05rWZyjo0wFEVK.rds"
cluster_ids_tsv = results_dir + "/nW05rWZyjo0wFEVK.tsv"
cluster_ids_rds

# + language="R" active=""
# cluster_ids = readRDS(
#     '/icgc/dkfzlsdf/analysis/hs_ontogeny/notebook-data/gNs4xcMJscaLLwlt/nW05rWZyjo0wFEVK.rds'
# )
# write.table(
#     cluster_ids,
#     file = "/icgc/dkfzlsdf/analysis/hs_ontogeny/notebook-data/gNs4xcMJscaLLwlt/nW05rWZyjo0wFEVK.tsv",
#     sep = "\t",
#     row.names = F,
#     col.names = T,
#     quote=F
# )
# -

cluster_ids_beta_values_dmr_coords = pd.read_csv(cluster_ids_tsv, sep="\t")

chrom_dtype = CategoricalDtype(
    np.sort(np.arange(1, 20).astype(str)),
    ordered=True,
)
chrom_dtype

cluster_ids_beta_values_dmr_coords["chr"] = (
    cluster_ids_beta_values_dmr_coords["chr"].astype(str).astype(chrom_dtype)
)

cluster_ids_beta_values_dmr_coords = (
    cluster_ids_beta_values_dmr_coords.rename(
        columns={"chr": "Chromosome", "start": "Start", "end": "End"}
    )
    .sort_values(["Chromosome", "Start", "End"])
    .reset_index(drop=True)
    .rename_axis(index="region_id")
)

cluster_ids_beta_values_dmr_coords

dmrs = cluster_ids_beta_values_dmr_coords[["Chromosome", "Start", "End"]]
dmrs

cluster_ids_df = cluster_ids_beta_values_dmr_coords[["cluster"]]
cluster_ids_df

# # Gene annotation

# ## Computation

gtf_fp = (
    "/icgc/dkfzlsdf/analysis/hs_ontogeny/databases/gene_annotations"
    "/gencode.vM19.annotation.no-prefix.gtf"
)

anno_res_paths_d = lib.run_gtfanno(
    granges_df=dmrs,
    gtf_fp=gtf_fp,
    output_trunk_path=results_dir + "/dmr-gene-annos",
    promoter_def=(-5000, 1000),
    distant_cis_regulatory_domain_def=(-50000, -5000),
    recompute=True,
)

# ## Paths to detailed annotations

anno_res_paths_d

for k, v in anno_res_paths_d.items():
    print(k)
    link_fn(v)

# ## Paths to merged annotations

# One row per DMR, one feature class per DMR (one of promoter, exon, etc.)

merged_annos = lib.merge_annos(
    gtfanno_result_fp=anno_res_paths_d["primary_annos_p"],
    grange_and_feature_ids=dmrs,
)
merged_annos.to_pickle(results_dir + "/merged-gene-annos.p")
merged_annos.to_csv(results_dir + "/merged-gene-annos.tsv", sep="\t", index=False)

merged_annos

print(results_dir + "/merged-gene-annos.p")
print(results_dir + "/merged-gene-annos.tsv")
link_fn(results_dir + "/merged-gene-annos.tsv")

# ## Genomic regions quick viz

# note about UTRs: i double checked the UTR calling; and I have set 5'-UTR precedence higher than promoters, so this is not an effect of covering 5'-UTRs with promoters. The picture is very similar to global hematopoiesis in general.

merged_annos['feat_class'].value_counts().sort_values().plot.bar()

# # Geneset enrichment

# ## Geneset inspection (readme!)

# **Note1**
# There are duplicates within the genesets; beware of this in other analyses too!

# **Note2**
# The cardinality of the genesets is in part quite high; it may be beneficial for the analysis or the interpretation of the analysis if this is repeated with more curated genelists. Such as it is, each cluster has a high number of gene overlaps, which may be more difficult to interprete. (See section at the end of the notebook)

# **Note3**
# A considerable number of genes occurs in multiple clusters; is this intended?

dc_genesets_csv = (
    "/home/kraemers/projects/dendritic_cells/local/genesets/concat_genesets.csv"
)

rosenbauer_genesets_df = pd.read_csv(dc_genesets_csv, sep="\t", header=0)

# ### Geneset contains duplicates (link to cleaned genesets at the end)

# The genelists contain multiple mentions of the same genes, this has to be corrected and should be relayed back to the Rosenbauer lab, since this could also affect other analyses, if these genelists would be used elsewherejf.
# Number of mentions of unique genes per main geneset (eg 20 genes in coeff4_down are mentioned twice)

rosenbauer_genesets_df.apply(
    lambda ser: ser.dropna().str.upper().value_counts().value_counts(), axis=0
)

# coeff4_down detailed

with pd.option_context("display.min_rows", 200):
    # pd.options.display.max_rows=200
    display(
        rosenbauer_genesets_df["coeff4_down"]
        .value_counts()
        .to_frame("frequency")
        .query("frequency > 1")
    )

# remove duplicates

rosenbauer_genesets_no_duplicates = rosenbauer_genesets_df.mask(
    rosenbauer_genesets_df.apply(
        lambda ser: ser.str.upper().duplicated(keep="first"), axis=0
    )
)

# verify

rosenbauer_genesets_no_duplicates.apply(
    lambda ser: ser.dropna().str.upper().value_counts().value_counts(), axis=0
)

rosenbauer_genesets_no_duplicates

rosenbauer_genesets_no_duplicates_tsv = (
    results_dir + "/rosenbauer-genesets_no-duplicates.tsv"
)
rosenbauer_genesets_no_duplicates.to_csv(
    rosenbauer_genesets_no_duplicates_tsv, sep="\t"
)

print(rosenbauer_genesets_no_duplicates_tsv)
link_fn(rosenbauer_genesets_no_duplicates_tsv)

# ### Genesets are large

# Number of unique genes

rosenbauer_genesets_no_duplicates.stack().nunique()

rosenbauer_genesets_no_duplicates.notnull().sum().sort_values().plot.bar()

# ### Genes occur in multiple clusters

# Example: just the non cell type specific clusters

n_cluster_per_gene = (
    rosenbauer_genesets_no_duplicates
    # filter for non cell type specific genesets
    .filter(regex=r"(down|up)$")
    .stack()
    .value_counts()
    .sort_values(ascending=False)
)

# Histogram of #clusters per gene

n_cluster_per_gene.value_counts().plot.bar()

with pd.option_context("display.min_rows", 50, "display.max_rows", 50):
    display(n_cluster_per_gene.head(30))

# ## Computation (Report with figures at the end of the section)

dmr_geneset_enrichments_output_dir = results_dir + "/dmr_geneset_enrichments_results"
dmr_geneset_enrichments_report_dir = (
    results_dir + "/dmr_geneset_enrichments_results/report"
)
# os.makedirs(dmr_geneset_enrichments_report_dir, exist_ok=True)

geneset_databases_d = lib.create_gmt_for_rosenbauer_genesets(
    rosenbauer_genesets_no_duplicates
)
geneset_databases_d

for fp in geneset_databases_d.values():
    assert os.path.exists(fp)

# currently need to clean output dir because all results there will be in report

# + active=""
# !rm -r {dmr_geneset_enrichments_output_dir}
# !rm -r {dmr_geneset_enrichments_report_dir}
# -

lib.run_geneset_enrichment_analysis(
    merged_gene_annos=merged_annos,
    cluster_ids=cluster_ids_df,
    geneset_databases_d=geneset_databases_d,
    output_dir=dmr_geneset_enrichments_output_dir,
    report_dir=dmr_geneset_enrichments_report_dir,
    max_pvalues=(0.05, 1),
    barcode_plot_args_d=dict(
        col_width_cm=0.5,
        row_height_cm=0.1,
        linewidth=0.5,
        vmin_quantile=0.05,
        divergent_cmap="RdYlGn_r",
        cluster_features=True,
    ),
    feature_annos=dmrs,
    n_cores=n_cores,
    recompute=True,
    filters=["promoter", "gene_regions", "all_annotated"],
    additional_formats=("pdf", "svg"),
    # vlim=(-5, 5),
)

# ## Overlap counts, overlapping genes

overlap_stats_pattern = (
    dmr_geneset_enrichments_output_dir + "/{anno_name}/{database}/overlap-stats.p"
)
cluster_overlap_stats_pattern = (
    dmr_geneset_enrichments_output_dir
    + "/{anno_name}/{database}/{clustering}.{filter}/cluster-overlap-stats.p"
)

# ### Table: dmrs vs hit in geneset (boolean)

hits_table = (overlap_stats_pattern[:-2] + "_hits.p").format(
    anno_name="gtfanno", database="rosenbauer_genesets_all_gmt"
)
hits_table_df = pd.read_pickle(hits_table)
display(hits_table_df.query("coeff1_up == 1").head())

hits_table_tsv = hits_table.replace("_hits.p", "_hits.tsv")
hits_table_df.to_csv(hits_table_tsv, sep="\t", index=True)
hits_table_tsv
link_fn(hits_table_tsv)

# ### Tables: clusters vs number of hits in cluster

# **Note** the caveats on the provided genesets (see above)

# Note that one gene is often hit by several DMRs, eg some DCRD, and hits in the promoter, introns etc.

cluster_sizes = cluster_ids_df["cluster"].value_counts()
cluster_sizes.sort_values()

cluster_overlap_counts_by_filter_tsv = (
    cluster_overlap_stats_pattern.format(
        anno_name="gtfanno",
        database="rosenbauer_genesets_all_gmt",
        clustering="cluster",
        filter="{filter}",
    )[:-2]
    + "_hits.tsv"
)
for filter_name in ["promoter", "gene_regions", "all_annotated"]:
    curr_cluster_overlap_counts_tsv = cluster_overlap_counts_by_filter_tsv.format(
        filter=filter_name
    )
    print(filter_name)
    print("Total # Hits")
    overlap_counts = pd.read_csv(curr_cluster_overlap_counts_tsv, sep="\t").set_axis(
        np.arange(1, 10)
    )
    display(overlap_counts)
    print("% Hits")
    display(
        overlap_counts.divide(
            cluster_sizes, axis=0
        ) * 100,
    )
    print(curr_cluster_overlap_counts_tsv)
    link_fn(curr_cluster_overlap_counts_tsv)

# ### Annotation: clusters vs marker genes found

# Genes found in clusters, filtered by gene body, promoter, or all regions

dc_genesets_csv = (
    "/home/kraemers/projects/dendritic_cells/local/genesets/concat_genesets.csv"
)
genesets_d_ser = {}
genesets_d_ser["rosenbauer_all"] = pd.read_csv(
    dc_genesets_csv, sep="\t", header=0
).T.agg(lambda ser: ser.str.upper().dropna().tolist(), axis=1)
genesets_d_ser["rosenbauer_coeff1and2"] = pd.read_csv(
    dc_genesets_csv,
    usecols=["coeff1_up", "coeff1_down", "coeff2_up", "coeff2_down"],
    sep="\t",
    header=0,
).T.agg(lambda ser: ser.str.upper().dropna().tolist(), axis=1)

cluster_marker_genes_df = lib.cluster_marker_genes_overview(
    genesets_d_ser=genesets_d_ser,
    primary_gene_annos=pd.read_pickle(anno_res_paths_d["primary_annos_p"]),
    granges_df=cluster_ids_beta_values_dmr_coords,
    cluster_ids=cluster_ids_df,
)
cluster_marker_genes_df.head()

cluster_marker_genes_tsv = results_dir + "/cluster-marker-genes.tsv"
cluster_marker_genes_df.to_csv(cluster_marker_genes_tsv, sep="\t")
cluster_marker_genes_tsv
link_fn(cluster_marker_genes_tsv)
