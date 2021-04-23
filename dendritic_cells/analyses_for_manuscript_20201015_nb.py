# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.7.0
#   kernelspec:
#     display_name: Python [conda env:mouse_hema_meth_py37_mamba_full] *
#     language: python
#     name: conda-env-mouse_hema_meth_py37_mamba_full-py
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
import seaborn as sns
from dendritic_cells.config import paper_context2
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


mpl.rcParams.update(paper_context2)

# %matplotlib inline

# ## Recompute flag

recompute = False

# # QC

# ## Flagstats

# +
extended_flagstats_pattern = (
    "/icgc/dkfzlsdf/project/mouse_hematopoiesis/sequencing"
    "/whole_genome_bisulfite_tagmentation_sequencing/view-by-pid"
    "/{pop}_{rep}/blood/paired/merged-alignment/qualitycontrol/merged/flagstats"
    "/blood_{pop}_{rep}_merged.mdup.bam_extendedFlagstats.txt"
)
original_flagstats_pattern = (
    "/icgc/dkfzlsdf/project/mouse_hematopoiesis/sequencing"
    "/whole_genome_bisulfite_tagmentation_sequencing/view-by-pid"
    "/{pop}_{rep}/blood/paired/merged-alignment/qualitycontrol/merged/flagstats"
    "/blood_{pop}_{rep}_merged.mdup.bam_flagstats.txt"
)
file_pattern = original_flagstats_pattern

flagstats_df_full = lib.collect_original_flagstats(
    file_pattern=original_flagstats_pattern
)
display(flagstats_df_full.head())
flagstats_df_full.reset_index().to_csv(results_dir + "/flagstats-full.csv")
flagstats_df_full.to_pickle(results_dir + "/flagstats-full.p")
print(results_dir + "/flagstats-full.csv")
link_fn(results_dir + "/flagstats-full.csv")
# -

flagstats_selected = (
    flagstats_df_full.loc[
        ["HSC", "MDP", "CDP", "cMoP", "Monocytes", "pDC", "cDC CD11b", "cDC CD8a"],
        "QC-passed reads",
    ]
    .unstack("Stat")
    .drop([("HSC", "1-1"), ("HSC", "4"), ("HSC", "5")], axis=0)
    .sort_index(ascending=False, axis=1)
)
display(flagstats_selected.head())
flagstats_selected.reset_index().to_csv(
    results_dir + "/flagstats-selected.csv", index=True
)
link_fn(results_dir + "/flagstats-selected.csv")
print(results_dir + "/flagstats-selected.csv")

# ## Metadata table

# +
from mouse_hema_meth.methylome.alignments_mcalls.meth_calling_paths import (
    ds1_metadata_table_tsv,
)

# compared to ds1, we
#  - restrict to the pops used in the papr
#  - drop one HSC replicate
#  - use monos-new instead of monos

pops = ["cdp", "cmop", "dc-cd11b", "dc-cd8a", "hsc", "mdp", "monos", "pdc"]
mcalls_metadata_table = (
    pd.read_csv(ds1_metadata_table_tsv, sep="\t")
    .set_index(["subject", "rep"])
    .loc[pops]
    .drop(("hsc", 4), axis=0)
)

mcalls_metadata_table.loc["monos", "bed_path"] = (
    mcalls_metadata_table.loc["monos", "bed_path"]
    .str.replace("monos", "monos-new")
    .to_numpy()
)
mcalls_metadata_table.loc["monos", "pickle_path"] = (
    mcalls_metadata_table.loc["monos", "pickle_path"]
    .str.replace("monos", "monos-new")
    .to_numpy()
)
# -

# ## Conversion rates

chh_meth = lib.gather_conversion_rate_table(mcalls_metadata_table)
display(chh_meth.head())
chh_meth.to_csv(results_dir + "/conversion-rates.csv")
print(results_dir + "/conversion-rates.csv")
link_fn(results_dir + "/conversion-rates.csv")

# ## Coverage

# computed in hema meth project
pops = ["cdp", "cmop", "dc-cd11b", "dc-cd8a", "hsc", "mdp", "monos", "pdc"]
pop_coverage_p = "/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/cohort_results/analyses/hierarchy/mcalling/qc/coverage/n-total_pop-level.p"
pop_ntotal = pd.read_pickle(pop_coverage_p)[pops]
pop_ntotal

# Restrict coverage to the three already published HSC replicates

n_total = pd.Series(0, index=pop_ntotal.index, dtype='i8')
for fp in mcalls_metadata_table.loc['hsc', 'pickle_path'].tolist():
    n_total += pd.read_pickle(fp)['n_total']
pop_ntotal['hsc'] = n_total

# Update monos to use the old monos

n_total = pd.Series(0, index=pop_ntotal.index, dtype='i8')
for fp in mcalls_metadata_table.loc['monos', 'bed_path'].tolist():
    print(fp)
    n_total += pd.read_csv(fp, sep='\t')['n_total']
    print(n_total.mean())
pop_ntotal['monos'] = n_total

coverage_mean_sd = pop_ntotal.describe().T
coverage_mean_sd

coverage_mean_sd.index.name = 'Population'
coverage_mean_sd_tsv = results_dir + '/coverage-mean-sd2.csv'
coverage_mean_sd.reset_index().to_csv(coverage_mean_sd_tsv, sep='\t')
print(coverage_mean_sd_tsv)
link_fn(coverage_mean_sd_tsv)

# ## Average methylation

(
    pop_quantiles_df,
    pop_quantiles_df_long,
    sampled_rep_meth_df,
    sampled_pop_meth_df,
    beta_value_pop_df,
    pop_beta_value_mean_df,
    rep_beta_value_mean_df,
) = lib.gather_meth_level_data(
    mcalls_metadata_table=mcalls_metadata_table.reset_index(),
    recompute=recompute,
    pop_quantiles_wide_p=results_dir + "/pop_beta-value_quantiles_wide.p",
    pop_quantiles_long_p=results_dir + "/pop_beta-value_quantiles_long.p",
    sampled_rep_betavalue_p=results_dir + "/pop_beta-value_sampled_rep.p",
    sampled_pop_betavalue_p=results_dir + "/pop_beta-value_sampled_pop.p",
    pop_betavalue_p=results_dir + "/pop_beta-values_all-cpgs.p",
    pop_betavalue_mean_df_p=results_dir + "/pop_beta-values_mean-sd.p",
    rep_beta_value_mean_p=results_dir + "/rep_beta-values_mean-sd.p",
)

rep_beta_value_mean_df_long = (
    rep_beta_value_mean_df.unstack()
    .rename_axis(["Population", "Replicate"], axis=0)
    .reset_index()
)
rep_beta_value_mean_df_long

rep_beta_value_mean_df_long.to_csv(results_dir + "/rep-beta-values_mean-sd.csv")
link_fn(results_dir + "/rep-beta-values_mean-sd.csv")

# # Retrieve and format input data

# Get cluster ids and dmr boundaries to TSV format, load into python, adjust column names and convert Chromosome to string categorical

dmr_calls_dir = (
    "/icgc/dkfzlsdf/analysis/hs_ontogeny/notebook-data/gNs4xcMJscaLLwlt/dmr-calls"
)
dmr_calls_dir

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

cluster_ids_beta_values_dmr_coords[
    [
        "hsc",
        "cdp",
        "cmop",
        "dc_cd11b",
        "dc_cd8a",
        "mono",
        "mdp",
        "pdc",
    ]
].isnull().sum()

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

# ## Meth levels in DMRs

beta_values = cluster_ids_beta_values_dmr_coords[
    [
        "hsc",
        "mdp",
        "cdp",
        "cmop",
        "mono",
        "pdc",
        "dc_cd11b",
        "dc_cd8a",
    ]
].dropna(how="any", axis=0)
zscores = ut.row_zscores(beta_values)
cluster_ids_no_nas = cluster_ids_df.loc[beta_values.index].copy()

beta_values.describe().reset_index().to_csv(results_dir + '/dmr-beta-value_agg-stats.tsv')
ut.dkfz_link(results_dir + '/dmr-beta-value_agg-stats.tsv')

fig, ax = plt.subplots(1, 1, figsize=(ut.cm(8), ut.cm(5)), constrained_layout=True, dpi=180)
sns.violinplot(data=beta_values.sample(1000), cut=0, inner='box', width=0.9)
ax.set_ylabel('DMR methylation')
ut.save_and_display(fig, png_path=results_dir + '/global-dmr-beta-values_violin.png')

fig, ax = plt.subplots(1, 1, figsize=(ut.cm(3), ut.cm(5)), constrained_layout=True, dpi=180)
sns.violinplot(data=beta_values[['pdc', 'dc_cd11b', 'dc_cd8a']].sample(1000), cut=0, inner='box', width=0.9)
ax.set_ylabel('DMR methylation')
ut.save_and_display(fig, png_path=results_dir + '/global-dmr-beta-values_violin_dcs-only.png')

fig, ax = plt.subplots(1, 1, figsize=(ut.cm(8), ut.cm(8)), constrained_layout=True, dpi=180)
sns.boxplot(data=beta_values.sample(1000), showfliers=False)
ut.save_and_display(fig, png_path=results_dir + '/global-dmr-beta-values_boxplot.png')

# # Characterize clustering

beta_values.isnull().sum()


# ## Heatmap

import clustering_tools as ct
import mouse_hema_meth.utils as ut

global_row_order_df = ct.hclust_within_clusters(
    data=zscores,
    cluster_ids_df=cluster_ids_no_nas,
    metric="euclidean",
    method="ward",
    n_jobs=9,
)

png_path = results_dir + "/heatmap_zscores_hclust-within-clusters.png"
lib.create_clustermaps_with_special_order(
    df=zscores,
    cluster_ids_df=cluster_ids_no_nas,
    png_path=png_path,
    cmap="RdBu_r",
    guide_title="Z-score(% Methylation)",
    global_row_order_df=global_row_order_df,
    figsize=(10 / 2.54, 10 / 2.54),
    n_per_cluster=500,
)
link_fn(png_path)
link_fn(png_path.replace(".png", ".pdf"))

png_path = results_dir + "/heatmap_beta-values_hclust-within-clusters.png"
lib.create_clustermaps_with_special_order(
    df=beta_values,
    cluster_ids_df=cluster_ids_no_nas,
    png_path=png_path,
    cmap="YlOrBr",
    guide_title="% Methylation",
    global_row_order_df=global_row_order_df,
    figsize=(10 / 2.54, 10 / 2.54),
    n_per_cluster=500,
)
link_fn(png_path)
link_fn(png_path.replace(".png", ".pdf"))

# ## Cluster sizes

fig, ax = plt.subplots(
    1, 1, constrained_layout=True, dpi=180, figsize=(5 / 2.54, 5 / 2.54)
)
cluster_sizes = cluster_ids_df["cluster"].value_counts()
display(cluster_sizes)
cluster_sizes.plot.bar(ax=ax)
ax.set(ylabel="Frequency", xlabel="Cluster ID")
fig.savefig(results_dir + "/cluster-sizes.pdf")
link_fn(results_dir + "/cluster-sizes.pdf")

# ## Gain/loss stratified counts

# not all dmrs have only one sign (deltas < 0.1 are not considered here)

masked_signs = np.sign(
    beta_values.subtract(beta_values["hsc"], axis=0).mask(lambda df: df.abs().le(0.1))
)

masked_signs.nunique(axis=1).ne(1).sum()

# we can use the max delta

max_delta_for_each_dmr = beta_values.subtract(beta_values["hsc"], axis=0).apply(
    lambda ser: ser.loc[ser.abs().idxmax()], axis=1
)
max_delta_for_each_dmr

gain_loss_classif = np.sign(max_delta_for_each_dmr).replace({1: "Gain", -1: "Loss"})
gain_loss_classif

gain_loss_classif.value_counts()

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
    recompute=False,
)

# ## Paths to detailed annotations

anno_res_paths_d

for k, v in anno_res_paths_d.items():
    print(k)
    link_fn(v)

# ## Paths to merged annotations

# One row per DMR, one feature class per DMR (one of promoter, exon, etc.)

if recompute:
    merged_annos = lib.merge_annos(
        gtfanno_result_fp=anno_res_paths_d["primary_annos_p"],
        grange_and_feature_ids=dmrs,
    )
    merged_annos.to_pickle(results_dir + "/merged-gene-annos.p")
    merged_annos.to_csv(results_dir + "/merged-gene-annos.tsv", sep="\t", index=False)
else:
    merged_annos = pd.read_pickle(results_dir + "/merged-gene-annos.p")

merged_annos

print(results_dir + "/merged-gene-annos.p")
print(results_dir + "/merged-gene-annos.tsv")
link_fn(results_dir + "/merged-gene-annos.tsv")

# ## Genomic regions quick viz

# note about UTRs: i double checked the UTR calling; and I have set 5'-UTR precedence higher than promoters, so this is not an effect of covering 5'-UTRs with promoters. The picture is very similar to global hematopoiesis in general.

merged_annos["feat_class"].value_counts().sort_values().plot.bar()

# # Concatenate all available annotations

all_annos = pd.concat(
    [cluster_ids_beta_values_dmr_coords, merged_annos[["gene_name", "feat_class"]]],
    axis=1,
)[
    [
        "Chromosome",
        "Start",
        "End",
        "cluster",
        "gene_name",
        "feat_class",
        "hsc",
        "cdp",
        "cmop",
        "dc_cd11b",
        "dc_cd8a",
        "mono",
        "mdp",
        "pdc",
    ]
]
all_annos_tsv = results_dir + "/all-annos.tsv"
all_annos.to_csv(all_annos_tsv, sep="\t", index=False, header=True)
all_annos.to_parquet(all_annos_tsv.replace(".tsv", ".parquet"))
all_annos.head()

print(all_annos_tsv)
link_fn(all_annos_tsv)

# ## Promoter DMR annos only

promoter_genes = all_annos.query('feat_class == "Promoter"').sort_values(
    ["cluster", "gene_name"]
)
promoter_genes_tsv = results_dir + "/promoter_genes.tsv"
promoter_genes.to_csv(promoter_genes_tsv, sep="\t", index=False, header=True)
promoter_genes.to_parquet(promoter_genes_tsv.replace(".tsv", ".parquet"))
promoter_genes.head()

print(promoter_genes_tsv)
link_fn(promoter_genes_tsv)

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
        overlap_counts.divide(cluster_sizes, axis=0) * 100,
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

# # Pairwise DMR analysis

# %load_ext rpy2.ipython

# +
import rpy2
import rpy2.robjects as ro
from rpy2.robjects.packages import importr

utils = importr("utils")
base = importr("base")
# -

import rpy2.robjects.pandas2ri as pandas2ri

pw_dmr_calls = pd.DataFrame(
    {"rds_path": dmr_calls_dir + "/" + pd.Series(os.listdir(dmr_calls_dir))}
)
pw_dmr_calls["pop"] = pw_dmr_calls["rds_path"].str.extract(
    r".*_dmrs_hsc_vs_([\w-]+)_0.01", expand=False
)
pw_dmr_calls

gain_loss_counts = pd.DataFrame(-1, columns=["Gain", "Loss"], index=pw_dmr_calls["pop"])
for _unused, row_ser in pw_dmr_calls.iterrows():
    # chr start end  length   nCG  meanMethy1  meanMethy2 diff.Methy  areaStat
    gain_loss_for_pop = (
        np.sign(
            pandas2ri.rpy2py(base.readRDS(row_ser["rds_path"])).eval(
                "meanMethy2 - meanMethy1"
            )
        )
        .value_counts()
        .sort_index()
        .set_axis(["Loss", "Gain"])
    )
    gain_loss_counts.loc[row_ser["pop"]] = gain_loss_for_pop
gain_loss_counts = gain_loss_counts.sort_values("Loss")

gain_loss_counts.head()

pw_counts_plot_df = (
    gain_loss_counts.stack()
    .reset_index()
    .set_axis(["Population", "Direction", "No. of DMRs"], axis=1)
)
pw_counts_plot_df

fig, ax = plt.subplots(
    1, 1, dpi=300, constrained_layout=True, figsize=(4 / 2.54, 3 / 2.54)
)
sns.barplot(
    x="Population",
    y="No. of DMRs",
    hue="Direction",
    palette={"Loss": "Blue", "Gain": "Red"},
    data=pw_counts_plot_df,
    order=gain_loss_counts.index,
)
ax.legend(loc="upper left", bbox_to_anchor=(1.05, 1))
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_yticks(np.arange(0, 40_000, 5_000))
ax.set_ylim(0, 32_000)
ax.set(xlabel="")
fig.savefig(results_dir + "gain-loss-pw-counts-barplot.pdf")
link_fn(results_dir + "gain-loss-pw-counts-barplot.pdf")
