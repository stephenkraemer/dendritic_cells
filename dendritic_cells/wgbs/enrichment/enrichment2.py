# python36_general
import json
import os
import pickle
import tempfile
from copy import copy

import matplotlib as mpl
from region_set_profiler import barcode_heatmap

mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from pybedtools import BedTool
from itertools import product
from pathlib import Path
from typing import Dict, Any, List
import figure_report as fr
import numpy as np
import pandas as pd
import region_set_profiler as rsp
import snakemake

from mouse_hematopoiesis.config.base_paths import (
    wgbs_cohort_results_dir, project_temp_dir)
from mouse_hematopoiesis.utils.plotting_tools import paper_context
import matplotlib as mpl


# Functions
# ######################################################################


def merge_annos(gtfanno_result_fp):
    """Merge annotation format with multiple rows per DMR

    Args:
        gtfanno_result_fp: path to the gtfanno result

    Returns:
        pd.DataFrame with columns ['gene_name', 'feat_class']
        - Multiple gene annotations per DMR are concatenated with a comma.
        - All gene names will be uppercased.
        - index: Chromosome, Start, End
    """

    gtfanno_result: pd.DataFrame = pd.read_pickle(gtfanno_result_fp)
    gene_annos = (gtfanno_result
                  # add gtfanno_uid to keep it in the index, even though it is redundant
                  .groupby(['Chromosome', 'Start', 'End', 'gtfanno_uid'])['gene_name']
                  .aggregate(lambda ser: ser.str.cat(sep=','))
                  ).str.upper()
    feat_class = (gtfanno_result
                  .groupby(['Chromosome', 'Start', 'End', 'gtfanno_uid'])['feat_class']
                  .first()
                  )
    gene_annos = pd.concat([gene_annos, feat_class], axis=1)
    assert (gene_annos.index.get_level_values('gtfanno_uid')
            == np.arange(gene_annos.shape[0])).all()
    gene_annos.index = gene_annos.index.droplevel('gtfanno_uid')
    return gene_annos


def process_cluster_overlap_stats(
        cluster_overlap_stats, max_pvalues, plot_args,
        barcode_plot_png_by_maxpvalue, cluster_overlap_stats_out_fp,
        test_arg_dict=None, sample=None, cores=1, plots_only=False):

    if not plots_only:
        if sample is not None and sample < cluster_overlap_stats.hits.shape[1]:
            random_sel = np.random.choice(cluster_overlap_stats.hits.shape[1],
                                          sample, replace=False)
            cluster_overlap_stats.hits = cluster_overlap_stats.hits.iloc[:, random_sel]

        print('Calculate test per feature')
        final_test_args = dict(simulate_pval=True, replicate=int(1e4),
                               workspace=1_000_000)
        if test_arg_dict is not None:
            final_test_args.update(test_arg_dict)
        cluster_overlap_stats.test_per_feature(method='hybrid', cores=cores,
                                               test_args=final_test_args)
        print('Calculate test per cluster per feature')
        cluster_overlap_stats.test_per_cluster_per_feature()

        with open(cluster_overlap_stats_out_fp, 'wb') as fout:
            pickle.dump(cluster_overlap_stats, fout)
        cluster_overlap_stats.hits.to_pickle(
                cluster_overlap_stats_out_fp[:-2] + '_hits.p')
        cluster_overlap_stats.hits.to_csv(
                cluster_overlap_stats_out_fp[:-2] + '_hits.tsv',
                sep='\t', index=False)

    for max_pvalue in max_pvalues:
        print('Create barcode figure', 'max_pvalue', max_pvalue)
        with mpl.rc_context(paper_context):
            fig = rsp.barcode_heatmap(cluster_overlap_stats, **plot_args,
                                      max_pvalue=max_pvalue)
            out_png = barcode_plot_png_by_maxpvalue.format(max_pvalue=max_pvalue)
            fig.set_dpi(90)
            fig.savefig(out_png)
            fig.savefig(out_png.replace('.png', '.pdf'))

# copy-pasted and adapted from mouse_hematopoiesis.wgbs.clustering.clustering1.geneset_enrichments
# run_geneset_overlap_stats
# added barcode_plot_args
# don't pass max_pvalue explicitely through config anymore, instead pass a list of max_pvalues explicitely
# add pdf output to run_geneset_overlap_stats
def run_geneset_overlap_stats(config: Dict[str, Any],
                              max_pvalues: List[float],
                              barcode_plot_args: Dict[str, any],
                              plots_only=False):
    """
    Args:
        barcode_plot_args: any argument except max_pvalue
        max_pvalues: one plot per max_pvalue in the list

    test function within this function is parallelized with config['cores'] (?)
    # query regions
    # annotation
    # database
      # overlap stats
    # clustering
    # filter / merge
      # cluster overlap stats with test results
      # plot1
      # plot2
    {root_output_dir}/{query_name}/{database}.{clustering}.{filter_or_merge}.{test}/file-type.suffix
    """


    overlap_stats_pattern = config['output_dir'] + '/{anno_name}/{database}/overlap-stats.p'
    cluster_overlap_stats_pattern = config['output_dir'] + '/{anno_name}/{database}/{clustering}.{filter}/cluster-overlap-stats.p'
    barcode_plot_png_pattern = config['output_dir'] + '/{anno_name}/{database}/{clustering}.{filter}/barcode-plot_maxpvalue-{max_pvalue}.png'

    for ((anno_name, gene_annos), (database_name, database_fp),
         (clustering_name, (clustering_column_name, clustering)), (cleaner_name, cleaner)) in product(
            config['annotations'].items(),
            config['databases'].items(),
            config['clusterings'].items(),
            config['filters'].items()):
        print('Anno', anno_name)
        print('Database', database_name)
        print('Partitioning:', clustering_name)
        print('Filter', cleaner_name)

        # A nested folder structure is created. All folders can be created at once
        # by creating the parent directory of the barcode plot png.
        overlap_stats_p = overlap_stats_pattern.format(
                anno_name=anno_name,
                database=database_name,
        )
        cluster_overlap_stats_p = cluster_overlap_stats_pattern.format(
                anno_name=anno_name,
                database=database_name,
                clustering=clustering_name,
                filter=cleaner_name)
        barcode_plot_png_by_maxpvalue = barcode_plot_png_pattern.format(
                anno_name=anno_name,
                database=database_name,
                clustering=clustering_name,
                filter=cleaner_name,
                max_pvalue='{max_pvalue}')
        Path(barcode_plot_png_by_maxpvalue).parent.mkdir(parents=True, exist_ok=True)

        # Prepare cluster ids as series with Granges index
        if isinstance(clustering, str):
            clustering = pd.read_pickle(clustering)
        if isinstance(gene_annos, str):
            gene_annos = pd.read_pickle(gene_annos)
        cluster_ids_ser = clustering[clustering_column_name]
        cluster_ids_ser.index = (cluster_ids_ser.reset_index()
                                 [['Chromosome', 'Start', 'End']].set_index(['Chromosome', 'Start', 'End'])
                                 .index)

        if not plots_only:
            print('Calculate OverlapStats')
            overlap_stats = rsp.GenesetOverlapStats(
                    annotations=gene_annos, genesets_fp=database_fp)
            overlap_stats.compute()
            with open(overlap_stats_p, 'wb') as fout:
                pickle.dump(overlap_stats, fout)
            overlap_stats.coverage_df.to_pickle(overlap_stats_p[:-2] + '_hits.p')
            overlap_stats.coverage_df.rename(columns={'Chromosome': '#Chromosome'}).to_csv(overlap_stats_p[:-2] + '_hits.bed', sep='\t')

            print('Calculate ClusterOverlapStats')
            # regions cleaner not yet implemented
            if isinstance(cleaner, pd.MultiIndex):
                cluster_overlap_stats = overlap_stats.aggregate(
                        cluster_ids_ser, index=cleaner)
            else:
                cluster_overlap_stats = overlap_stats.aggregate(cluster_ids_ser)
        else:
            with open(cluster_overlap_stats_p, 'rb') as fin:
                cluster_overlap_stats = pickle.load(fin)

        process_cluster_overlap_stats(
                cluster_overlap_stats,
                max_pvalues=max_pvalues,
                plot_args=barcode_plot_args,
                barcode_plot_png_by_maxpvalue=barcode_plot_png_by_maxpvalue,
                cluster_overlap_stats_out_fp=cluster_overlap_stats_p,
                test_arg_dict=config['test_arg_dict'],
                sample=config['sample'],
                plots_only=plots_only)





# Enrichments
# ######################################################################


# General config
# **********************************************************************
enrichment_output_dir = (wgbs_cohort_results_dir
                         + '/analyses/dendritic_cells/dmr_characterization/enrichments2')
geneset_enrichment_output_dir = enrichment_output_dir + '/gene-sets'
regionset_enrichment_output_dir = enrichment_output_dir + '/region-sets'
enrichment_report_dir = enrichment_output_dir + '/reports'

chromosomes = sorted([str(i) for i in range(1, 20)])

tmpdir = tempfile.TemporaryDirectory(dir=project_temp_dir)

# Shared objects
# **********************************************************************

gtfanno_result_fp = (
        wgbs_cohort_results_dir
        + '/analyses/dendritic_cells/dmr_characterization/annotation/gtfanno/gtfanno'
        + '_primary-annotations.p')
merged_annos = merge_annos(gtfanno_result_fp)

all_cluster_ids_p = '/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/cohort_results/analyses/dendritic_cells/dmr_characterization/clustering/cluster-ids.p'
all_cluster_ids = pd.read_pickle(all_cluster_ids_p)
all_cluster_ids.index = all_cluster_ids.index.droplevel(3)
assert merged_annos.index.equals(all_cluster_ids.index)


# Genesets
# **********************************************************************
promoter_index = merged_annos.query('feat_class == "Promoter"').index
gene_region_index = merged_annos.query('feat_class not in ["intergenic", "DCRD"]').index
all_annotated_index = merged_annos.query('feat_class != "intergenic"').index

clusterings = {k: (k, all_cluster_ids) for k in ['ds2_nopam']}

# Create gmt for dc genesets
dc_genesets_csv = '/home/kraemers/projects/dendritic_cells/local/genesets/concat_genesets.csv'
dc_genesets_gmt = '/home/kraemers/projects/dendritic_cells/local/genesets/concat_genesets.gmt'
dc_genesets_df = pd.read_csv(dc_genesets_csv, sep='\t', header=0)
gmt_df = dc_genesets_df.T.apply(lambda ser: ser.str.upper())
gmt_df.insert(0, 'description', 'NA')
with open(dc_genesets_gmt, 'wt') as fout:
    for i, row in gmt_df.reset_index().iterrows():
        fout.write(row.dropna().str.cat(sep='\t') + '\n')

barcode_plot_args = dict(clusters_as_rows=True,
                         col_width_cm=0.5, row_height_cm=0.1,
                         linewidth=0.5, vlim=(-5, 5), vmin_quantile=0.02)

config = {
    'annotations':   {'gtfanno': merged_annos['gene_name']},
    'databases':     {'mouse-hematopoiesis-marker-genes': dc_genesets_gmt},
    'clusterings':   clusterings,
    'filters':       {'promoter': promoter_index,
                      'gene_regions': gene_region_index,
                      'all_annotated': all_annotated_index,
                      },
    'output_dir':    geneset_enrichment_output_dir,
    'sample':        None,
    'cores':         24,
    'test_arg_dict': dict(replicate=int(1e4), workspace=1_000_000),
}
run_geneset_overlap_stats(config, max_pvalues=[1e-3, 1e-2, 1],
                          barcode_plot_args=barcode_plot_args,
                          plots_only=False)

# for plot fine-tuning
# cluster_overlap_stats_p = '/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/cohort_results/analyses/hierarchy/enrichments/full-hierarchy/clustering-v2/genesets/gtfanno/lab mouse-hematopoiesis marker genes/hierarchy-v2_29-clusters-merged.Promoter/cluster-overlap-stats.p'
# cluster_overlap_stats = pd.read_pickle(cluster_overlap_stats_p)
#
# with mpl.rc_context(paper_context):
#     fig = rsp.barcode_heatmap(cluster_overlap_stats, clusters_as_rows=True,
#                               col_width_cm=0.5, row_height_cm=0.1, max_pvalue=1e-2,
#                               linewidth=0.5)
#     fig.savefig('/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/cohort_results/analyses/hierarchy/enrichments/full-hierarchy/clustering-v2/genesets/gtfanno/lab mouse-hematopoiesis marker genes/hierarchy-v2_29-clusters-merged.Promoter/cluster-overlap-stats.temp.png')


# Create geneset report
# ----------------------------------------------------------------------
os.makedirs(enrichment_report_dir, exist_ok=True)
geneset_enrichment_patterns = {
    'gene set enrichment': (geneset_enrichment_output_dir
                            + '/gtfanno/{database}/{clustering}.{filter}/barcode-plot_maxpvalue-{max_pvalue}.png'),
}

section_cols = ['database', 'clustering', 'filter', 'max_pvalue']
metadata_table = fr.pattern_set_to_metadata_table(
        geneset_enrichment_patterns, wildcard_constraints = {'clustering': r'.+(?=\.)'})
metadata_table = metadata_table.sort_values(section_cols)
# metadata_table.sort_values(section_cols, inplace=True)
fr.copy_report_files_to_report_dir(metadata_table,
                                   root_dir=enrichment_output_dir,
                                   report_dir=enrichment_report_dir)
report_config = fr.convert_metadata_table_to_report_json(
        metadata_table, section_cols)
report_config.update(dict(toc_headings='h1, h2, h3, h4', autocollapse_depth='3'))
fr.Report({'marker_gene_enrichments': report_config}).generate(enrichment_report_dir)



# Region set enrichments
# **********************************************************************

df = pd.read_csv('/icgc/dkfzlsdf/analysis/hs_ontogeny/databases/enrichment_databases/lola_chipseq_2018-04-12/mm10/codex/index.txt', sep='\t', header=0)
df['name'] = df['antibody'] + '_' + df['cellType'] + '_' + np.arange(0, len(df)).astype(str)
codex_all_metadata_table_fp = '/icgc/dkfzlsdf/analysis/hs_ontogeny/databases/enrichment_databases/lola_chipseq_2018-04-12/mm10/codex/regions/codex_annotations_all.csv'
df['abspath'] = '/icgc/dkfzlsdf/analysis/hs_ontogeny/databases/enrichment_databases/lola_chipseq_2018-04-12/mm10/codex/regions/' + df.filename
assert os.path.exists(df.iloc[0]['abspath'])
df.to_csv(codex_all_metadata_table_fp, sep='\t', header=True, index=False)

# Run workflow
# ======================================================================

# cluster_ids = {k: (k, all_cluster_ids_p) for k in all_cluster_ids.columns}
# cluster_ids = {'merged_29-clusters': ('merged_29-clusters', all_cluster_ids_p)}
cluster_ids = {k: (k, all_cluster_ids_p) for k in ['ds2_nopam']}
metadata_tables = dict(
        encode='/icgc/dkfzlsdf/analysis/hs_ontogeny/databases/enrichment_databases/lola_chipseq_2018-04-12/mm10/encodeTFBSmm10/regions/encode_annotations.csv',
        codex_hematopoietic='/icgc/dkfzlsdf/analysis/hs_ontogeny/databases/enrichment_databases/lola_chipseq_2018-04-12/mm10/codex/regions/codex_annotations.csv',
        codex_all=codex_all_metadata_table_fp,
        homer='/icgc/dkfzlsdf/analysis/hs_ontogeny/databases/enrichment_databases/homer/homer_lola/homer_metadata-table.csv',
        cpg_islands='/icgc/dkfzlsdf/analysis/hs_ontogeny/databases/enrichment_databases/mm10-genome-annotation-par/mm10_genome_anntation/regions/cpg-islands_metadata-table.tsv',
)
tf_enrichment_config = {
    'tasks': dict(cluster_ids=cluster_ids,
                  metadata_tables=metadata_tables),
    'output_dir': regionset_enrichment_output_dir,
    'tmpdir': project_temp_dir,
    'chromosomes': chromosomes,
    'resource_size': 'L',
}
enrichment_config_path = tmpdir.name + '/tf-enrichment.json'
with open(enrichment_config_path, 'wt') as fout:
    json.dump(tf_enrichment_config, fout)

enrichment_jobscript = ('/home/kraemers/projects/mouse_hematopoiesis/src'
                        '/mouse_hematopoiesis/wgbs/clustering/clustering1/jobscript.sh')
enrichment_snakefile = ('/home/kraemers/projects/mouse_hematopoiesis/src'
                        '/mouse_hematopoiesis/wgbs/clustering/clustering1'
                        '/clustering-eval-enrichments.smk')
snakemake_args = dict(
        snakefile=enrichment_snakefile,
        configfile=enrichment_config_path,
        dryrun=False,
        cores = 24,
        # forcerun = ['create_default_visualizations'],
        # forcerun = ['compute_overlap_stats'],
        keepgoing=True,
        printreason=True,
)
# if run on cluster
snakemake_args.update(dict(
        latency_wait=60,
        nodes=2000,
        local_cores=24,
        jobscript=enrichment_jobscript,
        cluster="bsub -R rusage[mem={params.avg_mem}] -M {params.max_mem} -n {threads} -J {params.name} -W {params.walltime} -o /home/kraemers/temp/logs/",
))
snakemake_successful = snakemake.snakemake(**snakemake_args)
if not snakemake_successful:
    raise RuntimeError('Snakemake returned non-zero exit code')

# Filter codex_all

with open('/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/cohort_results/analyses/dendritic_cells/dmr_characterization/enrichments2/region-sets/ds2_nopam/codex_all/cluster-overlap-stats_codex_all.p', 'rb') as fin:
    cluster_overlap_stats = pickle.load(fin)

annos = pd.read_csv(codex_all_metadata_table_fp, sep='\t', header=0)
hematopoietic_pops = ['Multipotent myeloid progenitor',
                      'Haematopoietic progenitor',
                      'Leukemogenic',
                      'Mouse ErythroLeukaemic',
                      'Macrophages',
                      'T-Cells',
                      'B-Cells',
                      'Erythroid',
                      'Megakaryocyte Progenitors',
                      'AML pre-leukaemic',
                      'Myeloid progenitor cells',
                      'Plasmablasts',
                      'Megakaryocyte',
                      'Leukaemia',
                      'Hematopoietic Progenitor Cells',
                      'NK cells',
                      'Macrophage',
                      'Lymphoid cells',
                      'Dendritic',
                      'Primary murine bone marrow macrophage cells',
                      'Hematopoietic Stem Cells',
                      'Erythroid Progenitors',
                      'Erythroid progenitor',
                      'Myeloid Progenitors',
                      'Haemogenic endothelium',
                      'Haematopoietic progenitors',
                      'Megakaryoblastic cells',
                      'Haematopoietic precursor and progenitor cells',
                      'Haematopoietic stem and progenitor cells',
                      'Thymus',
                      'Mouse ErythroLeukaemic ',
                      ]

min_pvalue_df = cluster_overlap_stats.cluster_pvalues.min().to_frame('min_pvalue')

full_df = pd.concat([min_pvalue_df, annos.set_index('name')], axis=1).reset_index().rename(columns={"index": 'name'})
full_df = full_df.sort_values('min_pvalue', ascending=True)

# discard same cell type and antibody
unique_df = full_df.groupby(['cellType', 'antibody'], as_index=False).first()
# discard non-hematopoietic if there is a hematopoietic population
def filter_non_hematopoietic(group_df):
    is_hema = group_df.cellType.isin(hematopoietic_pops)
    if is_hema.any():
        return group_df.loc[is_hema, :]
    else:
        return group_df.iloc[[0], :]
single_non_hematopoietic = unique_df.groupby(['antibody'], group_keys=False).apply(filter_non_hematopoietic)

cluster_overlap_stats.cluster_pvalues = cluster_overlap_stats.cluster_pvalues.loc[:, single_non_hematopoietic['name']]
cluster_overlap_stats._odds_ratio = cluster_overlap_stats.log_odds_ratio.loc[ :, single_non_hematopoietic['name']]


codex_mix_dir = '/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/cohort_results/analyses/dendritic_cells/dmr_characterization/enrichments2/region-sets/ds2_nopam/codex_mix'
Path(codex_mix_dir).mkdir(parents=True, exist_ok=True)
trunk_path = codex_mix_dir + '/codex_mix'
plot_pattern_png = trunk_path + '_max-pvalue-{maxpvalue}_format-{format}_filter-on-{filteron}.png'

# cluster x features
with mpl.rc_context(paper_context):

    # print('2x2 cell filtered plots')
    for max_pvalue in [1, 1e-5, 1e-10]:
        print(max_pvalue)

        fig = rsp.barcode_heatmap(cluster_overlap_stats,
                                  clusters_as_rows=True,
                                  col_width_cm=0.5, row_height_cm=0.1,
                                  linewidth=0.5,
                                  max_pvalue=max_pvalue,
                                  # the default extend='both' leads to a matplotlib warning and the colorbar is missing
                                  cbar_args = dict(shrink=0.4, aspect=20),
                                  filter_on_per_feature_pvalue = False,
                                  vlim=(-5, 5),
                                  vmin_quantile=0.02,
                                  )

        out_png = plot_pattern_png.format(maxpvalue=max_pvalue,
                                          format='wide',
                                          filteron='2x2-cell')
        fig.savefig(out_png)
        fig.savefig(out_png.replace('.png', '.pdf'))

# rank

with mpl.rc_context(paper_context):

    max_pvalue = 1e-5

    fig = rsp.barcode_heatmap(cluster_overlap_stats,
                              clusters_as_rows=True,
                              col_width_cm=0.5,
                              row_height_cm=0.1,
                              linewidth=0.5,
                              max_pvalue=max_pvalue,
                              n_top_hits=50,
                              # the default extend='both' leads to a matplotlib warning and the colorbar is missing
                              cbar_args = dict(shrink=0.4, aspect=20),
                              filter_on_per_feature_pvalue = False,
                              vlim=(-5, 5),
                              vmin_quantile=0.02,
                              )

    out_png = plot_pattern_png.format(maxpvalue=f'{max_pvalue}-top-50',
                                      format='wide',
                                      filteron='2x2-cell')
    fig.savefig(out_png)
    fig.savefig(out_png.replace('.png', '.pdf'))

# rank and whitelist

whitelist = [
    'Aff4_Embryonic Stem Cell_394',
    'Ascl2_T-Cells_218',
    'Atf3_Primary murine bone marrow macrophage cells_232',
    'Atf4_Fibroblast_438',
    'Batf_Dendritic_221',
    'Batf_T-Cells_64',
    'Brd4_Leukaemia_216',
    'CTCF_B-Cells_159',
    'CTCF_Plasmablasts_160',
    # 'Cbfa2t2_Mouse ErythroLeukaemic _571',
    'Cbfa2t3_Mouse ErythroLeukaemic_40',
    'Cbfb_Megakaryoblastic cells_409',
    'Cbfb_T-Cells_417',
    'Cbx7_Haematopoietic progenitor_487',
    'Cbx8_Haematopoietic progenitor_488',
    'Cdk8_Embryonic stem cells_157',
    'Cdk9_Embryonic stem cells_158',
    'Cdx2_Motor neuron progenitors_529',
    'Cebpa_B-Cells_217',
    # 'Cebpa_Hematopoietic Progenitor Cells_171',
    'Cebpa_Hematopoietic Stem Cells_239',
    'Cebpa_Macrophage_207',
    # 'Cebpa_Macrophages_21',
    'Cebpb_Haemogenic endothelium_551',
    'Cebpb_Macrophages_22',
    # 'Cebpd_Adipocytes_28',
    'Chd2_Mouse ErythroLeukaemic_29',
    'Chd4_Embyonic stem cell_367',
    'Crebbp_Embryonic stem cells_214',
    'Ctcf_B-Cells_32',
    'Ctcf_Dendritic_444',
    'Ctcf_Erythroid_177',
    'Ctcf_Mouse ErythroLeukaemic_504',
    'Ctcf_Multipotent myeloid progenitor_16',
    'Ctcf_T-Cells_31',
    # 'Ctcf_Thymus_499',
    'Ctnnb1_Embryonic Stem Cell_152',
    'Ctr9_Embryonic Stem Cell_145',
    'Ddit3_Fibroblast_437',
    'Dpy30_Embryonic Stem Cell_354',
    'E2f1_Dendritic_450',
    'E2f4_Dendritic_453',
    'Ebf1_B-Cells_440',
    'Egr1_Dendritic_482',
    'Egr2_Dendritic_467',
    'Egr2_T-Cells_197',
    'Elf1_T-Cells_36',
    'Ell2_Embryonic Stem Cell_393',
    'Ell3_Embryonic stem cells_514',
    'Ep300_Megakaryocyte_166',
    'Ep300_Mouse ErythroLeukaemic_576',
    'Ep300_T-Cells_558',
    'Erg_Leukaemia_168',
    'Esrrb_Embryonic Stem Cell_250',
    'Ets1_Megakaryocyte Progenitors_41',
    'Ets1_T-Cells_43',
    'Ets2_Dendritic_457',
    'Etv6_T-Cells_68',
    'Ezh2_Embryonic Stem Cell_198',
    'Fli1_Haemogenic endothelium_549',
    'Fli1_Multipotent myeloid progenitor_304',
    'Fli1_T-Cells_357',
    'Fos_Mast_12',
    'Fosl2_T-Cells_70',
    'Foxa2_Pancreas Cells_176',
    'Foxo1_B-Cells_299',
    'Foxo1_T-Cells_167',
    'Foxp3_T-Cells_44',
    'Gata1_Erythroid_37',
    # 'Gata1_Erythroid Progenitors_382',
    # 'Gata1_Erythroid progenitor_434',
    'Gata1_Megakaryocyte_509',
    # 'Gata1_Megakaryocyte Progenitors_47',
    # 'Gata1_Mouse ErythroLeukaemic_45',
    # 'Gata1 (V205G)_Megakaryocyte Progenitors_436',
    'Gata2_Erythroid Progenitors_380',
    'Gata2_Haematopoietic progenitor_347',
    'Gata2_Megakaryocyte Progenitors_49',
    # 'Gata2_Multipotent myeloid progenitor_305',
    'Gata3_T-Cells_283',
    'Gata4_Cardiac muscle_313',
    'Gata6_Mouse Embryonic fibroblasts_237',
    'Gcgr_Adipocytes_51',
    'Gfi1_AML pre-leukaemic_50',
    'Gfi1_Lymphoid cells_208',
    # 'Gfi1b_Leukemogenic_18',
    'Gfi1b_Multipotent myeloid progenitor_306',
    'Hdac2_Embyonic stem cell_364',
    'Hif1a_Dendritic_462',
    'Hif1a_T-Cells_71',
    'Hnf1a_T-Cells_170',
    'Hoxa9_Hematopoietic Stem Cells_238',
    'Hoxb4_Haematopoietic precursor and progenitor cells_426',
    # 'Hoxb4_Haematopoietic stem and progenitor cells_428',
    'Hsf1_striatal cells_511',
    'Ikzf1_B-Cells_123',
    'Irf1_Dendritic_464',
    'Irf4_B-Cells_534',
    'Irf4_T-Cells_540',
    'Irf8_Dendritic_224',
    'Irf8_Myeloid progenitor cells_525',
    'Jarid2_Embryonic stem cell_182',
    # 'Jun_T-Cells_542',
    'JunB_T-Cells_230',
    'JunD_T-Cells_231',
    'Junb_Dendritic_222',
    'Junb_Macrophages_564',
    # 'Junb_T-Cells_541',
    # 'Jund_Mouse ErythroLeukaemic_565',
    # 'Jund_T-Cells_543',
    'Kdm1a_Embyonic stem cell_366',
    'Kdm2a_Embryonic stem cells_52',
    'Kdm2b_Embryonic Stem Cell_105',
    'Kdm4b_Embryonic Stem Cell_135',
    'Kdm4c_Embryonic Stem Cell_136',
    'Kdm6b_T-Cells_86',
    'Klf1_Erythroid_181',
    'Klf4_Embryonic Stem Cell_249',
    'Ldb1_Haematopoietic progenitor_345',
    # 'Ldb1_Mouse ErythroLeukaemic_567',
    'Ldb1_Mouse ErythroLeukaemic _566',
    'Lmo2_Multipotent myeloid progenitor_14',
    'Lyl1_Multipotent myeloid progenitor_307',
    'Maf_T-Cells_65',
    'Maff_Dendritic_468',
    'Mafk_Mouse ErythroLeukaemic_568',
    'Mapk8_Neuron_339',
    'Max_Mouse ErythroLeukaemic_569',
    'Mbd1_Embryonic Stem Cells_532',
    'Med1_B-Cells_122',
    'Med12_Embryonic Stem Cell_322',
    'Mef2a_Cardiac muscle_314',
    'Meis1_Multipotent myeloid progenitor_308',
    'Men1_T-Cells_225',
    'Mitf_Mast_10',
    'Mllt3_Leukaemia_379',
    'Mtf2_Embryonic Stem Cell_260',
    'Mxi1_Mouse ErythroLeukaemic_572',
    'Myb_Mouse ErythroLeukaemic_19',
    'Myb_Myeloid Progenitors_301',
    'Myc_Mouse ErythroLeukaemic_20',
    'Mycn_Embryonic Stem Cell_252',
    'Nanog_Embryoid bodies_164',
    'Ncoa3_Embryonic Stem Cell_498',
    'Nfe2_Erythroid_179',
    'Nfya_Neuron_342',
    'Nipbl_Embryonic Stem Cell_324',
    'Nkx2-5_Cardiac muscle_315',
    'Notch1_Leukaemia_385',
    'Nr5a2_Embryonic stem cells_266',
    'Pax5_B-Cells_513',
    'Phf19_Embryonic Stem Cell_374',
    'Polr2a_Primary murine bone marrow macrophage cells_233',
    'Pou2f2_B-Cells_574',
    'Pou5f1_Embryonic Stem Cell_144',
    'Pparg_Adipocytes_578',
    'Prdm14_Embryonic Stem Cell_337',
    'Prep_Embryonic Trunk_531',
    'Ptf1a_Pancreas Cells_172',
    'Rad21_B-Cells_349',
    'Rad21_Mouse ErythroLeukaemic_603',
    'Rag2_T-Cells_293',
    'Rbbp5_Embryonic Stem Cell_327',
    'Rbpj_Leukaemia_386',
    'Rbpjl_Pancreas Cells_175',
    'Rdbp_Mouse ErythroLeukaemic_573',
    'Rel_Dendritic_480',
    'Rela_Dendritic_474',
    'Rela_Macrophages_195',
    'Relb_Dendritic_472',
    'Rest_Embryonic stem cell_184',
    'Ring1B_Embryonic Stem Cell_124',
    'Ring1b_Embryonic Stem Cell_106',
    'Rnf2_Megakaryoblastic cells_410',
    'Rnf2_T-Cells_418',
    'Rorc_T-Cells_98',
    'Runx1_Haematopoietic progenitor_15',
    # 'Runx1_Haematopoietic progenitors_396',
    # 'Runx1_Haemogenic endothelium_548',
    # 'Runx1_Megakaryoblastic cells_408',
    'Runx1_Megakaryocyte_165',
    # 'Runx1_Multipotent myeloid progenitor_310',
    'Runx1_T-Cells_416',
    # 'Runx2_Osteoblast precursor_119',
    'Runx3_NK cells_200',
    'Runx3_T-Cells_190',
    'Rxra_Adipocytes_605',
    'STAT5_Fibroblast_433',
    'Setdb1_Embryonic Stem Cell_262',
    'Sin3a_Embryonic Stem Cell_335',
    'Smad1_Erythroid Progenitors_381',
    'Smad3_B-Cells_295',
    'Smc1a_Embryonic Stem Cell_319',
    'Smc3_Mouse ErythroLeukaemic_608',
    'Sox2_Embryonic Stem Cell_242',
    'Spi1_B-Cells_602',
    'Spi1_Dendritic_220',
    # 'Spi1_Erythroid Progenitors_297',
    'Spi1_Erythroid progenitor_298',
    'Spi1_Macrophages_589',
    'Spi1_Multipotent myeloid progenitor_309',
    # 'Spi1_Myeloid progenitor cells_121',
    'Spi1_T-Cells_398',
    'Srf_Cardiac muscle_316',
    'Stag1_Embryonic Stem Cells_486',
    'Stag2_Embryonic Stem Cells_485',
    'Stat1_Macrophages_423',
    'Stat1_T-Cells_560',
    'Stat3_Dendritic_360',
    'Stat3_T-Cells_355',
    'Stat4_T-Cells_302',
    'Stat5a_T-Cells_490',
    # 'Stat5a/b_T-Cells_356',
    'Stat5b_Dendritic_362',
    'Stat5b_T-Cells_491',
    'Stat6_Macrophages_617',
    'Stat6_T-Cells_212',
    'Supt5_Embryonic Stem Cell_270',
    # 'Supt5h_Embryonic Stem Cell_151',
    'Suz12_Embryonic Stem Cell_372',
    'Taf1_Embryonic Stem Cell_404',
    'Taf3_Embryonic Stem Cell_403',
    'Tal1_Erythroid_508',
    'Tal1_Haematopoietic progenitor_346',
    # 'Tal1_Haemogenic endothelium_552',
    # 'Tal1_Mouse ErythroLeukaemic _507',
    # 'Tal1_Multipotent myeloid progenitor_311',
    'Tbp_Embryonic Stem Cell_405',
    'Tbx21_T-Cells_419',
    'Tbx3_Embryoid bodies_162',
    'Tbx5_Cardiac muscle_317',
    'Tcf3_Multipotent myeloid progenitor_17',
    'Tcf7_Haematopoietic progenitors_395',
    'Tet1_Embryonic Stem Cell_334',
    'Tfcp2l1_Embryonic Stem Cell_245',
    'Tfe3_Embryonic Stem Cell_546',
    'Thap11_Embryonic stem cells_215',
    'Trp53_Embryonic Stem Cell_351',
    'Trp53S18P_Embryonic Stem Cell_353',
    'Wdr5_Embryonic Stem Cell_328',
    'Yy1_Embryonic Stem Cell_407',
    'Zfp143_Embryonic Stem cells_527',
    'Zfx_Embryonic Stem Cell_247',
    'p300_Embryonic Stem Cell_331',
    'pStat3_T-Cells_296'
]

clear_names = {
    'Macrophages': 'Macrophage',
    'Hematopoietic Stem Cells': 'HSC',
    'Leukaemia': 'Leukemia',
    'Erg_Leukaemia': 'Leukemia',
    # 'AML pre-leukaemic': 'Leukemic',
    'Megakaryoblastic cells': 'Megakaryoblast',
    'Multipotent myeloid progenitor': 'MPP',
}


def format_names(s):
    antibody, *pop, uid = s.split('_')
    pop = '_'.join(pop)
    pop = clear_names.get(pop, pop)
    return f'{antibody} | {pop}'

cluster_overlap_stats.cluster_pvalues = cluster_overlap_stats.cluster_pvalues.loc[:, whitelist].rename(columns=format_names)
cluster_overlap_stats._odds_ratio = cluster_overlap_stats.log_odds_ratio.loc[:, whitelist].rename(columns=format_names)


with mpl.rc_context(paper_context):

    max_pvalue = 1e-5

    fig = rsp.barcode_heatmap(cluster_overlap_stats,
                              clusters_as_rows=False,
                              col_width_cm=0.5,
                              row_height_cm=0.1,
                              linewidth=0.5,
                              max_pvalue=max_pvalue,
                              n_top_hits=50,
                              # the default extend='both' leads to a matplotlib warning and the colorbar is missing
                              cbar_args = dict(shrink=0.4, aspect=20),
                              filter_on_per_feature_pvalue = False,
                              vlim=(-5, 5),
                              vmin_quantile=0.02,
                              )

    out_png = plot_pattern_png.format(maxpvalue=f'{max_pvalue}-top-50-whitelist',
                                      format='wide',
                                      filteron='2x2-cell')
    fig.axes[0].set_xticklabels(fig.axes[0].get_xticklabels(), rotation=0)
    fig.savefig(out_png)
    fig.savefig(out_png.replace('.png', '.pdf'))


# Create report
# ======================================================================
enrichment_plot_pattern = (regionset_enrichment_output_dir +
                           '/{partitioning}/{database}/{database}_max-pvalue-{maxpvalue}_format-{format}_filter-on-{filteron}.png')
section_cols = ['partitioning', 'database', 'filteron', 'maxpvalue', 'format']
metadata_table = fr.pattern_to_metadata_table(enrichment_plot_pattern,
                                              field_constraints={
                                                  # 'database': 'ensembl_regulatory_regions',
                                              })
metadata_table.sort_values(section_cols, inplace=True)
fr.copy_report_files_to_report_dir(metadata_table, root_dir=enrichment_output_dir,
                                   report_dir=enrichment_report_dir)
report_config = fr.convert_metadata_table_to_report_json(
        metadata_table, section_cols)
report_config.update(dict(toc_headings='h1, h2, h3, h4, h5', autocollapse_depth='2'))
fr.Report({'region_set_enrichments': report_config}).generate(enrichment_report_dir)


# Genome annotations
# **********************************************************************

# Config
# ----------------------------------------------------------------------
# input
cluster_ids_p = '/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/cohort_results/analyses/hierarchy/clustering/full-hierarchy2-merged/cluster-ids.p'
annos_v2_p = '/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/cohort_results/analyses/hierarchy/annotation/hierarchy-dmrs/v2/hierarchy-dmrs-anno_primary-annotations.p'

# output
annotation_enrichment_output_dir = enrichment_output_dir + '/annotation-enrichment'
barcode_plot_png_pattern = (annotation_enrichment_output_dir
                            + '/{anno_name}/{clustering}/barcode-plot_maxpvalue-{max_pvalue}.png')
cluster_overlap_stats_obj_pattern = (annotation_enrichment_output_dir
                                     + '/{anno_name}/{clustering}/cluster-overlap-stats.p')

# prepare cluster overlap counts
# ----------------------------------------------------------------------
cluster_ids = pd.read_pickle(cluster_ids_p).iloc[:, 0]
cluster_ids.index = cluster_ids.index.droplevel(-1)
cluster_ids.name = 'cluster_id'
cluster_ids.head()

annos_df = pd.read_pickle(annos_v2_p).set_index('Chromosome Start End'.split())['feat_class']
annos_df_merged = (annos_df
                   .groupby(level=[0,1,2])
                   .first())
annos_df_merged.loc[annos_df_merged == "3'-UTR"] = 'Exon'
annos_df_merged = annos_df_merged.replace({'intron': 'Intron', 'intergenic': 'Intergenic', 'exon': 'Exon'})
assert annos_df_merged.index.equals(cluster_ids.index)
# annos_df_merged2 = (annos_df
#                    .groupby(level=[0,1,2])
#                    .agg(lambda ser: ser.iat[0]))
# assert annos_df_merged.equals(annos_df_merged2)

cluster_overlap_counts = pd.crosstab(index=cluster_ids, columns=annos_df_merged)
cluster_sizes = cluster_ids.value_counts().sort_index()
cluster_sizes.index.name = 'cluster_id'

# Genomic distribution plots
# ----------------------------------------------------------------------

genomic_region_anno_plot_dir = '/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/cohort_results/analyses/hierarchy/annotation/hierarchy-dmrs/v2/plots'
frequencies = cluster_overlap_counts.groupby('cluster_id').sum().apply(lambda ser: ser / ser.sum(), axis=1)
frequencies.columns.name = 'Genomic region'
frequencies.index.name = 'Cluster ID'
global_counts_df = cluster_overlap_counts.sum(axis=0).reset_index().set_axis(['Genomic region', 'Count'], axis=1, inplace=False).sort_values('Count', ascending=False)
freqs_long = frequencies.stack().to_frame('Fraction').reset_index()
freqs_long['Genomic region'] = pd.Categorical(freqs_long['Genomic region'], categories=global_counts_df['Genomic region'], ordered=True)

with mpl.rc_context(paper_context):
    fig, ax = plt.subplots(1, 1, figsize=(4/2.54, 4/2.54), constrained_layout=True, dpi=360)
    sns.barplot('Genomic region', 'Count', data=global_counts_df, color='gray', ax=ax)
    ax.set_xlabel('')
    ax.tick_params(axis='x', rotation=90)
    fig.savefig(genomic_region_anno_plot_dir + '/global-counts.png')
    fig.savefig(genomic_region_anno_plot_dir + '/global-counts.pdf')


plot_context = copy(paper_context)
plot_context['xtick.labelsize'] = 6
with mpl.rc_context(plot_context):
    g = sns.catplot(x='Cluster ID', y = 'Fraction', col='Genomic region', col_wrap=3,
                    data=freqs_long, kind='bar',
                    sharey=False, color='gray', height=4/2.54, aspect=1.7)
    g.set_titles("{col_name}")
    # g.fig.subplots_adjust(hspace=0.1, wspace=0.1)
    g.fig.set_dpi(360)
    # g.set_xticklabels(rotation=90, va='bottom')
    g.set_xticklabels(step=3, va='center')
    g.savefig(genomic_region_anno_plot_dir + '/frequencies_col-genomic-region.png')
    g.savefig(genomic_region_anno_plot_dir + '/frequencies_col-genomic-region.pdf')


with mpl.rc_context(paper_context):
    g = sns.catplot(x='Genomic region', y = 'Fraction', col='Cluster ID', col_wrap=3,
                    data=freqs_long, kind='bar',
                    sharey=True, color='gray', height=4/2.54, aspect=1)
    g.set_titles("Cluster {col_name}")
    g.fig.set_dpi(360)
    g.set_xticklabels(rotation=90)
    g.set(xlabel='')
    g.savefig(genomic_region_anno_plot_dir + '/frequencies_col-cluster-id.png')
    g.savefig(genomic_region_anno_plot_dir + '/frequencies_col-cluster-id.pdf')


genomic_region_tmpdir = tempfile.TemporaryDirectory(dir=project_temp_dir)
dmrs_bed = genomic_region_tmpdir.name + '/regions.bed'
(cluster_ids.index.to_frame()
 .set_axis('#Chromosome Start End'.split(), axis=1, inplace=False)
 .to_csv(dmrs_bed, sep='\t', header=True, index=False))
amit_cluster_id_files = [enhancer_cluster_ids_by_cluster_mm10_bed.format(clusterid=cid)
                         for cid in range(1, 10)]
# TODO: sort!
anno_bt = BedTool(dmrs_bed).annotate(files=amit_cluster_id_files,
                                     names=[str(i) for i in range(1, 10)],
                                     counts=True)
anno_bt_names = 'Chromosome Start End'.split() + [str(i) for i in range(1, 10)]
dmr_enhancer_overlap = pd.read_csv(anno_bt.fn, sep='\t', names=anno_bt_names, comment='#')
dmr_enhancer_overlap.Chromosome = dmr_enhancer_overlap.Chromosome.astype(str)
dmr_enhancer_overlap = dmr_enhancer_overlap.set_index('Chromosome Start End'.split()).sort_index()

dmr_enhancer_overlap['any_overlap'] = dmr_enhancer_overlap.sum(axis=1).gt(0)

# Total number of overlaps
dmr_enhancer_overlap['any_overlap'].sum()
# 36400

assert dmr_enhancer_overlap.index.equals(annos_df_merged.index)
dmr_enhancer_overlap_with_anno = pd.concat(
        [dmr_enhancer_overlap, annos_df_merged], axis=1)
feat_class_overlap_counts = (dmr_enhancer_overlap_with_anno
                             .groupby('feat_class')['any_overlap']
                             .agg([np.sum, lambda ser: ser.shape[0]])
                             .set_axis(['n_overlap', 'n_total'], axis=1, inplace=False)
                             )
feat_class_overlap_counts['fraction'] = feat_class_overlap_counts.eval('n_overlap / n_total')
feat_class_overlap_counts.to_csv(genomic_region_anno_plot_dir + '/amit-enhancer-overlap.tsv', sep='\t')


# Enrichment
# ----------------------------------------------------------------------
cluster_overlap_stats = rsp.ClusterOverlapStats(hits=cluster_overlap_counts,
                                                cluster_sizes=cluster_sizes)

clustering = 'hierarchy-v2_29-clusters-merged'
anno_name = 'hierarchy-anno-v2'
barcode_plot_png_fp = barcode_plot_png_pattern.format(
        anno_name=anno_name, clustering=clustering, max_pvalue='{max_pvalue}')
cluster_overlap_stats_obj_fp = cluster_overlap_stats_obj_pattern.format(
        anno_name=anno_name, clustering=clustering)
Path(cluster_overlap_stats_obj_fp).parent.mkdir(exist_ok=True, parents=True)

with mpl.rc_context(paper_context):
    process_cluster_overlap_stats(
            cluster_overlap_stats, max_pvalues=[1],
            plot_args = dict(clusters_as_rows=True,
                             col_width_cm=0.5, row_height_cm=0.2,
                             linewidth=0.5,
                             plot_stat = 'p-value',
                             vlim=(-3, 3),
                             vmin_quantile=0.05
                             ),
            cores=24,
            barcode_plot_png_by_maxpvalue=barcode_plot_png_fp,
            cluster_overlap_stats_out_fp=cluster_overlap_stats_obj_fp,
            test_arg_dict= dict(replicate=int(1e4), workspace=1_000_000),
            plots_only=True)

