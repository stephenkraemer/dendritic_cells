import json
from copy import deepcopy

import pandas as pd
import region_set_profiler as rsp

from dendritic_cells.wgbs.clustering_and_heatmaps import merged_dmrs_p, cluster_ids_p
from dendritic_cells.wgbs.annotation.bedtools_based_annotation import gencode_anno_p
from dendritic_cells.wgbs.enrichment.enrichment_plots import (
    masked_plot, norm_plot, plot_all_enriched_features, results_dir)
from dendritic_cells.base_paths import project_temp_dir


msigdb_gmts_dir = '/icgc/dkfzlsdf/analysis/hs_ontogeny/databases/region_set_profiler_databases/msigdb_gmts'
pathways_gmt = '/icgc/dkfzlsdf/analysis/hs_ontogeny/databases/region_set_profiler_databases/msigdb_gmts/canonical-pathways.gmt'
go_bp_gmt = '/icgc/dkfzlsdf/analysis/hs_ontogeny/databases/region_set_profiler_databases/msigdb_gmts/go_biological-process.gmt'



anno_df = pd.read_pickle(gencode_anno_p)
anno_df.set_index(['Chromosome', 'Start', 'End', 'region_id'], drop=False, inplace=True)

cluster_ids_df = pd.read_pickle(cluster_ids_p)
cluster_ids_ser = cluster_ids_df['ds2_nopam']

clustered_regions_anno = anno_df.loc[cluster_ids_ser.index, :]

# TODO: does gencode use the same gene symbols as my gmt files? how to check in code? should i switch to ensembl IDs?
gene_names_ser = clustered_regions_anno['gene_name'].astype(str)
# gene_names_ser.index = anno_df['region_id']
# gene_names_ser
compact_gene_anno = (gene_names_ser.groupby(gene_names_ser.index.names).agg(lambda ser: ','.join(ser))).to_frame('Gene')
assert cluster_ids_ser.index.equals(compact_gene_anno.index)

compact_gene_anno['Gene'] = compact_gene_anno['Gene'].str.upper()



genesets_fp = pathways_gmt
database = 'msigdb_pathways'
n_iter = 1e4

overlap_stats = rsp.GenesetOverlapStats(annotations=compact_gene_anno, genesets_fp=genesets_fp)
overlap_stats.compute()
cluster_stats = overlap_stats.aggregate(cluster_ids=cluster_ids_ser)

all_hits = cluster_stats.hits.copy(deep=True)
all_hits = all_hits.loc[:, ~all_hits.sum().eq(0)]

def meets_cochran(ser):
    expected = expected_freq(np.array([ser, cluster_stats2.cluster_sizes - ser]))
    emin = (np.round(expected) >= 1).all()
    perc_expected = ((expected > 5).sum() / expected.size) > 0.8
    return emin and perc_expected
meets_cochran_ser = all_hits.apply(meets_cochran)

chi_square_hits = all_hits.loc[:, meets_cochran_ser]
fisher_hits = all_hits.loc[:, ~meets_cochran_ser]

cluster_stats.hits = chi_square_hits
chi_square_pvalues = cluster_stats.test_for_enrichment('chi_square')

cluster_stats.hits = fisher_hits
test_args =  dict(simulate_pval=True, replicate=int(n_iter), workspace=200_000)
fisher_pvalues = cluster_stats.test_for_enrichment(
        'fisher', cores=24, test_args=test_args)

all_pvalues = pd.concat([chi_square_pvalues, fisher_pvalues], axis=0).sort_index()

is_significant = all_pvalues['qvalues'].lt(0.05)

cluster_stats.hits  = all_hits.loc[:, is_significant]

# overlap_counts_fp = f'/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/cohort_results/analyses/dendritic_cells/dmr_characterization/enrichments/coverage-stats/coverage-counts_{database}.p'
output_dir = results_dir / database
output_dir.mkdir(parents=True, exist_ok=True)

# TODO: automatically adjust for label length
masked_plot(cluster_stats, output_dir,
            col_width_cm = 2,
            row_height_cm = 1,
            quantile = 0.6,
            min_abs_log_odds = -1,
            linewidth = 0.3,
            param_name='quantile-0.6',
            cbar_size_cm = 2,
            row_labels_show=True,
            )

# norm_plot(
#         cluster_stats, output_dir,
#         col_width_cm = 2,
#         row_height_cm = 1,
#         quantile = 0.4,
#         min_abs_log_odds = 0,
#         linewidth = 0.3,
#         param_name='quantile-0.8_min-logodds-1',
#         cbar_size_cm = 2,
#         row_labels_show=True,
#         norm_plateau_height = 0.05
# )

# %%

import pickle
from pathlib import Path
from typing import Union

import pandas as pd
import numpy as np

import codaplot as co
from matplotlib import colors
import matplotlib.pyplot as plt
import region_set_profiler as rsp
from dendritic_cells.base_paths import wgbs_cohort_results_dir
from dendritic_cells.config import paper_context
import matplotlib as mpl

results_dir = Path(wgbs_cohort_results_dir).joinpath(
        'analyses/dendritic_cells/dmr_characterization/enrichments')

class MidpointNormalize(mpl.colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        mpl.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

# %%
def masked_plot(
        overlap_counts: Union[str, pd.DataFrame], output_dir, param_name,
        quantile, min_abs_log_odds,
        linewidth, col_width_cm, row_height_cm, cbar_size_cm,
        row_dendrogram=False, row_labels_show=True):

    print(param_name)

    if isinstance(overlap_counts, str):
        with open(overlap_counts, 'rb') as fin:
            overlap_counts = pickle.load(fin)
    assert isinstance(overlap_counts, rsp.ClusterOverlapStats)

    plot_data: pd.DataFrame = overlap_counts.log_odds_ratio
    # Drop NA because we want to cluster the features
    plot_data = plot_data.dropna(how='any', axis=1)


    # Set all log odds ratios below the quantile to 0
    # Possible improvement: separately for depletion and enrichment
    q = np.quantile(plot_data.abs(), quantile)
    plot_data = plot_data.where(plot_data.abs().gt(q), 0)
    plot_data = plot_data.loc[:, plot_data.abs().gt(min_abs_log_odds).any()]

    vmin = np.quantile(plot_data, 0.02)
    vmax = np.quantile(plot_data, 0.98)


    plot_data = plot_data.T
    width = plot_data.shape[1] * col_width_cm / 2.54
    height = plot_data.shape[0] * row_height_cm / 2.54
    print('width (cm)', width * 2.54, 'height (cm)', height * 2.54)
    shrink = cbar_size_cm / height

    print(plot_data.shape)

    with mpl.rc_context(paper_context):
        print('Clustered plot')
        cdg = (co.ClusteredDataGrid(main_df=plot_data)
               .cluster_rows(method='average', metric='cityblock'))
        gm = cdg.plot_grid(grid=[
            [
                co.Heatmap(df=plot_data,
                           cmap='RdBu_r',
                           vmin = vmin,
                           vmax = vmax,
                           norm = MidpointNormalize(vmin=vmin, vmax=vmax, midpoint=0),
                           row_labels_show=row_labels_show,
                           rasterized=False,
                           edgecolor='white',
                           linewidth=linewidth,
                           cbar_args = dict(shrink=shrink, aspect=20)
                           ),
            ]
        ],
                figsize=(width, height),
                height_ratios=[(1, 'rel')],
                row_dendrogram=row_dendrogram,
        )
        gm.create_or_update_figure()
        gm.fig.savefig(output_dir / f'all-significant_clustered_log-odds_masked-{param_name}.png')
        gm.fig.savefig(output_dir / f'all-significant_clustered_log-odds_masked-{param_name}.pdf')
# %%

# %%
def norm_plot(overlap_counts: Union[str, pd.DataFrame], output_dir, param_name,
              quantile, min_abs_log_odds,
              norm_plateau_height,
              linewidth, col_width_cm, row_height_cm, cbar_size_cm,
              row_dendrogram=False, row_labels_show=True):

    print(param_name)

    if isinstance(overlap_counts, str):
        with open(overlap_counts, 'rb') as fin:
            overlap_counts = pickle.load(fin)
    assert isinstance(overlap_counts, rsp.ClusterOverlapStats)

    class RangeNorm(colors.Normalize):
        def __init__(self, neg_q, pos_q,  vmin, vmax, height):
            self.vmin = vmin
            self.vmax = vmax
            self.height = height
            self.neg_q = neg_q
            self.pos_q = pos_q
            colors.Normalize.__init__(self, vmin, vmax, clip=True)
        def __call__(self, value, clip=None):
            # xp = [self.vmin, self.neg_q, 0, self.pos_q, self.vmax]
            # yp = [0, 0.5 - self.height/2, 0.5, 0.5 + self.height/2, 1]
            # return np.ma.masked_array(np.interp(value, xp, yp))
            xp = [self.vmin, self.vmax]
            yp = [0, 1]
            res = np.interp(value, xp, yp)
            midpoint = self.neg_q + (self.pos_q - self.neg_q) / 2
            # # res[(value > self.neg_q) & (value < self.pos_q)] = 0.5
            view1 = value[(value > self.neg_q) & (value < midpoint)]
            if view1.size != 0:
                cubic = ((view1 - midpoint) / np.abs(view1.min()))**4
                res[(value > self.neg_q) & (value < midpoint)] = 0.5 + -cubic * (0.5 - np.interp(self.neg_q, xp, yp))
            view2 = value[(value < self.pos_q) & (value > midpoint)]
            if view2.size != 0:
                cubic = ((view2 + midpoint) / np.abs(view2.max()))**4
                res[(value < self.pos_q) & (value > midpoint)] = 0.5 + -cubic * (0.5 - np.interp(self.pos_q, xp, yp))
            return np.ma.masked_array(res)

    plot_data: pd.DataFrame = overlap_counts.log_odds_ratio
    plot_data = plot_data.dropna(how='any', axis=1)

    # log_odds_flat = np.ravel(plot_data)
    # vmin, neg_sat, neg_quant = np.quantile(log_odds_flat[log_odds_flat < 0], [0.001, 0.2, 0.3])
    # pos_quant, pos_sat, vmax = np.quantile(log_odds_flat[log_odds_flat > 0], [0.7, 0.8, 0.999])
    # neg_quant = -0.7
    # pos_quant = +0.7
    vmin = np.quantile(plot_data, 0.01)
    vmax = np.quantile(plot_data, 0.99)
    q = np.quantile(plot_data.abs(), quantile)
    neg_quant = -q
    pos_quant = q

    plot_data = plot_data.loc[:, plot_data.abs().gt(min_abs_log_odds).any()]

    # plot_data = plot_data.where(plot_data.abs().gt(1), 0)

    norm = RangeNorm(neg_q=neg_quant, pos_q=pos_quant,
                     vmin=vmin, vmax=vmax, height=norm_plateau_height)

    # Plot the norm function
    x = np.linspace(-6, 6, 300)
    y = norm(x)
    fig, ax = plt.subplots(1, 1)
    ax.plot(x, y)
    x = vmin, neg_quant, 0, pos_quant, vmax
    y = [0, 0.5 - norm_plateau_height/2, 0.5, 0.5 + norm_plateau_height/2, 1]
    ax.scatter(x, y)
    fig.savefig(output_dir / 'test.png')

    plot_data = plot_data.T

    width = plot_data.shape[1] * col_width_cm / 2.54
    height = plot_data.shape[0] * row_height_cm / 2.54
    shrink = cbar_size_cm / height

    # with mpl.rc_context(paper_context):
    cdg = (co.ClusteredDataGrid(main_df=plot_data)
           .cluster_rows(method='average', metric='cityblock'))
    gm = cdg.plot_grid(grid=[
        [
            co.Heatmap(df=plot_data,
                       cmap='RdBu_r',
                       norm=norm,
                       row_labels_show=row_labels_show,
                       rasterized=False,
                       edgecolor='white',
                       linewidth=linewidth,
                       cbar_args = dict(
                               shrink=shrink,
                               aspect=20,
                       )
                       ),
        ]
    ],
            figsize=(width, height),
            height_ratios=[(1, 'rel')],
            row_dendrogram=row_dendrogram,
    )
    gm.create_or_update_figure()
    gm.fig.savefig(output_dir / f'all-significant_clustered_log-odds_norm-{param_name}.png')
    gm.fig.savefig(output_dir / f'all-significant_clustered_log-odds_norm-{param_name}.pdf')
# %%


# %%
def plot_all_enriched_features(col_width_cm, output_dir, overlap_counts, row_height_cm, sign_threshold, statistics, test_statistics_fp):

    print('Reading input data')

    if isinstance(overlap_counts, str):
        with open(overlap_counts, 'rb') as fin:
            overlap_counts = pickle.load(fin)
    assert isinstance(overlap_counts, rsp.ClusterOverlapStats)

    # test_statistics = pd.read_pickle(test_statistics_fp)

    for statistic in statistics:
        print(statistic)
        print('Prepare plot')
        if statistic == 'normalized_ratio':
            plot_data = overlap_counts.normalized_ratio
        elif statistic == 'log_odds_ratio':
            plot_data = overlap_counts.log_odds_ratio
        else:
            raise ValueError

        # is_significant = test_statistics.qvalues < sign_threshold
        # plot_data = plot_data.loc[:, is_significant].sort_index(axis=1)
        plot_data = plot_data.sort_index(axis=1)
        plot_data = plot_data.dropna(how='any', axis=1)
        if plot_data.lt(0).any(axis=None):
            cmap = 'RdYlGn_r'
        else:
            cmap = 'YlGnBu_r'
        width = plot_data.shape[1] * col_width_cm / 2.54
        height = plot_data.shape[0] * row_height_cm

        print('Clustered plot')
        cdg = (co.ClusteredDataGrid(main_df=plot_data)
               .cluster_cols(method='average', metric='cityblock'))
        gm = cdg.plot_grid(grid=[
            [
                co.Heatmap(df=plot_data,
                           cmap=cmap,
                           row_labels_show=True,
                           rasterized=True,
                           ),
            ]
        ],
                figsize=(width, height),
                height_ratios=[(1, 'rel')],
                col_dendrogram=True,
        )
        gm.create_or_update_figure()
        gm.fig.savefig(output_dir / f'all-significant_clustered_{statistic}.png')
        gm.fig.savefig(output_dir / f'all-significant_clustered_{statistic}.pdf')

        print('Alphabetic plot')
        cdg = (co.ClusteredDataGrid(main_df=plot_data))
        gm = cdg.plot_grid(grid=[
            [
                co.Heatmap(df=plot_data,
                           cmap=cmap,
                           row_labels_show=True,
                           rasterized=True,
                           ),
            ]
        ],
                figsize=(width, height),
                height_ratios=[(1, 'rel')],
        )
        gm.create_or_update_figure()
        gm.fig.savefig(output_dir / f'all-significant_alphabetic_{statistic}.png')
        gm.fig.savefig(output_dir / f'all-significant_alphabetic_{statistic}.pdf')
# %%




cluster_stats2 = deepcopy(cluster_stats)


import numpy as np
idx = np.random.choice(cluster_stats2.hits.shape[1], 48)
cluster_stats2.hits = cluster_stats.hits.iloc[:, idx]
from scipy.stats.contingency import expected_freq


def has_zero_expected(ser):
    expected = expected_freq(np.array([ser, cluster_stats2.cluster_sizes - ser]))
    return np.any(np.round(expected) == 0)

meets_cochran_ser

matrices = list(cluster_stats.hits.apply(lambda ser: np.array([ser, cluster_stats.cluster_sizes - ser]).tolist()))
from FisherExact import fisher_exact
pvalues = [fisher_exact(m, simulate_pval=True, replicate=10000, workspace=200000) for m in matrices[0:40]]
with open(project_temp_dir + '/matrices.json', 'w') as fin:
    json.dump(matrices, fin)


cluster_stats2.hits = cluster_stats.hits.loc[:, meets_cochran_ser]
matrices = list(cluster_stats2.hits.apply(lambda ser: np.array([ser, cluster_stats2.cluster_sizes - ser]).tolist()))
pvalues = [fisher_exact(m, simulate_pval=True, replicate=10000, workspace=200000) for m in matrices[0:100]]
pvalues = cluster_stats2.test_for_enrichment('fisher', cores=24, test_args=test_args)

has_zero_expected = cluster_stats.hits.apply(has_zero_expected)
cluster_stats2.hits = cluster_stats.hits.loc[:, ~has_zero_expected]

pvalues = cluster_stats2.test_for_enrichment('chi_square')


cluster_stats2.hits = cluster_stats.hits.loc[:, ~meets_cochran_ser & ~has_zero_expected]
test_args =  dict(simulate_pval=True, replicate=int(1e5), workspace=100_000_000, seed=123)
test_args =  dict(simulate_pval=False, replicate=int(1e5), workspace=1_000_000_000, seed=123)
pvalues2 = cluster_stats2.test_for_enrichment('fisher', cores=24, test_args=test_args)

test_stats = cluster_stats2.test_for_enrichment('fisher', cores=24,
                                                test_args=dict(
                                                        simulate_pval=False, workspace=100_000_000,
                                                        hybrid=True))
is_significant = test_stats['qvalues'].gt(1e-3)






