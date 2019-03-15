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

def main():
    databases = ['homer', 'codex', 'encode']
    sign_threshold_d = {
        'codex': 1e-10,
        'encode': 1e-10,
        'homer': 1e-3
    }
    col_width_cm = 0.3
    row_height_cm = 0.7
    statistics = ['normalized_ratio', 'log_odds_ratio']

    for database in databases:
        print('General plots for ', database)
        test_statistics_fp = f'/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/cohort_results/analyses/dendritic_cells/dmr_characterization/enrichments/tests/{database}_chi-square.p'
        overlap_counts_fp = f'/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/cohort_results/analyses/dendritic_cells/dmr_characterization/enrichments/coverage-stats/coverage-counts_{database}.p'
        sign_threshold = sign_threshold_d[database]
        output_dir = results_dir / database
        output_dir.mkdir(parents=True, exist_ok=True)
        plot_all_enriched_features(
                col_width_cm, output_dir, overlap_counts_fp,
                row_height_cm, sign_threshold, statistics, test_statistics_fp)

    # Codex
    # ==================================================================
    database = 'codex'
    print('Manual plots for ', database)
    test_statistics_fp = f'/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/cohort_results/analyses/dendritic_cells/dmr_characterization/enrichments/tests/{database}_chi-square.p'
    overlap_counts_fp = f'/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/cohort_results/analyses/dendritic_cells/dmr_characterization/enrichments/coverage-stats/coverage-counts_{database}.p'
    output_dir = results_dir / database
    output_dir.mkdir(parents=True, exist_ok=True)

    masked_plot(overlap_counts_fp, output_dir,
                col_width_cm = 0.6,
                row_height_cm = 0.07,
                quantile = 0.8,
                min_abs_log_odds = 1,
                linewidth = 0.3,
                param_name='quantile-0.8-pretty',
                cbar_size_cm = 2,
                row_labels_show=True,
                )

    masked_plot(overlap_counts_fp, output_dir,
                col_width_cm = 0.6,
                row_height_cm = 0.4,
                quantile = 0.8,
                min_abs_log_odds = 1,
                linewidth = 0.3,
                param_name='quantile-0.8',
                cbar_size_cm = 2,
                row_labels_show=True,
                )

    masked_plot(overlap_counts_fp, output_dir,
                col_width_cm = 0.6,
                row_height_cm = 0.4,
                quantile = 0.6,
                min_abs_log_odds = 0.6,
                linewidth = 0.3,
                param_name='quantile-0.6',
                cbar_size_cm = 2,
                row_labels_show=True,
                )

    norm_plot(
            overlap_counts_fp, output_dir,
            col_width_cm = 0.8,
            row_height_cm = 0.4,
            quantile = 0.7,
            min_abs_log_odds = 1.5,
            linewidth = 0.3,
            param_name='quantile-0.7_min-logodds-1.5',
            cbar_size_cm = 8,
            row_labels_show=True,
            norm_plateau_height = 0.05
    )


    # Encode
    # ==================================================================
    database = 'encode'
    print('Manual plots for ', database)
    test_statistics_fp = f'/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/cohort_results/analyses/dendritic_cells/dmr_characterization/enrichments/tests/{database}_chi-square.p'
    overlap_counts_fp = f'/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/cohort_results/analyses/dendritic_cells/dmr_characterization/enrichments/coverage-stats/coverage-counts_{database}.p'
    output_dir = results_dir / database
    output_dir.mkdir(parents=True, exist_ok=True)

    masked_plot(overlap_counts_fp, output_dir,
                col_width_cm = 0.6,
                row_height_cm = 0.4,
                quantile = 0.8,
                min_abs_log_odds = 1,
                linewidth = 0.3,
                param_name='quantile-0.8-min-logodds-1',
                cbar_size_cm = 2,
                row_labels_show=True,
                )

    """
    - irf8 and spi1 should be pdc specific - bug?
    """

    masked_plot(overlap_counts_fp, output_dir,
                col_width_cm = 0.6,
                row_height_cm = 0.4,
                quantile = 0.6,
                min_abs_log_odds = 0.5,
                linewidth = 0.3,
                param_name='quantile-0.6-min-logodds-0.5',
                cbar_size_cm = 2,
                row_labels_show=True,
                )

    norm_plot(
            overlap_counts_fp, output_dir,
            col_width_cm = 0.8,
            row_height_cm = 0.4,
            quantile = 0.6,
            min_abs_log_odds = 0.5,
            linewidth = 0.3,
            param_name='quantile-0.6_min-logodds-0.5',
            cbar_size_cm = 8,
            row_labels_show=True,
            norm_plateau_height = 0.05
    )



    # Homer
    # ==================================================================
    database = 'homer'
    print('Manual plots for ', database)
    test_statistics_fp = f'/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/cohort_results/analyses/dendritic_cells/dmr_characterization/enrichments/tests/{database}_chi-square.p'
    overlap_counts_fp = f'/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/cohort_results/analyses/dendritic_cells/dmr_characterization/enrichments/coverage-stats/coverage-counts_{database}.p'
    output_dir = results_dir / database
    output_dir.mkdir(parents=True, exist_ok=True)

    masked_plot(overlap_counts_fp, output_dir,
                col_width_cm = 0.8,
                row_height_cm = 0.4,
                quantile = 0.8,
                min_abs_log_odds = 1,
                linewidth = 0.3,
                param_name='quantile-0.8_min-logodds-1',
                cbar_size_cm = 2,
                row_labels_show=True,
                )

    masked_plot(overlap_counts_fp, output_dir,
                col_width_cm = 0.4,
                row_height_cm = 0.07,
                quantile = 0.8,
                min_abs_log_odds = 1,
                linewidth = 0.5,
                param_name='quantile-0.8_min-logodds-1_pretty',
                cbar_size_cm = 2,
                row_labels_show=False,
                row_dendrogram=True
                )

    masked_plot(overlap_counts_fp, output_dir,
                col_width_cm = 0.8,
                row_height_cm = 0.4,
                quantile = 0.6,
                min_abs_log_odds = 0.5,
                linewidth = 0.3,
                param_name='quantile-0.6_min-logodds-0.5',
                cbar_size_cm = 2,
                row_labels_show=True,
                )

    norm_plot(
            overlap_counts_fp, output_dir,
            col_width_cm = 0.8,
            row_height_cm = 0.4,
            quantile = 0.6,
            min_abs_log_odds = 0.5,
            linewidth = 0.3,
            param_name='quantile-0.6_min-logodds-0.5',
            cbar_size_cm = 8,
            row_labels_show=True,
            norm_plateau_height = 0.05
    )

    norm_plot(
            overlap_counts_fp, output_dir,
            col_width_cm = 0.8,
            row_height_cm = 0.4,
            quantile = 0.7,
            min_abs_log_odds = 1,
            linewidth = 0.3,
            param_name='quantile-0.7_min-logodds-1',
            cbar_size_cm = 8,
            row_labels_show=True,
            norm_plateau_height = 0.05
    )

    norm_plot(
            overlap_counts_fp, output_dir,
            col_width_cm = 0.8,
            row_height_cm = 0.4,
            quantile = 0.8,
            min_abs_log_odds = 0.8,
            linewidth = 0.3,
            param_name='quantile-0.8_min-logodds-0.8',
            cbar_size_cm = 8,
            row_labels_show=True,
            norm_plateau_height = 0.05
    )


    # MSIGDB
    # ==================================================================

    # hallmarks
    # ------------------------------------------------------------------
    database = 'msigdb_hallmarks'
    print('Manual plots for ', database)
    test_statistics_fp = f'/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/cohort_results/analyses/dendritic_cells/dmr_characterization/enrichments/tests/{database}_chi-square.p'
    overlap_counts_fp = f'/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/cohort_results/analyses/dendritic_cells/dmr_characterization/enrichments/coverage-stats/coverage-counts_{database}.p'
    output_dir = results_dir / database
    output_dir.mkdir(parents=True, exist_ok=True)

    # %%
    masked_plot(overlap_counts_fp, output_dir,
                col_width_cm = 1.5,
                row_height_cm = 0.4,
                quantile = 0.8,
                min_abs_log_odds = 0.8,
                linewidth = 0.3,
                param_name='quantile-0.8_min-logodds-1',
                cbar_size_cm = 2,
                row_labels_show=True,
                )


    norm_plot(
            overlap_counts_fp, output_dir,
            col_width_cm = 1.5,
            row_height_cm = 0.4,
            quantile = 0.8,
            min_abs_log_odds = 0.8,
            linewidth = 0.3,
            param_name='quantile-0.8_min-logodds-0.8',
            cbar_size_cm = 8,
            row_labels_show=True,
            norm_plateau_height = 0.05
    )
    # %%

    # oncogenic signatures
    # ------------------------------------------------------------------
    database = 'msigdb_oncogenic'
    print('Manual plots for ', database)
    test_statistics_fp = f'/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/cohort_results/analyses/dendritic_cells/dmr_characterization/enrichments/tests/{database}_chi-square.p'
    overlap_counts_fp = f'/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/cohort_results/analyses/dendritic_cells/dmr_characterization/enrichments/coverage-stats/coverage-counts_{database}.p'
    output_dir = results_dir / database
    output_dir.mkdir(parents=True, exist_ok=True)

    # %%
    masked_plot(overlap_counts_fp, output_dir,
                col_width_cm = 1.5,
                row_height_cm = 0.4,
                quantile = 0.8,
                min_abs_log_odds = 0.8,
                linewidth = 0.3,
                param_name='quantile-0.8_min-logodds-1',
                cbar_size_cm = 2,
                row_labels_show=True,
                )


    norm_plot(
            overlap_counts_fp, output_dir,
            col_width_cm = 1.5,
            row_height_cm = 0.4,
            quantile = 0.8,
            min_abs_log_odds = 0.8,
            linewidth = 0.3,
            param_name='quantile-0.8_min-logodds-0.8',
            cbar_size_cm = 8,
            row_labels_show=True,
            norm_plateau_height = 0.05
    )
    # %%

    # gene ontology
    # ==========================================================================














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


