"""Clustering and heatmaps

Hiearchical clusterin
"""

# %%
# %load_ext autoreload
# %autoreload 2
# import matplotlib
# matplotlib.use('Agg')
# %%

from pathlib import Path
from dendritic_cells.base_paths import project_temp_dir, wgbs_cohort_results_dir

RESULTS_DIR = Path(wgbs_cohort_results_dir).joinpath(
        'analyses/dendritic_cells/dmr_characterization')
merged_dmrs_p = RESULTS_DIR / 'merged-dmrs.p'
merged_dmrs_tsv = RESULTS_DIR / 'merged-dmrs.tsv'
merged_dmrs_bed = RESULTS_DIR / 'merged-dmrs.bed'
merged_dmrs_bed_no_header = RESULTS_DIR / 'merged-dmrs_no-header.bed'
meth_stats_flat_tsv = str(RESULTS_DIR / 'meth-stats.tsv')
meth_stats_p = RESULTS_DIR / 'meth-stats.p'
dist_mat_npy = RESULTS_DIR / 'dist-mat_z-score_euclidean.npy'
linkage_mat_npy = RESULTS_DIR / 'linkage-mat_z-score_euclidean_ward.npy'
cluster_ids_tsv = RESULTS_DIR / 'clustering/cluster-ids.tsv'
cluster_ids_p = RESULTS_DIR / 'clustering/cluster-ids.p'
cluster_ids_obj_p = RESULTS_DIR / 'clustering/cluster-ids_obj.p'
dmr_heatmaps_by_name_pdf = str(RESULTS_DIR / 'dmr-heatmaps/dmr-heatmap_sampling-{name}_cluster-cols-{clustercols}_pop-order-{poporder}.pdf')
dmr_heatmaps_by_name_png = dmr_heatmaps_by_name_pdf.replace('.pdf', '.png')
sampled_cluster_ids_shelve = str(RESULTS_DIR / f'dmr-heatmaps/sampled-clusterid_shelve')

RESULTS_DIR.mkdir(exist_ok=True, parents=True)
RESULTS_DIR.joinpath('clustering').mkdir(parents=False, exist_ok=True)

if __name__ == '__main__':

    import pickle
    import shelve
    from itertools import product
    from tempfile import TemporaryDirectory
    from typing import List

    import pandas as pd
    from pandas.api.types import CategoricalDtype
    import numpy as np

    from dendritic_cells.wgbs.dmr_calling import dc_dmr_metadata_table
    from mouse_hematopoiesis.wgbs.meth_calling2.meth_calling2 import calls_v1_metadata_table
    import methlevels as ml
    import subprocess
    import codaplot as co
    from scipy.spatial.distance import pdist
    from scipy.cluster.hierarchy import linkage


    import os
    import time
    from textwrap import dedent
    from typing import Optional, Dict

    import methlevels as ml

    from mouse_hematopoiesis.wgbs.meth_calling2.meth_calling2 import calls_v1_metadata_table
    from mouse_hematopoiesis.config.base_paths import (
        project_temp_dir, wgbs_analysis_dir, result_paths_report_dir)
    from mouse_hematopoiesis.utils import (
        get_output_overview, git_abbrev_hash_of_head2, json_read_file,
        write_output_paths_file)

    # %%

    RESULTS_DIR.joinpath('dmr-heatmaps').mkdir(parents=False, exist_ok=True)

    def main():

        chromosome_cat_dtype = CategoricalDtype(
                categories=sorted(str(i) for i in range(1, 20)), ordered=True)

        analysis_tmpdir_obj = TemporaryDirectory(dir=project_temp_dir)
        analysis_tmpdir = analysis_tmpdir_obj.name

        # get merged dmrs
        # ==================================================================

        merged_dmrs = merge_pw_dmr_calls(dc_dmr_metadata_table['dmr_bed'].tolist(),
                                         tmpdir=analysis_tmpdir,
                                         chromosome_cat_dtype=chromosome_cat_dtype)
        merged_dmrs.to_pickle(merged_dmrs_p)
        merged_dmrs.to_csv(merged_dmrs_tsv, sep='\t', header=True, index=False)
        (merged_dmrs.rename(columns={'Chromosome': '#Chromosome'})
         .to_csv(merged_dmrs_bed, sep='\t', header=True, index=False))
        merged_dmrs.to_csv(merged_dmrs_bed_no_header, sep='\t',
                           header=False, index=False)

        # get filtered meth stats
        # ==================================================================
        # merged_dmrs = pd.read_pickle(merged_dmrs_p)

        pop_order = [
            "hsc", "mpp1", "mpp2", "mpp34",
            "mdp", "cmop", "cdp",
            "dc-cd11b", "dc-cd8a", "pdc",
        ]
        calls_metadata_table = calls_v1_metadata_table.query(
                'subject in @pop_order')

        meth_stats = get_meth_stats(analysis_tmpdir, calls_metadata_table,
                                    merged_dmrs, pop_order)
        meth_stats.save_flat([meth_stats_flat_tsv])
        with open(meth_stats_p, 'wb') as fout:
            pickle.dump(meth_stats, fout, protocol=4)

        # Compute linkage matrix
        # ==================================================================
        # with open(meth_stats_p, 'rb') as fin:
        #     meth_stats = pickle.load(fin)
        dist_mat = pdist(meth_stats.stats['beta-value_zscores'], metric='euclidean')
        linkage_mat = linkage(dist_mat, method='ward')
        np.save(dist_mat_npy, dist_mat)
        np.save(linkage_mat_npy, linkage_mat)

        # Compute dynamic cutree partitioning
        # ==================================================================
        # with open(meth_stats_p, 'rb') as fin:
        #     meth_stats = pickle.load(fin)
        # dist_mat = np.load(dist_mat_npy)
        # linkage_mat = np.load(linkage_mat_npy)

        linkage_obj = co.Linkage(matrix=linkage_mat, dist_mat=dist_mat,
                                 index=meth_stats.counts.index)
        min_cluster_size = int(len(meth_stats.counts) * 0.01)
        # print('ds1_nopam')
        # linkage.dynamic_tree_cut('ds1_nopam', deepSplit=1,
        #                          minClusterSize=min_cluster_size, pamStage=False)
        print('ds2_nopam')
        linkage_obj.dynamic_tree_cut('ds2_nopam', deepSplit=2,
                                     minClusterSize=min_cluster_size, pamStage=False)
        # print('ds25_pam')
        # linkage.dynamic_tree_cut('ds25_pam', deepSplit=2.5,
        #                          minClusterSize=min_cluster_size*2, pamStage=True)
        linkage_obj.cluster_ids.df.to_csv(cluster_ids_tsv,
                                          sep='\t', header=True, index=True)
        linkage_obj.cluster_ids.df.to_pickle(cluster_ids_p)
        with open(cluster_ids_obj_p, 'wb') as fout:
            pickle.dump(linkage_obj, fout, protocol=4)

        # Get MethStats
        # ===================================================================

        # Functions
        # ----------------------------------------------------------------------

        # +
        def create_meth_stats(bed_calls: ml.BedCalls,
                              intervals: pd.DataFrame,
                              output_dir: str,
                              filter_args: Optional[Dict] = None,
                              prefix=None,
                              version='v1'):
            """Create complete methylation stats

            Creates MethStats and individual dataframes for
            - population and replicate level
            - element and interval level
            - auxiliary stats, such as AUC, n_cpg_relevant

            Args:
                bed_calls: determines the samples for which meth stats are calculated
                filter_args: if None, defaults to {'coverage': 30, 'delta': 0.1}.
                    Note that without any filtering, there may be undefined beta values,
                    which would make filtering in downstream processing code, e.g. for
                    clustering, necessary.

            Possible improvements:
            - option to activate defensive check that beta values etc. are not NA
            """

            print('Working on', prefix)
            print('Writing to', output_dir)

            if filter_args is None:
                filter_args = dict(coverage=30, n_cpg_min=3, min_delta=0.1)
            os.makedirs(output_dir, exist_ok=True)

            output_paths = _get_output_paths(output_dir, prefix, version)

            print('Calculating Replicate-level MethStats')
            t1 = time.time()
            meth_stats_replevel_v1 = (bed_calls
                                      .intersect(intervals, n_cores=24, elements=True,
                                                 additional_index_cols=['region_id'],
                                                 drop_additional_index_cols=True)
                                      # .aggregate_element_counts()
                                      )
            print(f'done after {t1 - time.time()}')
            meth_stats_replevel_v1 = meth_stats_replevel_v1.aggregate_element_counts()

            print('Aggregating to Subject-level')
            meth_stats_poplevel_v1 = meth_stats_replevel_v1.convert_to_populations()
            if filter_args:
                print('Filtering intervals')
                meth_stats_poplevel_v1 = meth_stats_poplevel_v1.filter_intervals(
                        **filter_args)
                print('Apply Population-Stats based QC filtering to Rep-level data')
                meth_stats_replevel_v1 = (meth_stats_replevel_v1
                                          .subset(meth_stats_poplevel_v1.counts.index))

            # This assertion is mainly intended to check the subsetting of the rep level data
            rep_level_region_ids = meth_stats_replevel_v1.counts.index.get_level_values('region_id')
            pop_level_region_ids = meth_stats_poplevel_v1.counts.index.get_level_values('region_id')
            assert rep_level_region_ids.equals(pop_level_region_ids)


            print('Calculate methylation statistics')
            (meth_stats_poplevel_v1
             .add_beta_value_stat()
             .add_zscores('beta-value')
             .add_deltas(stat_name='beta-value', root='hsc')
             # add AUC
             # add peak delta
             )


            print('Save Rep-level data')
            with open(output_paths['rep-level']['meth_stats_obj'], 'wb') as fout:
                pickle.dump(meth_stats_replevel_v1, fout)
            meth_stats_replevel_v1.save_flat_elements_df(
                    output_paths['rep-level']['elements_p'],
                    output_paths['rep-level']['elements_tsv'],
                    output_paths['rep-level']['elements_feather'],
            )
            meth_stats_replevel_v1.save_flat_intervals_df(
                    output_paths['rep-level']['intervals_p'],
                    output_paths['rep-level']['intervals_tsv'],
                    output_paths['rep-level']['intervals_feather'],
            )

            print('Save pop-level data')
            with open(output_paths['pop-level']['meth_stats_obj'], 'wb') as fout:
                pickle.dump(meth_stats_poplevel_v1, fout)
            meth_stats_poplevel_v1.save_flat_elements_df(
                    output_paths['pop-level']['elements_p'],
                    output_paths['pop-level']['elements_tsv'],
                    output_paths['pop-level']['elements_feather'],
            )
            meth_stats_poplevel_v1.save_flat_intervals_df(
                    output_paths['pop-level']['intervals_p'],
                    output_paths['pop-level']['intervals_tsv'],
                    output_paths['pop-level']['intervals_feather'],
            )

            return output_paths


        def _get_output_paths(output_dir, prefix, version='v1'):
            if prefix is not None:
                prefix += '_'
            output_paths = {
                'rep-level': {
                    'meth_stats_obj': output_dir + f'/{prefix}meth-stats-obj_rep-level_{version}.p',
                    'elements_tsv':   output_dir + f'/{prefix}element-meth-stats-df_rep-level_{version}.tsv',
                    'elements_p':     output_dir + f'/{prefix}element-meth-stats-df_rep-level_{version}.p',
                    'elements_feather':     output_dir + f'/{prefix}element-meth-stats-df_rep-level_{version}.feather',
                    'intervals_tsv':  output_dir + f'/{prefix}interval-meth-stats-df_rep-level_{version}.tsv',
                    'intervals_p':    output_dir + f'/{prefix}interval-meth-stats-df_rep-level_{version}.p',
                    'intervals_feather':  output_dir + f'/{prefix}interval-meth-stats-df_rep-level_{version}.feather',
                },
                'pop-level': {
                    'meth_stats_obj': output_dir + f'/{prefix}meth-stats-obj_pop-level_{version}.p',
                    'elements_tsv':   output_dir + f'/{prefix}element-meth-stats-df_pop-level_{version}.tsv',
                    'elements_p':     output_dir + f'/{prefix}element-meth-stats-df_pop-level_{version}.p',
                    'elements_feather':     output_dir + f'/{prefix}element-meth-stats-df_pop-level_{version}.feather',
                    'intervals_tsv':  output_dir + f'/{prefix}interval-meth-stats-df_pop-level_{version}.tsv',
                    'intervals_p':    output_dir + f'/{prefix}interval-meth-stats-df_pop-level_{version}.p',
                    'intervals_feather':    output_dir + f'/{prefix}interval-meth-stats-df_pop-level_{version}.feather',
                }
            }
            return output_paths

        # -


        # all pops methstats
        # ==============================================================================

        all_pops_methstats_output_dir = '/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/cohort_results/analyses/dendritic_cells/dmr_characterization/meth_in_regions/all_pops'
        Path(all_pops_methstats_output_dir).mkdir(parents=True, exist_ok=True)
        methstats_prefix = 'all-pops'
        methstats_version = 'v1'

        # +
        all_pops = [
            "hsc", "mpp1", "mpp2", "mpp34",
            "cmop", "monos",
            "mdp", "cdp",
            "dc-cd11b", "dc-cd8a", "pdc",
        ]

        all_pops_metadata_table_v2 = calls_v1_metadata_table.query('subject in @all_pops')
        all_pops_bed_calls_v2 = ml.BedCalls(metadata_table=all_pops_metadata_table_v2,
                                             tmpdir=project_temp_dir,
                                             pop_order=all_pops,
                                             beta_value_col=6,
                                             n_meth_col=7,
                                             n_total_col=8)

        all_cluster_ids_p = '/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/cohort_results/analyses/dendritic_cells/dmr_characterization/clustering/cluster-ids.p'
        all_cluster_ids = pd.read_pickle(all_cluster_ids_p)
        regions = all_cluster_ids.reset_index().drop("ds2_nopam", axis=1)

        create_meth_stats(bed_calls=all_pops_bed_calls_v2,
                          intervals=regions,
                          output_dir=all_pops_methstats_output_dir,
                          prefix=methstats_prefix,
                          version=methstats_version,
                          filter_args={})




        # Heatmap
        # ==================================================================
        # with open(meth_stats_p, 'rb') as fin:
        #     meth_stats = pickle.load(fin)
        # with open(cluster_ids_obj_p, 'rb') as fin:
        #     linkage_obj = pickle.load(fin)

        # Computing the subpart linkage currently is slow, > 3 min
        d = shelve.open(sampled_cluster_ids_shelve, protocol=4)
        sampled_cluster_ids_prop = linkage_obj.cluster_ids.sample_proportional(
                'ds2_nopam', n_total=5000, random_state=1, min_cluster_size=50)
        sampled_index_prop = sampled_cluster_ids_prop.df.index
        sampled_linkage_prop = linkage_obj.get_subpart(labels=sampled_index_prop)
        d['sampled_cluster_ids_prop'] = sampled_cluster_ids_prop
        d['sampled_index_prop'] = sampled_index_prop
        d['sampled_linkage_obj_prop'] = sampled_linkage_prop
        # sampled_cluster_ids_prop = d['sampled_cluster_ids_prop']
        # sampled_index_prop = d['sampled_index_prop']
        # sampled_linkage_prop = d['sampled_linkage_obj_prop']

        pop_orders = [
            ['cmop', 'mdp', 'cdp'],
            ['cmop', 'mdp', 'cdp',
             'dc-cd11b', 'dc-cd8a', 'pdc'],
            ['hsc', 'mpp1', 'mpp2', 'mpp34',
             'cmop', 'mdp', 'cdp',
             'dc-cd11b', 'dc-cd8a', 'pdc'
             ]
        ]

        for pop_order in pop_orders:
            print(pop_order)
            sampled_heatmap(meth_stats=meth_stats,
                            sampled_cluster_ids=sampled_cluster_ids_prop,
                            cluster_id_name='ds2_nopam',
                            sampled_index=sampled_index_prop,
                            sampled_linkage=sampled_linkage_prop,
                            name='prop-sampling',
                            linkage_obj=linkage_obj,
                            pop_order=pop_order,
                            out_pattern=dmr_heatmaps_by_name_pdf)



        sampled_cluster_ids_equal = linkage_obj.cluster_ids.sample_equal(
                'ds2_nopam', n_per_cluster=500, random_state=1, strict=True)
        sampled_index_equal = sampled_cluster_ids_equal.df.index
        sampled_linkage_equal = linkage_obj.get_subpart(labels=sampled_index_equal)
        d['sampled_cluster_ids_equal'] = sampled_cluster_ids_equal
        d['sampled_index_equal'] = sampled_index_equal
        d['sampled_linkage_obj_equal'] = sampled_linkage_equal
        sampled_cluster_ids_equal = d['sampled_cluster_ids_equal']
        sampled_index_equal = d['sampled_index_equal']
        sampled_linkage_equal = d['sampled_linkage_obj_equal']


        for pop_order in pop_orders:
            print(pop_order)
            sampled_heatmap(meth_stats=meth_stats,
                            sampled_cluster_ids=sampled_cluster_ids_equal,
                            cluster_id_name='ds2_nopam',
                            sampled_index=sampled_index_equal,
                            sampled_linkage=sampled_linkage_equal,
                            name='equal-sampling',
                            linkage_obj=linkage_obj,
                            pop_order=pop_order,
                            out_pattern=dmr_heatmaps_by_name_pdf)

        all_pops_methstats_p = '/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/cohort_results/analyses/dendritic_cells/dmr_characterization/meth_in_regions/all_pops/all-pops_meth-stats-obj_pop-level_v1.p'
        with open(all_pops_methstats_p, 'rb') as fin:
            all_pops_methstats = pickle.load(fin)

        assert all_pops_methstats.counts.index.equals(meth_stats.counts.index)

        out_pattern =  str(RESULTS_DIR / 'dmr-heatmaps/dmr-heatmap_sampling-{name}_cluster-cols-{clustercols}_pop-order-{poporder}_with-monos.pdf')
        sampled_heatmap(meth_stats=all_pops_methstats,
                        sampled_cluster_ids=sampled_cluster_ids_equal,
                        cluster_id_name='ds2_nopam',
                        sampled_index=sampled_index_equal,
                        sampled_linkage=sampled_linkage_equal,
                        name='equal-sampling',
                        linkage_obj=linkage_obj,
                        pop_order=all_pops,
                        out_pattern=out_pattern)





    def sampled_heatmap(meth_stats, sampled_cluster_ids: co.Linkage, cluster_id_name,
                        sampled_index, sampled_linkage, name, linkage_obj, pop_order,
                        out_pattern):

        zscores_df = meth_stats.stats['beta-value_zscores'].loc[sampled_index, pop_order]
        betas_df = meth_stats.stats['beta-value'].loc[sampled_index, pop_order]
        deltas_df = meth_stats.stats['beta-value_delta_hsc'].loc[sampled_index, pop_order]

        cdg_cols_clustered = co.ClusteredDataGrid(
                main_df=zscores_df,
                row_linkage=sampled_linkage.matrix)
        cdg_cols_clustered = cdg_cols_clustered.cluster_cols(
                method='complete', metric='cityblock')

        cdg_cols_ordered = co.ClusteredDataGrid(
                main_df=zscores_df,
                row_linkage=sampled_linkage.matrix)

        for cols_clustered, cdg in [(True, cdg_cols_clustered),
                                    (False, cdg_cols_ordered)]:
            print(cols_clustered)

            if name == 'equal_sampling':
                cluster_sizes = [co.ClusterSizePlot(
                        cluster_ids=linkage_obj.cluster_ids.df[cluster_id_name],
                        bar_height = 0.7, xlabel='Cluster size')]
            else:
                cluster_sizes = []
            # TODO: printing gm before create_or_update_figure is called:
            # -> AttributeError: 'GridManager' object has no attribute 'fig'
            gm = cdg.plot_grid(grid=[
                [
                    co.Heatmap(df=zscores_df,
                               cmap='RdBu_r', rasterized=True),
                    co.Heatmap(df=deltas_df,
                               cmap='RdBu_r', rasterized=True),
                    co.Heatmap(df=betas_df,
                               cmap='YlOrBr', rasterized=True),
                ] + cluster_sizes
            ],
                    figsize=(30 / 2.54, 15 / 2.54),
                    height_ratios=[(1, 'rel')],
                    row_annotation=sampled_cluster_ids.df,
                    row_anno_heatmap_args={'colors':      [(1, 1, 1), (.8, .8, .8)],
                                           'show_values': True},
                    row_anno_col_width=2.5 / 2.54,
                    col_dendrogram=cols_clustered,
            )
            gm.create_or_update_figure()
            out_pdf = out_pattern.format(
                    name=name,
                    clustercols=str(cols_clustered),
                    poporder = ','.join(pop_order)
            )
            gm.fig.savefig(out_pdf)
            gm.fig.savefig(out_pdf.replace('.pdf', '.png'))



    def get_meth_stats(analysis_tmpdir, calls_metadata_table, merged_dmrs, pop_order):
        bed_calls = ml.BedCalls(metadata_table=calls_metadata_table,
                                tmpdir=analysis_tmpdir,
                                pop_order=pop_order,
                                beta_value_col=6,
                                n_meth_col=7,
                                n_total_col=8)
        # TODO-here: save the calls, per element, per interval, per rep, per pop
        # TODO-methlevels: the different qc filter functions do not copy the element meth stats
        meth_stats = bed_calls.intersect(merged_dmrs, n_cores=24, elements=True)
        meth_stats = meth_stats.convert_to_populations().aggregate_element_counts()
        # the next stop looses the element info because of the filtering functions
        meth_stats = (meth_stats
                      .qc_filter(coverage=30, size=3, min_delta=0.1)
                      .add_beta_value_stat()
                      .add_zscores('beta-value')
                      .add_deltas(stat_name='beta-value', root='hsc')
                      )
        return meth_stats


    def merge_pw_dmr_calls(bed_fps: List[str], tmpdir: str,
                           chromosome_cat_dtype) -> pd.DataFrame:
        """Merge pairwise DMR calls (given ad BED files"""
        # restrict input types to avoid errors with caching
        assert isinstance(bed_fps, list)
        tmpdir = TemporaryDirectory(dir=tmpdir)
        output_fp = tmpdir.name + 'merged.bed'
        _ = subprocess.run(
                f"cat {' '.join(bed_fps)} | sort -k1,1 -k2,2n -k3,3n"
                f" | bedtools merge > {output_fp}",
                shell=True, executable='/bin/bash', check=True)
        merged_dmrs = pd.read_csv(output_fp, sep='\t', header=None,
                                  names=['Chromosome', 'Start', 'End'],
                                  dtype={'Chromosome': chromosome_cat_dtype,
                                         'Start': 'i8', 'End': 'i8'})
        # insert region id for CpG-level retrieval with BedCalls and use in MethStats
        merged_dmrs['region_id'] = np.arange(len(merged_dmrs), dtype='i8')
        # Assert: is sorted, has DMRs for every chromosome
        assert merged_dmrs['Chromosome'].is_monotonic_increasing
        assert merged_dmrs.groupby('Chromosome')['Start'].apply(lambda ser: ser.is_monotonic_increasing).all()
        assert merged_dmrs.groupby('Chromosome')['End'].apply(lambda ser: ser.is_monotonic_increasing).all()
        assert np.all(merged_dmrs['Chromosome'].unique() == chromosome_cat_dtype.categories)
        return merged_dmrs

    main()

