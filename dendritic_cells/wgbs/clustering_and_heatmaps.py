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
dmr_heatmaps_by_name_pdf = str(RESULTS_DIR / 'dmr-heatmaps/dmr-heatmap_{name}.pdf')
dmr_heatmaps_by_name_png = str(RESULTS_DIR / 'dmr-heatmaps/dmr-heatmap_{name}.png')
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


        # Heatmap
        # ==================================================================
        # with open(meth_stats_p, 'rb') as fin:
        #     meth_stats = pickle.load(fin)
        # with open(cluster_ids_obj_p, 'rb') as fin:
        #     linkage_obj = pickle.load(fin)
        d = shelve.open(sampled_cluster_ids_shelve, protocol=4)

        # Computing the subpart linkage currently is slow, > 3 min
        sampled_cluster_ids_prop = linkage_obj.cluster_ids.sample_proportional(
                'ds2_nopam', n_total=5000, random_state=1, min_cluster_size=50)
        sampled_index_prop = sampled_cluster_ids_prop.df.index
        sampled_linkage_prop = linkage_obj.get_subpart(labels=sampled_index_prop)
        d['sampled_cluster_ids_prop'] = sampled_cluster_ids_prop
        d['sampled_index_prop'] = sampled_index_prop
        d['sampled_linkage_obj_prop'] = sampled_linkage_prop

        sampled_heatmap(meth_stats, sampled_cluster_ids_prop, 'ds2_nopam',
                        sampled_index_prop, sampled_linkage_prop, 'prop-sampling',
                        linkage_obj)



        sampled_cluster_ids_equal = linkage_obj.cluster_ids.sample_equal(
                'ds2_nopam', n_per_cluster=500, random_state=1, strict=True)
        sampled_index_equal = sampled_cluster_ids_equal.df.index
        sampled_linkage_equal = linkage_obj.get_subpart(labels=sampled_index_equal)
        d['sampled_cluster_ids_prop'] = sampled_cluster_ids_prop
        d['sampled_index_prop'] = sampled_index_prop
        d['sampled_linkage_obj_prop'] = sampled_linkage_prop
        sampled_heatmap(meth_stats, sampled_cluster_ids_equal, 'ds2_nopam',
                        sampled_index_equal, sampled_linkage_equal, 'equal_sampling',
                        linkage_obj)


    def sampled_heatmap(meth_stats, sampled_cluster_ids: co.Linkage, cluster_id_name,
                        sampled_index, sampled_linkage, name, linkage_obj):
        cdg = (co.ClusteredDataGrid(main_df=meth_stats.stats['beta-value_zscores'].loc[sampled_index],
                                    row_linkage=sampled_linkage.matrix)
               .cluster_cols(method='complete', metric='cityblock'))

        if name == 'equal_sampling':
            cluster_sizes = [co.ClusterSizePlot(cluster_ids=linkage_obj.cluster_ids.df[cluster_id_name],
                                                bar_height = 0.7, xlabel='Cluster size')]
        else:
            cluster_sizes = []
        # TODO: printing gm before create_or_update_figure is called:
        # -> AttributeError: 'GridManager' object has no attribute 'fig'
        gm = cdg.plot_grid(grid=[
            [
                co.Heatmap(df=meth_stats.stats['beta-value_zscores'].loc[sampled_index],
                           cmap='RdBu_r', rasterized=True),
                co.Heatmap(df=meth_stats.stats['beta-value_delta_hsc'].loc[sampled_index],
                           cmap='RdBu_r', rasterized=True),
                co.Heatmap(df=meth_stats.stats['beta-value'].loc[sampled_index],
                           cmap='YlOrBr', rasterized=True),
            ] + cluster_sizes
        ],
                figsize=(30 / 2.54, 15 / 2.54),
                height_ratios=[(1, 'rel')],
                row_annotation=sampled_cluster_ids.df,
                row_anno_heatmap_args={'colors':      [(1, 1, 1), (.8, .8, .8)],
                                       'show_values': True},
                row_anno_col_width=2.5 / 2.54,
                col_dendrogram=True,
        )
        gm.create_or_update_figure()
        gm.fig.savefig(dmr_heatmaps_by_name_pdf.format(name=name))
        gm.fig.savefig(dmr_heatmaps_by_name_png.format(name=name))



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

