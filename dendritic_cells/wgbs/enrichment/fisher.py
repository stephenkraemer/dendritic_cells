import pickle

with open(snakemake.input[0], 'rb') as fin:
    cluster_overlap_stats = pickle.load(fin)

test_stats_df = cluster_overlap_stats.test_for_enrichment(
        method='fisher', test_args=dict(
                simulate_pval=True, replicate=1e6, workspace=300,
                seed=1))
test_stats_df.to_pickle(snakemake.output[0])
