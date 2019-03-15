import pickle
import pandas as pd
import region_set_profiler as rsp

cluster_ids_df = pd.read_pickle(snakemake.input.cluster_ids)
regions_df = cluster_ids_df.index.to_frame().reset_index(drop=True).drop('region_id', axis=1)
coverage_stats = rsp.OverlapStats(
        regions=regions_df,
        tmpdir=snakemake.params.analysis_tmpdir,
        chromosomes=snakemake.params.chromosomes,
        metadata_table=snakemake.input.metadata_table)
coverage_stats.compute(snakemake.threads)
with open(snakemake.output[0], 'wb') as fout:
    pickle.dump(coverage_stats, fout, protocol=4)

