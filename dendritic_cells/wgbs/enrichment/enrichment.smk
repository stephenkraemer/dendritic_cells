"""

snakemake(
snakefile='/home/kraemers/projects/dendritic_cells/dendritic_cells/wgbs/enrichment/enrichment.smk',
cores=24,
nodes=2000,
dryrun=False,
forcerun=['compute_overlap_stats'],
latency_wait=60,
jobscript='/home/kraemers/projects/dendritic_cells/dendritic_cells/wgbs/enrichment/jobscript.sh',
cluster="bsub -R rusage[mem={params.avg_mem}] -M {params.max_mem} -n {threads} -J {params.name} -W {params.walltime} -o /home/kraemers/temp/logs/",
)

"""
import pickle
import pandas as pd
from glob import glob
from pathlib import Path
from tempfile import TemporaryDirectory
import region_set_profiler as rsp
from dendritic_cells.base_paths import project_temp_dir, wgbs_cohort_results_dir

cluster_ids = '/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/cohort_results/analyses/dendritic_cells/dmr_characterization/clustering/cluster-ids.p'
codex_bed_files_dir = '/icgc/dkfzlsdf/analysis/hs_ontogeny/databases/enrichment_databases/lola_chipseq_2018-04-12/mm10/codex/regions'
encode_bed_files_dir = '/icgc/dkfzlsdf/analysis/hs_ontogeny/databases/enrichment_databases/lola_chipseq_2018-04-12/mm10/encodeTFBSmm10/regions'
bed_files_dict = dict(
    codex = glob(codex_bed_files_dir + '/*.bed'),
    encode = glob(codex_bed_files_dir + '/*.bed'),
    homer = glob('/icgc/dkfzlsdf/analysis/hs_ontogeny/databases/enrichment_databases/homer/homer_lola/*/regions/*.bed'),
    msigdb_oncogenic = glob('/icgc/dkfzlsdf/analysis/hs_ontogeny/databases/region_set_profiler_databases/msigdb/oncogenic_signatures'),
    msigdb_hallmarks = glob('/icgc/dkfzlsdf/analysis/hs_ontogeny/databases/region_set_profiler_databases/msigdb/hallmarks'),
)
metadata_table_dict = dict(
    encode = encode_bed_files_dir + '/encode_annotations.csv',
    codex = codex_bed_files_dir + '/codex_annotations.csv',
    msigdb_oncogenic = '/icgc/dkfzlsdf/analysis/hs_ontogeny/databases/region_set_profiler_databases/msigdb/oncogenic-signatures_metadata-table.tsv',
    msigdb_hallmarks = '/icgc/dkfzlsdf/analysis/hs_ontogeny/databases/region_set_profiler_databases/msigdb/hallmarks_metadata-table.tsv',
    )


results_dir = Path(wgbs_cohort_results_dir).joinpath(
    'analyses/dendritic_cells/dmr_characterization/enrichments')
analysis_tempdir_obj = TemporaryDirectory(dir=project_temp_dir)
analysis_tmpdir = analysis_tempdir_obj.name
chromosomes = sorted([str(i) for i in range(1, 20)])

overlap_stats_by_db = str(results_dir / 'coverage-stats/coverage-stats_{database}.p')
cluster_overlap_stats_by_db = str(results_dir / 'coverage-stats/coverage-counts_{database}.p')
test_statistics_by_db_test = str(results_dir / 'tests/{database}_{test}.p')

# databases = ['codex', 'encode', 'homer']
databases = ['msigdb_oncogenic', 'msigdb_hallmarks']
rule all:
    input:
        expand(overlap_stats_by_db, database=databases),
        expand(cluster_overlap_stats_by_db, database=databases),
        # expand(test_statistics_by_db_test, database=databases, test=['chi-square']),


rule compute_overlap_stats:
    input:
        # Required because the cluster ids define the regions where
        # enrichment is allowed
        cluster_ids=cluster_ids,
        # TODO: this is not correct - it takes all bed files, even if they are not
        # read, because they are not in the metadata table
        bed_files=lambda wildcards: bed_files_dict[wildcards.database],
        metadata_table = lambda wildcards: metadata_table_dict[wildcards.database],
    params:
        analysis_tmpdir=analysis_tmpdir,
        chromosomes=chromosomes,
        walltime='00:10',
        max_mem=32000,
        avg_mem=16000,
        name='coverage-stats_{database}',
    threads: 24
    output:
        overlap_stats_by_db,
    script:
        'get_overlap_stats.py'

localrules: compute_cluster_overlap_stats
rule compute_cluster_overlap_stats:
    input:
        cluster_ids=cluster_ids,
        overlap_stats=overlap_stats_by_db,
    output:
        cluster_overlap_stats_by_db,
    run:
        with open(input.overlap_stats, 'rb') as fin:
            overlap_stats = pickle.load(fin)
        cluster_ids_df = pd.read_pickle(input.cluster_ids)
        cluster_counts = overlap_stats.aggregate(
            cluster_ids=cluster_ids_df.iloc[:, 0].reset_index('region_id', drop=True),
            min_counts=20)
        with open(output[0], 'wb') as fout:
          pickle.dump(cluster_counts, fout)


localrules: chi_square_test
rule chi_square_test:
    input:
        cluster_overlap_stats_by_db,
    output:
        expand(test_statistics_by_db_test, database='{database}', test='chi-square'),
    run:
        with open(input[0], 'rb') as fin:
            cluster_overlap_stats = pickle.load(fin)
        cluster_overlap_stats.hits += 5
        cluster_overlap_stats.cluster_sizes += 5
        test_stats_df = cluster_overlap_stats.test_for_enrichment(method='chi_square')
        test_stats_df.to_pickle(output[0])


# rule generalized_fisher_test:
#     input:
#         cluster_overlap_stats_by_db,
#     output:
#         expand(test_statistics_by_db_test, database='{database}', test='gen-fisher'),
#     run:
