"""Pairwise DMR calls using DSS"""
# %%
import json
import os
import tempfile
from pathlib import Path
from pkg_resources import resource_filename

import dss_workflow

# TODO: add meth. calling code to this project (extract from mouse_hematopoiesis umbrella project)
import mouse_hematopoiesis.wgbs.meth_calling2.meth_calling2 as mcalling
from dendritic_cells.base_paths import (
    wgbs_cohort_results_dir, project_temp_dir, snakemake_conda_prefix,
    wgbs_metadata_dir,
)
# %%

dc_calls_config = {
    'chromosomes': sorted([str(i) for i in range(1, 20)]),
    'output_dir': str(Path(wgbs_cohort_results_dir).joinpath(
            f'pairwise-dmr-calls/DSS/{mcalling.calls_v1_str}')),
    'metadata_table': mcalling.calls_v1_paths['metadata_table'],
    'group_column': 'subject',
    'sample_id_column': 'sample_id',
    'bed_by_chrom_column': 'bed_by_chrom_path',
    'comparisons': [
        ('hsc', 'mpp1'),
        ('hsc', 'mpp2'),
        ('hsc', 'mpp34'),
        ('hsc', 'mdp'),
        ('hsc', 'cdp'),
        ('hsc', 'pdc'),
        ('hsc', 'dc-cd11b'),
        ('hsc', 'dc-cd8a'),
        ('hsc', 'cmop'),
    ],
    'parameters': [{
        'pvalue': 0.01,
        'delta': 0.1,
        'minlen': 50,
        'minCG': 2,
        'merge_dist': 50,
        'pct_sign': 0.5,
        'smooth': False,
        'smooth_span': 500,
    }],
    'create_dmr_coverage_bedgraph': False,
    'cpg_index_file': mcalling.calls_v1_paths['cg_index_path'],
}

dc_dmr_metadata_table = dss_workflow.create_dmr_metadata_table(dc_calls_config)
dc_dmr_metadata_table_fp = os.path.join(wgbs_metadata_dir, 'final_dmr_calls.tsv')

if __name__ == '__main__':
    import snakemake

    with tempfile.TemporaryDirectory(dir=project_temp_dir) as tmpdir:
        config_fp = Path(tmpdir) / 'config.json'
        with open(config_fp, 'wt') as fout:
            json.dump(dc_calls_config, fout)
        dss_snakefile = dss_workflow.get_snakefile_path()
        snakemake_succeeded = snakemake.snakemake(
                snakefile=dss_snakefile,
                configfile=config_fp,
                nodes=1000,
                latency_wait=180,
                jobscript=resource_filename('dendritic_cells', 'snakemake_jobscript.sh'),
                cluster=("bsub -R rusage[mem={params.avg_mem}] -M {params.max_mem} "
                         "-n {threads} -J {params.name} -W {params.walltime} "
                         "-o /home/kraemers/temp/logs/"),
                conda_prefix=snakemake_conda_prefix,
                use_conda=True,
                forceall=True,
                dryrun=False)
        if not snakemake_succeeded:
            raise RuntimeError('Execution of DSS workflow failed')

    Path(dc_dmr_metadata_table_fp).parent.mkdir(exist_ok=True, parents=True)
    dc_dmr_metadata_table.to_csv(dc_dmr_metadata_table_fp, sep='\t', header=True,
                                 index=False)


