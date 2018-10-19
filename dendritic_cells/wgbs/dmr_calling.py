"""Pairwise DMR calls using DSS"""
import json
import tempfile
from subprocess import run
from pathlib import Path

from dmr_calling import create_dmr_metadata_table

# TODO: add meth. calling code to this project (extract from mouse_hematopoiesis umbrella project)
import mouse_hematopoiesis.wgbs.meth_calling2.meth_calling2 as mcalling
from dendritic_cells.base_paths import (
    wgbs_cohort_results_dir, project_temp_dir
)

dc_calls_config = {
    'chromosomes': sorted([str(i) for i in range(1, 2)]),
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

dc_dmr_metadata_table = create_dmr_metadata_table(dc_calls_config)

if __name__ == '__main__':

    with tempfile.TemporaryDirectory(dir=project_temp_dir) as tmpdir:
        config_fp = Path(tmpdir) / 'config.json'
        with open(config_fp, 'wt') as fout:
            json.dump(dc_calls_config, fout)
        run(f"""snakemake \
                --snakefile /home/kraemers/projects/dmr_calling/dmr-calling.snakefile \
                --configfile {config_fp} \
                --jobs 1000 \
                --latency-wait 180 \
                --jobscript /home/kraemers/projects/mouse_hematopoiesis/src/mouse_hematopoiesis/wgbs/dmr_calling2/jobscript.sh \
                --cluster "bsub -R rusage[mem={{params.avg_mem}}] -M {{params.max_mem}} -n {{threads}} -J {{params.name}} -W {{params.walltime}} -o /home/kraemers/temp/logs/" \
                """, shell=True, check=True)


