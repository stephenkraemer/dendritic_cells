from pathlib import Path
from dendritic_cells.wgbs.clustering_and_heatmaps import merged_dmrs_bed_no_header
from dendritic_cells.base_paths import project_temp_dir, wgbs_cohort_results_dir

results_dir = Path(wgbs_cohort_results_dir).joinpath(
        'analyses/dendritic_cells/dmr_characterization/annotation/test')
results_dir.mkdir(parents=True, exist_ok=True)

prefix = 'merged-dmrs_no-header'
final_results = results_dir.joinpath(f'{prefix}_finalhits.txt')
allhits = results_dir.joinpath(f'{prefix}_allhits.txt')
anno_df_p = results_dir.joinpath(f'{prefix}_anno.p')
anno_df_tsv = results_dir.joinpath(f'{prefix}_anno.tsv')

if __name__ == '__main__':

    import os
    import json
    import subprocess
    from tempfile import TemporaryDirectory
    import pandas as pd
    import numpy as np

    gencode19_gtf = "/icgc/dkfzlsdf/analysis/hs_ontogeny/databases/gene_annotations/gencode.vM19.annotation.appris-principal.gtf"

    uropa_config = {
        "queries": [
            {
                "feature": "Selenocysteine",
                "show.attributes": ["gene_name", "gene_id", "transcript_id"],
                "internals": "any",
                "distance": 0,
                "feature.anchor": "center",
            },
            {
                "feature": "start_codon",
                "show.attributes": ["gene_name", "gene_id", "transcript_id"],
                "internals": "any",
                "distance": 0,
                "feature.anchor": "center",
            },
            {
                "feature": "stop_codon",
                "show.attributes": ["gene_name", "gene_id", "transcript_id"],
                "internals": "any",
                "distance": 0,
                "feature.anchor": "center",
            },
            {
                'name': 'proximal_promoter',
                "feature": "transcript",
                "distance": [5000, 1000],
                "feature.anchor": "start",
                "show.attributes": ["gene_name", "gene_id", "transcript_id"],
                "internals": "none",
            },
            {
                'name': 'distal_promoter',
                "feature": "transcript",
                "distance": [50000, 1000],
                "feature.anchor": "start",
                "show.attributes": ["gene_name", "gene_id", "transcript_id"],
                "internals": "none",
            },
            {
                "feature": "UTR",
                "show.attributes": ["gene_name", "gene_id", "transcript_id"],
                "internals": "center",
                "distance": 0,
                "feature.anchor": "center",
            },
            {
                "feature": "exon",
                "show.attributes": ["gene_name", "gene_id", "transcript_id"],
                "internals": "center",
                "distance": 0,
                "feature.anchor": "center",
            },
            {
                "feature": "CDS",
                "show.attributes": ["gene_name", "gene_id", "transcript_id"],
                "internals": "center",
                "distance": 0,
                "feature.anchor": "center",
            },
            # intron
            {
                'name': 'intron',
                "feature": "transcript",
                "show.attributes": ["gene_name", "gene_id", "transcript_id"],
                "internals": "any",
                "distance": 0,
                "feature.anchor": "center",
            },
        ],
        "priority": "True",
        "gtf": gencode19_gtf,
        "bed": str(merged_dmrs_bed_no_header),
    }

    query_names = []
    for q in uropa_config['queries']:
        try:
            query_names.append(q.pop('name'))
        except KeyError:
            query_names.append(q['feature'])


    tmpdir = TemporaryDirectory(dir=project_temp_dir)
    tmpdir_path = Path(tmpdir.name)
    uropa_config_json = tmpdir_path.joinpath('uropa-config.json')
    with open(uropa_config_json, 'w') as fout:
        json.dump(uropa_config, fout)

    # uropa places results in the current working directory
    os.chdir(results_dir)
    subprocess.run(f'uropa -i {uropa_config_json} -p {prefix}'
                   ' -r -t 24', shell=True, executable='/bin/bash', check=True)

    final_results_df = pd.read_csv(final_results, sep='\t', header=0, index_col=0)
    anno_df = final_results_df.sort_values(
            ['peak_chr', 'peak_start', 'peak_end'])
    query_is_na = anno_df['query'].eq(','.join(np.arange(len(query_names)).astype(str)))
    anno_df['Region class'] = anno_df['query'].map({str(i): query_name
                                                    for i, query_name in enumerate(query_names)})
    anno_df.loc[query_is_na, 'Region class'] = 'NA'
    anno_df['Region class'].unique()
    anno_df.to_pickle(anno_df_p)
    anno_df.to_csv(anno_df_tsv, sep='\t', header=True, index=False)
