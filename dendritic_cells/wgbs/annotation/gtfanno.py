import gtfanno
import pandas as pd
from tempfile import TemporaryDirectory
from mouse_hematopoiesis.config.base_paths import project_temp_dir, wgbs_cohort_results_dir
from pathlib import Path



cluster_ids = pd.read_pickle('/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/cohort_results/analyses/dendritic_cells/dmr_characterization/clustering/cluster-ids.p')

output_trunk_path = (wgbs_cohort_results_dir
                     + '/analyses/dendritic_cells/dmr_characterization/annotation/gtfanno/gtfanno')
Path(output_trunk_path).parent.mkdir(parents=True, exist_ok=True)

with TemporaryDirectory(dir=project_temp_dir) as tmpdir:
    cluster_ids.to_csv(tmpdir + '/query.bed', sep='\t', header=False, index=True)
    gtfanno.annotate(query_bed=tmpdir + '/query.bed',
                     gtf_fp=("/icgc/dkfzlsdf/analysis/hs_ontogeny/databases/gene_annotations"
                             "/gencode.vM19.annotation.no-prefix.gtf"),
                     trunk_path=output_trunk_path,
                     tmpdir=project_temp_dir,
                     promoter = (-5000, 1000),
                     distant_cis_regulatory_domain = (-50000, -5000),
                     )
