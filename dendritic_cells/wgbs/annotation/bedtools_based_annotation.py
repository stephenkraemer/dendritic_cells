from subprocess import run
from pathlib import Path
from dendritic_cells.wgbs.clustering_and_heatmaps import merged_dmrs_bed, merged_dmrs_p
from dendritic_cells.base_paths import project_temp_dir, wgbs_cohort_results_dir

results_dir = Path(wgbs_cohort_results_dir).joinpath(
        'analyses/dendritic_cells/dmr_characterization/annotation/bedtools-based')
results_dir.mkdir(parents=True, exist_ok=True)
gencode_anno_p = results_dir.joinpath('merged-dmrs_gene-annotation_gencode_vM19.p')
gencode_anno_tsv = results_dir.joinpath('merged-dmrs_gene-annotation_gencode_vM19.tsv')

def main():
    tmpdir_obj = TemporaryDirectory(dir=project_temp_dir)
    tmpdir_path = Path(tmpdir_obj.name)

    gencode_dir = Path("/icgc/dkfzlsdf/analysis/hs_ontogeny/databases/gene_annotations")

    gencode19_gtf_appris_principal = (
        "/icgc/dkfzlsdf/analysis/hs_ontogeny/databases/gene_annotations"
        "/gencode.vM19.annotation.appris-principal.no-prefix.gtf")
    gencode19_gtf = (
        "/icgc/dkfzlsdf/analysis/hs_ontogeny/databases/gene_annotations"
        "/gencode.vM19.annotation.no-prefix.gtf")

    all_tss_area_bed = gencode_dir.joinpath('all-tss_slop-100000-5000.bed')

    strand_dtype = CategoricalDtype(['+', '-'], ordered=True)

    # %% compute tss area intersects, ~ 30s
    gencode_df = pd.read_csv(
            gencode19_gtf, sep='\t', header=None, comment='#',
            names=['feat_chrom', 'source', 'feature', 'Start', 'End',
                   'score', 'feat_strand', 'frame', 'attribute'],
            dtype={'feat_chrom': str, 'Start': 'i8', 'End': 'i8', 'feat_strand': strand_dtype})

    tss_n_upstream = 100_000
    tss_n_downstream = 5000
    transcripts = gencode_df.query('feature == "transcript"').copy()
    on_plus_strand = transcripts['feat_strand'] == '+'
    transcripts['TSS'] = -1
    transcripts['feat_start'] = -1
    transcripts['feat_end'] = -1
    transcripts['feat_class'] = 'TSS_area'
    transcripts = expand_gtf_attributes(transcripts)

    # custom slop
    transcripts.loc[on_plus_strand, 'TSS'] = transcripts.loc[on_plus_strand, 'Start']
    transcripts.loc[~on_plus_strand, 'TSS'] = transcripts.loc[~on_plus_strand, 'End']
    transcripts.loc[on_plus_strand, 'feat_start'] = transcripts.loc[on_plus_strand, 'Start'] - tss_n_upstream
    transcripts.loc[on_plus_strand, 'feat_end'] = transcripts.loc[on_plus_strand, 'Start'] + tss_n_downstream
    transcripts.loc[~on_plus_strand, 'feat_start'] = transcripts.loc[~on_plus_strand, 'End'] - tss_n_downstream
    transcripts.loc[~on_plus_strand, 'feat_end'] = transcripts.loc[~on_plus_strand, 'End'] + tss_n_upstream
    transcripts = transcripts.sort_values(['feat_chrom', 'feat_start', 'feat_end', 'TSS'])
    transcripts.loc[transcripts['feat_start'].lt(0), 'feat_start'] = 0
    transcripts_cols = ['feat_chrom', 'feat_start', 'feat_end', 'TSS', 'feat_strand',
                        'feat_class', 'gene_name', 'gene_id', 'transcript_id', 'appris_principal_score']
    transcripts[transcripts_cols].to_csv(all_tss_area_bed, sep='\t', header=False, index=False)

    all_tss_area_bt = BedTool(str(all_tss_area_bed))
    merged_dmrs_bt = BedTool(str(merged_dmrs_bed))
    tss_intersect_bt = merged_dmrs_bt.intersect(all_tss_area_bt, wa=True, wb=True)
    tss_intersect_df = pd.read_csv(tss_intersect_bt.fn, sep='\t',
                                   names=['Chromosome', 'Start', 'End', 'region_id']
                                         + transcripts_cols)
    tss_intersect_df['perc_feature'] = np.nan
    tss_intersect_df['perc_region'] = np.nan
    tss_intersect_df['distance'] = -1e8
    tss_intersect_df['center'] = tss_intersect_df.eval('Start + (End - Start)/2')
    tss_intersect_df['feat_center'] = np.nan
    tss_intersect_df['has_center'] = False
    tss_intersect_df['distance'] = tss_intersect_df.eval('center - TSS')
    assert tss_intersect_df['distance'].ne(-1e8).all()
    # tss_intersect_df.loc[tss_intersect_df.eval('Start <= TSS <= End'), 'distance'] = 0
    # tss_intersect_df.loc[tss_intersect_df.eval('End < TSS'), 'distance'] = tss_intersect_df.eval('End - TSS')
    # tss_intersect_df.loc[tss_intersect_df.eval('Start > TSS'), 'distance'] = tss_intersect_df.eval('Start - TSS')


    full_cols = ['Chromosome', 'Start', 'End', 'region_id', 'center',
                 'feat_class', 'perc_feature', 'perc_region', 'distance', 'has_center',
                 'gene_name','gene_id', 'transcript_id', 'appris_principal_score',
                 'feat_chrom', 'feat_start', 'feat_end', 'feat_center', 'feat_strand']
    tss_intersect_df_full = tss_intersect_df[full_cols]




    # %% compute exon, intron overlap, ~45s
    transcript_parts = gencode_df.loc[~gencode_df['feature'].isin(['gene', 'start_codon', 'stop_codon']), :]
    transcript_parts_fp = tmpdir_path.joinpath('transcript_parths.gtf')
    transcript_parts.to_csv(transcript_parts_fp, sep='\t', header=False, index=False)

    transcript_parts_bt = BedTool(str(transcript_parts_fp))
    transcript_parts_anno = merged_dmrs_bt.intersect(transcript_parts_bt, wa=True, wb=True)
    transcript_parts_anno.head()

    transcript_parts_df = pd.read_csv(transcript_parts_anno.fn, sep='\t', header=None)
    transcript_parts_df.columns = ['Chromosome', 'Start', 'End', 'region_id'] + [
        'feat_chrom', 'source', 'feat_class', 'feat_start', 'feat_end',
        'score', 'feat_strand', 'frame', 'attribute']
    start = transcript_parts_df.eval('Start - feat_start').where(lambda ser: ser.gt(0), 0)
    feat_size = transcript_parts_df.eval('feat_end - feat_start')
    end = transcript_parts_df.eval('End - feat_start').where(lambda ser: ser.lt(feat_size), feat_size)
    overlap_size = end - start
    region_size = transcript_parts_df.eval('End - Start')

    transcript_parts_df['center'] = transcript_parts_df.eval('Start + (End - Start)/2')
    transcript_parts_df['feat_center'] = transcript_parts_df.eval('feat_start + (feat_end - feat_start)/2')
    transcript_parts_df['distance'] = transcript_parts_df.eval('center - feat_center')
    transcript_parts_df['has_center'] = transcript_parts_df['distance'].lt(feat_size / 2)

    transcript_parts_df['perc_feature'] = overlap_size / feat_size
    transcript_parts_df['perc_region'] = overlap_size / region_size
    transcript_parts_df['distance'] = np.nan

    transcript_parts_df = expand_gtf_attributes(transcript_parts_df)

    transcript_parts_df_full = transcript_parts_df[full_cols]


    # %% classify into proximal and distal cis regulatory regions
    promoter_anno = tss_intersect_df_full.copy()
    is_proximal_promoter = promoter_anno.eval('-5000 <= distance <= 1000')
    is_distant_cis_regulatory_domain = promoter_anno.eval('-20000 <= distance < -5000')
    promoter_anno['feat_class'] = np.nan
    promoter_anno.loc[is_proximal_promoter, 'feat_class'] = 'Promoter'
    promoter_anno.loc[is_proximal_promoter, 'has_center'] = True
    promoter_anno.loc[is_distant_cis_regulatory_domain, 'feat_class'] = 'UCRD'
    promoter_anno.loc[is_distant_cis_regulatory_domain, 'has_center'] = True
    promoter_anno = promoter_anno.loc[~promoter_anno['feat_class'].isna(), :]

    # %% concatenate and type casts
    full_annos = (pd.concat([promoter_anno, transcript_parts_df_full], axis=0)
                  .sort_values(['Chromosome', 'Start', 'End']))

    precedence = pd.Series([
        'start_codon',
        'stop_codon',
        'Promoter',
        'UTR',
        'exon',
        'CDS',
        'UCRD',
        'transcript',
    ])

    # %% Filter according to precedence
    def filter_annotations(group_df):
        highest_class = precedence.iloc[precedence.isin(group_df['feat_class']).idxmax()]
        class_df = group_df.query('feat_class == @highest_class').sort_values(['appris_principal_score', 'perc_region'])
        # TODO: sort by appris score, then by overlap
        if highest_class == 'transcript':
            class_df['feat_class'] = 'intron'
        if class_df['gene_name'].nunique() == 1:
            return class_df.iloc[[0], :]
        else:
            return class_df.groupby('gene_name', as_index=False).nth(0)

    center_annos = full_annos.loc[full_annos['has_center'], :]

    # filter, takes ~
    # could maybe be sped up by removing more introns, perhaps with better intron annotation?
    t1 = time.time()
    filtered_annos = center_annos.groupby('region_id', group_keys=False).apply(filter_annotations)
    print(time.time() - t1)
    # cores = 24
    # filtered_annos_l = Parallel(cores)(delayed(filter_annotations)(group_df) for unused_name, group_df in grouped)
    # filtered_annos = pd.concat(filtered_annos_l, axis=0)

    filtered_annos.to_pickle(results_dir.joinpath('filtered-annos_no-intergenic.p'))
    # filtered_annos = pd.read_pickle(results_dir.joinpath('filtered-annos_no-intergenic.p'))

    ids_annotated_regions = filtered_annos['region_id'].unique()
    merged_dmrs_df = pd.read_pickle(merged_dmrs_p)
    intergenic_regions = merged_dmrs_df.loc[~merged_dmrs_df['region_id'].isin(ids_annotated_regions), :].copy()
    intergenic_regions['feat_class'] = 'intergenic'

    # error: chromosome dtypes are different
    all_regions_annotated = (pd.concat([filtered_annos, intergenic_regions],
                                       sort=False, axis=0))
    all_regions_annotated['Chromosome'] = all_regions_annotated['Chromosome'].astype(str)
    all_regions_annotated.sort_values(['Chromosome', 'Start', 'End'], inplace=True)
    assert (all_regions_annotated['region_id'].unique() == np.arange(53231)).all()
    all_regions_annotated['region_id'].value_counts().value_counts()
    all_regions_annotated['feat_class'].value_counts()

    all_regions_annotated.to_pickle(gencode_anno_p)
    all_regions_annotated.to_csv(gencode_anno_tsv, sep='\t', header=True, index=False)

    filtered_annos['feat_class'].value_counts()

    gencode_df_w_attributes = expand_gtf_attributes(gencode_df)
    principal_transcripts = gencode_df_w_attributes.query('appris_principal_score > 0 and feature == "transcript"').copy()
    # get TSS
    tss_on_plus_strand = principal_transcripts['feat_strand'].eq('+')
    principal_transcripts.loc[tss_on_plus_strand, 'End'] = principal_transcripts.loc[tss_on_plus_strand, 'Start'] + 1
    principal_transcripts.loc[~tss_on_plus_strand, 'Start'] = principal_transcripts.loc[~tss_on_plus_strand, 'End'] - 1
    principal_transcripts = principal_transcripts.sort_values(['feat_chrom', 'Start', 'End'])
    principal_transcripts_fp = tmpdir_path / 'principal-transcripts.gtf'
    principal_transcripts.iloc[:, 0:9].to_csv(principal_transcripts_fp, sep='\t', header=False, index=False)


    # pybedtools.featurefuncs.TSS has a bug

    gtf_princ_tss_bt = BedTool(str(principal_transcripts_fp))
    closest_tss_bt = merged_dmrs_bt.closest(gtf_princ_tss_bt, D='b', fu=True, t='first')
    closest_tss_df = pd.read_csv(closest_tss_bt.fn, sep='\t', header=None)
    distances = closest_tss_df.iloc[:, -1]
    distances = distances.loc[(distances < 100_000) & (distances > -100_000)]

    import matplotlib as mpl
    from matplotlib.axes import Axes # for autocompletion in pycharm
    from matplotlib.figure import Figure  # for autocompletion in pycharm
    mpl.use('Agg') # import before pyplot import!
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    import seaborn as sns

    fig, ax = plt.subplots(1, 1)
    sns.distplot(distances.values, bins=1000, kde=False, ax=ax)
    fig.savefig(results_dir.joinpath('tss-distance-dist.png'))
    fig.savefig(results_dir.joinpath('tss-distance-dist.pdf'))





    # merge missing regions

    # if everything is from the same gene, return the hit with the highest overlap (or from the canonical transcript)

    # df1 = pd.DataFrame({'a': pd.Categorical(['1', '2'])})
    # df2 = pd.DataFrame({'a': pd.Categorical(['2', '3'])})
    # df1.dtypes
    # df2.dtypes
    # res = pd.concat([df1, df2], axis=0)
    # res.dtypes
    # res['a']

    # full_annos_typed = full_annos.assign(
    #         Chromosome=pd.Categorical(full_annos['Chromosome'].astype(str), ordered=True),
    #         feat_class=pd.Categorical(full_annos['feat_class'], categories=precedence, ordered=True),
    #         gene_name=pd.Categorical(full_annos['gene_name'], ordered=True),
    #         gene_id=pd.Categorical(full_annos['gene_id'], ordered=True),
    #         transcript_id=pd.Categorical(full_annos['transcript_id'], ordered=True),
    #         feat_strand=full_annos['feat_strand'].astype(strand_dtype)
    # )
    #
    #

def expand_gtf_attributes(df):
    df['gene_id'] = df['attribute'].str.extract('gene_id "(.*?)";')
    df['transcript_id'] = df['attribute'].str.extract('transcript_id "(.*?)";')
    df['gene_name'] = df['attribute'].str.extract('gene_name "(.*?)";')
    df['appris_principal_score'] = df['attribute'].str.extract('tag "appris_principal_(\d)";').astype(float)
    df['appris_principal_score'] = df['appris_principal_score'].fillna(0)
    return df


if __name__ == '__main__':
    import os
    import json
    import subprocess
    from tempfile import TemporaryDirectory
    import pandas as pd
    from pandas.api.types import CategoricalDtype
    import numpy as np
    from pybedtools import BedTool
    from joblib import Parallel, delayed

    main()














