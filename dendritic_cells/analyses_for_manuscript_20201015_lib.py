"""
gtfanno version: b81595d
"""

import os
import pickle
import tempfile
from itertools import product
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import figure_report as fr
import gtfanno as ga
import matplotlib as mpl
import matplotlib.pyplot as plt
import mouse_hema_meth.paths as mhpaths
import mouse_hema_meth.shared_vars as mhvars
import mouse_hema_meth.styling as mhstyle
import numpy as np
import pandas as pd
import region_set_profiler as rsp
from typing_extensions import Literal

mpl.use("Agg")


PROJECT_TEMPDIR = "/icgc/dkfzlsdf/analysis/hs_ontogeny/temp"


def run_gtfanno(
    granges_df: pd.DataFrame,
    gtf_fp: str,
    output_trunk_path: str,
    distant_cis_regulatory_domain_def: Tuple[int, int],
    promoter_def: Tuple[int, int],
    recompute: bool,
) -> Dict[str, str]:
    """Run gtfanno

    Parameters
    ----------
    granges_df
        Chromosome Start End
        Chromosome should be Categorical string, in alphabetical sorting order
        not sure whether this is a hard requirement, but it is the convention
        I have always used so far

    Returns
    -------
        output paths dict
    """
    out_paths_d = dict(
        primary_annos_bed=output_trunk_path + "_primary-annotations.bed",
        primary_annos_p=output_trunk_path + "_primary-annotations.p",
        all_annos_bed=output_trunk_path + "_all-annotations.bed",
        all_annos_p=output_trunk_path + "_all-annotations.p",
    )
    if recompute:
        print("preparing query bed")
        temp_dir_obj = tempfile.TemporaryDirectory(dir=PROJECT_TEMPDIR)
        temp_dir_name = temp_dir_obj.name
        query_bed = temp_dir_name + "/query.bed"
        granges_df.to_csv(query_bed, sep="\t", header=False, index=False)
        print("Calling gtfanno.annotate")
        ga.annotate(
            query_bed=query_bed,
            gtf_fp=gtf_fp,
            trunk_path=output_trunk_path,
            tmpdir=temp_dir_name,
            promoter=promoter_def,
            distant_cis_regulatory_domain=distant_cis_regulatory_domain_def,
        )
    return out_paths_d


def merge_annos(
    gtfanno_result_fp: str, grange_and_feature_ids: pd.DataFrame
) -> pd.DataFrame:
    """Merge annotation format with multiple rows per DMR

    Parameters
    ----------
        gtfanno_result_fp
            path to the gtfanno result (*_primary-annotations.p)
        grange_and_feature_ids
            region_index (eg named region_id, cpg_id) // Chromosome Start End
            the region_index is the desired index, in the desired sorting order
            for the gtfanno result
            currently gtfanno drops original index, need this to merge it back


    Returns
    -------
        pd.DataFrame with columns ['Chromosome', 'Start', 'End', 'gene_name', 'feat_class'], and region_id index
        - Multiple gene annotations per DMR are concatenated with a comma.
    """

    assert grange_and_feature_ids.index.name is not None

    gtfanno_result: pd.DataFrame = pd.read_pickle(gtfanno_result_fp)

    gene_annos = (
        gtfanno_result
        # add gtfanno_uid to keep it in the index, even though it is redundant
        .groupby(["Chromosome", "Start", "End", "gtfanno_uid"], observed=True)[
            "gene_name"
        ].aggregate(lambda ser: ser.str.cat(sep=","))
    )

    # Assert that there is only one feature class per query region
    assert (
        gtfanno_result.groupby(
            ["Chromosome", "Start", "End", "gtfanno_uid"], observed=True
        )["feat_class"]
        .nunique()
        .eq(1)
        .all()
    )

    feat_class = gtfanno_result.groupby(
        ["Chromosome", "Start", "End", "gtfanno_uid"], observed=True
    )["feat_class"].first()

    gene_annos = pd.concat([gene_annos, feat_class], axis=1)
    # noinspection PyUnresolvedReferences
    assert (
        gene_annos.index.get_level_values("gtfanno_uid")
        == np.arange(gene_annos.shape[0])
    ).all()
    gene_annos.index = gene_annos.index.droplevel("gtfanno_uid")

    # merge and set feature_id column as index (eg region_id)
    gene_annos_with_region_id_idx = (
        pd.merge(
            gene_annos,
            grange_and_feature_ids.reset_index(),
            on=["Chromosome", "Start", "End"],
            how="inner",
        )
        .set_index(grange_and_feature_ids.index.name)
        .sort_values(["Chromosome", "Start", "End"])
    )

    # assert that the resulting gene annos are sorted
    pd.testing.assert_frame_equal(
        gene_annos_with_region_id_idx[["Chromosome", "Start", "End"]],
        grange_and_feature_ids,
    )

    # return gene_annos_with_region_id_idx.drop(["Chromosome", "Start", "End"], axis=1)
    return gene_annos_with_region_id_idx


def process_cluster_overlap_stats(
    cluster_overlap_stats,
    max_pvalues,
    plot_args,
    barcode_plot_png_by_maxpvalue,
    cluster_overlap_stats_out_fp,
    plotting_context,
    test_arg_dict=None,
    sample=None,
    cores=1,
    plots_only=False,
    additional_formats=None,
):
    if additional_formats is None:
        additional_formats = tuple()

    if not plots_only:
        if sample is not None and sample < cluster_overlap_stats.hits.shape[1]:
            random_sel = np.random.choice(
                cluster_overlap_stats.hits.shape[1], sample, replace=False
            )
            cluster_overlap_stats.hits = cluster_overlap_stats.hits.iloc[:, random_sel]

        # print("Calculate test per feature")
        final_test_args = dict(
            simulate_pval=True, replicate=int(1e4), workspace=1_000_000
        )
        if test_arg_dict is not None:
            final_test_args.update(test_arg_dict)
        cluster_overlap_stats.test_per_feature(
            method="hybrid", cores=cores, test_args=final_test_args
        )
        # print("Calculate test per cluster per feature")
        cluster_overlap_stats.test_per_cluster_per_feature()

        with open(cluster_overlap_stats_out_fp, "wb") as fout:
            pickle.dump(cluster_overlap_stats, fout)

        cluster_overlap_stats.hits.to_pickle(
            cluster_overlap_stats_out_fp[:-2] + "_hits.p"
        )
        cluster_overlap_stats.hits.to_csv(
            cluster_overlap_stats_out_fp[:-2] + "_hits.tsv", sep="\t", index=False
        )

        cluster_overlap_stats.cluster_pvalues.to_csv(
            cluster_overlap_stats_out_fp[:-2] + "_element-pvalues.tsv",
            sep="\t",
            index=False,
        )

        cluster_overlap_stats.log_odds_ratio.to_csv(
            cluster_overlap_stats_out_fp[:-2] + "_log-odds.tsv", sep="\t", index=False
        )

    cluster_overlap_stats: rsp.ClusterOverlapStats
    for max_pvalue in max_pvalues:
        # print("Create barcode figure", "max_pvalue", max_pvalue)
        cluster_overlap_stats_filtered = cluster_overlap_stats.filter(
            "cluster_pvalues", max_pvalue
        )
        if not cluster_overlap_stats_filtered.cluster_pvalues.empty:
            with mpl.rc_context(plotting_context):
                fig = rsp.barcode_heatmap(
                    cluster_overlap_stats_filtered,
                    **plot_args,
                )
                out_png = barcode_plot_png_by_maxpvalue.format(max_pvalue=max_pvalue)
                fig.set_dpi(90)
                fig.savefig(out_png)
                if "pdf" in additional_formats:
                    fig.savefig(out_png.replace(".png", ".pdf"))
                if "svg" in additional_formats:
                    fig.savefig(out_png.replace(".png", ".svg"))
            plt.close()
        else:
            print(
                f"WARNING: not features left for pvalue {max_pvalue} threshold. No plot created"
            )


# copy-pasted from mouse_hematopoiesis.wgbs.clustering.clustering1.merged_hierarchy_v2_enrichments_lib
def run_geneset_overlap_stats(
    config: Dict[str, Any],
    max_pvalues: List[float],
    plotting_context: Dict,
    barcode_plot_args: Dict[str, any],
    plots_only=False,
    additional_formats=None,
):
    """
    Args:
        barcode_plot_args: any argument except max_pvalue
        max_pvalues: one plot per max_pvalue in the list

    test function within this function is parallelized with config['cores'] (?)
    # query regions
    # annotation
    # database
      # overlap stats
    # clustering
    # filter / merge
      # cluster overlap stats with test results
      # plot1
      # plot2
    {root_output_dir}/{query_name}/{database}.{clustering}.{filter_or_merge}.{test}/file-type.suffix
    """

    overlap_stats_pattern = (
        config["output_dir"] + "/{anno_name}/{database}/overlap-stats.p"
    )
    cluster_overlap_stats_pattern = (
        config["output_dir"]
        + "/{anno_name}/{database}/{clustering}.{filter}/cluster-overlap-stats.p"
    )
    barcode_plot_png_pattern = (
        config["output_dir"]
        + "/{anno_name}/{database}/{clustering}.{filter}/barcode-plot_maxpvalue-{max_pvalue}.png"
    )

    for (
        (anno_name, gene_annos),
        (database_name, database_fp),
        (clustering_name, (clustering_column_name, clustering)),
        (cleaner_name, cleaner),
    ) in product(
        config["annotations"].items(),
        config["databases"].items(),
        config["clusterings"].items(),
        config["filters"].items(),
    ):
        # print("Anno", anno_name)
        # print("Database", database_name)
        # print("Partitioning:", clustering_name)
        # print("Filter", cleaner_name)

        # A nested folder structure is created. All folders can be created at once
        # by creating the parent directory of the barcode plot png.
        overlap_stats_p = overlap_stats_pattern.format(
            anno_name=anno_name, database=database_name
        )
        cluster_overlap_stats_p = cluster_overlap_stats_pattern.format(
            anno_name=anno_name,
            database=database_name,
            clustering=clustering_name,
            filter=cleaner_name,
        )
        barcode_plot_png_by_maxpvalue = barcode_plot_png_pattern.format(
            anno_name=anno_name,
            database=database_name,
            clustering=clustering_name,
            filter=cleaner_name,
            max_pvalue="{max_pvalue}",
        )
        Path(barcode_plot_png_by_maxpvalue).parent.mkdir(parents=True, exist_ok=True)

        # Prepare cluster ids as series with Granges index
        if isinstance(clustering, str):
            clustering = pd.read_pickle(clustering)
        if isinstance(gene_annos, str):
            gene_annos = pd.read_pickle(gene_annos)
        cluster_ids_ser = clustering[clustering_column_name]
        cluster_ids_ser.index = (
            cluster_ids_ser.reset_index()[["Chromosome", "Start", "End"]]
            .set_index(["Chromosome", "Start", "End"])
            .index
        )

        if not plots_only:
            # print("Calculate OverlapStats")
            overlap_stats = rsp.GenesetOverlapStats(
                annotations=gene_annos, genesets_fp=database_fp
            )
            overlap_stats.compute()
            with open(overlap_stats_p, "wb") as fout:
                pickle.dump(overlap_stats, fout)
            overlap_stats.coverage_df.to_pickle(overlap_stats_p[:-2] + "_hits.p")
            overlap_stats.coverage_df.rename(
                columns={"Chromosome": "#Chromosome"}
            ).to_csv(overlap_stats_p[:-2] + "_hits.bed", sep="\t")

            # print("Calculate ClusterOverlapStats")
            # regions cleaner not yet implemented
            if isinstance(cleaner, pd.MultiIndex):
                cluster_overlap_stats = overlap_stats.aggregate(
                    cluster_ids_ser, index=cleaner
                )
            else:
                cluster_overlap_stats = overlap_stats.aggregate(cluster_ids_ser)
        else:
            with open(cluster_overlap_stats_p, "rb") as fin:
                cluster_overlap_stats = pickle.load(fin)

        process_cluster_overlap_stats(
            cluster_overlap_stats,
            plotting_context=plotting_context,
            max_pvalues=max_pvalues,
            plot_args=barcode_plot_args,
            barcode_plot_png_by_maxpvalue=barcode_plot_png_by_maxpvalue,
            cluster_overlap_stats_out_fp=cluster_overlap_stats_p,
            test_arg_dict=config["test_arg_dict"],
            sample=config["sample"],
            plots_only=plots_only,
            cores=config["cores"],
            additional_formats=additional_formats,
        )


def run_geneset_enrichment_analysis(
    merged_gene_annos: pd.DataFrame,
    feature_annos: pd.DataFrame,
    cluster_ids: pd.DataFrame,
    geneset_databases_d: Dict[str, str],
    output_dir: str,
    report_dir: str,
    max_pvalues: Tuple[float, ...],
    barcode_plot_args_d: Dict[str, Any],
    recompute: bool,
    filters: List[Literal["promoter", "gene_regions", "all_annotated"]],
    additional_formats=Optional[Tuple[str, ...]],
    n_cores=24,
):
    """

    Parameters
    ----------
    merged_gene_annos
    cluster_ids
    geneset_databases_d
        geneset name -> path to GMT file
    max_pvalues
    barcode_plot_args_d
        passed to
    feature_annos
        Chromosome Start End, historical expectation in code: ['Chromosome', 'Start', 'End'] index, need this DF to create that index

    Implementation notes
    --------------------

    *Workflow*
    - currently no separation of plotting and enrichment computation, this needs to be changed

    *Report*
    - the report contains all results on disk matching the expected pattern, ie the report is not filtered
      by the results of the current run
    - to clean up the report, delete the output_dir and rerun the desired enrichments


    *Performance*
    in order to improve performance, you could
    - provide flag arg for switching off per feature test, these tests are often not used
      (in favor of per feature, per cluster tests)
    - provide flag arg for not saving so many individual results to disk
      - while screening many partitionings, it may be overkill to write out
        many intermediate results to disk, as this workflow does


    """

    if recompute:
        # replace integer index with ['Chromosome', 'Start', 'End'] index, expected for historical reasons
        pd.testing.assert_index_equal(merged_gene_annos.index, cluster_ids.index)
        pd.testing.assert_index_equal(merged_gene_annos.index, feature_annos.index)
        merged_gene_annos = merged_gene_annos.copy()
        cluster_ids = cluster_ids.copy()
        grange_index = feature_annos.set_index(mhvars.grange_col_names).index
        merged_gene_annos.index = grange_index
        cluster_ids.index = grange_index

        merged_gene_annos["gene_name"] = merged_gene_annos["gene_name"].str.upper()

        # Compute filters for geneset enrichment, ie to do enrichment only against promoters, gene regions etc.
        promoter_index = merged_gene_annos.query('feat_class == "Promoter"').index
        gene_region_index = merged_gene_annos.query(
            'feat_class not in ["intergenic", "DCRD"]'
        ).index
        all_annotated_index = merged_gene_annos.query(
            'feat_class != "intergenic"'
        ).index

        # ### Create config

        clusterings = {k: (k, cluster_ids) for k in cluster_ids.columns}

        config = {
            "annotations": {"gtfanno": merged_gene_annos["gene_name"]},
            "databases": geneset_databases_d,
            "clusterings": clusterings,
            "filters": pd.Series(
                {
                    "promoter": promoter_index,
                    "gene_regions": gene_region_index,
                    "all_annotated": all_annotated_index,
                }
            )
            .loc[filters]
            .to_dict(),
            "output_dir": output_dir,
            "sample": None,
            "cores": n_cores,
            "test_arg_dict": dict(replicate=int(1e4), workspace=1_000_000),
        }

        run_geneset_overlap_stats(
            config,
            max_pvalues=max_pvalues,
            barcode_plot_args=barcode_plot_args_d,
            plots_only=False,
            plotting_context=mhstyle.paper_context,
            additional_formats=additional_formats,
        )

        # ### Create geneset report

        print("Generate report")
        geneset_enrichment_patterns = {
            "geneset_enrichment": (
                output_dir
                + "/gtfanno/{database}/{clustering}.{filter}/barcode-plot_maxpvalue-{max_pvalue}.png"
            )
        }
        os.makedirs(report_dir, exist_ok=True)

        section_cols = ["database", "clustering", "filter", "max_pvalue"]
        metadata_table = fr.pattern_set_to_metadata_table(
            geneset_enrichment_patterns,
            wildcard_constraints={"clustering": r".+(?=\.)"},
        )
        metadata_table = metadata_table.sort_values(section_cols)
        # metadata_table.sort_values(section_cols, inplace=True)
        fr.copy_report_files_to_report_dir(
            metadata_table,
            root_dir=output_dir,
            report_dir=report_dir,
        )
        report_config = fr.convert_metadata_table_to_report_json(
            metadata_table, section_cols
        )
        report_config.update(
            dict(toc_headings="h1, h2, h3, h4", autocollapse_depth="3")
        )
        fr.Report({"marker_gene_enrichments": report_config}).generate(report_dir)

    report_path = report_dir + "/marker_gene_enrichments.html"
    print("Report: ", report_path)
    mhpaths.link_fn(report_path)

    # for plot fine-tuning
    # cluster_overlap_stats_p = '/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/cohort_results/analyses/hierarchy/enrichments/full-hierarchy/clustering-v2/genesets/gtfanno/lab mouse-hematopoiesis marker genes/hierarchy-v2_29-clusters-merged.Promoter/cluster-overlap-stats.p'
    # cluster_overlap_stats = pd.read_pickle(cluster_overlap_stats_p)
    #
    # with mpl.rc_context(paper_context):
    #     fig = rsp.barcode_heatmap(cluster_overlap_stats, clusters_as_rows=True,
    #                               col_width_cm=0.5, row_height_cm=0.1, max_pvalue=1e-2,
    #                               linewidth=0.5)
    #     fig.savefig('/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/cohort_results/analyses/hierarchy/enrichments/full-hierarchy/clustering-v2/genesets/gtfanno/lab mouse-hematopoiesis marker genes/hierarchy-v2_29-clusters-merged.Promoter/cluster-overlap-stats.temp.png')

    #
    # # ### Important result paths
    #
    #
    # dm(
    #     f"""
    # Results are under: {os.path.dirname(geneset_enrichment_patterns["geneset_enrichment"])}
    #
    # fields:
    # - database: one of {config['databases'].keys()}
    # - filter: eg only {config['filters'].keys()}
    #
    # files:
    # - cluster-overlap-stats.p: full ClusterOverlapStats object
    # - log-odds.tsv: log odds ratios for every cluster/feature combination
    # - element-pvalues.tsv: p-values for every cluster/feature log odds ratio, from 2x2 fisher tests
    # - barcode-plot_maxpvalue-{{max_pvalue}}.png
    #     - max_pvalue: displayed items must have pvalue below this threshold
    # """
    # )

    """alternative from d3a

# ### Geneset enrichments figure report

# %%capture
geneset_enrichment_patterns = {
    "gene set enrichment": (
        geneset_enrichment_output_dir
        + "/gtfanno/{database}/{clustering}.{filter}/barcode-plot_maxpvalue-{max_pvalue}.png"
    )
}
section_cols = ["database", "clustering", "filter", "max_pvalue"]
metadata_table = fr.pattern_set_to_metadata_table(
    geneset_enrichment_patterns, wildcard_constraints={"clustering": r".+(?=\.)"}
)
metadata_table = metadata_table.sort_values(section_cols)
fr.copy_report_files_to_report_dir(
    metadata_table, root_dir=enrichment_output_dir, report_dir=enrichment_report_dir
)
report_config = fr.convert_metadata_table_to_report_json(metadata_table, section_cols)
report_config.update(dict(toc_headings="h1, h2, h3, h4", autocollapse_depth="3"))
fr.Report({"msigdb_gene_enrichments": report_config}).generate(enrichment_report_dir)

# #### Link to report (on http server)

dkfz_link(enrichment_report_dir + "/msigdb_gene_enrichments.html")
    """


def create_gmt_for_rosenbauer_genesets(rosenbauer_genesets_df):
    gmts_d = dict(
        rosenbauer_genesets_all_gmt=(
            "/home/kraemers/projects/dendritic_cells/local/genesets/rosenbauer_all.gmt"
        ),
        rosenbauer_genesets_coeff12_only_gmt=(
            "/home/kraemers/projects/dendritic_cells/local/genesets/rosenbauer_coeff12.gmt"
        ),
    )
    gmt_df = rosenbauer_genesets_df.T.apply(lambda ser: ser.str.upper())
    gmt_df.insert(0, "description", "NA")

    with open(gmts_d["rosenbauer_genesets_all_gmt"], "wt") as fout:
        for i, row in gmt_df.reset_index().iterrows():
            fout.write(row.dropna().str.cat(sep="\t") + "\n")

    with open(gmts_d["rosenbauer_genesets_coeff12_only_gmt"], "wt") as fout:
        for i, row in (
            gmt_df.loc[["coeff1_up", "coeff1_down", "coeff2_up", "coeff2_down"]]
            .reset_index()
            .iterrows()
        ):
            fout.write(row.dropna().str.cat(sep="\t") + "\n")

    return gmts_d


def cluster_marker_genes_overview(
    genesets_d_ser: Dict[str, pd.Series],
    primary_gene_annos: pd.DataFrame,
    granges_df: pd.DataFrame,
    cluster_ids: pd.DataFrame,
):
    """Annotate clusters with geneset genes

    Parameters
    ----------
    genesets_d_ser
        {'geneset_database_name': pd.Series({'geneset1': ['geneA', 'geneB'], 'geneset2': ...})
    primary_gene_annos
        as returned by gtfanno
        this version of the function assumes that the dataframe does not yet contain a 'region_id' column automatically and merges this in manually
    granges_df
        region_id // Chromosome Start End
        used for addion region_id to primary_gene_annos; index must match cluster_ids
    cluster_ids
        region_id // partitioning_name1 ...

    Returns
    -------
    pd.DataFrame
        int_index // geneset_database_name 	genomic_regions 	geneset_name 	cluster_id 	gene_name
    """

    pd.testing.assert_index_equal(granges_df.index, cluster_ids.index)
    primary_gene_annos = primary_gene_annos.copy().reset_index(drop=True)
    primary_gene_annos_with_region_id_column = pd.merge(
        primary_gene_annos,
        granges_df[["Chromosome", "Start", "End"]].reset_index(),
        on=["Chromosome", "Start", "End"],
        how="left",
    )
    pd.testing.assert_frame_equal(
        primary_gene_annos[["Chromosome", "Start", "End"]],
        primary_gene_annos_with_region_id_column[["Chromosome", "Start", "End"]],
    )
    assert primary_gene_annos_with_region_id_column["region_id"].notnull().all()
    primary_gene_annos = primary_gene_annos_with_region_id_column.set_index("region_id")
    primary_gene_annos["gene_name"] = primary_gene_annos["gene_name"].str.upper()

    assert "DCRD" in primary_gene_annos["feat_class"].values
    genomic_region_filters = dict(
        promoter=primary_gene_annos["feat_class"] == "Promoter",
        gene_region=primary_gene_annos["feat_class"].isin(["intergenic", "DCRD"]),
        all_annotated=primary_gene_annos["feat_class"] != "intergenic",
    )

    res_d = {}
    for geneset_name, geneset_ser in genesets_d_ser.items():
        for cell_type_name, curr_genes_l in geneset_ser.items():
            for (
                genomic_region_name,
                genomic_region_bool_idx,
            ) in genomic_region_filters.items():
                is_in_geneset_bool_idx = primary_gene_annos["gene_name"].isin(
                    pd.Series(curr_genes_l).str.upper().to_numpy()
                )
                gene_anno_ser_filtered = primary_gene_annos.loc[
                    genomic_region_bool_idx & is_in_geneset_bool_idx,
                    "gene_name",
                ]
                for part_name, cluster_id_ser in cluster_ids.iteritems():
                    res_d[
                        (geneset_name, genomic_region_name, part_name, cell_type_name)
                    ] = (
                        gene_anno_ser_filtered.groupby(cluster_id_ser)
                        .agg(lambda ser: list(ser.sort_values()))
                        .reset_index()
                        .set_axis(["cluster_id", "gene_name"], axis=1)
                        .set_index("cluster_id")
                    )

    res_df = pd.concat(
        res_d,
        axis=0,
        names=[
            "geneset_database_name",
            "genomic_regions",
            "partitioning_name",
            "geneset_name",
        ],
    ).reset_index()

    # clean up unused categories, not necessary if the full looping above is actually used
    res_df = res_df.drop('partitioning_name', axis=1)

    return res_df
