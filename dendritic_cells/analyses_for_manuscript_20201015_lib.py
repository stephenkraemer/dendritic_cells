"""
gtfanno version: b81595d
"""

import os
import pickle
import re
import tempfile
from itertools import product
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import methlevels as ml
import codaplot as co
import figure_report as fr
import gtfanno as ga
import matplotlib as mpl
import matplotlib.pyplot as plt
import mouse_hema_meth.paths as mhpaths
import mouse_hema_meth.shared_vars as mhvars
import mouse_hema_meth.styling as mhstyle
import mouse_hema_meth.utils as ut
import numpy as np
# pyright workaround
import pandas
import pandas.testing
import pandas.api
pd = pandas
import region_set_profiler as rsp
import seaborn as sns
from matplotlib.figure import Figure
from typing_extensions import Literal


import seaborn as sns
import scipy.stats
from itertools import combinations, product

# isort: off
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import StrVector, IntVector

from functools import partial
import rpy2.ipython.html

rpy2.ipython.html.init_printing()
rpy2.ipython.html.html_rdataframe = partial(
    rpy2.ipython.html.html_rdataframe, table_class="docutils"
)


from rpy2.robjects.conversion import localconverter
from rpy2.robjects import pandas2ri, numpy2ri

import mouse_hema_meth.rpy2_utils as rpy2_utils

numpy2ri.activate()
pandas2ri.activate()

import rpy2.robjects.lib.ggplot2 as gg
from rpy2.ipython.ggplot import image_png


# def image_png2(p, figsize):
#     display(image_png(p, figsize[0] * 72, figsize[1] * 72))


import mouse_hema_meth.rpy2_styling as mh_rpy2_styling

# lemon = importr("lemon")

NULL = ro.NULL
import rpy2.rinterface as ri

NA = ri.NA_Logical

# isort: on

mpl.use("Agg")


PROJECT_TEMPDIR = "/icgc/dkfzlsdf/analysis/hs_ontogeny/temp"

print('reloaded analyses lib')


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
        'geneset_database_name' -> pd.Series({'geneset1': ['geneA', 'geneB'], 'geneset2': ...})
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
            int_index // geneset_database_name genomic_regions geneset_name cluster_id gene_name
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
    res_df = res_df.drop("partitioning_name", axis=1)

    return res_df


def create_clustermaps_with_special_order(
    df: pd.DataFrame,
    cluster_ids_df: pd.DataFrame,
    png_path: str,
    cmap: str,
    guide_title: str,
    global_row_order_df: Optional[pd.DataFrame] = None,
    figsize=(10, 10),
    n_per_cluster=None,
) -> None:
    """Create special-order clustermaps for all given partitionings

    Uses clustermap_special order, see there for details


    Parameters
    ----------
    df
    cluster_ids_df
        one or more partitionings (one column per partitioning)
    global_row_order_df
        as returned by hclust_within_clusters
        index: same as df (consecutive integer index starting at 0)
        columns: same as cluster_ids_df (partitioning names)
        values: position of feature in globally ordered data (according to partitioning + hclust)
    png_path
    """

    pd.testing.assert_index_equal(df.index, cluster_ids_df.index, check_names=False)
    if global_row_order_df is not None:
        pd.testing.assert_index_equal(
            df.index, global_row_order_df.index, check_names=False
        )

    for part_name, cluster_id_ser in cluster_ids_df.iteritems():
        if global_row_order_df is None:
            global_row_order_ser = None
        else:
            global_row_order_ser = global_row_order_df[part_name]

        fig = clustermap_special_order(
            df,
            cluster_id_ser,
            global_row_order_ser,
            figsize=figsize,
            n_per_cluster=n_per_cluster,
            cmap=cmap,
            guide_title=guide_title,
        )
        fig.savefig(png_path)
        fig.savefig(png_path.replace(".png", ".pdf"))


def clustermap_special_order(
    df,
    cluster_id_ser,
    global_row_order_ser,
    n_per_cluster: int = 200,
    figsize=None,
    rasterized=True,
    guide_title="guide title",
    cmap="RdBu_r",
):
    """Clustermap respecting special row_order argument

    Row order is taken from global_row_order_ser if defined (which could for example come
    from hierarchical clustering, within a primary partitioning)
    Otherwise, row order is created based on the cluster_ids_df partitionings, and row
    order within each cluster is random.

    Parameters
    ----------
    df
    cluster_id_ser
    global_row_order_ser
        as returned by hclust_within_clusters
        index: same as df (consecutive integer index starting at 0)
        values: position of feature in globally ordered data (according to partitioning + hclust)
    n_per_cluster
        number of features per cluster shown in heatmap

    Returns
    -------


    """

    pd.testing.assert_index_equal(df.index, cluster_id_ser.index, check_names=False)
    if global_row_order_ser is not None:
        pd.testing.assert_index_equal(
            df.index, global_row_order_ser.index, check_names=False
        )

    # sample max n_per_cluster features per cluster
    sampled_cluster_ids = (
        cluster_id_ser.groupby(cluster_id_ser, group_keys=False)
        .apply(
            lambda ser: ser.sample(n_per_cluster)
            if n_per_cluster < ser.shape[0]
            else ser
        )
        .sort_index()
    )  # sort index just for faster retrieval during indexing operations
    if global_row_order_ser is None:
        row_order = sampled_cluster_ids.reset_index(drop=True).sort_values().index
    else:
        # get the row_order for the heatmap from the global row order
        row_order = (
            # get global order position for selected elements
            global_row_order_ser.loc[sampled_cluster_ids.index]
            # reset index to eventually get integer indices as required by row_order
            .reset_index(drop=True)
            # sort on order positions
            .sort_values()
            # get original integer index values in the now correctly ordered data
            # this can be used as integer row order argument for eg co.cross_plot
            .index
        )
    fig = nice_clustermap(
        stat_df=df.loc[sampled_cluster_ids.index],
        cluster_id_ser=sampled_cluster_ids,
        row_order=row_order,
        figsize=figsize,
        rasterized=rasterized,
        cmap=cmap,
        guide_title=guide_title,
    )
    return fig


def nice_clustermap(
    stat_df: pd.DataFrame,
    cluster_id_ser: pd.Series,
    row_order: Union[pd.Series, np.ndarray],
    colors: Dict = None,
    figsize=(20 / 2.54, 20 / 2.54),
    rasterized=True,
    cmap="magma",
    guide_title="% sign cpgs",
) -> Figure:
    """Nice clustermap for given cluster_ids and integer index row_order

    Parameters
    ----------
    stat_df
    cluster_id_ser
    row_order
    colors
    figsize

    Returns
    -------
    Figure

    """

    # create dict cluster_id -> color if necessary
    if colors is None:
        n_clusters = cluster_id_ser.nunique()
        colors = dict(
            zip(np.unique(cluster_id_ser), sns.color_palette("Set1", n_clusters))
        )

    array_to_figure_res, plot_arr = co.cross_plot(
        center_plots=[
            dict(
                df=stat_df,
                cmap=cmap,
                guide_title=guide_title,
                yticklabels=False,
                xticklabels=True,
                rasterized=rasterized,
                guide_args=dict(shrink=0.4, aspect=4),
            )
        ],
        # center_margin_ticklabels=True,
        pads_around_center=[(0.2 / 2.54, "abs")],
        figsize=figsize,
        # constrained_layout=True,
        # layout_pads=dict(wspace=0, hspace=0, h_pad=0, w_pad=0),
        left_plots=[
            dict(
                _func=co.frame_groups,
                direction="y",
                colors=colors,
                linewidth=1,
                add_labels=True,
                labels=None,
                label_colors=None,
                label_groups_kwargs=dict(rotation=0),
            )
        ],
        left_col_sizes=[(0.5 / 2.54, "abs")],
        row_order=row_order,
        # col_linkage=False,
        row_spacing_group_ids=cluster_id_ser,
        # # col_spacing_group_ids=col_clusters,
        row_spacer_sizes=0.01,
        # # col_spacer_sizes=0.05,
        # default_plotting_func_kwargs=dict(guide_args=dict(shrink=0.3, aspect=8)),
    )
    return array_to_figure_res["fig"]


def collect_original_flagstats(file_pattern):

    flagstats_paths_df = ut.get_files_df(file_pattern)
    res = {}
    for pop_rep_idx, path in flagstats_paths_df.set_index(["pop", "rep"])[
        "path"
    ].iteritems():
        with open(path) as fin:
            flagstats_data_ser = pd.Series(
                fin.readlines(),
                index=[
                    "No. of reads",
                    "No. of duplicates",
                    "No. mapped",
                    "No. paired",
                    "No. read1",
                    "No. read2",
                    "No. of proper pairs",
                    "No. both mapped",
                    "No. of Singletons",
                    "No. diffchrom",
                    "No. diffchrom (mapQ >=5)",
                ],
            )
        flagstats_counts = (
            flagstats_data_ser.str.extract(r"^(\d+) \+ (\d+) ")
            .astype("i4")
            .set_axis(["QC-passed reads", "QC-failed reads"], axis=1)
            .T
        )

        flagstats_perc = flagstats_counts.divide(
            flagstats_counts["No. of reads"], axis=0
        ).rename(columns=lambda s: re.sub(r"No. (of )?", "% ", s))

        flagstats_df_both = pd.concat([flagstats_counts, flagstats_perc], axis=1)
        try:
            res[
                mhstyle.nice_pop_names_d[pop_rep_idx[0]], pop_rep_idx[1]
            ] = flagstats_df_both
        except KeyError:
            # not a hema meth pop with a nice name, skip
            continue

    res_df = (
        pd.concat(res)
        .rename_axis(["Population", "Replicate", "Read QC status"], axis=0)
        .rename_axis("Stat", axis=1)
        .stack()
        .unstack("Read QC status")
        .sort_index()
    )

    return res_df


def gather_conversion_rate_table(mcalls_metadata_table):
    """Gather CG and CH methylation stats for chromosome 1 and MT chromosome

    methylctools based

    Parameters
    ----------
    mcalls_metadata_table
        should be prefiltered to desired pops and reps, index: (pop, rep)
    """

    metrics_file_by_sampleid_chrom = "/icgc/dkfzlsdf/project/mouse_hematopoiesis/sequencing/whole_genome_bisulfite_tagmentation_sequencing/view-by-pid/{sample_id}/blood/paired/merged-alignment/methylation/merged/methylationCallingMetrics/blood_{sample_id}_merged.mdup.bam.{chrom}.metrics.csv"

    res = {}
    for pop_rep_idx, row_ser in mcalls_metadata_table.iterrows():
        for curr_chrom in ["1", "MT"]:
            with open(
                metrics_file_by_sampleid_chrom.format(
                    chrom=curr_chrom, sample_id=row_ser["sample_id"]
                )
            ) as f:
                lines = f.readlines()
                res[
                    mhstyle.nice_pop_names_d[pop_rep_idx[0]], pop_rep_idx[1], curr_chrom
                ] = (pd.Series(lines[2:4]).str.strip().str.split("\t", expand=True))

    res_df = (
        (
            pd.concat(res)
            .reset_index()
            .set_axis(
                [
                    "Population",
                    "Replicate",
                    "Chromosome",
                    "Unused",
                    "Chromosome_b",
                    "Motif",
                    "N methylated",
                    "N unmethylated",
                    "Beta value",
                ],
                axis=1,
            )
        )
        .drop("Unused", axis=1)
        .astype({"N methylated": "i4", "N unmethylated": "i4", "Beta value": "f8"})
    )
    res_df["% Methylation"] = (
        res_df.eval("`N methylated` / (`N methylated` + `N unmethylated`)") * 100
    )
    assert (
        res_df["% Methylation"].round(2) == (res_df["Beta value"] * 100).round(2)
    ).all()

    res_df["Est. conversion rate"] = 100 - res_df["% Methylation"]

    pd.testing.assert_series_equal(
        res_df["Chromosome"], res_df["Chromosome_b"], check_names=False
    )

    return (
        res_df.drop("Chromosome_b", axis=1)
        .sort_values(["Motif", "Chromosome", "Population", "Replicate"])
        .query('Motif != "CG" and Chromosome == "1"')
    )


def gather_meth_level_data(
    mcalls_metadata_table: pd.DataFrame,
    recompute,
    pop_quantiles_wide_p,
    pop_quantiles_long_p,
    sampled_rep_betavalue_p,
    sampled_pop_betavalue_p,
    pop_betavalue_p,
    pop_betavalue_mean_df_p,
    rep_beta_value_mean_p,
):

    # Get total number of CpGs for preallocating some results
    n_cpgs = pd.read_pickle(
        mcalls_metadata_table.query('subject == "hsc"').iloc[0]["pickle_path"]
    ).shape[0]

    # this currently assumes that the meth calls file only contains autosomes,
    # and that the qc analysis is also only performed on autosomes
    # let's put a reminder in case the calls change in the future
    assert np.all(
        np.sort(
            pd.read_pickle(
                mcalls_metadata_table.query('subject == "hsc"').iloc[0]["pickle_path"]
            )["#chrom"]
            .unique()
            .astype(str)
        )
        == np.sort(np.arange(1, 20).astype(str))
    )

    if recompute:

        # coverage quantiles to be queried
        # noinspection PyUnresolvedReferences
        queried_quantiles = np.sort(
            np.linspace(0, 1, 11).tolist() + [0.25, 0.75, 0.025, 0.975]
        )

        # random integer index to sample CpGs from each replicate
        random_int_idx = np.sort(
            np.random.RandomState(1).choice(n_cpgs, size=100_000, replace=False)
        )

        # ~ 7 min
        # compute coverage per population by aggregating replicates
        # then compute some stats
        pop_quantile_ser_d = {}
        sampled_rep_meth_d = {}
        sampled_pop_meth_d = {}
        pop_beta_value_d = {}
        rep_mean_d = {}
        for pop, group_df in mcalls_metadata_table.groupby("subject"):
            print(pop)
            # preallocate pop coverage series for genome-wide CpGs
            nmeth_pop = pd.Series(0, index=pd.RangeIndex(n_cpgs), dtype="i4")
            ntotal_pop = pd.Series(0, index=pd.RangeIndex(n_cpgs), dtype="i4")
            # Iterate over replicates, compute total pop coverage
            for _unused, row_ser in group_df.iterrows():
                if row_ser["subject"] == "monos":
                    # for monos-new, no pickle
                    curr_calls = pd.read_csv(row_ser["bed_path"], sep="\t").rename(
                        columns={"#chrom": "Chromosome", "start": "Start", "end": "End"}
                    )
                else:
                    curr_calls = pd.read_pickle(row_ser.pickle_path)
                nmeth_rep = curr_calls["n_meth"]
                ntotal_rep = curr_calls["n_total"]
                beta_rep = curr_calls["beta_value"]
                rep_mean_d[row_ser.subject, row_ser.rep] = beta_rep.agg(
                    ["mean", "std", "var", "median"]
                )
                nmeth_pop += nmeth_rep
                ntotal_pop += ntotal_rep
                sampled_rep_meth_d[(row_ser.subject, row_ser.rep)] = beta_rep.iloc[
                    random_int_idx
                ].reset_index(drop=True)

            beta_values_pop = nmeth_pop / ntotal_pop
            # Coverage quantiles for pop
            pop_quantile_ser_d[pop] = beta_values_pop.quantile(queried_quantiles)
            sampled_pop_meth_d[pop] = beta_values_pop.iloc[random_int_idx].reset_index(
                drop=True
            )
            # full coverage series
            pop_beta_value_d[pop] = beta_values_pop

        pop_quantiles_df = pd.DataFrame(pop_quantile_ser_d)
        pop_quantiles_df.to_pickle(pop_quantiles_wide_p)
        pop_quantiles_df.to_csv(ut.tsv(pop_quantiles_wide_p), sep="\t")

        pop_quantiles_df_long = (
            pop_quantiles_df.rename(columns=mhstyle.nice_pop_names_d)
            .rename_axis("Population", inplace=False, axis=1)
            .loc[[0.1, 0.25, 0.5, 0.75, 0.9]]
            .stack()
            .unstack(0)
            .reset_index()
        )
        pop_quantiles_df_long.columns = pop_quantiles_df_long.columns.astype(str)
        pop_quantiles_df_long.to_pickle(pop_quantiles_long_p)
        pop_quantiles_df_long.to_csv(ut.tsv(pop_quantiles_long_p), sep="\t")

        sampled_rep_meth_df = pd.DataFrame(sampled_rep_meth_d)
        sampled_rep_meth_df.columns.names = ["Population", "Replicate"]
        sampled_rep_meth_df.to_pickle(sampled_rep_betavalue_p)
        sampled_rep_meth_df.to_csv(ut.tsv(sampled_rep_betavalue_p), sep="\t")

        sampled_pop_meth_df = pd.DataFrame(sampled_pop_meth_d)
        sampled_pop_meth_df.to_pickle(sampled_pop_betavalue_p)
        sampled_pop_meth_df.to_csv(ut.tsv(sampled_pop_betavalue_p), sep="\t")

        beta_value_pop_df = pd.DataFrame(pop_beta_value_d)
        beta_value_pop_df.to_pickle(pop_betavalue_p)

        pop_beta_value_mean_df = beta_value_pop_df.agg(["mean", "std"])
        pop_beta_value_mean_df.to_pickle(pop_betavalue_mean_df_p)

        rep_beta_value_mean_df = pd.concat(rep_mean_d)
        rep_beta_value_mean_df.to_pickle(rep_beta_value_mean_p)
    else:
        pop_quantiles_df = pd.read_pickle(pop_quantiles_wide_p)
        pop_quantiles_df_long = pd.read_pickle(pop_quantiles_long_p)
        sampled_rep_meth_df = pd.read_pickle(sampled_rep_betavalue_p)
        sampled_pop_meth_df = pd.read_pickle(sampled_pop_betavalue_p)
        beta_value_pop_df = pd.read_pickle(pop_betavalue_p)
        pop_beta_value_mean_df = pd.read_pickle(pop_betavalue_mean_df_p)
        rep_beta_value_mean_df = pd.read_pickle(rep_beta_value_mean_p)

    return (
        pop_quantiles_df,
        pop_quantiles_df_long,
        sampled_rep_meth_df,
        sampled_pop_meth_df,
        beta_value_pop_df,
        pop_beta_value_mean_df,
        rep_beta_value_mean_df,
    )




def get_rep_meth_levels(mcalls_metadata_table, dmrs, n_cores) -> ml.MethStats:

    pop_order = [
        "hsc",
        "mdp",
        "cdp",
        "cmop",
        "monos",
        "pdc",
        "dc-cd11b",
        "dc-cd8a",
    ]

    bed_calls_ds1 = ml.BedCalls(
        metadata_table=mcalls_metadata_table.reset_index(),
        tmpdir=mhpaths.project_temp_dir,
        #     pop_order=mcalls_metadata_table['sample_id'],
        pop_order=pop_order,
        beta_value_col=6,
        n_meth_col=7,
        n_total_col=8,
    )

    meth_stats_rep = bed_calls_ds1.aggregate(
        intervals_df=dmrs,
        n_cores=n_cores,
        subjects=None,
        additional_index_cols=None,
    )



def get_meth_stats_arrays(meth_stats_rep):

    beta_values = meth_stats_rep.counts.loc[:, pd.IndexSlice[:, :, "beta_value"]]
    beta_values = beta_values.droplevel(-1, axis=1).reset_index(drop=True)

    q1_median_q3_rep_wide = (
        beta_values.describe(percentiles=[0.01, 0.25, 0.5, 0.75, 0.99])
        .T.reset_index()
        .rename(
            columns={
                "Subject": "Population",
                "25%": "Q1",
                "50%": "median",
                "75%": "Q3",
                "1%": "min1",
                "99%": "max99",
            }
        )
    )
    q1_median_q3_rep_wide

    q1_median_q3_rep_long = (
        beta_values.describe()
        .T.reset_index()
        .rename(
            columns={"Subject": "Population", "25%": "Q1", "50%": "median", "75%": "Q3", 'mean': 'mean'}
        )[["Population", "Replicate", "Q1", "median", "Q3", 'mean']]
        .set_index(["Population", "Replicate"])
        .stack()
        .to_frame()
        .reset_index()
        .set_axis(["Population", "Replicate", "stat", "value"], axis=1)
    )
    q1_median_q3_rep_long

    return q1_median_q3_rep_wide, q1_median_q3_rep_long

def rest():
    df = q1_median_q3_rep_wide
    pops = ["pdc", "dc-cd11b", "dc-cd8a"]

    stats_l = []
    for stat, (popa, popb) in product(["Q1", "median", "Q3"], product(pops, pops)):
        print(stat, popa, popb)

        popa = "hsc"
        popb = "pdc"
        stat = "median"

        mw_u, pvalue = scipy.stats.mannwhitneyu(
            [0.8, 0.81, 0.79],
            [0.4, 0.39, 0.41],
            # df.query("Population == @popa")[stat].to_numpy(),
            # df.query("Population == @popb")[stat].to_numpy(),
            use_continuity=True,
            alternative="two-sided",
        )
        pvalue

        stats_l.append([stat, popa, popb, mw_u, pvalue])
    stats_df = pd.DataFrame(stats_l).set_axis(
        ["stat", "popA", "popB", "U", "pvalue"], axis=1
    )

    kruskal_format_means = pd.pivot(
        q1_median_q3_rep_wide.query("Population in @pops"),
        index="Population",
        columns="Replicate",
        values="mean",
    )

    import scikit_posthocs

    stat, p_value = scipy.stats.kruskal(
        *[kruskal_format_means.loc[pop].to_numpy() for pop in pops],
    )

    dunn_res_df = scikit_posthocs.posthoc_dunn(
        kruskal_format_means.to_numpy(),
        p_adjust='fdr_bh',
        sort=True,
    )

    stat, pvalue = scipy.stats.f_oneway(
        *[kruskal_format_means.loc[pop].to_numpy() for pop in pops],
    )


    import statsmodels


    df = kruskal_format_means.stack().reset_index()

    kruskal_format_means


    res = statsmodels.stats.multicomp.pairwise_tukeyhsd(
        df[0],
        df['Population'].to_numpy(),
        alpha=0.05)

    res.pvalues
    res.summary()


    # wilcox.test(c(0.8, 0.79, 0.81), c(0.4, 0.39, 0.41), paired=F, exact=F)

    plot_pops = ["pdc", "dc-cd8a", "dc-cd11b"]

    results_dir = "/icgc/dkfzlsdf/analysis/hs_ontogeny/notebook-data/gNs4xcMJscaLLwlt"
    point_plot_quartiles_png = results_dir + "/point-plot-quartiles.png"

    q1_median_q3_rep_wide

    ggplot_data = (
        q1_median_q3_rep_long.query("Population in @plot_pops")
        .sort_values(
            "value",
            ascending=False,
        )
        .groupby(["Population", "stat"])
        .apply(lambda df: df.assign(group_order=np.arange(1, df.shape[0] + 1)))
    )

    g = (
        gg.ggplot(ggplot_data)
        + gg.aes_string(x="Population", y="value", group="group_order", color="stat")
        + gg.geom_point(position=gg.position_dodge(width=0.5), size=1)
        + mh_rpy2_styling.gg_paper_theme
        + gg.labs(y='Methylation (%)', x='')
    )
    a = 3

    rpy2_utils.image_png2(g, (ut.cm(6), ut.cm(6)))

    ut.save_and_display(
        g,
        png_path=point_plot_quartiles_png,
        # additional_formats=tuple(),
        height=ut.cm(6),
        width=ut.cm(6),
    )

    q1_median_q3_rep_wide

    g = (
        gg.ggplot(
            q1_median_q3_rep_wide.query("Population in @plot_pops").assign(
                sample=lambda df: df["Population"].astype(str)
                + df["Replicate"].astype(str)
            )
        )
        + gg.geom_boxplot(
            gg.aes_string(
                x="Population",
                fill="Population",
                group="sample",
                lower="Q1",
                upper="Q3",
                middle="median",
                ymin="min1",
                ymax="max99",
                # position=gg.position_dodge(width=0.5),
            ),
            stat="identity",
        )
        # + mh_rpy2_styling.gg_paper_theme
        + gg.theme(axis_text_x=gg.element_text(angle=90, hjust=1))
        + gg.scale_fill_brewer(guide=False)
    )
    a = 3
    ut.save_and_display(
        g,
        png_path=point_plot_quartiles_png,
        additional_formats=tuple(),
        height=ut.cm(6),
        width=ut.cm(7),
    )
    # image_png2(g, (ut.cm(12), ut.cm(12)))

    beta_values.loc[:, ("hsc", "1")]


def _format_stats_df(df, index_name):

    new_column_names = (
        "Population" if df.columns.nlevels == 1 else ["Population", "Replicate"]
    )
    return (
        df.reset_index(drop=True)
        .pipe(_stringify_columns)
        .rename_axis(index=index_name, columns=new_column_names)
        .rename(columns=mhstyle.nice_pop_names_d, level=0)
    )


def _stringify_columns(df):
    if df.columns.nlevels == 1:
        return df.set_axis(df.columns.astype(str), axis=1, inplace=False)
    else:
        idx = pd.MultiIndex.from_arrays(
            [df.columns.get_level_values(0).astype(str), df.columns.get_level_values(1)]
        )
        df.columns = idx
        return df
