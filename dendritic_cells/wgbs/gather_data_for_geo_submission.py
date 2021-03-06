# executed with mouse_hema_meth_py37_mamba_full

import jobsub
import subprocess
import pandas as pd
import mouse_hema_meth.utils as ut
import os
from pathlib import Path
import dendritic_cells.wgbs.gather_data_for_geo_submission_lib as lib


def main():

    PROJECT_TEMPDIR = "/icgc/dkfzlsdf/analysis/hs_ontogeny/temp3"
    os.makedirs(PROJECT_TEMPDIR, exist_ok=True)

    results_dir = "/icgc/dkfzlsdf/analysis/hs_ontogeny/notebook-data/NlIosmWpSIbw7MA8"
    os.makedirs(results_dir, exist_ok=True)

    # get metadata table for FASTQ files
    # ==================================

    fastq_pattern = (
        "/icgc/dkfzlsdf/project/mouse_hematopoiesis/sequencing/whole_genome_bisulfite_tagmentation_sequencing"
        "/view-by-pid/{sample}/blood/paired"
        "/{run}/sequence/{filename}.fastq.gz"
    )
    samples = [
        "mdp_1",
        "mdp_2",
        "monos-new_1",
        "monos-new_2",
        "monos-new_3",
        "cmop_1",
        "cmop_2",
        "cdp_1",
        "cdp_2",
        "cdp_3",
        "cdp_4",
        "pdc_1",
        "pdc_2",
        "pdc_3",
        "dc-cd11b_1",
        "dc-cd11b_2",
        "dc-cd11b_3",
        "dc-cd8a_1",
        "dc-cd8a_2",
        "dc-cd8a_3",
    ]
    fastq_metadata_table = ut.get_files_df(fastq_pattern)

    # add full_filename: {run_id}_{fastq_name}.fastq.gz
    fastq_metadata_table["full_filename"] = (
        fastq_metadata_table[["run", "filename"]].apply(
            lambda ser: ser.str.cat(sep="_"), axis=1
        )
        + ".fastq.gz"
    )

    # Restrict to and order by samples
    fastq_metadata_table = fastq_metadata_table.set_index("sample").loc[samples]

    # Create TSV 1
    # ============

    # sample1\tfastq11\tfastq12\t...
    # sample2\tfastq21\tfastq22\t...
    # This TSV can be copy-pasted into the samples section of the geo metadata sheet
    # it is correct that there are a lot (24) fastq files:
    # 4 libraries * 3 lanes * 2 reads = 24 fastqs per replicate
    full_filenames_ser = (
        fastq_metadata_table.groupby("sample")["full_filename"]
        .agg(lambda ser: ser.str.cat(sep="\t"))
        .loc[samples]
    )
    full_filenames_ser = (
        full_filenames_ser.index.to_series() + "\t" + full_filenames_ser
    )
    Path(results_dir + "/hBnx1AMysAAVQWsi.tsv").write_text(
        full_filenames_ser.str.cat(sep="\n")
    )
    ut.dkfz_link(results_dir + "/hBnx1AMysAAVQWsi.tsv")

    # Get FASTQ md5sums
    # =================



    lsf_setup_commands = """\
            source ~/bash_profile
            conda activate mouse_hema_meth_py37_mamba_full
    """

    jobs = []
    for sample, row_ser in fastq_metadata_table.iterrows():
        print(sample)
        job = jobsub.LsfJob(
            func=lib.get_md5sum2,
            lsf_resources=jobsub.LsfResources(
                average_memory_gb=1,
                max_memory_gb=1,
                walltime="01:00",
                logdir="/home/kraemers/temp/logs",
                name=f"md5sum_{row_ser['full_filename']}",
                cores=1,
            ),
            kwargs=dict(fastq_fp=row_ser["path"], sample=row_ser["full_filename"]),
        )
        job.submit(
            tmpdir=PROJECT_TEMPDIR, setup_commands=lsf_setup_commands, dryrun=False
        )
        jobs.append(job)
    jobsub.wait_for_jobs(jobs)

    results = [job.get_result(delete_temp_file=False) for job in jobs]
    full_filenames_fastq_md5_tu, md5res_str = zip(*results)
    md5sum_lines_split = [line.strip().split() for line in md5res_str]
    md5_res_df = pd.DataFrame(md5sum_lines_split, columns=['md5sum', 'abspath'])
    md5_res_df['full_filename'] = full_filenames_fastq_md5_tu
    assert md5_res_df.notnull().all().all()
    md5_res_df[['full_filename', 'md5sum']].to_csv(
        results_dir + '/abcdeaPLoB3fQorpulFFH.tsv', sep='\t', index=False
    )
    ut.dkfz_link(results_dir + '/abcdeaPLoB3fQorpulFFH.tsv')

    # initally had to use local parallalization, because the cluster had a problem
    # this took ~ 3h
    # results = Parallel(n_jobs=24)(
    #     delayed(lib.get_md5sum2)(row_ser["path"], row_ser["full_filename"])
    #     for sample, row_ser in fastq_metadata_table.iterrows()
    # )
    # NOTE: there is an error in this code - it does not add the correct full_filename to the output table, which needs to be in the metadata sheet
    # also, in the end i just wrote out the df to csv, this can now be simplified to DataFrame.to_csv
    # full_filenames, md5sum_lines = zip(*results)
    # md5sum_lines_split = [line.strip().split()[::-1] for line in md5sum_lines]
    # # this dataframe needs a full_filename column
    # raw_data_df = pd.DataFrame(md5sum_lines_split, columns=["filename", "md5sum"])
    # Path(results_dir + "/4fL3D7henDSKzAtJ.tsv").write_text(
    #     raw_data_df.apply(lambda ser: ser.str.cat(sep="\t"), axis=1).str.cat(sep="\n")
    # )
    # ut.dkfz_link(results_dir + "/4fL3D7henDSKzAtJ.tsv")


    # Get md5sums for mcalls
    # ======================

    # Get mcalls metadata table
    mcalls_file_pattern = (
        "/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/results_per_pid"
        "/v1_bistro-0.2.0_odcf-alignment/{sample}/meth/meth_calls"
        "/mcalls_{sample}_CG_chrom-merged_strands-merged.bed.gz"
    )
    mcalls_metadata_table = ut.get_files_df(mcalls_file_pattern)
    mcalls_metadata_table = mcalls_metadata_table.set_index("sample").loc[samples]
    mcalls_metadata_table["filename"] = mcalls_metadata_table["path"].str.extract(
        r".*/(.*.bed.gz)$"
    )

    # Get md5sums - quite fast, ~ 1min
    md5sum_results = []
    for sample, row_ser in mcalls_metadata_table.iterrows():
        print(sample)
        md5sum_results.append(
            (
                sample,
                subprocess.run(
                    ["md5sum", row_ser["path"]],
                    check=True,
                    capture_output=True,
                    encoding="utf-8",
                ).stdout,
            )
        )

    # Create tsv: filename, md5sum
    mcalls_samples, mcalls_paths, mcalls_md5sums = zip(
        *[(sample, *res.strip().split()[::-1]) for sample, res in md5sum_results]
    )
    mcalls_metadata_table.loc[mcalls_samples, "md5sum"] = mcalls_md5sums
    Path(results_dir + "/8LbKdGx7hhpQvjyw.tsv").write_text(
        (
            mcalls_metadata_table[["filename", "md5sum"]]
            .apply(lambda ser: ser.str.cat(sep="\t"), axis=1)
            .str.cat(sep="\n")
        )
    )
    ut.dkfz_link(results_dir + "/8LbKdGx7hhpQvjyw.tsv")

    # PE experiment info
    # ==================
    fastq_metadata_table["read"] = fastq_metadata_table["path"].str.extract(
        r".*_R(\d).fastq.gz$"
    )
    assert fastq_metadata_table['read'].notnull().all()
    fastq_metadata_table["stem"] = fastq_metadata_table["path"].str.extract(
        r"(.*)_R\d.fastq.gz$"
    )
    assert fastq_metadata_table['stem'].notnull().all()

    def agg_pe_filenames(ser):
        assert ser.shape[0] == 2
        return ser.str.cat(sep="\t")

    pe_info_df = (
        fastq_metadata_table.sort_values("read")
        .groupby("stem")["full_filename"]
        .agg(agg_pe_filenames)
    )

    Path(results_dir + "/XL72U2utUGgpYnLU.tsv").write_text(pe_info_df.str.cat(sep="\n"))
    ut.dkfz_link(results_dir + "/XL72U2utUGgpYnLU.tsv")

    # Collect all data in one folder for upload to GEO
    # ================================================

    # geo_upload_folder = results_dir + "/OV3eavXF6PkzOpYz/dc-methylomes"
    geo_upload_folder = results_dir + "/kgBLAl9d1HSGDSCI/dc-methylomes"
    os.makedirs(geo_upload_folder, exist_ok=True)

    # FASTQs
    for sample, row_ser in fastq_metadata_table.loc[[
        'monos-new_1',
        'monos-new_2',
        'monos-new_3',
    ]].iterrows():
        link_name = geo_upload_folder + "/" + row_ser["full_filename"]
        try:
            os.remove(link_name)
        except FileNotFoundError:
            pass
        os.symlink(row_ser["path"], link_name)
    assert fastq_metadata_table["full_filename"].duplicated().sum() == 0

    # Mcalls
    for idx, row_ser in mcalls_metadata_table.iterrows():
        curr_filename = row_ser['filename'].replace('monos-new', 'monos')
        link_name = geo_upload_folder + "/" + curr_filename
        try:
            os.remove(link_name)
        except FileNotFoundError:
            pass
        os.symlink(row_ser["path"], link_name)


