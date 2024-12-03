"""Functions for calling all the querying and generating the output data"""

from collections import defaultdict
import concurrent
import os

import dxpy

from .io import make_cnv_session, read_excluded_regions_to_df
from .query import make_url


def generate_single_sample_output(
    sample,
    sample_data,
    url_duration,
    ex_intervals_url,
    bed_url,
    expiry_date,
    job_output,
) -> dict:
    """
    Generates all URLs and session file for a given sample

    Parameters
    ----------
    sample : str
        sample name
    sample_data : dict
        file IDs of sample related files
    url_duration : int
        URL duration in seconds
    ex_intervals_url : str
        URL of excluded intervals file
    bed_url : str
        URL of bed file
    expiry_date : str
        date of URL expiration
    job_output : str
        output folder set for job

    Returns
    -------
    outputs : dict
        dict of sample outputs (i.e. URLs/dataframes)
    """
    dx_project = os.environ.get("DX_PROJECT_CONTEXT_ID")

    outputs = defaultdict(lambda: defaultdict(dict))
    outputs["sample"] = sample

    bam_url = make_url(sample_data["Alignment BAM"], dx_project, url_duration)
    bai_url = make_url(sample_data["Alignment BAI"], dx_project, url_duration)

    outputs["bam_url"] = bam_url
    outputs["bai_url"] = bai_url

    # Loop over clinical indications
    for clin_ind in sample_data["clinical_indications"]:
        snv_files = sample_data["clinical_indications"][clin_ind].get("SNV")
        # If we have SNV reports, for each report make URLs and append to
        # SNV value
        if snv_files:
            outputs["clinical_indications"][clin_ind]["SNV"] = []

            for snv_file in snv_files:
                coverage_url = make_url(
                    snv_file["Coverage report"], dx_project, url_duration
                )
                snv_url = make_url(
                    snv_file["SNV variant report"], dx_project, url_duration
                )

                # Read in summary text file if file exists, else leave
                # as None
                if snv_file["Summary text file"]:
                    snv_file["Summary text file"] = (
                        dxpy.open_dxfile(
                            snv_file["Summary text file"],
                            project=dx_project,
                            mode="r",
                        )
                        .read()
                        .replace("Clinical report summary:\n", "")
                    )

                outputs["clinical_indications"][clin_ind]["SNV"].append(
                    {
                        "SNV count": snv_file["SNV count"],
                        "coverage_url": (
                            f'=HYPERLINK("{coverage_url}", "{coverage_url}")'
                        ),
                        "coverage_summary": snv_file["Summary text file"],
                        "snv_url": f'=HYPERLINK("{snv_url}", "{snv_url}")',
                    }
                )

        # If we have CNV reports, for each report make URLs/session files
        # and append to CNV value
        cnv_files = sample_data["clinical_indications"][clin_ind].get("CNV")
        if cnv_files:
            cnv_bed = make_url(
                sample_data["CNV visualisation"], dx_project, url_duration
            )
            outputs["clinical_indications"][clin_ind]["CNV"] = []

            for cnv_file in cnv_files:
                session_file_end = cnv_file["Session file name"]
                cnv_url = make_url(
                    cnv_file["CNV variant report"], dx_project, url_duration
                )
                cnv_seg = make_url(
                    cnv_file["CNV calls for IGV"], dx_project, url_duration
                )
                cnv_session = make_cnv_session(
                    sample,
                    session_file_end,
                    bam_url,
                    bai_url,
                    cnv_bed,
                    cnv_seg,
                    ex_intervals_url,
                    bed_url,
                    job_output,
                    expiry_date,
                )
                cnv_session_url = make_url(
                    cnv_session, dx_project, url_duration
                )
                cnv_excluded_regions_df = read_excluded_regions_to_df(
                    file_id=cnv_file["CNV excluded regions"],
                    project=dx_project,
                )

                outputs["clinical_indications"][clin_ind]["CNV"].append(
                    {
                        "CNV count": cnv_file["CNV count"],
                        "cnv_bed": cnv_bed,
                        "cnv_seg": cnv_seg,
                        "cnv_session_fileid": cnv_session,
                        "cnv_url": f'=HYPERLINK("{cnv_url}", "{cnv_url}")',
                        "cnv_session_url": (
                            f'=HYPERLINK("{cnv_session_url}",'
                            f' "{cnv_session_url}")'
                        ),
                        "cnv_excluded_regions_df": cnv_excluded_regions_df,
                    }
                )

    return outputs


def generate_all_sample_outputs(
    file_data,
    url_duration,
    ex_intervals_url,
    bed_file_url,
    expiry_date,
    job_output,
) -> dict:
    """
    Wrapper to call generate_single_sample_output in parallel, returning
    a chongus dict of all sample data to generate the final output xlsx
    """
    all_sample_outputs = {}

    with concurrent.futures.ThreadPoolExecutor(max_workers=32) as executor:
        # submit jobs mapping each id to describe call
        concurrent_jobs = {
            executor.submit(
                generate_single_sample_output,
                sample,
                sample_data,
                url_duration,
                ex_intervals_url,
                bed_file_url,
                expiry_date,
                job_output,
            ): sample
            for sample, sample_data in file_data.items()
        }

        for future in concurrent.futures.as_completed(concurrent_jobs):
            # access returned output as each is returned in any order
            try:
                data = future.result()
                all_sample_outputs[data["sample"]] = data
            except Exception as exc:
                # catch any errors that might get raised during querying
                raise RuntimeError(
                    f"Error getting data for {concurrent_jobs[future]}: {exc}"
                ) from exc

    return all_sample_outputs
