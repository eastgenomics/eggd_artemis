"""Functions for calling all the querying and generating the output data"""

from collections import defaultdict
import concurrent
import os
from typing import Tuple

import dxpy
from mergedeep import merge

from .io import make_cnv_session, read_excluded_regions_to_df
from .query import (
    get_cnv_file_ids,
    get_cnv_call_details,
    find_snv_files,
    make_url,
)
from .util_functions import get_excluded_intervals


def gather_snv_data(snv_paths, dx_project) -> dict:
    """
    Finds all xlsx reports and the associated job input files using
    query.find_snv_files from the given SNV path(s)
    Verify that the number of processed SNV reports matches the number of
    `.xlsx` report files discovered in the path.

    Args:
        snv_paths (str): comma separated string of SNV path(s) to search
        dx_project (str) : DNAnexus project ID to search

    Returns:
        snv_data (dict): dict of all sample SNV file data
    Raises:
        RuntimeError: if the number of SNV reports found does not match
            the expected number based on unique sample names.
    """
    snv_data = {}

    for path in snv_paths.split(","):
        print(f"Gathering reports from: {path}")
        snv_reports = list(
            dxpy.bindings.search.find_data_objects(
                name="*xlsx",
                name_mode="glob",
                folder=path,
                project=dx_project,
                describe=True,
            )
        )

        no_reports_expected = len(snv_reports)

        snv_files = find_snv_files(snv_reports)
        merge(snv_data, snv_files)
        # Count total reports across all samples and clinical indications
        total_reports_processed = 0
        for sample, sample_data in snv_files.items():
            for clin_ind in sample_data.get("clinical_indications", {}):
                snv_reports_list = sample_data["clinical_indications"][clin_ind].get("SNV", [])
                if snv_reports_list:
                    total_reports_processed += len(snv_reports_list)

        print(f"Found {no_reports_expected} SNV reports")

        if total_reports_processed != no_reports_expected:
            raise RuntimeError(
                f"Found {no_reports_expected} SNV report files, "
                f"but processed {total_reports_processed} reports in the data structure."
            )
        print(f"Successfully processed {total_reports_processed} SNV reports from {path}")

    return snv_data


def gather_cnv_data(cnv_paths, dx_project, url_duration) -> Tuple[dict, str]:
    """
    Finds all xlsx reports and the associated job input files from the
    given CNV path(s).
    Check if the number of reports found matches the expected number of unique
    sample names based on the CNV report names.


    Args:
        cnv_paths (str): comma separated string of CNV path(s) to search
        dx_project (str): DNAnexus project ID to search
        url_duration (int): maximum duration of generated download URL

    Returns:
        cnv_data (dict): dict of all sample CNV file data
        ex_intervals_url (str): download URL for the excluded intervals file
    Raises:
        RuntimeError: if the number of CNV reports found does not match
            the expected number based on unique sample names.
    """
    cnv_data = {}

    for path in cnv_paths.split(","):
        print(f"Gathering reports from: {path}")
        cnv_reports = list(
            dxpy.bindings.search.find_data_objects(
                name="*xlsx",
                name_mode="glob",
                folder=path,
                project=dx_project,
                describe=True,
            )
        )

        no_reports_expected = len(cnv_reports)

        gcnv_job_info = get_cnv_call_details(cnv_reports)
        cnv_files = get_cnv_file_ids(cnv_reports, gcnv_job_info)

        # Count total reports across all samples and clinical indications
        total_reports_processed = 0
        for sample, sample_data in cnv_files.items():
            for clin_ind in sample_data.get("clinical_indications", {}):
                cnv_reports_list = sample_data["clinical_indications"][clin_ind].get("CNV", [])
                if cnv_reports_list:
                    total_reports_processed += len(cnv_reports_list)

        print(f"Found {no_reports_expected} CNV reports")

        if total_reports_processed != no_reports_expected:
            raise RuntimeError(
                f"Found {no_reports_expected} CNV report files, "
                f"but processed {total_reports_processed} reports in the data structure."
            )
        print(f"Successfully processed {total_reports_processed} CNV reports from {path}")
        # Get excluded intervals file
        excluded_intervals = get_excluded_intervals(gcnv_job_info)
        ex_intervals_url = make_url(
            excluded_intervals, dx_project, url_duration
        )

        merge(cnv_data, cnv_files)

    return cnv_data, ex_intervals_url


def generate_single_sample_output(
    sample,
    sample_data,
    url_duration,
    ex_intervals_url,
    bed_url,
    expiry_date,
    job_output,
    build,
    select_tracks,
) -> dict:
    """
    Generates all URLs and session file for a given sample

    Args:
        sample (str): sample name
        sample_data (dict): file IDs of sample related files
        url_duration (int): URL duration in seconds
        ex_intervals_url (str): URL of excluded intervals file
        bed_url (str): URL of bed file
        expiry_date (str): date of URL expiration
        job_output (str): output folder set for job
        build (int): genome build to add reference files to the session file for
        select_tracks (str): comma separated string of IGV reference tracks to select


    Returns:
    outputs (dict): dict of sample outputs (i.e. URLs/dataframes)
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
                    build,
                    select_tracks,
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
    build,
    select_tracks,
) -> dict:
    """
    Wrapper to call generate_single_sample_output in parallel, returning
    a chongus dict of all sample data to generate the final output xlsx

    Args:
        sample (str): sample name
        sample_data (dict): file IDs of sample related files
        url_duration (int): URL duration in seconds
        ex_intervals_url (str): URL of excluded intervals file
        bed_url (str): URL of bed file
        expiry_date (str): date of URL expiration
        job_output (str): output folder set for job
        build (int): genome build to add reference files to the session file for
        select_tracks (str): comma separated string of IGV reference tracks to select

    Returns:
        outputs (dict): dict of sample outputs (i.e. URLs/dataframes)
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
                build,
                select_tracks,
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
