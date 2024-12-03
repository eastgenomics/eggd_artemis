"""Main app entrypoint"""

import concurrent
import os
import dxpy
import datetime
import logging
from collections import defaultdict

from glob import glob
import pip

pip.main(["install", "--no-index", "--no-deps", *glob("packages/*")])


from mergedeep import merge


from utils.io import (
    make_cnv_session,
    read_excluded_regions_to_df,
    write_output_file,
)
from utils.query import (
    get_cnv_call_details,
    get_cnv_file_ids,
    get_multiqc_report,
    make_url,
    find_snv_files,
)
from utils.util_functions import (
    get_excluded_intervals,
    remove_unnecessary_outputs,
)


def generate_sample_outputs(
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


@dxpy.entry_point("main")
def main(
    url_duration,
    lock_cells=True,
    snv_path=None,
    cnv_path=None,
    bed_file=None,
    qc_status=None,
    multiqc_report=None,
):
    # Set up logging
    logger = logging.getLogger(__name__)
    logger.addHandler(dxpy.DXLogHandler())
    logger.propagate = False
    logger.setLevel(logging.DEBUG)

    # Get the project ID
    DX_PROJECT = os.environ.get("DX_PROJECT_CONTEXT_ID")
    print(DX_PROJECT)

    # Set the environment context to allow upload
    dxpy.set_workspace_id(DX_PROJECT)

    # Get output folder set for this job
    job_output = dxpy.bindings.dxjob.DXJob(
        os.environ.get("DX_JOB_ID")
    ).describe()["folder"]

    # Make output folder for job
    dxpy.api.project_new_folder(
        DX_PROJECT, input_params={"folder": job_output, "parents": True}
    )

    # Get name of project for output naming
    project_name = dxpy.describe(DX_PROJECT)["name"]
    project_name = "_".join(project_name.split("_")[1:-1])

    # Gather required SNV files if SNV path is provided
    if snv_path:
        logger.info("Gathering Small variant files")

        snv_data = {}

        for path in snv_path.split(","):
            print(f"Gathering reports from: {path}")
            snv_reports = list(
                dxpy.bindings.search.find_data_objects(
                    name="*xlsx",
                    name_mode="glob",
                    folder=path,
                    project=DX_PROJECT,
                    describe=True,
                )
            )

            # Get SNV ids
            snv_files = find_snv_files(snv_reports)
            merge(snv_data, snv_files)

        print(f"Size of snv data dict: {len(snv_data.keys())}")

        # If multiqc report not given, search for multiqc report
        if not multiqc_report:
            print(
                "No MultiQC report given, searching instead "
                f"with path {snv_path}"
            )
            multiqc_report = get_multiqc_report(
                snv_path.split(",")[0], DX_PROJECT
            )

    # Gather required CNV files if CNV path is provided
    if cnv_path:
        logger.info("Gathering CNV files")

        # Create folder for session files
        dxpy.api.project_new_folder(
            DX_PROJECT,
            input_params={
                "folder": f"{job_output}/igv_sessions",
                "parents": True,
            },
        )

        cnv_data = {}

        for path in cnv_path.split(","):
            print(f"Gathering reports from: {path}")
            cnv_reports = list(
                dxpy.bindings.search.find_data_objects(
                    name="*xlsx",
                    name_mode="glob",
                    folder=path,
                    project=DX_PROJECT,
                    describe=True,
                )
            )

            gcnv_job_info = get_cnv_call_details(cnv_reports)
            cnv_files = get_cnv_file_ids(cnv_reports, gcnv_job_info)

            # Get excluded intervals file
            excluded_intervals = get_excluded_intervals(gcnv_job_info)
            ex_intervals_url = make_url(
                excluded_intervals, DX_PROJECT, url_duration
            )

            merge(cnv_data, cnv_files)

        print(f"Size of cnv data dict: {len(cnv_data.keys())}")

        # If multiqc report not given or found earlier, search for it
        if not multiqc_report:
            print(
                "No MultiQC report given, searching instead "
                f"with path {cnv_path}"
            )
            multiqc_report = get_multiqc_report(
                cnv_path.split(",")[0], DX_PROJECT
            )
    else:
        ex_intervals_url = ""

    logger.info("Making URLs for additional files")

    # If a bed file is provided, add to a link to the output
    if bed_file:
        bed_file_url = make_url(
            bed_file, "project-Fkb6Gkj433GVVvj73J7x8KbV", url_duration
        )
    else:
        # Setting as empty to avoid session errors
        bed_file_url = ""

    # If a QC Status xlsx is provided, add to a link to the output
    if qc_status:
        qc_status_url = make_url(qc_status, DX_PROJECT, url_duration)
    else:
        qc_status_url = "No QC status file provided"

    file_data = {}

    if snv_path and cnv_path:
        merge(file_data, snv_data, cnv_data)
    elif snv_path:
        file_data = snv_data
    elif cnv_path:
        file_data = cnv_data
    else:
        logger.debug("No paths given, exiting...")
        exit(1)

    today = datetime.datetime.now().strftime("%Y-%m-%d")
    expiry_date = (
        datetime.datetime.strptime(
            datetime.datetime.now().strftime("%Y-%m-%d"), "%Y-%m-%d"
        )
        + datetime.timedelta(seconds=url_duration)
    ).strftime("%Y-%m-%d")

    # generate all outputs for each sample
    logger.info("Generating per sample outputs")
    all_sample_outputs = {}

    with concurrent.futures.ThreadPoolExecutor(max_workers=32) as executor:
        # submit jobs mapping each id to describe call
        concurrent_jobs = {
            executor.submit(
                generate_sample_outputs,
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
                print(
                    f"Error getting data for {concurrent_jobs[future]}: {exc}"
                )

    multiqc_url = make_url(multiqc_report, DX_PROJECT, url_duration)

    # Remove download URLs for reports with no variants in and remove
    # excluded regions dataframe if no excluded regions
    all_sample_outputs = remove_unnecessary_outputs(
        all_sample_outputs, snv_reports=True, cnv_reports=True
    )

    write_output_file(
        all_sample_outputs,
        today,
        expiry_date,
        multiqc_url,
        qc_status_url,
        project_name,
        lock_cells,
    )

    # Upload output to the platform
    output = {}
    url_file = dxpy.upload_local_file(
        f"{project_name}_{today}.xlsx", folder=job_output, tags=[expiry_date]
    )
    output["url_file"] = dxpy.dxlink(url_file)

    # Find session files to link to output
    """
    Example all_sample_outputs format we're looping over:
    'sample-name': {
        'sample': 'sample-name',
        'bam_url': 'url_for_bam',
        'bai_url':  'url_for_bai',
        'clinical_indications': {
            'R141.1': {
                'SNV': [{
                    'SNV count': '7',
                    'coverage_url: '=HYPERLINK("hyperlink", "hyperlink"),
                    'coverage_summary' : 'coverage summary text',
                    'snv_url': '=HYPERLINK("hyperlink", "hyperlink"),
                }],
                'CNV': [{
                    'CNV count': '1',
                    'cnv_bed': 'bed_url',
                    'cnv_seg': 'seg_url',
                    'cnv_session_fileid': 'file-XYZ',
                    'cnv_url': '=HYPERLINK("hyperlink", "hyperlink"),
                    'cnv_session_url': '=HYPERLINK("hyperlink", "hyperlink"),
                    'cnv_excluded_regions_df': pd.DataFrame of excluded regions file
                }]
            }
        }
    }
    """
    session_files = []
    for sample, sample_info in all_sample_outputs.items():
        for clin_ind in sample_info.get("clinical_indications"):
            if sample_info["clinical_indications"][clin_ind].get("CNV"):
                for cnv_item in sample_info["clinical_indications"][clin_ind][
                    "CNV"
                ]:
                    session_file_id = cnv_item["cnv_session_fileid"]
                    session_files.append(session_file_id)

    print(f"Found {len(session_files)} session files to link to job output")
    if session_files:
        output["session_files"] = [dxpy.dxlink(item) for item in session_files]

    return output


dxpy.run()
