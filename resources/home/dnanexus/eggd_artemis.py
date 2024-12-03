"""Main app entrypoint"""

import os
import dxpy
import datetime
import logging

from glob import glob
import pip

pip.main(["install", "--no-index", "--no-deps", *glob("packages/*")])

from mergedeep import merge

from utils.generate import generate_all_sample_outputs
from utils.io import write_output_file
from utils.query import (
    get_cnv_call_details,
    get_cnv_file_ids,
    get_multiqc_report,
    make_url,
    find_snv_files,
)
from utils.util_functions import (
    add_session_file_ids_to_job_output,
    get_excluded_intervals,
    remove_unnecessary_outputs,
)


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

    logger.info("Generating per sample outputs")
    all_sample_outputs = generate_all_sample_outputs(
        file_data=file_data,
        url_duration=url_duration,
        ex_intervals_url=ex_intervals_url,
        bed_file_url=bed_file_url,
        expiry_date=expiry_date,
        job_output=job_output,
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

    output = add_session_file_ids_to_job_output(
        all_sample_outputs=all_sample_outputs, job_output=output
    )

    return output


dxpy.run()
