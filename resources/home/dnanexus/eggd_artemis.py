"""Main app entrypoint for the app"""

from glob import glob
import dxpy
import datetime
import logging
import pip

pip.main(["install", "--no-index", "--no-deps", *glob("packages/*")])

from mergedeep import merge

from utils.generate import (
    gather_cnv_data,
    gather_snv_data,
    generate_all_sample_outputs,
)
from utils.io import write_output_file
from utils.query import get_multiqc_report, make_url
from utils.util_functions import (
    add_session_file_ids_to_job_output,
    initialise_project,
    remove_unnecessary_outputs,
)


@dxpy.entry_point("main")
def main(
    url_duration,
    build,
    select_tracks=None,
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

    project_name, project_id, job_output_folder = initialise_project()

    snv_data = {}
    cnv_data = {}

    # Gather required SNV files if SNV path is provided
    if snv_path:
        logger.info("Gathering Small variant files")
        snv_data = gather_snv_data(snv_paths=snv_path, dx_project=project_id)

        # If multiqc report not given, search for multiqc report
        if not multiqc_report:
            print(
                "No MultiQC report given, searching instead "
                f"with path {snv_path}"
            )
            multiqc_report = get_multiqc_report(
                snv_path.split(",")[0], project_id
            )

    # Gather required CNV files if CNV path is provided
    if cnv_path:
        logger.info("Gathering CNV files")

        # Create folder for session files
        dxpy.api.project_new_folder(
            project_id,
            input_params={
                "folder": f"{job_output_folder}/igv_sessions",
                "parents": True,
            },
        )

        cnv_data, ex_intervals_url = gather_cnv_data(
            cnv_paths=cnv_path,
            dx_project=project_id,
            url_duration=url_duration,
        )

        # If multiqc report not given or found earlier, search for it
        if not multiqc_report:
            print(
                "No MultiQC report given, searching instead "
                f"with path {cnv_path}"
            )
            multiqc_report = get_multiqc_report(
                cnv_path.split(",")[0], project_id
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
        qc_status_url = make_url(qc_status, project_id, url_duration)
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
        job_output=job_output_folder,
        build=build,
        select_tracks=select_tracks,
    )

    multiqc_url = make_url(multiqc_report, project_id, url_duration)

    # Remove download URLs for reports with no variants in and remove
    # excluded regions dataframe if no excluded regions
    all_sample_outputs = remove_unnecessary_outputs(
        all_sample_outputs, snv_reports=True, cnv_reports=True
    )

    output_xlsx_file = write_output_file(
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
        filename=output_xlsx_file, folder=job_output_folder, tags=[expiry_date]
    )
    output["url_file"] = dxpy.dxlink(url_file)

    output = add_session_file_ids_to_job_output(
        all_sample_outputs=all_sample_outputs, job_output=output
    )

    return output


dxpy.run()
