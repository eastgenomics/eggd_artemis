"""General utility functions"""

import os
from typing import Tuple, Union

import dxpy


def add_session_file_ids_to_job_output(all_sample_outputs, job_output) -> dict:
    """
    Search through the all samples output dict for any session file URLs
    generated and add them to the job output to link to the job output spec

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

    Parameters
    ----------
    all_sample_outputs : dict
        dict of all sample output data
    job_output : dict
        output mapping dict of the job output

    Returns
    -------
    dict
        job output dict with added session file IDs
    """
    session_files = []

    for _, sample_info in all_sample_outputs.items():
        for clin_ind in sample_info.get("clinical_indications"):
            if sample_info["clinical_indications"][clin_ind].get("CNV"):
                for cnv_item in sample_info["clinical_indications"][clin_ind][
                    "CNV"
                ]:
                    session_file_id = cnv_item["cnv_session_fileid"]
                    session_files.append(session_file_id)

    print(f"Found {len(session_files)} session files to link to job output")

    if session_files:
        job_output["session_files"] = [
            dxpy.dxlink(item) for item in session_files
        ]

    return job_output


def get_excluded_intervals(gcnv_output_dict) -> Union[str, None]:
    """Get the excluded regions file from the gcnv dictionary

    Args:
        gcnv_output_dict (dict): dictionary of gcnv i/o files

    Returns:
        excluded_file (str | None): file id of the excluded regions file
            or None if excluded_intervals.bed not found in gCNV output
    """
    for k, v in gcnv_output_dict.items():
        if k.endswith("excluded_intervals.bed"):
            return v
    return None


def initialise_project() -> Tuple[str, str, str]:
    """
    Set required project data, get the project name and destination for
    downstream naming

    Returns:
        project_name (str): name of DNAnexus project
        project_id (str): ID of DNAnexus project
        job_output (str): destination folder set for the job
    """
    project_id = os.environ.get("DX_PROJECT_CONTEXT_ID")

    # Set the environment context to allow upload
    dxpy.set_workspace_id(project_id)

    # Get output folder set for this job
    job_output = dxpy.bindings.dxjob.DXJob(
        os.environ.get("DX_JOB_ID")
    ).describe()["folder"]

    # Make output folder for job
    dxpy.api.project_new_folder(
        project_id, input_params={"folder": job_output, "parents": True}
    )

    # Get name of project for output naming
    project_name = dxpy.describe(project_id)["name"]
    project_name = "_".join(project_name.split("_")[1:-1])

    return project_name, project_id, job_output


def filter_reference_tracks(select_tracks, reference_tracks) -> list:
    """
    Filter the reference tracks from the defaults.py URL tracks list by
    those provided to the input select tracks

    Args:
        select_tracks (str): comma separated string of IGV reference
            tracks to select
        reference_tracks (list): list of reference track URL dicts

    Returns:
        reference_tracks (list): filtered list of reference tracks

    Raises:
        ValueError: raised if invalid track names provided to select_tracks
    """
    select_tracks = [x.strip().lower() for x in select_tracks.split(",")]
    available_reference_tracks = [x["name"].lower() for x in reference_tracks]

    invalid_tracks = [
        x for x in select_tracks if x not in available_reference_tracks
    ]
    if invalid_tracks:
        raise ValueError(
            "Invalid track names provided to select from build URL"
            f" tracks: {invalid_tracks}"
        )

    return [x for x in reference_tracks if x["name"].lower() in select_tracks]


def check_session_track_order(session_tracks) -> list:
    """
    Ensures that each track for the IGV session has the `order` key set,
    this ensures the order of tracks displaying is consistent.

    For any that are not defined, we will set it to the index within the
    track list or igv.js will shuffle them around.
    """
    orders_set = [x.get("order") for x in session_tracks]

    if all(orders_set):
        return session_tracks

    print("setting order")

    max_set_order = max([x for x in orders_set if x])
    ordered_tracks = []

    for track in session_tracks:
        if not track.get("order"):
            print("unordered track")
            max_set_order += 1

            print(track)
            print(max_set_order)
            track["order"] = max_set_order

        ordered_tracks.append(track)

    return ordered_tracks


def set_order_map(snv_only=False) -> dict:
    """Set the order of the session depending on input

    Args:
        snv_only (bool, optional): If True, returns snv order map

    Returns:
        order_map (dict): order of inputs to the session template
    """

    if snv_only == True:

        order_map = {
            "6": {"url": "", "indexURL": ""},
            "7": {"url": "", "name": "vcf"},
        }

    else:

        order_map = {
            "6": {"url": "", "name": "CNV-bed"},
            "7": {"url": "", "indexURL": ""},
            "8": {"url": "", "name": "variant-bed"},
            "9": {"url": "", "name": "excluded-regions"},
            "10": {"url": "", "name": "CNV-targets-bed"},
        }

    return order_map


def remove_unnecessary_outputs(
    all_sample_outputs, snv_reports, cnv_reports
) -> dict:
    """
    Remove links to download Excel reports if there are no variants and remove
    excluded regions dataframe if there are no excluded regions

    Parameters
    ----------
    all_sample_outputs : dict
        generated URLs, dataframes and file metadata for each sample
    snv_reports : boolean
        if True, SNV report links will be removed if no SNVs pass filtering
    cnv_reports : boolean
        if True, CNV report links will be removed if no CNVs detected

    Returns
    -------
    all_sample_outputs: dict
        generated URLs, dataframes and file metadata for each sample with Excel
          download URLs removed when no variants are present and excluded
          region dataframe removed when there are no excluded regions.
    """
    for sample, sample_data in all_sample_outputs.items():
        for clin_ind in sample_data["clinical_indications"]:
            file_data = sample_data["clinical_indications"][clin_ind]

            # If variant count was never present then replace with None
            # so we can easily not include this field in final .xlsx
            # If snv_reports option True, replace SNV report URL with text
            # if variant count is zero.
            if "SNV" in file_data:
                for file in file_data["SNV"]:
                    if file.get("SNV count") == "Unknown":
                        file["SNV count"] = None
                    if snv_reports:
                        if file.get("SNV count") == "0":
                            file["snv_url"] = "No SNVs post-filtering"

            # Do same for CNV reports, remove CNV report link if
            # cnv_reports option True
            if "CNV" in file_data:
                for file in file_data["CNV"]:
                    if file.get("CNV count") == "Unknown":
                        file["CNV count"] = None
                    if cnv_reports:
                        if file.get("CNV count") == "0":
                            file["cnv_url"] = "No CNVs detected"
                    # Column names/header are read into the first row of df,
                    # therefore if df length is < 2 it contains no data
                    if len(file.get("cnv_excluded_regions_df")) < 2:
                        file["cnv_excluded_regions_df"] = (
                            "No CNV excluded regions"
                        )
    return all_sample_outputs
