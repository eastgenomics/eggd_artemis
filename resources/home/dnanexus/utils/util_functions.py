"""General utility functions"""


def get_excluded_intervals(gcnv_output_dict):
    """Get the excluded regions file from the gcnv dictionary

    Args:
        gcnv_output_dict (dict): dictionary of gcnv i/o files

    Returns:
        excluded_file (str): file id of the excluded regions file
    """
    for k, v in gcnv_output_dict.items():
        if k.endswith("excluded_intervals.bed"):
            return v
    return None


def set_order_map(snv_only=False):
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


def remove_unnecessary_outputs(all_sample_outputs, snv_reports, cnv_reports):
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
