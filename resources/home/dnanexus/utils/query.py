"""dxpy querying functions"""

from collections import defaultdict
import concurrent
from typing import Union

import dxpy


def get_multiqc_report(path_to_reports, project) -> Union[str, None]:
    """
    Find appropriate multiqc file

    Args:
        snv_path (string): path to snv reports

    Returns:
        multiqc_file (str): file id of the multiqc file or None if not found
    """
    # Get path to single from reports path
    single = f"/output/{path_to_reports.split('/')[2]}"

    # Find MultiQC jobs in the project
    multiqc_reports = list(
        dxpy.bindings.search.find_jobs(
            name_mode="glob",
            name="*MultiQC*",
            state="done",
            project=project,
            describe=True,
        )
    )

    # In case there's more than one MultiQC job
    # Find the one that has the exists in the single path in question
    multiqc = None
    for report in multiqc_reports:
        if report["describe"]["folder"].startswith(single):
            multiqc = report["describe"]["output"]["multiqc_html_report"]

    return multiqc


def make_url(file_id, project, url_duration) -> str:
    """
    Given a file id create a download url

    Args:
        file_id (string/dict): dxpy file id or dnanexus file link
        project (string): id of project
        url_duration (int, optional): URL duration in seconds.

    Returns:
        file_url (string): Download url of requested file
    """
    # Bind dxpy file object to pass to make_download_url command
    file_info = dxpy.bindings.dxfile.DXFile(dxid=file_id, project=project)

    # Extract the file name to allow it to be used in the url
    file_name = file_info.describe()["name"]

    # Duration currently defaults to 28 days unless provided as input
    file_url = file_info.get_download_url(
        duration=url_duration,
        preauthenticated=True,
        filename=file_name,
        project=project,
    )[0]

    # Without unsetting and clearing the workstation environment
    # the output urls have the local hostname which doesn't work
    # Find hostname
    workstation_hostname = file_url.split("/")[2]
    dx_hostname = "dl.ec1.dnanex.us"

    # Replace with eu specific hostname and update protocol
    file_url = file_url.replace(workstation_hostname, dx_hostname)
    file_url = file_url.replace("http", "https")

    return file_url


def get_cnv_call_details(reports) -> dict:
    """
    Gets the job id of the cnv call job and stores input and output
    files in a dictionary for later use

    Args:
        reports (list): list of cnv report objects from dxpy search

    Returns:
        calling_files (dict): dictionary storing all input and
            output information
    """
    # Get the job id of the generate workbook job
    gen_xlsx_job = reports[0]["describe"]["createdBy"]["job"]

    # Find the reports workflow analysis id
    reports_analysis = dxpy.bindings.dxjob.DXJob(dxid=gen_xlsx_job).describe()[
        "parentAnalysis"
    ]

    # Find the input vcf id
    try:
        vcf_id = dxpy.bindings.dxanalysis.DXAnalysis(
            dxid=reports_analysis
        ).describe()["input"]["stage-cnv_vep.vcf"]
    except KeyError:
        vcf_id = dxpy.bindings.dxanalysis.DXAnalysis(
            dxid=reports_analysis
        ).describe()["input"]["stage-GFYvJF04qq8VKgq34j30pZZ3.vcf"]

    # Find the cnv call job id
    cnv_call_job = dxpy.describe(vcf_id)["createdBy"]["job"]

    # Store all file ids and names in a dictionary
    calling_files = {}

    # Find the output files of the cnv call job
    cnv_details = dxpy.bindings.dxjob.DXJob(dxid=cnv_call_job).describe()
    gcnv_output = cnv_details["output"]["result_files"]
    gcnv_input = cnv_details["input"]["bambais"]

    # get and store name and file ID of all input and output files
    with concurrent.futures.ThreadPoolExecutor(max_workers=32) as executor:
        # submit jobs mapping each id to describe call
        concurrent_jobs = {
            executor.submit(
                dxpy.describe, file, fields={"name": True, "id": True}
            ): file
            for file in gcnv_input + gcnv_output
        }

        for future in concurrent.futures.as_completed(concurrent_jobs):
            # access returned output as each is returned in any order
            try:
                data = future.result()
                calling_files[data["name"]] = data["id"]
            except Exception as exc:
                # catch any errors that might get raised during querying
                print(
                    f"Error getting data for {concurrent_jobs[future]}: {exc}"
                )

    print(f"Found {len(calling_files.keys())} gCNV input / output files")

    return calling_files


def get_cnv_file_ids(reports, gcnv_dict) -> dict:
    """
    Gather all related file ids for cnv files

    Args:
        reports (list): list of cnv report objects from dxpy search
        gcnv_dict(dict): dictionary output of gcnv call i/o files

    Returns:
        cnv_data (dict): Nested dictionary of files with sample name
            as key and nested dict of CNV files by clinical indication

    Example:
    {
        'X123456-GM1234567-23NGCEN23-1234-F-99347387: {
            'sample': 'X123456-GM1234567-23NGCEN23-1234-F-99347387',
            'Alignment BAI': {'$dnanexus_link: 'file-XX'},
            'Alignment BAI': {'$dnanexus_link: 'file-XY'},
            'clinical_indications': {
                'R414.1_APC associated Polyposis_G': {
                    'CNV': [
                        {
                            'CNV variant report': 'file-ZX',
                            'CNV calls for IGV': {'$dnanexus_link': 'file-TU'},
                            'CNV count': 0,
                            'CNV excluded regions': 'file-CD',
                            'Session file name': 'R414.1_CNV_1'
                        }
                    ],
                'R208.1_Inherited breast cancer and ovarian cancer_P': {
                    'CNV': [
                        {
                            'CNV variant report': 'file-OP',
                            'CNV calls for IGV': {'$dnanexus_link': 'file-LM'},
                            'CNV count': 0,
                            'CNV excluded regions': 'file-AB',
                            'Session file name': 'R208.1_CNV_1'
                        }
                    ]
                }
            },
            'CNV visualisation': 'file-ABC'
        }
    }
    """

    def _find(report):
        """
        Find files for single sample

        Args:
            report (dict): dx describe return for single sample

        Returns:
            data (dict): files found for sample
        """
        report_name = report["describe"]["name"]
        sample = report_name.split("_")[0]
        # Get end of CNV report name for naming IGV session file later
        # because if multiple panels exist for 1 sample we don't want to name
        # the separate session files the same
        cnv_file_ending = report_name.split("_", 1)[1].replace(".xlsx", "")

        # Get file 'details' of the CNV xlsx report so that we can query
        # the clinical indication and variant count. If file has no 'details'
        # set default return to "Unknown"
        file_details = dxpy.DXFile(report["id"]).get_details()
        clin_ind = file_details.get("clinical_indication", "Unknown")
        cnv_variant_count = file_details.get("variants", "Unknown")

        gen_xlsx_job = report["describe"]["createdBy"]["job"]

        excluded_regions_id = dxpy.bindings.dxjob.DXJob(
            dxid=gen_xlsx_job
        ).describe()["input"]["additional_files"][0]["$dnanexus_link"]

        # Find the reports workflow analysis id
        reports_analysis = dxpy.bindings.dxjob.DXJob(
            dxid=gen_xlsx_job
        ).describe()["parentAnalysis"]
        reports_details = dxpy.bindings.dxanalysis.DXAnalysis(
            dxid=reports_analysis
        ).describe()

        # Get the CNV report and seg file using the analysis id
        try:
            cnv_workbook_id = reports_details["output"][
                "stage-cnv_generate_workbook.xlsx_report"
            ]
        except KeyError:
            cnv_workbook_id = reports_details["output"][
                "stage-GFfYY9j4qq8ZxpFpP8zKG7G0.xlsx_report"
            ]
        try:
            seg_id = reports_details["output"][
                "stage-cnv_additional_tasks.seg_file"
            ]
        except KeyError:
            seg_id = reports_details["output"][
                "stage-GG2z5yQ4qq8vb2xp4pB8XByz.seg_file"
            ]

        # gCNV job has all input bams and all outputs go through info
        # saved in the dictionary and find the sample specific files required
        bam = bai = gcnv_bed = ""
        for k, v in gcnv_dict.items():
            if k.startswith(f"{sample}"):
                if k.endswith(".bam"):
                    bam = v
                elif k.endswith(".bai"):
                    bai = v
                elif k.endswith(".bed.gz"):
                    gcnv_bed = v

        # Store in dictionary to return
        data = {
            "sample": sample,
            "Alignment BAM": bam,
            "Alignment BAI": bai,
            "CNV visualisation": gcnv_bed,
            "clinical_indications": {
                clin_ind: {
                    "CNV": [
                        {
                            "CNV variant report": cnv_workbook_id,
                            "CNV count": str(cnv_variant_count),
                            "CNV excluded regions": excluded_regions_id,
                            "CNV calls for IGV": seg_id,
                            "Session file name": cnv_file_ending,
                        }
                    ]
                }
            },
        }
        return data

    cnv_data = defaultdict(
        lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    )

    with concurrent.futures.ThreadPoolExecutor(max_workers=32) as executor:
        # submit jobs mapping each id to describe call
        concurrent_jobs = {
            executor.submit(_find, report): report for report in reports
        }

        for future in concurrent.futures.as_completed(concurrent_jobs):
            # access returned output as each is returned in any order
            try:
                data = future.result()

                # Get info from the returned dict (will only ever be one
                # clinical indication / CNV file dict)
                sample = data["sample"]
                clin_ind = list(data["clinical_indications"])[0]
                cnv_info = data["clinical_indications"][clin_ind]["CNV"][0]

                # If sample and clinical indication already exists, append CNV
                # file info to existing list.
                # If sample doesn't exist yet just add the data as is
                cnv_data[sample]["sample"] = sample
                cnv_data[sample]["Alignment BAM"] = data["Alignment BAM"]
                cnv_data[sample]["Alignment BAI"] = data["Alignment BAI"]
                cnv_data[sample]["CNV visualisation"] = data[
                    "CNV visualisation"
                ]
                cnv_data[sample]["clinical_indications"][clin_ind][
                    "CNV"
                ].append(cnv_info)

            except Exception as exc:
                # catch any errors that might get raised during querying
                print(
                    f"Error getting data for {concurrent_jobs[future]}: {exc}"
                )
                # Propagate the exception to the caller so doesn't silently fail.
                raise exc

    return cnv_data


def find_snv_files(reports, build) -> dict:
    """
    Gather files related to SNV reports

    Args:
        reports (list): List of SNV report dxpy describe dicts
        build (int): Genome build used for SNV calling, e.g. 37 or 38.

    Returns:
        snv_data (dict): Nested dictionary of files with sample name
            as key and list of SNV files within each clinical indication

    Example:
    {
        'X123456-GM1234567-23NGCEN23-1234-F-99347387: {
            'sample': 'X123456-GM1234567-23NGCEN23-1234-F-99347387',
            'Alignment BAI': {'$dnanexus_link: 'file-XX'},
            'Alignment BAI': {'$dnanexus_link: 'file-XY'},
            'clinical_indications': {
                'R414.1_APC associated Polyposis_G': {
                    'SNV': [
                        {
                            'SNV variant report': 'file-ZX',
                            'Coverage report': {'$dnanexus_link': 'file-GH'},
                            'Summary text file': 'file-NM',
                            'SNV count': 1
                        }
                    ],
                'R208.1_Inherited breast cancer and ovarian cancer_P': {
                    'SNV': [
                        {
                            'SNV variant report': 'file-GH',
                            'Coverage report': {'$dnanexus_link': 'file-GH'},
                            'Summary text file': 'file-JK',
                            'SNV count': 3
                        }
                    ]
                }
            }
        }
    }
    """

    def _find(report, build):
        """
        Find files for single sample

        Args:
            report (dict): dx describe return for single sample
            build (int): Genome build used for SNV calling, e.g. 37 or 38.

        Returns:
            data (dict): files found for sample#
        """
        # Get sample name
        sample = report["describe"]["name"].split("_")[0]
        # Get file 'details' of the SNV xlsx report so that we can query
        # the clinical indication and variant count. If file has no 'details'
        # set default return to "Unknown"
        file_details = dxpy.DXFile(report["id"]).get_details()
        clinical_indication = file_details.get(
            "clinical_indication", "Unknown"
        )
        snv_variant_count = file_details.get("included", "Unknown")

        # Get the job id that created the report
        job_id = report["describe"]["createdBy"]["job"]

        # Get the workflow id that included the job
        report_parent_analysis = dxpy.bindings.dxjob.DXJob(dxid=job_id).describe()[
            "parentAnalysis"
        ]

        # Get the vcf file id and athena coverage file id
        report_parent_details = dxpy.bindings.dxanalysis.DXAnalysis(
            dxid=report_parent_analysis
        ).describe()
        try:
            vcf_file = report_parent_details["input"]["stage-rpt_vep.vcf"]
        except KeyError:
            vcf_file = report_parent_details["input"][
                "stage-G9Q0jzQ4vyJ3x37X4KBKXZ5v.vcf"
            ]
        try:
            coverage_report = report_parent_details["output"][
                "stage-rpt_athena.report"
            ]
        except KeyError:
            coverage_report = report_parent_details["output"][
                "stage-Fyq5z18433GfYZbp3vX1KqjB.report"
            ]

        summary_text = (
            report_parent_details.get("output", {})
            .get("stage-rpt_athena.summary_text", {})
            .get("$dnanexus_link")
        )
        if not summary_text:
            print(
                "No summary .txt file found in output of eggd_athena stage"
                f" for SNV reports workflow ({report_parent_analysis})"
            )
        # Logic for extracting bam and bai files
        mappings_bam = mappings_bai = None
        # Extract the additional regions calling job id from the vcf metadata
        if build == 38:
            # For genome build 38, the additional calling job is stored in the
            # vcf file metadata
            additional_calling_job_id = dxpy.describe(vcf_file)["createdBy"][
                "job"
            ]
            if not additional_calling_job_id:
                print(
                    "No additional calling job id found in vcf file metadata. "
                    "Assuming the vcf was created by a sentieon job."
                )
            else:
                parent_vcf_job_details = dxpy.bindings.dxjob.DXJob(
                    dxid=additional_calling_job_id
                ).describe()
        else:
            # If the vcf file was created by a sentieon job, we can find the
            # sentieon job id from the vcf file metadata
            sentieon_job_id = dxpy.describe(vcf_file).get(
                "createdBy", {}
            ).get("job", None)
            if not sentieon_job_id:
                print("No sentieon job id found in vcf file metadata.")
            else:
                parent_vcf_job_details = dxpy.bindings.dxjob.DXJob(
                    dxid=sentieon_job_id
                ).describe()

        parent_dias_single_analysis = (
            parent_vcf_job_details.get("parentAnalysis", None)
            )
        dias_single_analysis_details = dxpy.bindings.dxanalysis.DXAnalysis(
                dxid=parent_dias_single_analysis
            ).describe()

        if dias_single_analysis_details:
            # Get bam & bai job id from sention job metadata
            try:
                mappings_bam_stage = dias_single_analysis_details["output"]["stage-sentieon_dnaseq.mappings_bam"]
                mappings_bam = mappings_bam_stage.get("$dnanexus_link", None)
                mappings_bai_stage = dias_single_analysis_details["output"]["stage-sentieon_dnaseq.mappings_bam_bai"]
                mappings_bai = mappings_bai_stage.get("$dnanexus_link", None)
            except KeyError as err:
                print(
                    "No mappings bam or bai found in output of sentieon_dnaseq stage"
                    f" for dias single workflow ({parent_dias_single_analysis})"
                )
                raise err
        else:
            # If no parent analysis found
            print("No parent analysis found for dias single workflow.")

        # Store in dictionary to return
        data = {
            "sample": sample,
            "Alignment BAM": mappings_bam,
            "Alignment BAI": mappings_bai,
            "clinical_indications": {
                clinical_indication: {
                    "SNV": [
                        {
                            "SNV variant report": report["describe"]["id"],
                            "Coverage report": coverage_report,
                            "Summary text file": summary_text,
                            "SNV count": str(snv_variant_count),
                        }
                    ]
                }
            },
        }

        return data

    # Create a nested defaultdict to store SNV data
    snv_data = defaultdict(
        lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    )

    with concurrent.futures.ThreadPoolExecutor(max_workers=32) as executor:
        # submit jobs mapping each id to describe call
        concurrent_jobs = {
            executor.submit(_find, report, build): report for report in reports
        }

        for future in concurrent.futures.as_completed(concurrent_jobs):
            # access returned output as each is returned in any order
            try:
                data = future.result()
                # Get info from the returned dict (only ever 1 clinical
                # indication and one set of SNV files)
                sample = data["sample"]
                clin_ind = list(data["clinical_indications"])[0]
                snv_files = data["clinical_indications"][clin_ind]["SNV"][0]

                # Add in sample name, BAM and BAI info then append the SNV
                # file info to the relevant clinical indication
                snv_data[sample]["sample"] = data["sample"]
                snv_data[sample]["Alignment BAM"] = data["Alignment BAM"]
                snv_data[sample]["Alignment BAI"] = data["Alignment BAI"]
                snv_data[sample]["clinical_indications"][clin_ind][
                    "SNV"
                ].append(snv_files)

            except Exception as exc:
                # catch any errors that might get raised during querying
                print(
                    f"Error getting data for {concurrent_jobs[future]}: {exc}"
                )
                # Propagate the exception to the caller so doesn't silently fail.
                raise exc

    return snv_data
