#!/usr/bin/env python
# eggd_artemis
import concurrent
import os
import pip
import re
import dxpy
import json
import datetime
import logging
from copy import deepcopy

# Install required packages
for package in os.listdir("/home/dnanexus/packages"):
    print(f"Installing {package}")
    pip.main(["install", "--no-index", "--no-deps", f"packages/{package}"])

from mergedeep import merge
from openpyxl.styles import DEFAULT_FONT, Font
import pandas as pd


def find_snv_files(reports):
    """Gather files related to SNV reports

    Args:
        reports (list): List of SNV report dxpy describe dicts

    Returns:
        snv_data (dict): Nested dictionary of files with sample name
            as key and a list of files as values with key labels
    """
    def _find(report):
        """
        Find files for single sample

        Args:
            report (dict): dx describe return for single sample

        Returns:
            data (dict): files found for sample
        """
        # Get sample name
        sample = report['describe']['name'].split("_")[0]

        # Get the 'details' of the SNV xlsx report
        # If file has no details this will be empty response i.e. '{}' so
        # calling get for the include_variants key will result in None
        snv_variant_count = dxpy.bindings.dxdataobject_functions.get_details(
            report['id']
        ).get('include_variants')

        # Get the job id that created the report
        job_id = report['describe']['createdBy']['job']

        # Get the workflow id that included the job
        parent_analysis = dxpy.bindings.dxjob.DXJob(
            dxid=job_id).describe()["parentAnalysis"]

        # Get the vcf file id and athena coverage file id
        parent_details = dxpy.bindings.dxanalysis.DXAnalysis(
            dxid=parent_analysis).describe()
        try:
            vcf_file = parent_details["input"]["stage-rpt_vep.vcf"]
        except KeyError:
            vcf_file = parent_details["input"]["stage-G9Q0jzQ4vyJ3x37X4KBKXZ5v.vcf"]
        try:
            coverage_report = parent_details["output"]["stage-rpt_athena.report"]
        except KeyError:
            coverage_report = parent_details["output"]["stage-Fyq5z18433GfYZbp3vX1KqjB.report"]

        # Extract the sention job id from the vcf metadata
        sention_job_id=dxpy.describe(vcf_file)["createdBy"]["job"]
        sentieon_details = dxpy.bindings.dxjob.DXJob(dxid=sention_job_id).describe()

        # Get bam & bai job id from sention job metadata
        mappings_bam = sentieon_details["output"]["mappings_bam"]
        mappings_bai = sentieon_details["output"]["mappings_bam_bai"]

        # Store in dictionary to return
        data = {
            "sample": sample,
            "SNV variant report": report['describe']['id'],
            "SNV count": snv_variant_count,
            "Coverage report": coverage_report,
            'Alignment BAM': mappings_bam,
            'Alignment BAI': mappings_bai
        }

        return data

    snv_data = {}

    with concurrent.futures.ThreadPoolExecutor(max_workers=32) as executor:
        # submit jobs mapping each id to describe call
        concurrent_jobs = {
            executor.submit(_find, report) for report in reports
        }

        for future in concurrent.futures.as_completed(concurrent_jobs):
            # access returned output as each is returned in any order
            try:
                data = future.result()
                snv_data[data['sample']] = data
            except Exception as exc:
                # catch any errors that might get raised during querying
                print(f"Error getting data for {concurrent_jobs[future]}: {exc}")

    return snv_data


def get_cnv_call_details(reports):
    """ Gets the job id of the cnv call job and stores input and output
        files in a dictionary for later use

    Args:
        reports (list): list of cnv report objects from dxpy search

    Returns:
        calling_files (str): dictionary storing all input and output information
    """
    # Get the job id of the generate workbook job
    gen_xlsx_job = reports[0]["describe"]["createdBy"]["job"]

    # Find the reports workflow analysis id
    reports_analysis = dxpy.bindings.dxjob.DXJob(
        dxid=gen_xlsx_job).describe()["parentAnalysis"]

    # Find the input vcf id
    try:
        vcf_id = dxpy.bindings.dxanalysis.DXAnalysis(
            dxid=reports_analysis).describe()["input"]["stage-cnv_vep.vcf"]
    except KeyError:
        vcf_id = dxpy.bindings.dxanalysis.DXAnalysis(
            dxid=reports_analysis).describe()["input"]["stage-GFYvJF04qq8VKgq34j30pZZ3.vcf"]

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
            ) for file in gcnv_input + gcnv_output
        }

        for future in concurrent.futures.as_completed(concurrent_jobs):
            # access returned output as each is returned in any order
            try:
                data = future.result()
                calling_files[data['name']] = data['id']
            except Exception as exc:
                # catch any errors that might get raised during querying
                print(f"Error getting data for {concurrent_jobs[future]}: {exc}")

    print(f"Found {len(calling_files.keys())} gCNV input / output files")

    return calling_files


def get_cnv_file_ids(reports, gcnv_dict):
    """ Gather all related file ids for cnv files

    Args:
        reports (list): list of cnv report objects from dxpy search
        gcnv_dict(dict): dictionary output of gcnv call i/o files

    Returns:
        cnv_data (dict): Nested dictionary of files with sample name
            as key and a list of files as values with key labels
    """
    def _find(report):
        """
        Find files for single sample

        Args:
            report (dict): dx describe return for single sample

        Returns:
            data (dict): files found for sample
        """
        sample = report['describe']['name'].split("_")[0]

        # Get the 'details' of the CNV xlsx report
        # If file has no details this will be empty response i.e. '{}' so
        # calling get for the include_variants key will result in None
        cnv_variant_counts = dxpy.bindings.dxdataobject_functions.get_details(
            report['id']
        ).get('include_variants')

        gen_xlsx_job = report["describe"]["createdBy"]["job"]

        # Find the reports workflow analysis id
        reports_analysis = dxpy.bindings.dxjob.DXJob(
            dxid=gen_xlsx_job).describe()["parentAnalysis"]
        reports_details = dxpy.bindings.dxanalysis.DXAnalysis(
            dxid=reports_analysis).describe()

        # Get the CNV report and seg file using the analysis id
        try:
            cnv_workbook_id = reports_details["output"]["stage-cnv_generate_workbook.xlsx_report"]
        except KeyError:
            cnv_workbook_id = reports_details["output"]["stage-GFfYY9j4qq8ZxpFpP8zKG7G0.xlsx_report"]
        try:
            seg_id = reports_details["output"]["stage-cnv_additional_tasks.seg_file"]
        except KeyError:
            seg_id = reports_details["output"]["stage-GG2z5yQ4qq8vb2xp4pB8XByz.seg_file"]

        # gCNV job has all input bams and all outputs go through info
        # saved in the dictionary and find the sample specific files required
        for k, v in gcnv_dict.items():
            if k.startswith(f'{sample}'):
                if k.endswith(".bam"):
                    bam = v
                elif k.endswith(".bai"):
                    bai = v
                elif k.endswith(".bed.gz"):
                    gcnv_bed = v

        data = {
            'sample': sample,
            'Alignment BAM': bam,
            'Alignment BAI': bai,
            'CNV variant report': cnv_workbook_id,
            'CNV count': cnv_variant_counts,
            'CNV visualisation': gcnv_bed,
            'CNV calls for IGV': seg_id
        }

        return data

    cnv_data = {}

    with concurrent.futures.ThreadPoolExecutor(max_workers=32) as executor:
        # submit jobs mapping each id to describe call
        concurrent_jobs = {
            executor.submit(_find, report) for report in reports
        }

        for future in concurrent.futures.as_completed(concurrent_jobs):
            # access returned output as each is returned in any order
            try:
                data = future.result()
                cnv_data[data['sample']] = data
            except Exception as exc:
                # catch any errors that might get raised during querying
                print(f"Error getting data for {concurrent_jobs[future]}: {exc}")

    return cnv_data


def get_excluded_intervals(gcnv_output_dict):
    """ Get the excluded regions file from the gcnv dictionary

    Args:
        gcnv_output_dict (dict): dictionary of gcnv i/o files

    Returns:
        excluded_file (str): file id of the excluded regions file
    """
    for k, v in gcnv_output_dict.items():
        if k.endswith("excluded_intervals.bed"):
            return v
    return None


def get_multiqc_report(path_to_reports, project):
    """ Find appropriate multiqc file

    Args:
        snv_path (string): path to snv reports

    Returns:
        multiqc_file (str): file id of the multiqc file
    """
    # Get path to single from reports path
    single=f"/output/{path_to_reports.split('/')[2]}"

    # Find MultiQC jobs in the project
    multiqc_reports=list(dxpy.bindings.search.find_jobs(
            name_mode='glob',
            name="*MultiQC*",
            state="done",
            project=project,
            describe=True))

    # In case there's more than one MultiQC job
    # Find the one that has the exists in the single path in question
    for report in multiqc_reports:
        if report['describe']["folder"].startswith(single):
            multiqc = report["describe"]["output"]['multiqc_html_report']

    return multiqc


def make_url(file_id, project, url_duration):
    """ Given a file id create a download url

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
            duration=url_duration, preauthenticated=True,
            filename=file_name, project=project)[0]

    # Without unsetting and clearing the workstation environment
    # the output urls have the local hostname which doesn't work
    # Find hostname
    workstation_hostname = file_url.split('/')[2]
    dx_hostname = 'dl.ec1.dnanex.us'

    # Replace with eu specific hostname and update protocol
    file_url = file_url.replace(workstation_hostname, dx_hostname)
    file_url = file_url.replace('http', 'https')

    return file_url


def set_order_map(snv_only=False):
    """ Set the order of the session depending on input

    Args:
        snv_only (bool, optional): If True, returns snv order map

    Returns:
        order_map (dict): order of inputs to the session template
    """

    if snv_only == True:

        order_map = {
            '6': {
                'url': '',
                'indexURL': ''
            },
            '7': {
                'url': '',
                'name': 'vcf'
            }
        }

    else:

        order_map = {
            '6': {
                'url': '',
                'name': 'CNV-bed'
            },
            '7': {
                'url': '',
                'indexURL': ''
            },
            '8': {
                'url': '',
                'name': 'variant-bed'
            },
            '9': {
                'url': '',
                'name': 'excluded-regions'
            },
            '10': {
                'url': '',
                'name': 'CNV-targets-bed'
            }

        }

    return order_map


def make_cnv_session(
    sample, bam_url, bai_url, bed_url,
    seg_url, excluded_url, targets_url, job_output, expiry_date):
    """ Create a session file for IGV

    Args:
        sessions_list (list): list of session files to upload in the end
        sample (string): sample name
        bam_url (string): URL of the bam file
        bai_url (string): URL of the bai file
        bed_url (string): URL of the bed file
        seg_url (string): URL of the seg file
        excluded_url (string): URL of the excluded regions file
        targets_url (string): URL of the targets file
        job_output (str) : output folder set for job
        expiry_date (string): date of expiry of file created

    Returns:
        session_file (string): IGV session file URL
    """
    order_map = set_order_map()

    # Copy template to avoid overwriting
    with open('/home/dnanexus/cnv-template.json') as fh:
        template = json.load(fh)

    template_copy = deepcopy(template)

    template_copy['tracks'] = []

    # Order 6 - Visualisation bed for CNV calls
    order_map['6']['url'] = bed_url

    # Order 7 - Patient bam and index
    order_map['7']['name'] = sample
    order_map['7']['url'] = bam_url
    order_map['7']['indexURL'] = bai_url

    # Order 8 - Variant seg file
    order_map['8']['url'] = seg_url

    # Order 9 - Excluded regions from CNV calling bed
    order_map['9']['url'] = excluded_url

    # Order 10 - Target regions bed
    order_map['10']['url'] = targets_url

    # Map urls to template
    for track in template['tracks']:
        if track['order'] in [6, 7, 8, 9, 10]:
            order = str(track['order'])
            track['url'] = order_map[order]['url']

        if track['order'] == 7:
            track['indexURL'] = order_map[order]['indexURL']

        template_copy['tracks'].append(track)


    output_name = f"{sample}_igv.json"

    with open(output_name, "w") as outfile:
        json.dump(template_copy, outfile, indent=4)

    session_file = dxpy.upload_local_file(
        output_name,
        folder=f'{job_output}/igv_sessions',
        tags=[expiry_date],
        wait_on_close=True)

    session_file_id = session_file.get_id()

    return session_file_id


def generate_sample_urls(
    sample, sample_data, url_duration, ex_intervals_url,
    bed_url, expiry_date, job_output) -> dict:
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
    dict
        dict of sample URLs
    """
    dx_project = os.environ.get("DX_PROJECT_CONTEXT_ID")

    urls = {}
    urls['sample'] = sample

    bam_url = make_url(sample_data['Alignment BAM'], dx_project, url_duration)
    bai_url = make_url(sample_data['Alignment BAI'], dx_project, url_duration)

    urls['bam_url'] = bam_url
    urls['bai_url'] = bai_url

    if "SNV variant report" in sample_data:
        coverage_url = make_url(
            sample_data['Coverage report'], dx_project, url_duration)
        snv_url = make_url(
            sample_data['SNV variant report'], dx_project, url_duration)

        urls['coverage_url'] = (
            f'=HYPERLINK("{coverage_url}", "{coverage_url}")'
        )
        urls['snv_url'] = f'=HYPERLINK("{snv_url}", "{snv_url}")'

    if 'CNV variant report' in sample_data:
        cnv_url = make_url(
            sample_data['CNV variant report'], dx_project, url_duration)

        cnv_bed = make_url(
            sample_data['CNV visualisation'], dx_project, url_duration)
        cnv_seg = make_url(
            sample_data['CNV calls for IGV'], dx_project, url_duration)

        cnv_session = make_cnv_session(
            sample, bam_url, bai_url, cnv_bed, cnv_seg, ex_intervals_url,
            bed_url, job_output, expiry_date)

        cnv_session_url = make_url(cnv_session, dx_project, url_duration)

        urls['cnv_bed'] = cnv_bed
        urls['cnv_seg'] = cnv_seg
        urls['cnv_session_fileid'] = cnv_session
        urls['cnv_url'] = f'=HYPERLINK("{cnv_url}", "{cnv_url}")'
        urls['cnv_session_url'] = (
            f'=HYPERLINK("{cnv_session_url}", "{cnv_session_url}")'
        )

    return urls


def write_output_file(
    sample_urls, file_dict, today, expiry_date, multiqc_url, qc_url, project_name):
    """
    Writes output xlsx file with all download URLs

    Parameters
    ----------
    sample_urls : dict
        generated URLs for each sample
    file_dict : dict
        dict created earlier with variant counts
    today : str
        today date
    expiry_date : str
        date of URL experation
    multiqc_url : str
        download URL for multiQC report
    qc_url : str
        download URL for QC report
    project_name: str
        project name

    Outputs
    -------
    xlsx file
    """
    print("Writing output file")

    multiqc_url = f'=HYPERLINK("{multiqc_url}", "{multiqc_url}")'
    if qc_url.startswith('http'):
        qc_url = f'=HYPERLINK("{qc_url}", "{qc_url}")'

    df = pd.DataFrame(columns=['a', 'b'])
    df = df.append({'a': 'Run:', 'b': project_name}, ignore_index=True)
    df = df.append({}, ignore_index=True)
    df = df.append({'a': 'Date Created:', 'b': today}, ignore_index=True)
    df = df.append({'a': 'Expiry Date:', 'b': expiry_date}, ignore_index=True)
    df = df.append({}, ignore_index=True)
    df = df.append({'a': 'Run Level Files'}, ignore_index=True)
    df = df.append({'a': 'MultiQC report', 'b': multiqc_url}, ignore_index=True)
    df = df.append({'a': 'QC Status Report', 'b': qc_url}, ignore_index=True)
    df = df.append({}, ignore_index=True)
    df = df.append({}, ignore_index=True)
    df = df.append({'a': 'Per Sample Files'}, ignore_index=True)
    df = df.append({}, ignore_index=True)

    sample_urls_plus_metadata = merge(sample_urls, file_dict)

    sample_order = sorted(sample_urls_plus_metadata.keys())

    for sample in sample_order:
        urls = sample_urls_plus_metadata.get(sample)

        df = df.append({'a': sample}, ignore_index=True)

        url_fields = {
            'coverage_url': 'Coverage report:',
            'SNV count': 'Number of SNVs post-filtering',
            'snv_url': 'Small variant report',
            'CNV count': 'Number of CNVs',
            'cnv_url': 'CNV variant report',
            'cnv_session_url': 'CNV IGV Session:',
            'bam_url': 'Alignment BAM',
            'bai_url': 'Alignment BAI'
        }

        for field, label in url_fields.items():
            if urls.get(field):
                if field == 'SNV count':
                    # If SNV count exists and 0 variants pass filtering
                    # Remove the SNV report URL altogether and replace with
                    # NMD text
                    if int(urls.get(field)) == 0:
                        urls['snv_url'] = 'No SNVs passed filtering'

                if field == 'CNV count':
                    # If CNV count exists and 0 variants pass filtering
                    # Remove the CNV report URL altogether and replace with
                    # NMD text
                    if (int(urls.get(field)) == 0):
                        urls['cnv_url'] = 'No CNVs detected'

                if field == 'bam_url':
                    df = df.append({}, ignore_index=True)

                df = df.append(
                    {'a': label, 'b': urls.get(field)}, ignore_index=True
                )

        df = df.append({}, ignore_index=True)
        df = df.append({}, ignore_index=True)


    writer = pd.ExcelWriter(f'{project_name}_{today}.xlsx', engine='openpyxl')
    df.to_excel(writer, index=False, header=False)

    # set column widths
    sheet = writer.sheets['Sheet1']
    sheet.column_dimensions['A'].width = 20
    sheet.column_dimensions['B'].width = 145

    sheet['A1'].font = Font(bold=True, name=DEFAULT_FONT.name)
    sheet['A6'].font = Font(bold=True, name=DEFAULT_FONT.name)
    sheet['A11'].font = Font(bold=True, name=DEFAULT_FONT.name)

    # make sample IDs bold
    for cell in sheet.iter_rows(max_col=1):
        if re.match(r'X[\d]+', cell[0].value):
            sheet[cell[0].coordinate].font = Font(
                bold=True, name=DEFAULT_FONT.name)
        elif re.match(r'[a-zA-Z0-9]+-[a-zA-Z0-9]+-[a-zA-Z0-9]+-[a-zA-Z0-9]+-[MFU]-[a-zA-Z0-9]+', cell[0].value):
            sheet[cell[0].coordinate].font = Font(
                bold=True, name=DEFAULT_FONT.name)

    # make hyperlinks blue
    for cell in sheet.iter_rows(min_col=2, max_col=2):
        if 'HYPERLINK' in str(cell[0].value):
            sheet[cell[0].coordinate].font = Font(
                color='00007f', name=DEFAULT_FONT.name)

    writer.book.save(f'{project_name}_{today}.xlsx')

    print(f"Output written to {project_name}_{today}.xlsx")


@dxpy.entry_point('main')
def main(url_duration, snv_path=None, cnv_path=None, bed_file=None, qc_status=None):

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
        os.environ.get('DX_JOB_ID')).describe()['folder']

    # Make output folder for job
    dxpy.api.project_new_folder(
        DX_PROJECT,
        input_params={
            "folder": job_output,
            "parents": True})

    # Get name of project for output naming
    project_name = dxpy.describe(DX_PROJECT)['name']
    project_name = '_'.join(project_name.split('_')[1:-1])

    # Gather required SNV files if SNV path is provided
    multiqc = None
    if snv_path:
        logger.info("Gathering Small variant files")

        snv_data = {}

        for path in snv_path.split(','):
            print(f'Gathering reports from: {path}')
            snv_reports = list(dxpy.bindings.search.find_data_objects(
                name="*xlsx",
                name_mode='glob',
                folder=path,
                project=DX_PROJECT,
                describe=True))

            # Get SNV ids
            snv_files = find_snv_files(snv_reports)
            merge(snv_data, snv_files)

        print(f"Size of snv data dict: {len(snv_data.keys())}")

        # Get multiqc report
        multiqc = get_multiqc_report(snv_path.split(',')[0], DX_PROJECT)

    # Gather required CNV files if CNV path is provided
    if cnv_path:
        logger.info("Gathering CNV  files")

        # Create folder for session files
        dxpy.api.project_new_folder(
            DX_PROJECT,
            input_params={
                "folder": f"{job_output}/igv_sessions",
                "parents": True})

        cnv_data = {}

        for path in cnv_path.split(','):
            print(f'Gathering reports from: {path}')
            cnv_reports = list(dxpy.bindings.search.find_data_objects(
                name="*xlsx",
                name_mode='glob',
                folder=path,
                project=DX_PROJECT,
                describe=True))

            gcnv_job_info = get_cnv_call_details(cnv_reports)
            cnv_files = get_cnv_file_ids(cnv_reports, gcnv_job_info)

            # Get excluded intervals file
            excluded_intervals = get_excluded_intervals(gcnv_job_info)
            ex_intervals_url = make_url(
                excluded_intervals, DX_PROJECT, url_duration)

            merge(cnv_data, cnv_files)

        print(f"Size of cnv data dict: {len(cnv_data.keys())}")

        # Get multiqc report
        if not multiqc:
            multiqc = get_multiqc_report(cnv_path.split(',')[0], DX_PROJECT)
    else:
        ex_intervals_url = ''

    logger.info("Making URLs for additional files")

    # If a bed file is provided, add to a link to the output
    if bed_file:
        bed_file_url = make_url(
            bed_file, 'project-Fkb6Gkj433GVVvj73J7x8KbV', url_duration)
    else:
        # Setting as empty to avoid session errors
        bed_file_url = ''

    # If a QC Status xlsx is provided, add to a link to the output
    if qc_status:
        qc_status_url = make_url(qc_status, DX_PROJECT, url_duration)
    else:
        qc_status_url = 'No QC status file provided'

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
    expiry_date=(
        datetime.datetime.strptime(datetime.datetime.now().strftime('%Y-%m-%d'),
        '%Y-%m-%d') + datetime.timedelta(seconds=url_duration)
    ).strftime('%Y-%m-%d')


    # generate all urls for each sample
    logger.info("Generating per sample URLs")
    all_sample_urls = {}

    with concurrent.futures.ThreadPoolExecutor(max_workers=32) as executor:
        # submit jobs mapping each id to describe call
        concurrent_jobs = {
            executor.submit(
                generate_sample_urls, sample, sample_data, url_duration,
                ex_intervals_url, bed_file_url, expiry_date, job_output
            ) for sample, sample_data in file_data.items()
        }

        for future in concurrent.futures.as_completed(concurrent_jobs):
            # access returned output as each is returned in any order
            try:
                data = future.result()
                all_sample_urls[data['sample']] = data
            except Exception as exc:
                # catch any errors that might get raised during querying
                print(f"Error getting data for {concurrent_jobs[future]}: {exc}")

    multiqc_url = make_url(multiqc, DX_PROJECT, url_duration)

    write_output_file(
        all_sample_urls, file_data, today, expiry_date, multiqc_url,
        qc_status_url, project_name
    )

    # Upload output to the platform
    output = {}
    url_file = dxpy.upload_local_file(
        f'{project_name}_{today}.xlsx',
        folder=job_output,
        tags=[expiry_date]
    )
    output["url_file"] = dxpy.dxlink(url_file)

    session_files = [
        x.get('cnv_session_fileid') for x in all_sample_urls.values()
        if x.get('cnv_session_fileid')
    ]
    print(f"Found {len(session_files)} session files to link to job output")
    if session_files:
        output["session_files"] = [dxpy.dxlink(item) for item in session_files]


    return output

dxpy.run()
