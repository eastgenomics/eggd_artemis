#!/usr/bin/env python
# eggd_artemis

import os
import pip
import dxpy
import json
import datetime
import logging
from copy import deepcopy

# Install merge dict package
for package in os.listdir("/home/dnanexus/packages"):
    print(f"Installing {package}")
    pip.main(["install","--no-index","--no-deps",f"packages/{package}"])

from mergedeep import merge


def find_snv_files(reports):
    """Gather files related to SNV reports

    Args:
        reports (list): List of SNV report dxpy describe dicts

    Returns:
        snv_data (dict): Nested dictionary of files with sample name
            as key and a list of files as values with key labels
    """
    snv_data = {}

    for report in reports:

        # Get sample name
        sample = report['describe']['name'].split("_")[0]

        # Get the job id that created the report
        job_id = report['describe']['createdBy']['job']

        # Get the workflow id that included the job
        parent_analysis = dxpy.bindings.dxjob.DXJob(
            dxid=job_id).describe()["parentAnalysis"]

        # Get the vcf file id
        vcf_file = dxpy.bindings.dxanalysis.DXAnalysis(
            dxid=parent_analysis).describe()["input"]["stage-G9Q0jzQ4vyJ3x37X4KBKXZ5v.vcf"]

        # Get the athena coverage file id
        coverage_report = dxpy.bindings.dxanalysis.DXAnalysis(
            dxid=parent_analysis).describe()["output"]["stage-Fyq5z18433GfYZbp3vX1KqjB.report"]

        # Extract the sention job id from the vcf metadata
        sention_job_id=dxpy.describe(vcf_file)["createdBy"]["job"]

        # Get bam & bai job id from sention job metadata
        mappings_bam = dxpy.bindings.dxjob.DXJob(
            dxid=sention_job_id).describe()["output"]["mappings_bam"]

        mappings_bai = dxpy.bindings.dxjob.DXJob(
            dxid=sention_job_id).describe()["output"]["mappings_bam_bai"]

        # Store in dictionary to return
        snv_data[sample] = {
            "SNV variant report": report['describe']['id'],
            "Coverage report": coverage_report,
            'Alignment BAM': mappings_bam,
            'Alignment BAI': mappings_bai
        }

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
    vcf_id = dxpy.bindings.dxanalysis.DXAnalysis(
        dxid=reports_analysis).describe()["input"]["stage-GFYvJF04qq8VKgq34j30pZZ3.vcf"]

    # Find the cnv call job id
    cnv_call_job = dxpy.describe(vcf_id)["createdBy"]["job"]

    # Store all file ids and names in a dictionary
    calling_files={}

    # Find the output files of the cnv call job
    gcnv_output=(dxpy.bindings.dxjob.DXJob(dxid=cnv_call_job).describe()["output"]["result_files"])
    gcnv_input=(dxpy.bindings.dxjob.DXJob(dxid=cnv_call_job).describe()["input"]["bambais"])

    # Store the input files in the dictionary with sample name as key
    # and file-id as value
    for file in gcnv_output:
        file_details = dxpy.describe(file,fields={"name":True,"id":True})
        calling_files[file_details['name']] = file_details['id']

    # Store the output files in the dictionary
    for file in gcnv_input:
        file_details = dxpy.describe(file,fields={"name":True,"id":True})
        calling_files[file_details['name']] = file_details['id']

    return calling_files

def get_cnv_file_ids(reports,gcnv_dict):
    """ Gather all related file ids for cnv files

    Args:
        reports (list): list of cnv report objects from dxpy search
        gcnv_dict(dict): dictionary output of gcnv call i/o files

    Returns:
        cnv_data (dict): Nested dictionary of files with sample name
            as key and a list of files as values with key labels
    """
    cnv_data = {}

    for report in reports:
        sample = report['describe']['name'].split("_")[0]
        gen_xlsx_job = report["describe"]["createdBy"]["job"]

        # Find the reports workflow analysis id
        reports_analysis = dxpy.bindings.dxjob.DXJob(
            dxid=gen_xlsx_job).describe()["parentAnalysis"]
        # Get the CNV report and seg file using the analysis id
        cnv_workbook_id = dxpy.bindings.dxanalysis.DXAnalysis(
            dxid=reports_analysis).describe()["output"]["stage-GFfYY9j4qq8ZxpFpP8zKG7G0.xlsx_report"]
        seg_id = dxpy.bindings.dxanalysis.DXAnalysis(
            dxid=reports_analysis).describe()["output"]["stage-GG2z5yQ4qq8vb2xp4pB8XByz.seg_file"]

        # gCNV job has all input bams and all outputs go through info
        # saved in the dictionary and find the sample specific files required
        for k,v in gcnv_dict.items():
            if k.startswith(f'{sample}-'):
                if k.endswith(".bam"):
                    bam = v
                elif k.endswith(".bai"):
                    bai = v
                elif k.endswith(".bed.gz"):
                    gcnv_bed = v

        cnv_data[sample] = {
            'Alignment BAM': bam,
            'Alignment BAI': bai,
            'CNV variant report': cnv_workbook_id,
            'CNV visualisation': gcnv_bed,
            'CNV calls for IGV': seg_id
        }


    return cnv_data

def get_excluded_intervals(gcnv_output_dict):
    """ Get the excluded regions file from the gcnv dictionary

    Args:
        gcnv_output_dict (dict): dictionary of gcnv i/o files

    Returns:
        excluded_file (str): file id of the excluded regions file
    """
    for k,v in gcnv_output_dict.items():
        if k.endswith("excluded_intervals.bed"):
            excluded_file = v

    return excluded_file

def get_multiqc_report(path_to_reports,project):
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
        file_info = dxpy.bindings.dxfile.DXFile(
                dxid=file_id, project=project)

        # Extract the file name to allow it to be used in the url
        file_name = file_info.describe()["name"]

        # Duration currently defaults to 28 days unless provided as input
        file_url = file_info.get_download_url(
                duration=url_duration, preauthenticated=True,
                project=project, filename=file_name)[0]

        # Without unsetting and clearing the workstation environment
        # the output urls have the local hostname which doesn't work
        # Find hostname
        workstation_hostname = file_url.split('/')[2]
        dx_hostname = 'dl.ec1.dnanex.us'

        # Replace with eu specific hostname and update protocol
        file_url = file_url.replace(workstation_hostname,dx_hostname)
        file_url = file_url.replace('http','https')

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
    sessions_list, sample, bam_url, bai_url, bed_url,
    seg_url, excluded_url, targets_url, DX_PROJECT, expiry_date):
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
        DX_PROJECT (string): id of project
        expiry_date (string): date of expiry of file created

    Returns:
        session_file (string): IGV session file URL
        sessions_list (list): input list with current session appended
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

    # Get current job id
    DX_JOB_ID = os.environ.get("DX_JOB_ID")

    # Get output folder set for this job
    job_output = dxpy.bindings.dxjob.DXJob(DX_JOB_ID).describe()['folder']
    output_folder = f'{job_output}/igv_sessions'

    # Set the environment context to allow upload
    dxpy.set_workspace_id(DX_PROJECT)

    # Create folder if it doesn't exist
    dxpy.api.project_new_folder(
        DX_PROJECT,
        input_params={
            "folder": output_folder,
            "parents": True})

    session_file = dxpy.upload_local_file(
        output_name,
        folder=output_folder,
        tags=[expiry_date],
        wait_on_close=True)

    # Append session file to output list
    sessions_list.append(session_file)

    session_file_id = session_file.get_id()

    return session_file_id, sessions_list


@dxpy.entry_point('main')
def main(url_duration, snv_path=None, cnv_path=None,bed_file=None,qc_status=None):

    # Set up logging
    logger = logging.getLogger(__name__)
    logger.addHandler(dxpy.DXLogHandler())
    logger.propagate = False
    logger.setLevel(logging.DEBUG)

    # Get the project ID
    DX_PROJECT = os.environ.get("DX_PROJECT_CONTEXT_ID")
    print(DX_PROJECT)

    # Get name of project for output naming
    DX_PROJECT_NAME = dxpy.describe(DX_PROJECT)['name']

    # Gather required SNV files if SNV path is provided
    if snv_path is not None:
        logger.info("Gathering Small variant files")

        snv_data = {}

        for path in snv_path.split(','):
            print (f'Gathering reports from: {path}')
            snv_reports = list(dxpy.bindings.search.find_data_objects(
                name="*xlsx",
                name_mode='glob',
                folder=path,
                project=DX_PROJECT,
                describe=True))

            # Get SNV ids
            snv_files = find_snv_files(snv_reports)
            merge (snv_data,snv_files)
        # Get multiqc report
        multiqc = get_multiqc_report(snv_path.split(',')[0],DX_PROJECT)

    # Gather required CNV files if CNV path is provided
    if cnv_path is not None:
        logger.info("Gathering CNV  files")

        cnv_data = {}

        for path in cnv_path.split(','):
            print (f'Gathering reports from: {path}')
            cnv_reports = list(dxpy.bindings.search.find_data_objects(
                name="*xlsx",
                name_mode='glob',
                folder=path,
                project=DX_PROJECT,
                describe=True))

            gcnv_job_info = get_cnv_call_details(cnv_reports)
            cnv_files = get_cnv_file_ids(cnv_reports,gcnv_job_info)

            # Get excluded intervals file
            excluded_intervals = get_excluded_intervals(gcnv_job_info)
            ex_intervals_url = make_url(excluded_intervals, DX_PROJECT, url_duration)

            merge(cnv_data,cnv_files)
        # Get multiqc report
        multiqc = get_multiqc_report(cnv_path.split(',')[0],DX_PROJECT)

    logger.info("Making URLs for additional files")
    # If a bed file is provided, add to a link to the output
    if bed_file:
        bed_file_url = make_url(bed_file, 'project-Fkb6Gkj433GVVvj73J7x8KbV',url_duration)
    else:
        # Setting as empty to avoid session errors
        bed_file_url = ''

    # If a QC Status xlsx is provided, add to a link to the output
    if qc_status:
        qc_status_url = make_url(qc_status, DX_PROJECT,url_duration)
    else:
        qc_status_url = 'No QC status file provided'

    data = {}

    if snv_path and cnv_path:
        merge(data, snv_data, cnv_data)
        session_files = []

    elif snv_path:
        data = snv_data
    elif cnv_path:
        data = cnv_data
        session_files = []
    else:
        logger.debug("No paths given, exiting...")
        exit(1)


    # Write output file
    # Note2self - add version number to output file name??
    today = datetime.datetime.now().strftime("%y%m%d")
    output_name = f"{DX_PROJECT_NAME[4:]}_{today}.tsv"

    # Set timestamps
    now = datetime.datetime.now().strftime("%Y-%m-%d")
    days2expiry = int(url_duration/86400)
    expiry = datetime.timedelta(days=days2expiry)

    expiry_date=(
        datetime.datetime.strptime(datetime.datetime.now().strftime('%Y-%m-%d'),
         '%Y-%m-%d') + datetime.timedelta(seconds=url_duration)).strftime('%Y-%m-%d')


    # Write file
    logger.info("Writing output file")
    with open(output_name, 'w') as f:
            # Write run specific details at top of file
            f.write(f"Run:\t{DX_PROJECT_NAME}\n\n")
            f.write(f"Date Created:\t{str(datetime.datetime.now().strftime('%Y-%m-%d'))}\n")
            f.write(f"Expiry Date:\t{str(expiry_date)}\n\n")
            f.write("Run level files\n")
            f.write(f"MultiQC report\t{make_url(multiqc, DX_PROJECT, url_duration)}\n")
            f.write(f"QC Status report\t{qc_status_url}\n\n")
            f.write('Per Sample files\n\n')

            # For each sample, write out the available sample URLs
            for sample,details in data.items():

                f.write(f"\nSample ID:\t{sample}\n")

                if "SNV variant report" in details:
                    f.write(f"Coverage report:\t{make_url(details['Coverage report'], DX_PROJECT, url_duration)}\n")
                    f.write(f"Small variant report:\t{make_url(details['SNV variant report'], DX_PROJECT, url_duration)}\n")

                if 'CNV variant report' in details:
                    f.write(f"CNV variant report:\t{make_url(details['CNV variant report'], DX_PROJECT, url_duration)}\n\n")

                bam = make_url(details['Alignment BAM'], DX_PROJECT, url_duration)
                bai = make_url(details['Alignment BAI'] ,DX_PROJECT, url_duration)

                if 'CNV variant report' in details:

                    cnv_bed = make_url(details['CNV visualisation'], DX_PROJECT, url_duration)
                    cnv_seg = make_url(details['CNV calls for IGV'], DX_PROJECT, url_duration)

                    cnv_session, session_files = make_cnv_session(
                        session_files, sample, bam, bai, cnv_bed, cnv_seg,
                        ex_intervals_url, bed_file_url, DX_PROJECT, expiry_date)

                    cnv_session_url=make_url(cnv_session, DX_PROJECT, url_duration)

                    f.write(f"CNV IGV Session:\t{cnv_session_url}\n\n")

                elif "SNV variant report" in details:
                    f.write(f"Alignment BAM:\t{bam}\n")
                    f.write(f"Alignment BAI:\t{bai}\n")

    # Upload output to the platform
    output = {}
    url_file = dxpy.upload_local_file(output_name,tags=[expiry_date])
    output["url_file"] = dxpy.dxlink(url_file)

    if session_files != []:
        output["session_files"] = [dxpy.dxlink(item) for item in session_files]


    return output

dxpy.run()
