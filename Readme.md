<!-- dx-header -->
# eggd_artemis (DNAnexus Platform App)


## What does this app do?

Gathers required files and creates file containing urls to allow file download within CUH & the LGLs.

## What inputs are required for this app to run?
* `snv_path`[optional]: List of absolute paths to folder containing the variant reports for small variants
* `cnv_path`[optional]: List of absolute paths to folder containing the variant reports for CNVs
* `url_duration`[optional]: Time (in seconds) until the generated links expire. (Default = 4838400)
* `bed_file`[optional]: Static capture bed file
* `qc_status`[optional]: Input of xlsx file containing agreed QC status for samples in the run
* `multiqc_report`[optional]: Input of the MultiQC report - if not provided, will search for a MultiQC job in the project to find this
* `lock_cells`[optional]: Determines whether to protect any populated cells in the output .xlsx from editing (Default=True)

## How does this app work?

Given a list of paths the app finds the appropriate files needed to share with the scientists by backwards searching through the file metadata for the jobs they were created from and extracting desired files. Current files outputted are:
* Run level file URLs:
  * MultiQC Report `html`
  * QC Status Report `xlsx`
* For small variant analysis, DNAnexus URLs for:
  * Variant Report `xlsx`
  * Coverage Report `html`
  * Alignment `bam`
  * Alignment index `bai`
* For copy number variant analysis, DNAnexus URLs for:
  * CNV variant report `xlsx`
  * CNV calls for IGV `seg`
  * CNV visualisation `bed`
  * CNV target regions `bed`
  * Alignment `bam`
  * Alignment index `bai`

## What does this app output
* `url file`: File contained URLs for the files described above

## Notes
* This app works specifically with the output of the dias pipeline, specifically expecting the path inputs to be in specific format to identify the MultiQC reports.
* Please note although both paths are optional, if none are given the app will fail

## This app was created by East Genomics GLH
