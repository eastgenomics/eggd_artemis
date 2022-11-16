<!-- dx-header -->
# eggd_artemis (DNAnexus Platform App)


## What does this app do?

Gathers required files and creates file containing urls to allow file download within CUH & the LGLs.

## What inputs are required for this app to run?
* `snv_path`: List of absolute paths to folder containing the variant reports for small variants
* `cnv_path` : List of absolute paths to folder containing the variant reports for CNVs
* `url_duration`: Time (in seconds) until the generated links expire
* `bed_file`: Static capture bed file


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


## This app was created by East Genomics GLH
