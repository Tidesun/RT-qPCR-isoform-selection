# Introduction

This is a script to select RT-qPCR primers given a gene isoform annotation. <br>
The package include 3 steps: <br>
1. Selecting unique regions as the target sequences of RT-qPCR experiment for each isoform that is not shared with other isoforms.
2. Trying to design primer for each target sequence by primer3 (via Primer3-py the python version of primer3)
3. Use Magic-BLAST to check the specificity of the RT-qPCR product against the input annotation and NONCODE database(http://www.noncode.org)

# Installation
```
conda env create --file=environment.yml
```
# Usage
```
usage: main.py [-h] -gtf GTF_ANNOTATION_PATH -genome REFERENCE_GENOME_FASTA -o
               OUTPUT_DIR [--conda_env_name CONDA_ENV_NAME] [-t THREADS]
               [--min_region_length MIN_REGION_LENGTH]
               [--NONCODE_fasta_url NONCODE_FASTA_URL]
               [--only_identify_unique_region ONLY_IDENTIFY_UNIQUE_REGION]
               [--isoform_list_fpath ISOFORM_LIST_FPATH]

RT-qPCR-isoform-selction

optional arguments:
  -h, --help            show this help message and exit

required named arguments:
  -gtf GTF_ANNOTATION_PATH, --gtf_annotation_path GTF_ANNOTATION_PATH
                        The path of annotation file
  -genome REFERENCE_GENOME_FASTA, --reference_genome_fasta REFERENCE_GENOME_FASTA
                        The path of reference genome
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        The path of output directory

optional named arguments:
  --conda_env_name CONDA_ENV_NAME
                        Conda env name for the
                        tool.[default:qPCR_isoform_selection]
  -t THREADS, --threads THREADS
                        Number of threads.[default:1]
  --min_region_length MIN_REGION_LENGTH
                        Min region length.[default:95]
  --min_intron_length MIN_INTRON_LENGTH
                        Min intron length.[default:1000]
  --num_primers_returns NUM_PRIMERS_RETURNS
                        Maxium number of primers designed for each unique
                        region of transcript.[default:1000]
  --NONCODE_fasta_url NONCODE_FASTA_URL
                        NONCDOE reference fasta.[default:http://www.noncode.or
                        g/datadownload/NONCODEv6_human.fa.gz]
  --only_identify_unique_region ONLY_IDENTIFY_UNIQUE_REGION
                        Whether only identify unique region of isoforms and do
                        not find primers.[default:False]
  --isoform_list_fpath ISOFORM_LIST_FPATH
                        Path of isoform_list_fpath.[default:None]
```
# Input files
## --isoform_list_fpath
- Refer to `examples/isoform_list.txt`. A one column tsv file with isoform IDs consistent with the annotation file. <br>
# Output files
## unique_region.tsv
- SEQUENCE_ID: Gene-isoform-region_id <br>
- SEQUENCE_TEMPLATE: region sequence<br>
- region_type: junction or exon<br>
- Chr <br>
- Strand<br>
- 5_end_and_3_end_distance_difference: abs(distance_to_5_end - distance_to_3_end). Lower value indicates it is located far from 5' and 3' end. For primer design, the target region should be as far from 5' and 3' end as possible.<br>
- exon_start: The start position of all the exons of the region. 1 base indexed. <br>
- exon_end: The end position of the all the exon of the region. 1 base indexed. <br>
- intron_length: The total length of intron region.
- exon_length: The total length of the exon region.
## unique_region_with_primers.tsv
- Same format as unique region. Only show the regions with valid primers.<br>
## all_primers_overview.tsv
- The primers' info of all isoforms. Refer to https://primer3.org/manual.html for the meaning of the fields. <br>
## designs/[ISOFORM_ID]/
- *.gtf: primer annotation
- primer_overview.tsv: same format as all_primers_overview.tsv but for only one isoform.
