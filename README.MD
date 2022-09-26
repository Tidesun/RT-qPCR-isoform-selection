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
  --NONCODE_fasta_url NONCODE_FASTA_URL
                        NONCDOE reference fasta.[default:http://www.noncode.or
                        g/datadownload/NONCODEv6_human.fa.gz]
  --only_identify_unique_region ONLY_IDENTIFY_UNIQUE_REGION
                        Whether only identify unique region of isoforms and do
                        not find primers.[default:False]
```