import argparse
from pathlib import Path

from prepare_target_sequence import prepare_target_sequence_dict
from design_primers import design_primers
def parse_arguments():
    """
    Parse the arguments
    """
    parser = argparse.ArgumentParser(description="RT-qPCR-isoform-selction",add_help=True)
    required = parser.add_argument_group('required named arguments')
    required.add_argument('-gtf','--gtf_annotation_path', type=str, help="The path of annotation file",required=True)
    required.add_argument('-genome','--reference_genome_fasta', type=str, help="The path of reference genome",required=True)
    required.add_argument('-o','--output_dir', type=str, help="The path of output directory",required=True)
    optional = parser.add_argument_group('optional named arguments')
    optional.add_argument('--conda_env_name', type=str, help="Conda env name for the tool.[default:qPCR_isoform_selection]",required=False,default='qPCR_isoform_selection')
    optional.add_argument('-t','--threads', type=int, help="Number of threads.[default:1]",required=False,default=1)
    optional.add_argument('--min_region_length', type=int, help="Min region length.[default:95]",required=False,default=95)
    optional.add_argument('--NONCODE_fasta_url', type=str, help="NONCDOE reference fasta.[default:http://www.noncode.org/datadownload/NONCODEv6_human.fa.gz]",required=False,default='http://www.noncode.org/datadownload/NONCODEv6_human.fa.gz')
    optional.add_argument('--only_identify_unique_region', type=str, help="Whether only identify unique region of isoforms and do not find primers.[default:False]",required=False,default='False')
    optional.add_argument('--isoform_list_fpath', type=str, help="Path of isoform_list_fpath.[default:None]",required=False,default=None)
    args = parser.parse_args()
    return args
import time
def main():
    args = parse_arguments()
    print('\n'.join(f'{k}={v}' for k, v in vars(args).items()))
    ref_file_path = args.gtf_annotation_path
    threads = args.threads
    output_dir = args.output_dir
    reference_genome_path = args.reference_genome_fasta
    lower_region_length = args.min_region_length
    conda_env_name = args.conda_env_name
    NONCODE_fasta_url = args.NONCODE_fasta_url
    isoform_list_fpath = args.isoform_list_fpath
    Path(f'{output_dir}/temp').mkdir(exist_ok=True,parents=True)
    Path(f'{output_dir}/temp/blast_res/').mkdir(exist_ok=True,parents=True)
    Path(f'{output_dir}/temp/target_sequences/').mkdir(exist_ok=True,parents=True)
    st = time.time()
    prepare_target_sequence_dict(ref_file_path,reference_genome_path,lower_region_length,output_dir, isoform_list_fpath,threads)
    duration = time.time() - st
    duration_minutes = duration/60
    print(f'Prepare target sequence done in {duration} seconds / {duration_minutes} minutes!',flush=True)
    if args.only_identify_unique_region == 'False':
        st = time.time()
        design_primers(ref_file_path,reference_genome_path,output_dir,conda_env_name,NONCODE_fasta_url,threads)
        duration = time.time() - st
        duration_minutes = duration/60
        print(f'Design primers done in {duration} seconds / {duration_minutes} minutes!',flush=True)
main()