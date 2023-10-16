import subprocess
from pathlib import Path
import multiprocessing as mp
import pickle

from Bio.Seq import Seq
import pandas as pd
from primer3 import bindings
import config

from GTFBasics import GTFFile
def design_primer_for_isoform(isoform_target_sequence_dict,allow_exon=True):
    all_primer_design = {}
    all_raw_primer_design_result = {}
    if allow_exon:
        region_types = ['junction','exon']
    else:
        region_types = ['junction']
    for region_type in region_types:
        if len(all_primer_design) != 0:
            break
        for region_dict,info in zip(isoform_target_sequence_dict[region_type],isoform_target_sequence_dict[f'{region_type}_info']):
            try:
                res = bindings.designPrimers(
                region_dict,
                {
                    'PRIMER_OPT_SIZE': 20,
                    'PRIMER_MIN_SIZE': 18,
                    'PRIMER_MAX_SIZE': 25,
                    'PRIMER_MIN_GC':40,
                    'PRIMER_OPT_GC_PERCENT':50,
                    'PRIMER_MAX_GC':60,
                    'PRIMER_NUM_RETURN':config.num_primers_returns,
                    'PRIMER_SECONDARY_STRUCTURE_ALIGNMENT':1,
                    'PRIMER_PRODUCT_SIZE_RANGE': [[80,250]],
                })
            except Exception as e:
                print(e)
                print(len(region_dict['SEQUENCE_TEMPLATE']))
                continue
            primer_design = {}
            for key,seq in res.items():
                if 'SEQUENCE' in key:
                    primer_id = key.split('_')[2]
                    direction = key.split('_')[1]
                    if primer_id not in primer_design:
                        primer_design[primer_id] = {}
                    primer_design[primer_id][direction] = seq
                    if seq in info['seq'].replace('^',''):
                        primer_design[primer_id][direction+'_start'] = info['seq'].replace('^','').index(seq)
                    else:
                        primer_design[primer_id][direction+'_start'] = info['seq'].replace('^','').index(str(Seq(seq).reverse_complement()))
                    primer_design[primer_id][direction+'_len'] = res[f'PRIMER_{direction}_{primer_id}'][1]
                    primer_design[primer_id][direction+'_PENALTY'] = res[f'PRIMER_{direction}_{primer_id}_PENALTY']
                    primer_design[primer_id][direction+'_TM'] = res[f'PRIMER_{direction}_{primer_id}_TM']
                    primer_design[primer_id][direction+'_GC_PERCENT'] = res[f'PRIMER_{direction}_{primer_id}_GC_PERCENT']
                    primer_design[primer_id][direction+'_SELF_ANY_TH'] = res[f'PRIMER_{direction}_{primer_id}_SELF_ANY_TH']
                    primer_design[primer_id][direction+'_SELF_END_TH'] = res[f'PRIMER_{direction}_{primer_id}_SELF_END_TH']
                    primer_design[primer_id][direction+'_END_STABILITY'] = res[f'PRIMER_{direction}_{primer_id}_END_STABILITY']
                    primer_design[primer_id][direction+'_HAIRPIN_TH'] = res[f'PRIMER_{direction}_{primer_id}_HAIRPIN_TH']
                    primer_design[primer_id]['PAIR_COMPL_ANY_TH'] = res[f'PRIMER_PAIR_{primer_id}_COMPL_ANY_TH']
                    primer_design[primer_id]['PAIR_COMPL_END_TH'] = res[f'PRIMER_PAIR_{primer_id}_COMPL_END_TH']
                    primer_design[primer_id]['PAIR_PRODUCT_SIZE'] = res[f'PRIMER_PAIR_{primer_id}_PRODUCT_SIZE']
                    primer_design[primer_id]['PAIR_PENALTY'] = res[f'PRIMER_PAIR_{primer_id}_PENALTY']
                    primer_design[primer_id]['target_sequence'] = info['seq']
            invalid_primer_id_set = set()
            if region_type == 'junction':
                strand = info['strand']
                assert strand in ['+','-']
                exon_starts,exon_ends = info['exon_start'],info['exon_end']
                exon_lens = [e-s+1 for s,e in zip(exon_starts,exon_ends)]
               
                
                for primer_id,design in primer_design.copy().items():
                    is_invalid_primer_design = False
                    if strand == '+':
                        design_LEFT_start_genome = design['LEFT_start']
                        design_LEFT_len = design['LEFT_len']
                        design_RIGHT_start_genome = design['RIGHT_start']
                        design_RIGHT_len = design['RIGHT_len']
                    elif strand == '-':
                        design_LEFT_start_genome = design['RIGHT_start']
                        design_LEFT_len = design['RIGHT_len']
                        design_RIGHT_start_genome = design['LEFT_start']
                        design_RIGHT_len = design['LEFT_len']

                    left_remaing_flanking = design_LEFT_start_genome
                    for exon_index in range(0,len(exon_lens)):
                        exon_len = exon_lens[exon_index]
                        if left_remaing_flanking < exon_len:
                            left_start = exon_starts[exon_index] + left_remaing_flanking
                            if left_start + design_LEFT_len - 1 >= exon_ends[exon_index]:
                                # del primer_design[primer_id]
                                is_invalid_primer_design = True
                            else:
                                left_end = left_start + design_LEFT_len - 1
                            break
                        else:
                            left_remaing_flanking -= exon_len
                    if is_invalid_primer_design:
                        # print(primer_id)
                        invalid_primer_id_set.add(primer_id)
                        continue 

                    right_remaing_flanking = design_RIGHT_start_genome
                    for exon_index in range(0,len(exon_lens)):
                        exon_len = exon_lens[exon_index]
                        if right_remaing_flanking < exon_len:
                            right_start = exon_starts[exon_index] + right_remaing_flanking
                            if right_start + design_RIGHT_len - 1 >= exon_ends[exon_index]:
                                # del primer_design[primer_id]
                                is_invalid_primer_design = True
                            else:
                                right_end = right_start + design_RIGHT_len - 1
                            break
                        else:
                            right_remaing_flanking -= exon_len
                    if is_invalid_primer_design:
                        # print(primer_id)
                        invalid_primer_id_set.add(primer_id)
                        continue 
                    # if left and right primer in the same exon; dropped
                    for exon_start,exon_end in zip(exon_starts,exon_ends):
                        if (left_start >= exon_start) & (left_end <= exon_end) & (right_start >= exon_start) & (right_end <= exon_end):
                            # del primer_design[primer_id]
                            is_invalid_primer_design = False
                            break
                    if is_invalid_primer_design:
                        # print(primer_id)
                        invalid_primer_id_set.add(primer_id)
                        continue 
                    else:
                        if strand == '+':
                            primer_design[primer_id]['LEFT_start_genome'] = left_start
                            primer_design[primer_id]['LEFT_end_genome'] = left_end
                            primer_design[primer_id]['RIGHT_start_genome'] = right_start
                            primer_design[primer_id]['RIGHT_end_genome'] = right_end
                        elif strand == '-':
                            primer_design[primer_id]['LEFT_start_genome'] = right_start
                            primer_design[primer_id]['LEFT_end_genome'] = right_end
                            primer_design[primer_id]['RIGHT_start_genome'] = left_start
                            primer_design[primer_id]['RIGHT_end_genome'] = left_end
                    # if design['LEFT_start'] >=  junc_pos or design['RIGHT_start'] - design['RIGHT_len'] +1 <= junc_pos:
                    #     # if junc point outside the product
                    #     # if ((junc_pos < design['LEFT_start'] + design['LEFT_len'] - 1) & (junc_pos > design['LEFT_start'])) or  ((junc_pos > design['RIGHT_start'] - design['RIGHT_len'] + 1) & (junc_pos < design['RIGHT_start'])):
                    #     #     # if junc point inside the primer
                    #     #     pass
                    #     # else:
                    #         # if junc point not inside the product nor the primer
                    #         del primer_design[primer_id]
                    #         continue
                    # primer_design[primer_id]['product_exons_len_difference'] = abs(abs(junc_pos - design['LEFT_start']) - abs(junc_pos - (design['RIGHT_start']-design['RIGHT_len']+1)))
            if 'PRIMER_PAIR_4_COMPL_END_STUCT' in res:
                print('Possible secondary structure.')
            if len(primer_design) == 0:
                continue
            region_id = region_dict['SEQUENCE_ID']
            for primer_id,design in primer_design.items():
                if primer_id in invalid_primer_id_set:
                    continue
                if region_id not in all_primer_design:
                    all_primer_design[region_id] = {}
                all_primer_design[region_id][primer_id] = design
            # all_primer_design[region_id] = primer_design
            all_raw_primer_design_result[region_id] = res
    return all_primer_design,all_raw_primer_design_result
def build_blast_library(ref_file_path,reference_genome_path,output_dir,conda_env_name,NONCODE_url):
    
    gtf_file = GTFFile(ref_file_path,reference_genome_path)
    ref_transcriptome = f'{output_dir}/temp/reference.fa'
    gtf_file.write_fa(ref_transcriptome)
    completed = subprocess.run(f'wget {NONCODE_url} -O {output_dir}/temp/NONCODE.fa.gz && zcat {output_dir}/temp/NONCODE.fa.gz > {output_dir}/temp/NONCODE.fa', shell=True,stdout=subprocess.DEVNULL)
    completed = subprocess.run(f'conda run -n {conda_env_name} makeblastdb -in {output_dir}/temp/reference.fa -dbtype nucl -parse_seqids -out {output_dir}/temp/reference', shell=True,stdout=subprocess.PIPE,executable = '/bin/bash')
    completed = subprocess.run(f'conda run -n {conda_env_name} makeblastdb -in {output_dir}/temp/NONCODE.fa -dbtype nucl -parse_seqids -out {output_dir}/temp/NONCODE', shell=True,stdout=subprocess.PIPE,executable = '/bin/bash')
def check_for_blast(all_primer_design,isoform,output_dir,conda_env_name,isoform_set=None,allow_shared_region=False):
    is_isoform_primers_available = False
    with open(f'{output_dir}/temp/blast_res/{isoform}.fa','w') as f:
        for region,primer_design in all_primer_design.items():
            for primer_id,primer in primer_design.items():
                f.write(f'>{region}_{primer_id}\n')
                primer_seq_1 = Seq(primer['LEFT']).upper()
                primer_seq_2 = Seq(primer['RIGHT']).reverse_complement().upper()
                f.write(f'{primer_seq_1}NNNNNNNNNN{primer_seq_2}\n')
                is_isoform_primers_available = True
    if not is_isoform_primers_available:
        return {}
    completed = subprocess.run(f'conda run -n {conda_env_name} magicblast -query {output_dir}/temp/blast_res/{isoform}.fa -db {output_dir}/temp/reference -reftype transcriptome -num_threads 1 -outfmt tabular > {output_dir}/temp/blast_res/{isoform}.blast.tsv', shell=True,stdout=subprocess.PIPE,executable = '/bin/bash')
    completed = subprocess.run(f'conda run -n {conda_env_name} magicblast -query {output_dir}/temp/blast_res/{isoform}.fa -db {output_dir}/temp/NONCODE -reftype transcriptome -num_threads 1 -outfmt tabular > {output_dir}/temp/blast_res/{isoform}.NONCODE_blast.tsv', shell=True,stdout=subprocess.PIPE,executable = '/bin/bash')
    all_problem_design_set = set()
    with open(f'{output_dir}/temp/blast_res/{isoform}.blast.tsv','r') as f:
        f.readline()
        f.readline()
        f.readline()
        for line in f:
            fields = line.split('\t')
            if len(fields) < 2:
                continue
            if not allow_shared_region:
                if fields[1] not in fields[0]:
                    all_problem_design_set.add(fields[0])
            # else:
            #     if fields[1] not in isoform_set:
            #         all_problem_design_set.add(fields[0])
    with open(f'{output_dir}/temp/blast_res/{isoform}.NONCODE_blast.tsv','r') as f:
        f.readline()
        f.readline()
        f.readline()
        for line in f:
            fields = line.split('\t')
            if len(fields) < 2:
                continue
            if fields[1] != '-':
                all_problem_design_set.add(fields[0])
    good_primer_design = {}
    for region,primer_design in all_primer_design.items():
        for primer_id,primer in primer_design.items():
            if f'{region}_{primer_id}' not in all_problem_design_set:
                region_id = region.split('_')[-1]
                if region_id not in good_primer_design:
                    good_primer_design[region_id] = {}
                primer['primer_id'] = primer_id
                primer['region_id'] = region_id
                primer['region'] = region
                good_primer_design[region_id][primer_id] = primer
#     Path(f'{output_dir}/temp/blast_res/{isoform}.blast.tsv').unlink()
#     Path(f'{output_dir}/temp/blast_res/{isoform}.NONCODE_blast.tsv').unlink()
#     Path(f'{output_dir}/temp/blast_res/{isoform}.fa').unlink()
    return good_primer_design
def get_primer_gtf_lines(good_primer_design,region_info,isoform,region_type,gname):
    primer_pos_info = []
    for primer_id in good_primer_design:
        assert good_primer_design[primer_id]['region_id'] == str(region_info['region_id'])
        primer_info = good_primer_design[primer_id]
        # positive strand
        left_span_junction = False
        right_span_junction = False
        assert region_info['strand'] in ['+','-']
        if region_info['strand'] == '+':
            if region_type == 'junction':
                left_primer_start = primer_info['LEFT_start_genome']
                left_primer_end = primer_info['LEFT_end_genome']
                right_primer_start = primer_info['RIGHT_start_genome']
                right_primer_end = primer_info['RIGHT_end_genome']
            # if region_type == 'junction':
            #     left_primer_start = primer_info['LEFT_start'] + region_info['exon_0_start']
            #     exon_0_len = region_info['exon_0_end'] - region_info['exon_0_start'] + 1
            #     # span the splice site
            #     if left_primer_start + primer_info['LEFT_len'] - 1 > region_info['exon_0_end']:
            #         left_span_junction = True
            #         length_on_exon_0 = region_info['exon_0_end'] - left_primer_start + 1
            #         assert length_on_exon_0 > 0
            #         left_primer_end = region_info['exon_1_start'] + (primer_info['LEFT_len'] - length_on_exon_0) - 1
            #     # not span the splice site
            #     else:
            #         left_primer_end = left_primer_start + primer_info['LEFT_len'] - 1
            #     exon_1_len = region_info['exon_1_end'] - region_info['exon_1_start'] + 1
            #     right_primer_end = (primer_info['RIGHT_start'] + 1 - exon_0_len) + region_info['exon_1_start'] - 1
            #     # span the splice site
            #     if right_primer_end - primer_info['RIGHT_len'] + 1 < region_info['exon_1_start']:
            #         right_span_junction = True
            #         length_on_exon_1 = right_primer_end - region_info['exon_1_start'] + 1
            #         assert length_on_exon_1 > 0
            #         right_primer_start = region_info['exon_0_end'] - (primer_info['RIGHT_len'] - length_on_exon_1) + 1
            #     # not span the splice site
            #     else:
            #         right_primer_start = right_primer_end - primer_info['RIGHT_len'] + 1
            elif region_type == 'exon':
                left_primer_start = primer_info['LEFT_start'] + region_info['exon_start']
                left_primer_end = left_primer_start + primer_info['LEFT_len'] - 1
                right_primer_end = primer_info['RIGHT_start'] + region_info['exon_start']
                right_primer_start = right_primer_end - primer_info['RIGHT_len'] + 1
        # negative strand
        else:
            if region_type == 'junction':
                left_primer_start = primer_info['LEFT_start_genome']
                left_primer_end = primer_info['LEFT_end_genome']
                right_primer_start = primer_info['RIGHT_start_genome']
                right_primer_end = primer_info['RIGHT_end_genome']

                # left_primer_end = region_info['exon_1_end'] - primer_info['LEFT_start'] 
                # exon_1_len = region_info['exon_1_end'] - region_info['exon_1_start'] + 1
                # # span the splice site
                # if left_primer_end - primer_info['LEFT_len'] + 1 < region_info['exon_1_start']:
                #     left_span_junction = True
                #     length_on_exon_1 = left_primer_end - region_info['exon_1_start'] + 1
                #     assert length_on_exon_1 > 0
                #     left_primer_start = region_info['exon_0_end'] - (primer_info['LEFT_len'] - length_on_exon_1) + 1
                # # not span the splice site
                # else:
                #     left_primer_start = left_primer_end - primer_info['LEFT_len'] + 1
                # exon_0_len = region_info['exon_0_end'] - region_info['exon_0_start'] + 1
                
                # right_primer_start = region_info['exon_0_end'] - (primer_info['RIGHT_start'] + 1 - exon_1_len) + 1
                # # span the splice site
                # if right_primer_start + primer_info['RIGHT_len'] - 1 > region_info['exon_0_end']:
                #     right_span_junction = True
                #     length_on_exon_0 =  region_info['exon_0_end'] - right_primer_start + 1
                #     if length_on_exon_0 <= 0:
                #         print(region_info)
                #         print(primer_info)
                #     assert length_on_exon_0 > 0
                #     right_primer_end = region_info['exon_1_start'] + (primer_info['RIGHT_len'] - length_on_exon_0) - 1
                # # not span the splice site
                # else:
                #     right_primer_end = right_primer_start + primer_info['RIGHT_len'] - 1
            elif region_type == 'exon':
                left_primer_end = region_info['exon_end'] - primer_info['LEFT_start']
                left_primer_start = left_primer_end - primer_info['LEFT_len'] + 1
                right_primer_start = region_info['exon_end'] - primer_info['RIGHT_start'] 
                right_primer_end = right_primer_start + primer_info['RIGHT_len'] - 1
        if region_info['strand'] == '+':
            reverse_strand = '-'
        else:
            reverse_strand = '+'
        if region_type == 'junction':
            primer_pos_info.append([region_info['chr'],'Primer3','exon',str(left_primer_start),str(left_primer_end),'.',region_info['strand'],'.',f'primer_id "{primer_id}"; direction "LEFT"; target_type "JUNCTION"; gene_id "{gname}"; transcript_id "{isoform}_{primer_id};'])
            primer_pos_info.append([region_info['chr'],'Primer3','exon',str(right_primer_start),str(right_primer_end),'.',reverse_strand,'.',f'primer_id "{primer_id}"; direction "RIGHT"; target_type "JUNCTION"; gene_id "{gname}"; transcript_id "{isoform}_{primer_id};'])
            # if left_span_junction:
            #     assert left_primer_start <= region_info['exon_0_end']
            #     assert region_info['exon_1_start'] <= left_primer_end
            #     primer_pos_info.append([region_info['chr'],'Primer3','exon',str(left_primer_start),str(region_info['exon_0_end']),'.',region_info['strand'],'.',f'primer_id "{primer_id}"; direction "LEFT"; target_type "JUNCTION"; gene_id "{gname}"; transcript_id "{isoform}_{primer_id};'])
            #     primer_pos_info.append([region_info['chr'],'Primer3','exon',str(region_info['exon_1_start']),str(left_primer_end),'.',region_info['strand'],'.',f'primer_id "{primer_id}"; direction "LEFT"; target_type "JUNCTION" gene_id "{gname}"; transcript_id "{isoform}_{primer_id};'])
            # else:
            #     assert left_primer_start <= left_primer_end
            #     primer_pos_info.append([region_info['chr'],'Primer3','exon',str(left_primer_start),str(left_primer_end),'.',region_info['strand'],'.',f'primer_id "{primer_id}"; direction "LEFT"; target_type "JUNCTION"; gene_id "{gname}"; transcript_id "{isoform}_{primer_id};'])
            # if right_span_junction:
            #     assert right_primer_start <= region_info['exon_0_end']
            #     assert region_info['exon_1_start'] <= right_primer_end
            #     primer_pos_info.append([region_info['chr'],'Primer3','exon',str(right_primer_start),str(region_info['exon_0_end']),'.',reverse_strand,'.',f'primer_id "{primer_id}"; direction "RIGHT"; target_type "JUNCTION"; gene_id "{gname}"; transcript_id "{isoform}_{primer_id};'])
            #     primer_pos_info.append([region_info['chr'],'Primer3','exon',str(region_info['exon_1_start']),str(right_primer_end),'.',reverse_strand,'.',f'primer_id "{primer_id}"; direction "RIGHT"; target_type "JUNCTION"; gene_id "{gname}"; transcript_id "{isoform}_{primer_id};'])
            # else:
            #     assert right_primer_start <= right_primer_end
            #     primer_pos_info.append([region_info['chr'],'Primer3','exon',str(right_primer_start),str(right_primer_end),'.',reverse_strand,'.',f'primer_id "{primer_id}"; direction "RIGHT"; target_type "JUNCTION"; gene_id "{gname}"; transcript_id "{isoform}_{primer_id};'])
        elif region_type == 'exon':
            assert right_primer_end - right_primer_start + 1 == primer_info['RIGHT_len']
            assert left_primer_end - left_primer_start + 1 == primer_info['LEFT_len']
            primer_pos_info.append([region_info['chr'],'Primer3','exon',str(left_primer_start),str(left_primer_end),'.',region_info['strand'],'.',f'primer_id "{primer_id}"; direction "LEFT"; target_type "EXON"; gene_id "{gname}"; transcript_id "{isoform}_{primer_id};'])
            primer_pos_info.append([region_info['chr'],'Primer3','exon',str(right_primer_start),str(right_primer_end),'.',reverse_strand,'.',f'primer_id "{primer_id}"; direction "RIGHT"; target_type "EXON"; gene_id "{gname}"; transcript_id "{isoform}_{primer_id};'])
        
    return ['\t'.join(fields) for fields in primer_pos_info]
def check_primer_3_and_blast_single_thread(worker_id,output_dir,conda_env_name):
    if not Path(f'{output_dir}/temp/target_sequences/{worker_id}').exists():
        return None
    with open(f'{output_dir}/temp/target_sequences/{worker_id}','rb') as f:
        all_isoform_target_sequence_dict,all_isoform_gname_dict = pickle.load(f)
    all_good_primer_design = {}
    all_gtf_lines = {}
    for isoform,isoform_target_sequence_dict in all_isoform_target_sequence_dict.items():
        # use primer 3 to design primer
        all_primer_design,all_raw_primer_design_result = design_primer_for_isoform(isoform_target_sequence_dict)
        # check by blast
        good_primer_design = check_for_blast(all_primer_design,isoform,output_dir,conda_env_name,None,False)
        if len(good_primer_design) == 0 :
            continue 
        all_good_primer_design[isoform] = good_primer_design
        for region_id,region_dict in good_primer_design.items():
            # bad practice. search the isoform_target_sequence_dict to find the region
            region_info = None
            region_type = None
            for region in isoform_target_sequence_dict['junction_info']:
                if str(region['region_id'])  == str(region_id):
                    region_info = region
                    region_type = 'junction'
                    break
            if region_info is None and region_type is None:
                for region in isoform_target_sequence_dict['exon_info']:
                    if str(region['region_id'])  == str(region_id):
                        region_info = region
                        region_type = 'exon'
                        break
            gname = all_isoform_gname_dict[isoform]
            region_gtf_lines = get_primer_gtf_lines(good_primer_design[region_id],region_info,isoform,region_type,gname)
            if isoform not in all_gtf_lines:
                all_gtf_lines[isoform] = {}
            if f'{region_type}_region_{region_id}' not in all_gtf_lines[isoform]:
                all_gtf_lines[isoform][f'{region_type}_region_{region_id}'] = {}
            all_gtf_lines[f'{isoform}'][f'{region_type}_region_{region_id}'] = region_gtf_lines
    for isoform in all_gtf_lines:
        for region in all_gtf_lines[isoform]:
            Path(f'{output_dir}/designs/{isoform}/').mkdir(parents=True,exist_ok=True)
            with open(f'{output_dir}/designs/{isoform}/{region}.gtf','w') as f:
                for line in all_gtf_lines[isoform][region]:
                    f.write(line)
                    f.write('\n')
    list_of_gene_primer_design_df = []
    for isoform in all_good_primer_design:
        all_excel_lines = []
        for region in all_good_primer_design[isoform]:
            for primer_id in all_good_primer_design[isoform][region]:
                region_name = all_good_primer_design[isoform][region][primer_id]['region']
                fields = ['Isoform','SEQUENCE_ID','Primer_id']
                line = [isoform,region_name,primer_id]
                for key in all_good_primer_design[isoform][region][primer_id]:
                    if key in ['primer_id','region_id','region']:
                        continue
                    fields.append(key)
                    line.append(all_good_primer_design[isoform][region][primer_id][key])
                all_excel_lines.append(line)
        gene_primer_design_df = pd.DataFrame(all_excel_lines)
        gene_primer_design_df.columns = fields
        gene_primer_design_df.to_csv(f'{output_dir}/designs/{isoform}/primer_overview.tsv',sep='\t',index=False)
        list_of_gene_primer_design_df.append(gene_primer_design_df)
    return list_of_gene_primer_design_df
def callback_error(result):
    print('ERR:', result,flush=True)
def check_primer_3_and_blast(ref_file_path,reference_genome_path,output_dir,conda_env_name,NONCODE_fasta_url,threads):
    build_blast_library(ref_file_path,reference_genome_path,output_dir,conda_env_name,NONCODE_fasta_url)
    pool = mp.Pool(threads)
    futures = []
    for i in range(threads):
        futures.append(pool.apply_async(check_primer_3_and_blast_single_thread,(i,output_dir,conda_env_name,),error_callback=callback_error))

    list_of_all_gene_primer_design_df = []
    for future in futures:
        list_of_gene_primer_design_df = future.get()
        if list_of_gene_primer_design_df is not None:
            list_of_all_gene_primer_design_df += list_of_gene_primer_design_df
    pool.close()
    pool.join()
    if len(list_of_all_gene_primer_design_df) != 0:
        all_gene_primer_design_df = pd.concat(list_of_all_gene_primer_design_df)
        return all_gene_primer_design_df
    else:
        return None
def design_primers(ref_file_path,reference_genome_path,output_dir,conda_env_name,NONCODE_fasta_url,threads):

    all_gene_primer_design_df = check_primer_3_and_blast(ref_file_path,reference_genome_path,output_dir,conda_env_name,NONCODE_fasta_url,threads)
    unique_region_df = pd.read_csv(f'{output_dir}/unique_region.tsv',sep='\t').set_index('SEQUENCE_ID')
    # unique_region_df['exon_1_start'] = unique_region_df['exon_1_start'].astype('Int64')
    # unique_region_df['exon_1_end'] = unique_region_df['exon_1_end'].astype('Int64')
    if all_gene_primer_design_df is not None:
        unique_region_with_primers_df = unique_region_df.loc[set(all_gene_primer_design_df['SEQUENCE_ID'])]
        unique_region_with_primers_df.to_csv(f'{output_dir}/unique_region_with_primers.tsv',sep='\t',na_rep='nan')
        all_gene_primer_design_df.to_csv(f'{output_dir}/all_primers_overview.tsv',sep='\t',index=False)
