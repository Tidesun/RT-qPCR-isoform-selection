import bisect
import re
import pickle

from Bio.Seq import Seq
import pysam
from intervaltree import IntervalTree
import numpy as np
import pandas as pd

from parse_annotation import parse_annotation
def get_sequence(start,end,chr,strand,reference_fasta):
    sequence = Seq(reference_fasta.fetch(chr,start - 1,end)).upper()
    if strand == -1:
        sequence = sequence.reverse_complement()
    return str(sequence)
def get_gene_interval_tree_dict(gene_points_dict,gene_exons_dict):
    gene_interval_tree_dict = dict()
    for chr_name in gene_points_dict:
        gene_interval_tree_dict[chr_name] = IntervalTree()
        for gene_name in gene_exons_dict[chr_name]:
            for [start_pos,end_pos,_] in gene_exons_dict[chr_name][gene_name]:
                # Interval tree exclude end position
                gene_interval_tree_dict[chr_name].addi(start_pos, end_pos+1, gene_name)
    return gene_interval_tree_dict
def get_gene_strand_dict(ref_file_path):
    gene_strand_dict = {}
    with open(ref_file_path,'r') as f:
        for line in f:
            if line.lstrip()[0] == "#":
                continue
            fields = line.split('\t')
            gene_id = re.findall('gene_id "([^"]*)"', fields[8])[0]
            if gene_id not in gene_strand_dict:
                gene_strand_dict[gene_id] = fields[6]
    return gene_strand_dict
def get_isoform_exons_dict(raw_isoform_exons_dict):
    isoform_exons_dict = {}
    for rname in raw_isoform_exons_dict:
        for gname in raw_isoform_exons_dict[rname]:
            for isoform in raw_isoform_exons_dict[rname][gname]:
                isoform_exons_dict[isoform] = raw_isoform_exons_dict[rname][gname][isoform]
    return isoform_exons_dict
def get_target_sequence(rname,gname,chr_name,reference_fasta,lower_region_length,gene_regions_dict,gene_points_dict,gene_isoforms_dict,gene_interval_tree_dict,gene_strand_dict,shared_region_isoform=None,unique_region=True,shared_all_region=False):
    point_dict = {}
    for pos,point in gene_points_dict[rname][gname].items():
        point_dict[f'P{point}'] = pos
    all_isoforms = gene_isoforms_dict[rname][gname]
    all_isoform_target_sequence_dict = {}
    region_id = 0

    for (region,isoform_set) in gene_regions_dict[rname][gname].items():
        target_sequence = None
        if unique_region:
            if len(isoform_set) != 1:
                continue
            isoform = list(isoform_set)[0]
        elif not shared_all_region:
            if shared_region_isoform not in isoform_set:
                continue
            isoform = shared_region_isoform
        elif shared_all_region:
            if len(all_isoforms) != len(isoform_set):
                continue
            isoform = shared_region_isoform
        if ':' in region and '-' in region:
            if region.count(':') > 2:
                continue
            exons = region.split('-')
            [start_pt_0,end_pt_0] = exons[0].split(':')
            [start_pt_1,end_pt_1] = exons[1].split(':')
            start_pos_0,end_pos_0 = point_dict[start_pt_0],point_dict[end_pt_0]
            start_pos_1,end_pos_1 = point_dict[start_pt_1],point_dict[end_pt_1]
            if (end_pos_1 - start_pos_1 + 1) + (end_pos_0 - start_pos_0 + 1)  >= lower_region_length:
                if len(gene_interval_tree_dict[rname].overlap(start_pos_0 + 1,end_pos_0)) == 1:
                    if len(gene_interval_tree_dict[rname].overlap(start_pos_1 + 1,end_pos_1)) == 1:
                    
                        seq_dict = {}
                        region_id += 1
                        seq_dict['SEQUENCE_ID'] = f'{gname}_{isoform}_{region_id}'
                        target_sequence_0 = get_sequence(start_pos_0,end_pos_0,chr_name,gene_strand_dict[gname],reference_fasta)
                        target_sequence_1 = get_sequence(start_pos_1,end_pos_1,chr_name,gene_strand_dict[gname],reference_fasta)
                        if gene_strand_dict[gname] == 1:
                            seq_dict['SEQUENCE_TEMPLATE'] = target_sequence_0 + target_sequence_1
                            seq_dict['SEQUENCE_OVERLAP_JUNCTION_LIST'] = len(target_sequence_0)
                            target_sequence = target_sequence_0 +'^'+target_sequence_1
                        else:
                            seq_dict['SEQUENCE_TEMPLATE'] =  target_sequence_1 + target_sequence_0
                            seq_dict['SEQUENCE_OVERLAP_JUNCTION_LIST'] = len(target_sequence_1)
                            target_sequence = target_sequence_1 +'^'+target_sequence_0
                        if isoform not in all_isoform_target_sequence_dict:
                            all_isoform_target_sequence_dict[isoform] = {'junction':[],'junction_info':[],'exon':[],'exon_info':[]}
                        all_isoform_target_sequence_dict[isoform]['junction'].append(seq_dict)
                        strand = '+' if gene_strand_dict[gname] == 1 else '-'
                        all_isoform_target_sequence_dict[isoform]['junction_info'].append(\
                        {'region_id':region_id,'seq':target_sequence,'chr':chr_name,'strand':strand,'exon_0_start':start_pos_0,'exon_0_end':end_pos_0,'exon_1_start':start_pos_1,'exon_1_end':end_pos_1})
        elif ':' in region and not '-' in region:
            [start_pt,end_pt] = region.split(':')
            start_pos,end_pos = point_dict[start_pt],point_dict[end_pt]
            if end_pos - start_pos + 1 >= lower_region_length:
                if len(gene_interval_tree_dict[rname].overlap(start_pos + 1,end_pos)) == 1:
                    seq_dict = {}
                    region_id += 1
                    seq_dict['SEQUENCE_ID'] = f'{gname}_{isoform}_{region_id}'
                    target_sequence = get_sequence(start_pos,end_pos,chr_name,gene_strand_dict[gname],reference_fasta)
                    seq_dict['SEQUENCE_TEMPLATE'] = target_sequence
                    if isoform not in all_isoform_target_sequence_dict:
                        all_isoform_target_sequence_dict[isoform] = {'junction':[],'junction_info':[],'exon':[],'exon_info':[]}
                    all_isoform_target_sequence_dict[isoform]['exon'].append(seq_dict)
                    strand = '+' if gene_strand_dict[gname] == 1 else '-'
                    all_isoform_target_sequence_dict[isoform]['exon_info'].append(\
                        {'region_id':region_id,'seq':target_sequence,'chr':chr_name,'strand':strand,'exon_start':start_pos,'exon_end':end_pos})
    return all_isoform_target_sequence_dict 
import numpy as np
def get_distance_to_two_ends(region_start,region_end,isoform,isoform_exons_dict):
    start_index = bisect.bisect_left(isoform_exons_dict[isoform]['start_pos'],region_start)
    end_index = bisect.bisect_right(isoform_exons_dict[isoform]['end_pos'],region_end)
    start_distance = 0
    end_distance = 0
    for exon_index in range(0,start_index):
        exon_start = isoform_exons_dict[isoform]['start_pos'][exon_index]
        exon_end =  isoform_exons_dict[isoform]['end_pos'][exon_index]
        if exon_end < region_start:
            start_distance += exon_end - exon_start + 1
        else:
            start_distance += region_start - exon_start + 1
    for exon_index in range(end_index,len(isoform_exons_dict[isoform]['start_pos'])):
        exon_start = isoform_exons_dict[isoform]['start_pos'][exon_index]
        exon_end =  isoform_exons_dict[isoform]['end_pos'][exon_index]
        if exon_start > region_end:
            end_distance += exon_end - exon_start + 1
        else:
            end_distance += exon_end - region_end  + 1
    return np.abs(start_distance-end_distance)
def get_distance_to_five_and_three_ends(all_isoform_target_sequence_dict,isoform_exons_dict):
    '''
    We want to first use the region that are far from 5' and 3' ends (region located in the center of isoform).
    '''
    all_regions_diff_dict = {}
    for isoform,isoform_target_sequence_dict in all_isoform_target_sequence_dict.items():
        all_regions_diff_dict[isoform] = {}
        for region_dict in isoform_target_sequence_dict['junction_info']:
            region_diff = get_distance_to_two_ends(region_dict['exon_0_end'],region_dict['exon_1_start'],isoform,isoform_exons_dict)
#             all_regions_diff.append([str(region_dict['region_id']),region_diff,'junction'])
            all_regions_diff_dict[isoform][str(region_dict['region_id'])] = region_diff
        for region_dict in isoform_target_sequence_dict['exon_info']:
            region_diff = get_distance_to_two_ends(region_dict['exon_start'],region_dict['exon_end'],isoform,isoform_exons_dict)
            all_regions_diff_dict[isoform][str(region_dict['region_id'])] = region_diff
    return all_regions_diff_dict
def output_target_sequence_info(all_isoform_target_sequence_dict,isoform_exons_dict,output_dir):
    all_regions_diff_dict = get_distance_to_five_and_three_ends(all_isoform_target_sequence_dict,isoform_exons_dict)
    rows = []
    for isoform in all_isoform_target_sequence_dict:
        for design_dict,info_dict in zip(all_isoform_target_sequence_dict[isoform]['exon'],all_isoform_target_sequence_dict[isoform]['exon_info']):
            region_id = design_dict['SEQUENCE_ID'].split('_')[-1]
            row = [design_dict['SEQUENCE_ID'],design_dict['SEQUENCE_TEMPLATE'],'exon',info_dict['chr'],info_dict['strand'],all_regions_diff_dict[isoform][region_id],info_dict['exon_start'],info_dict['exon_end']]
            rows.append(row)
        for design_dict,info_dict in zip(all_isoform_target_sequence_dict[isoform]['junction'],all_isoform_target_sequence_dict[isoform]['junction_info']):
            region_id = design_dict['SEQUENCE_ID'].split('_')[-1]
            row = [design_dict['SEQUENCE_ID'],design_dict['SEQUENCE_TEMPLATE'],'junction',info_dict['chr'],info_dict['strand'],all_regions_diff_dict[isoform][region_id],info_dict['exon_0_start'],info_dict['exon_0_end'],info_dict['exon_1_start'],info_dict['exon_1_end']]
            rows.append(row)
    df = pd.DataFrame(rows)
    df.columns=['SEQUENCE_ID','SEQUENCE_TEMPLATE','region_type','Chr','Strand','5_end_and_3_end_distance_difference','exon_0_start','exon_0_end','exon_1_start','exon_1_end']
    df['exon_1_start'] = df['exon_1_start'].astype('Int64')
    df['exon_1_end'] = df['exon_1_end'].astype('Int64')
    df = df.set_index('SEQUENCE_ID')
    df.to_csv(f'{output_dir}/unique_region.tsv',sep='\t',na_rep='nan')
    return df
def read_isoform_list(isoform_list_fpath,all_isoform_target_sequence_dict):
    isoform_set = set()
    with open(isoform_list_fpath,'r') as f:
        for line in f:
            isoform_id = line.strip()
            if isoform_id in all_isoform_target_sequence_dict:
                isoform_set.add(isoform_id)
    print(str(len(isoform_set))+' isoforms identified in the isoform list file and have unique region!',flush=True)
    return isoform_set
def prepare_target_sequence_dict(ref_file_path,reference_genome_path, lower_region_length,output_dir,isoform_list_fpath,threads):
    reference_fasta = pysam.FastaFile(reference_genome_path)
    [gene_exons_dict, gene_points_dict, gene_isoforms_dict, genes_regions_len_dict,
        _, gene_regions_dict, gene_isoforms_length_dict,raw_isoform_exons_dict,raw_gene_exons_dict,same_structure_isoform_dict] = parse_annotation(ref_file_path, int((threads+1)//2))
    gene_strand_dict = get_gene_strand_dict(ref_file_path)
    gene_interval_tree_dict = get_gene_interval_tree_dict(gene_points_dict,gene_exons_dict)
    # get all isoform target sequence
    all_isoform_target_sequence_dict = {}
    all_isoform_gname_dict = {}
    for rname in gene_regions_dict:
        chr_name = rname
        if rname not in reference_fasta.references:
            if 'chr'+rname in reference_fasta.references:
                chr_name = 'chr'+rname
            else:
                continue
        for gname in gene_regions_dict[rname]:
            if gname in gene_strand_dict:
                isoform_target_sequence_dict = get_target_sequence(rname,gname,chr_name,reference_fasta,lower_region_length,gene_regions_dict,gene_points_dict,gene_isoforms_dict,gene_interval_tree_dict,gene_strand_dict)
                all_isoform_target_sequence_dict.update(isoform_target_sequence_dict)
                for isoform in all_isoform_target_sequence_dict:
                    all_isoform_gname_dict[isoform] = gname
    isoform_exons_dict = get_isoform_exons_dict(raw_isoform_exons_dict)
    target_sequence_info_df = output_target_sequence_info(all_isoform_target_sequence_dict,isoform_exons_dict,output_dir)
    print(str(len(all_isoform_target_sequence_dict))+' isoforms have unique region!',flush=True)
    if isoform_list_fpath is not None:
        selected_isoform_set = read_isoform_list(isoform_list_fpath,all_isoform_target_sequence_dict)
    # debug
    # backup = all_isoform_target_sequence_dict.copy()
    # all_isoform_target_sequence_dict = {}
    # for isoform in backup:
    #     all_isoform_target_sequence_dict[isoform] = backup[isoform]
    #     if len(all_isoform_target_sequence_dict) > 1000:
    #         break
    # split the isoform_target_sequence_dict by threads
    chunksize, extra = divmod(len(selected_isoform_set), threads)
    if extra:
        chunksize += 1
    print(f'{chunksize} isoforms assigned to each thread for primer design.',flush=True)
    worker_id = 0
    all_isoform_target_sequence_dict_single_thread = {}
    for isoform in all_isoform_target_sequence_dict:
        if isoform_list_fpath is not None:
            if isoform not in selected_isoform_set:
                continue
        all_isoform_target_sequence_dict_single_thread[isoform] = all_isoform_target_sequence_dict[isoform]
        if len(all_isoform_target_sequence_dict_single_thread) >= chunksize:
            with open(f'{output_dir}/temp/target_sequences/{worker_id}','wb') as f:
                pickle.dump((all_isoform_target_sequence_dict_single_thread,all_isoform_gname_dict),f)
            all_isoform_target_sequence_dict_single_thread = {}
            worker_id += 1
    if len(all_isoform_target_sequence_dict_single_thread) >= 0:
        with open(f'{output_dir}/temp/target_sequences/{worker_id}','wb') as f:
            pickle.dump((all_isoform_target_sequence_dict_single_thread,all_isoform_gname_dict),f)
        all_isoform_target_sequence_dict_single_thread = {}
        worker_id += 1
