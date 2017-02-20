#!/usr/bin/python
# -*- coding: UTF-8 -*-
# Introduction: This script is used to prepare and run rHGT detection automatically
# Created by galaxy on 2016/10/19 11:19

import os
import re
import shutil
import argparse
import numpy as np
from Bio.Emboss.Applications import NeedleCommandline
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool
from collections import defaultdict


def load_strains_label(strain_file):
    """
    This function is used to load strains information.
    :param strain_file: the path of the file contains the ids and names of strains.
    :return: a python dict
    """
    strain_dict = defaultdict()
    strain_list = []
    with open(strain_file, 'r') as f:
        for a_line in f.readlines():
            a_list = a_line.strip().split('\t')
            strain_name = a_list[1].split(' ')[2]
            strain_id = a_list[2]
            strain_dict[strain_id] = strain_name
            strain_list.append(strain_name)
    return strain_dict


def each_gene_needle_run(pair_gene_dir, tmp_gene_converted_dir, pair_gene_alignment_dir, gene, strain_dict):
    """
    This function is used to call needle program to do pair-wised sequence alignment
    :param pair_gene_dir: each homologous gene directory
    :param tmp_gene_converted_dir: used to put some temporary files and will be deleted in the end
    :param pair_gene_alignment_dir: each homologous gene pair-wised alignment directory
    :param gene: each gene id
    :param strain_dict: inherit from load_strains_label function with strain information
    :return: the alignment result of each gene
    """
    tmp_gene_fasta = os.path.join(pair_gene_dir, gene + '.fasta')
    converted_records = []
    re_pattern = re.compile(r'fig\|(\d+\.\d+)\.peg\.\d+\s(.*)')
    in_pattern = re.compile(r'Identity.*\((\d+\.\d+)%\)')
    annotation = ''
    for record in SeqIO.parse(tmp_gene_fasta, 'fasta'):
        m = re.search(re_pattern, record.description)
        strain_id = m.group(1)
        annotation = m.group(2)
        record.id = strain_dict[strain_id]
        final_record = SeqRecord(record.seq, record.id, description='')
        converted_records.append(final_record)
    the_strain_fasta = os.path.join(tmp_gene_converted_dir, 'a.fasta')
    other_strain_fasta = os.path.join(tmp_gene_converted_dir, 'b.fasta')
    SeqIO.write(converted_records[0], the_strain_fasta, 'fasta')
    SeqIO.write(converted_records[1], other_strain_fasta, 'fasta')
    result_file = os.path.join(pair_gene_alignment_dir, "{0}.txt".format(gene))
    needle_cline = NeedleCommandline()
    needle_cline.asequence = the_strain_fasta
    needle_cline.bsequence = other_strain_fasta
    needle_cline.gapopen = 10
    needle_cline.gapextend = 0.5
    needle_cline.outfile = result_file
    os.popen(str(needle_cline))
    os.remove(the_strain_fasta)
    os.remove(other_strain_fasta)
    gene_alignment_result = ''
    with open(result_file, 'r') as f:
        for a_line in f.readlines():
            if 'Identity' in a_line:
                m = re.search(in_pattern, a_line.strip())
                identity = m.group(1)
                gene_alignment_result = '{0}\t{1}\t{2}\n'.format(gene, identity, annotation)
    return gene_alignment_result


def each_strain_pair_run(strain_pair, all_genes_dir, result_dir, strain_dict, strain_results_dir):
    """
    This function is used to call the functions above to do each strain pair genes alignment.
    :param strain_pair: each strain pair name, such as A_B or B_A
    :param all_genes_dir: a directory with all genes of each strain pair
    :param result_dir: result directory
    :param strain_dict: inherit from load_strains_label function with strain information
    :param strain_results_dir: a directory with all strain results
    :return: ......
    """
    gene_list = []
    strain_pair_genes_dir = os.path.join(all_genes_dir, strain_pair)
    for rs, ds, fs in os.walk(strain_pair_genes_dir):
        for f in fs:
            fname = os.path.splitext(f)[0]
            gene_list.append(fname)
    strain_pair_result_dir = os.path.join(result_dir, strain_pair)
    if os.path.exists(strain_pair_result_dir):
        shutil.rmtree(strain_pair_result_dir)
    os.makedirs(strain_pair_result_dir)
    tmp_gene_converted_dir = os.path.join(strain_pair_result_dir, 'tmp')
    if os.path.exists(tmp_gene_converted_dir):
        shutil.rmtree(tmp_gene_converted_dir)
    os.makedirs(tmp_gene_converted_dir)
    pair_results = 'Gene\tIdentity\tAnnotation\n'
    for gene in gene_list:
        gene_alignment_result = each_gene_needle_run(strain_pair_genes_dir, tmp_gene_converted_dir,
                                                     strain_pair_result_dir, gene, strain_dict)
        pair_results += gene_alignment_result
    shutil.rmtree(tmp_gene_converted_dir)
    pair_result_collection_file = os.path.join(strain_results_dir, '{0}.txt'.format(strain_pair))
    with open(pair_result_collection_file, 'w') as f:
        f.write(pair_results)


def first_part(base_path):
    gbk_folder = os.path.join(base_path, 'gbk')
    fasta_folder = os.path.join(base_path, 'genome_fasta')
    if not os.path.exists(fasta_folder):
        os.makedirs(fasta_folder)
    i = 0
    for root, dirs, files in os.walk(gbk_folder):
        for each_file in files:
            strain_id = os.path.splitext(each_file)[0]
            gbk = os.path.join(gbk_folder, each_file)
            fasta = os.path.join(fasta_folder, strain_id + '.fasta')
            SeqIO.convert(gbk, 'genbank', fasta, 'fasta')
            i += 1
    print('{0} genomes fasta format files have been converted.'.format(str(i)))
    ani_folder = os.path.join(base_path, 'ANI_out')
    cmd = "average_nucleotide_identity.py -i {0} -o {1} -m ANIm -g".format(fasta_folder, ani_folder)
    os.popen(cmd)
    return


def second_part(base_path, processes):
    strain_pair_dir = os.path.join(base_path, 'all_strain_pairs')
    output_dir = os.path.join(base_path, 'all_strain_pairs_genes_alignment')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    strain_label_file = os.path.join(base_path, 'strains.txt')
    p = Pool(processes)
    all_strain_dict = load_strains_label(strain_label_file)
    strain_results_collection_dir = os.path.join(base_path, 'all_strain_pairs_results_collection')
    if not os.path.exists(strain_results_collection_dir):
        os.makedirs(strain_results_collection_dir)
    pairs = []
    for root, dirs, files in os.walk(strain_pair_dir):
        for each_dir in dirs:
            pairs.append(each_dir)
    for each_pair in pairs:
        p.apply_async(each_strain_pair_run, args=(each_pair, strain_pair_dir, output_dir,
                                                  all_strain_dict, strain_results_collection_dir))
    p.close()
    p.join()
    message = 'All alignments have been finished.'
    return message


def third_part(base_path, ani_out_dir):
    strain_ani_matrix_file = os.path.join(ani_out_dir, 'ANIm_percentage_identity.tab')
    matrix_strain_dict = defaultdict()
    tmp_matrix_lines = ''
    with open(strain_ani_matrix_file, 'r') as f1:
        strain_list = f1.readlines()[0].strip().split('\t')[1:]
        for each_strain in strain_list:
            strain_index = strain_list.index(each_strain)
            matrix_strain_dict[each_strain] = int(strain_index)
        for each_line in f1.readlines()[1:]:
            tmp_list = each_line.strip().split('\t')
            tmp_line = ''
            for i in tmp_list[1:]:
                i = float('%.3f' % i) * 100
                tmp_line += '{0}\t'.format(str(i))
            tmp_matrix_lines += tmp_line.strip('\t') + '\n'
    tmp_matrix_file = os.path.join(base_path, 'tmp.txt')
    with open(tmp_matrix_file, 'w') as f2:
        f2.write(tmp_matrix_lines)
    ani_matrix = np.loadtxt(tmp_matrix_file, delimiter='\t')
    strain_results_collection_dir = os.path.join(base_path, 'all_strain_pairs_results_collection')
    all_result_dir = os.path.join(base_path, 'all_strain_pairs_results')
    if not os.path.exists(all_result_dir):
        os.makedirs(all_result_dir)
    strain_pair_pictures_dir = os.path.join(base_path, 'all_strain_pairs_pictures')
    if not os.path.exists(strain_pair_pictures_dir):
        os.makedirs(strain_pair_pictures_dir)
    result_dict = defaultdict(list)
    for root, dirs, files in os.walk(strain_results_collection_dir):
        for each_file in files:
            file_name = os.path.splitext(each_file)[0]
            file_path = os.path.join(strain_results_collection_dir, each_file)
            pairs = file_name.split('_')
            if pairs[0] not in result_dict:
                result_dict[pairs[0]] = []
            pair_name = str(pairs[0]) + '~' + str(pairs[1])
            pair_ani = ani_matrix[matrix_strain_dict[pairs[0]]][matrix_strain_dict[pairs[1]]]
            if pair_ani < 95:
                with open(file_path) as f2:
                    for line in f2.readlines()[1:]:
                        result_line = '{0} ({1})\t{2}'.format(pair_name, pair_ani, line)
                        result_dict[pairs[0]].append(result_line)
    r_script_1 = os.path.join(base_path, 'draw_identity_histograms.R')
    header_line = 'Pair\tGene\tIdentity\tAnnotation\n'
    for each_strain, results in result_dict.items():
        strain_result_file = os.path.join(all_result_dir, each_strain + '.txt')
        with open(strain_result_file, 'w') as f3:
            f3.write(header_line)
            for each_result in results:
                f3.write(each_result)
        each_result_picture = os.path.join(strain_pair_pictures_dir, '{0}_results.pdf'.format(each_strain))
        r_cmd_1 = "Rscript {0} {1} {2}".format(r_script_1, strain_result_file, each_result_picture)
        os.popen(r_cmd_1)
    os.remove(tmp_matrix_file)
    message = 'All alignments pictures have been drew and placed in {0}'.format(strain_pair_pictures_dir)
    return message


def fourth_part(base_path):
    print('Processing HGT detection...')
    results_collection_dir_name = 'all_strain_pairs_results_collection'
    r_script_2 = os.path.join(base_path, 'rHGT_alpha.R')
    param_min = 50.0
    param_max = 98.5
    r_cmd_2 = "Rscript {0} {1} {2} {3} {4}".format(r_script_2, base_path, results_collection_dir_name,
                                                   str(param_min), str(param_max))
    os.popen(r_cmd_2)


def auto_run():
    return


def main_process(part):
    if part == 1:
        pass
    elif part == 2:
        pass
    elif part == 3:
        pass
    print('----------------------Finish----------------------')


if __name__ == '__main__':
    newParser = argparse.ArgumentParser()
    newParser.add_argument("-n", type=int, dest="threads", help="How many threads will be used?")
    newParser.add_argument("-p", type=int, dest="part", help="Which part will be run?")
    print('----------------------Start----------------------')
    program_path = os.getcwd()
    args = newParser.parse_args()
    run_part = args.part


