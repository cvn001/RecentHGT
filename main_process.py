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


if __name__ == '__main__':
    newParser = argparse.ArgumentParser()
    newParser.add_argument("-n", "--threads", type=int, dest="threads", help="How many threads will be used?")
    newParser.add_argument("-p", "--part", type=int, dest="part", help="Which part?")
    print('----------------------Start----------------------')
    program_path = os.getcwd()
    args = newParser.parse_args()
    part = args.part
    if part == 1:
        strain_pair_dir = os.path.join(program_path, 'all_strain_pairs')
        output_dir = os.path.join(program_path, 'all_strain_pairs_genes_alignment')
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        strain_label_file = os.path.join(program_path, 'strains.txt')
        processes = args.threads
        p = Pool(processes)
        all_strain_dict = load_strains_label(strain_label_file)
        strain_results_collection_dir = os.path.join(program_path, 'all_strain_pairs_results_collection')
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
    elif part == 2:
        strain_ANI_matrix = os.path.join(program_path, 'ANI_matrix.txt')
        ani_matrix = np.loadtxt(strain_ANI_matrix, delimiter='\t')
        matrix_strain_id_file = os.path.join(program_path, 'matrix_strain_id.txt')
        matrix_strain_dict = defaultdict()
        strain_name_dict = defaultdict()
        with open(matrix_strain_id_file, 'r') as mx:
            for each_line in mx.readlines():
                mx_list = each_line.strip().split('\t')
                s_list = mx_list[1].split('_')
                matrix_strain_dict[s_list[1]] = int(mx_list[0])
                strain_name_dict[s_list[1]] = mx_list[1]
        strain_results_collection_dir = os.path.join(program_path, 'all_strain_pairs_results_collection')
        all_result_dir = os.path.join(program_path, 'all_strain_pairs_results')
        if not os.path.exists(all_result_dir):
            os.makedirs(all_result_dir)
        strain_pair_pictures_dir = os.path.join(program_path, 'all_strain_pairs_pictures')
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
                pair_name = strain_name_dict[str(pairs[0])] + '~' + strain_name_dict[str(pairs[1])]
                pair_ANI = ani_matrix[matrix_strain_dict[pairs[0]]][matrix_strain_dict[pairs[1]]]
                if pair_ANI < 85:
                    with open(file_path) as f0:
                        for line in f0.readlines()[1:]:
                            result_line = '{0} ({1})\t{2}'.format(pair_name, pair_ANI, line)
                            result_dict[pairs[0]].append(result_line)
        R_script_1 = os.path.join(program_path, 'draw_identity_histograms.R')
        header_line = 'Pair\tGene\tIdentity\tAnnotation\n'
        for each_strain, results in result_dict.items():
            strain_result_file = os.path.join(all_result_dir, each_strain + '.txt')
            with open(strain_result_file, 'w') as f1:
                f1.write(header_line)
                for each_result in results:
                    f1.write(each_result)
            each_result_picture = os.path.join(strain_pair_pictures_dir, '{0}_results.pdf'.format(each_strain))
            R_cmd_1 = "Rscript {0} {1} {2}".format(R_script_1, strain_result_file, each_result_picture)
            os.popen(R_cmd_1)
    elif part == 3:
        print('Processing HGT detection...')
        results_collection_dir_name = 'all_strain_pairs_results_collection'
        R_script_2 = os.path.join(program_path, 'rHGT_alpha.R')
        param_min = 50.0
        param_max = 98.5
        R_cmd_2 = "Rscript {0} {1} {2} {3} {4}".format(R_script_2, program_path,
                                                       results_collection_dir_name,
                                                       str(param_min), str(param_max))
        os.popen(R_cmd_2)
    print('----------------------Finish----------------------')
