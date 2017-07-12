#!/usr/bin/python
# -*- coding: UTF-8 -*-
# Introduction: simulation
# Created by galaxy on 2017/6/21 20:30

import os
import re
import random
import shutil
import subprocess
import numpy as np
from multiprocessing import Pool, cpu_count
from collections import defaultdict, OrderedDict
from Bio.Emboss.Applications import NeedleCommandline
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def random_without_same(mi, ma, num):
    temp = list(range(mi, ma))
    random.shuffle(temp)
    random_list = temp[0:num]
    return random_list


def random_aa_seq(length):
    all_aa = 'ACDEFGHIKLMNPQRSTVWY'
    result_seq = ''
    for i in range(length):
        result_seq += random.choice(all_aa)
    return result_seq


def generate_coding_seqs(length=0, number=0):
    standard_codon_table = {'M': ['ATG'],
                            'W': ['TGG'],
                            'K': ['AAA', 'AAG'],
                            'Q': ['CAA', 'CAG'],
                            'Y': ['TAT', 'TAC'],
                            'N': ['AAT', 'AAC'],
                            'E': ['GAA', 'GAG'],
                            'H': ['CAT', 'CAC'],
                            'D': ['GAT', 'GAC'],
                            'F': ['TTT', 'TTC'],
                            'C': ['TGT', 'TGC'],
                            'STOP': ['TAA', 'TAG', 'TGA'],
                            'I': ['ATT', 'ATC', 'ATA'],
                            'A': ['GCT', 'GCC', 'GCA', 'GCG'],
                            'P': ['CCT', 'CCC', 'CCA', 'CCG'],
                            'T': ['ACT', 'ACC', 'ACA', 'ACG'],
                            'G': ['GGT', 'GGC', 'GGA', 'GGG'],
                            'V': ['GTT', 'GTC', 'GTA', 'GTG'],
                            'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
                            'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
                            'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC']}
    all_seqs_dict = defaultdict()
    for i in range(number):
        a = int(length / 3 - 2)
        random_aa = random_aa_seq(a)
        random_dna = 'ATG'
        for j in random_aa:
            random_dna += random.choice(standard_codon_table[j])
        random_dna += random.choice(standard_codon_table['STOP'])
        dna_list = []
        for k in random_dna:
            dna_list.append(k)
        all_seqs_dict[i] = dna_list
    return all_seqs_dict


def generate_mutation_rate(shape, scale, number=1):
    mutation_rate_dict = OrderedDict()
    s = np.random.gamma(shape, scale, number)
    for i in range(number):
        mutation_rate_dict[i] = s[i] * (10 ** -7) * 1000
    return mutation_rate_dict


def each_generation(last_seq_dict, mutation_rate_dict):
    dna = 'ATCG'
    new_seq_dict = defaultdict()
    for id, seq_list in last_seq_dict.items():
        each_mutation_rate = mutation_rate_dict[id]
        random_index = random.random()
        if random_index < each_mutation_rate:
            position = random.randint(3, len(seq_list) - 4)
            # other_dna = dna.replace(seq_list[position], '')
            seq_list[position] = random.choice(dna)
        new_seq_dict[id] = seq_list
    return new_seq_dict


def simulate_generations(first_seq_dict, generation_number, mutation_rate_dict, out_dir, id, item):
    if mutation_rate_dict is None:
        mutation_rate_dict = {0: 1.0}
    new_seq_dict = each_generation(first_seq_dict, mutation_rate_dict)
    for i in range(generation_number):
        if i % 10000 == 0:
            print('Seq_{0}: {1}'.format(str(id), str(int(i / 10000))))
        new_seq_dict = each_generation(new_seq_dict, mutation_rate_dict)
    if item == 0:
        out_file = os.path.join(out_dir, 'Seq_{0}.fasta'.format(str(id)))
    else:
        out_file = os.path.join(out_dir, 'Seq_{0}_HGT.fasta'.format(str(id)))
    with open(out_file, 'w') as f1:
        for seq_id in range(len(new_seq_dict)):
            seq_line = ''.join(str(e) for e in new_seq_dict[seq_id])
            each_line = '>{0}\n{1}\n'.format(str(seq_id), seq_line)
            f1.write(each_line)
    mutation_rate_file = os.path.join(out_dir, 'mt_rate_{0}.txt'.format(str(id)))
    with open(mutation_rate_file, 'w') as f1:
        for m, n in mutation_rate_dict.items():
            result_line = '{0}\t{1}\n'.format(str(m), str(n))
            f1.write(result_line)


def run_in_batch(seq_length=0, seq_number=0, generation_number=0, shape=1.0, scale=1.0, out_dir=''):
    origin_seq_dict = generate_coding_seqs(seq_length, seq_number)
    mutation_rate_dict = generate_mutation_rate(shape, scale, generation_number)
    p = Pool(2)
    for each_id in range(2):
        p.apply_async(simulate_generations, args=(origin_seq_dict, generation_number,
                                                  mutation_rate_dict, out_dir, each_id, 0))
    p.close()
    p.join()


def run_hgt(seq_dict, hgt_generation, out_dir, id):
    print('run HGT {0}'.format(id))
    mutation_rate_file = os.path.join(out_dir, 'mt_rate_{0}.txt'.format(str(id)))
    mutation_rate_dict = OrderedDict()
    with open(mutation_rate_file, 'r') as f1:
        for each_line in f1.readlines():
            a_list = each_line.strip().split('\t')
            mutation_rate_dict[int(a_list[0])] = float(a_list[1])
    simulate_generations(seq_dict, hgt_generation, mutation_rate_dict, out_dir, id, item=1)


def simulate_hgt_event(generation_number, hgt_generation, hgt_number, seq_length, seq_number, shape, scale, out_dir):
    first_run_generations = generation_number - hgt_generation
    # 先进行正常的迭代
    run_in_batch(seq_length, seq_number, first_run_generations, shape, scale, out_dir)
    # 读取生成的中间状态序列，保存为两个字典
    dict_1 = OrderedDict()
    tmp_1_file = os.path.join(out_dir, 'Seq_0.fasta')
    for record in SeqIO.parse(tmp_1_file, 'fasta'):
        tmp_1_list = []
        for m in str(record.seq):
            tmp_1_list.append(m)
        dict_1[int(str(record.id))] = tmp_1_list
    dict_2 = OrderedDict()
    tmp_2_file = os.path.join(out_dir, 'Seq_1.fasta')
    for record in SeqIO.parse(tmp_2_file, 'fasta'):
        tmp_2_list = []
        for n in str(record.seq):
            tmp_2_list.append(n)
        dict_2[int(str(record.id))] = tmp_2_list
    # 生成随机数目的转移基因
    hgt_seq_id_list = random_without_same(0, seq_number, hgt_number)
    all_seq_id_file = os.path.join(out_dir, 'all_gene_category.txt')
    with open(all_seq_id_file, 'w') as f1:
        for each_id in range(seq_number):
            if each_id in hgt_seq_id_list:
                result_line = '{0}\tHGT\n'.format(str(each_id))
            else:
                result_line = '{0}\tVGT\n'.format(str(each_id))
            f1.write(result_line)
    # 将这些转移基因的序列从基因组A转移到基因组B
    for j in hgt_seq_id_list:
        dict_2[j] = dict_1[j]
    p = Pool(2)
    for i in '01':
        if i == '0':
            p.apply_async(run_hgt, args=(dict_1, hgt_generation, out_dir, i))
        else:
            p.apply_async(run_hgt, args=(dict_2, hgt_generation, out_dir, i))
    p.close()
    p.join()


def each_seq_align(each_id, record_list, pair_aln_dir):
    tmp_a_seq = os.path.join(pair_aln_dir, '{0}_a.fasta'.format(str(each_id)))
    SeqIO.write(record_list[0], tmp_a_seq, 'fasta')
    tmp_b_seq = os.path.join(pair_aln_dir, '{0}_b.fasta'.format(str(each_id)))
    SeqIO.write(record_list[1], tmp_b_seq, 'fasta')
    result_file = os.path.join(pair_aln_dir, "{0}.txt".format(str(each_id)))
    needle_cline = NeedleCommandline()
    needle_cline.asequence = tmp_a_seq
    needle_cline.bsequence = tmp_b_seq
    needle_cline.gapopen = 10
    needle_cline.gapextend = 0.5
    needle_cline.outfile = result_file
    devnull = open(os.devnull, 'w')
    try:
        subprocess.call(str(needle_cline), shell=True, stdout=devnull, stderr=devnull)
    except OSError:
        sys.exit(1)
    os.remove(tmp_a_seq)
    os.remove(tmp_b_seq)
    in_pattern = re.compile(r'Identity.*\((\d+\.\d+)%\)')
    gene_alignment_result = ''
    with open(result_file, 'r') as f1:
        for a_line in f1.readlines():
            if 'Identity' in a_line:
                m = re.search(in_pattern, a_line.strip())
                similarity = m.group(1)
                gene_alignment_result = '{0}\t{1}\n'.format(str(each_id), str(similarity))
    os.remove(result_file)
    with open(result_file, 'w') as f2:
        f2.write(gene_alignment_result)


def pairwise_seq_align(out_dir, seq_number, item):
    all_seq_dict = OrderedDict()
    if item == 0:
        item_line = ''
    else:
        item_line = '_HGT'
    seq_0 = os.path.join(out_dir, 'Seq_0{0}.fasta'.format(item_line))
    seq_1 = os.path.join(out_dir, 'Seq_1{0}.fasta'.format(item_line))
    print('Using Seq_0{0}.fasta Seq_1{0}.fasta'.format(item_line))
    for each_seq in [seq_0, seq_1]:
        for record in SeqIO.parse(each_seq, 'fasta'):
            final_record = SeqRecord(record.seq, record.id, description='')
            all_seq_dict.setdefault(str(record.id), []).append(final_record)
    pair_aln_dir = os.path.join(out_dir, 'pairwise_alignment')
    if not os.path.exists(pair_aln_dir):
        os.makedirs(pair_aln_dir)
    p = Pool(cpu_count())
    for each_id, record_list in all_seq_dict.items():
        p.apply_async(each_seq_align, args=(each_id, record_list, pair_aln_dir))
    p.close()
    p.join()
    result_dict = OrderedDict()
    for i in range(seq_number):
        f_path = os.path.join(pair_aln_dir, '{0}.txt'.format(str(i)))
        with open(f_path, 'r') as f1:
            for each_line in f1.readlines():
                a_list = each_line.strip().split('\t')
                result_dict[str(i)] = a_list[1]
    final_result_file = os.path.join(out_dir, 'pairwise_alignment_result.txt')
    all_seq_id_file = os.path.join(out_dir, 'all_gene_category.txt')
    all_seq_dict = defaultdict()
    with open(all_seq_id_file, 'r') as f2:
        for each_line in f2.readlines():
            b_list = each_line.strip().split('\t')
            all_seq_dict[b_list[0]] = b_list[1]
    with open(final_result_file, 'w') as f3:
        header = 'SeqID\tCategory\tSimilarity\n'
        f3.write(header)
        for m, n in result_dict.items():
            result_line = '{0}\t{1}\t{2}\n'.format(str(m), all_seq_dict[m], n)
            f3.write(result_line)
    shutil.rmtree(pair_aln_dir)


def draw_distribution(out_dir, generations_number, hgt_generation, hgt_number):
    final_result_file = os.path.join(out_dir, 'pairwise_alignment_result.txt')
    r_script = os.path.join(my_path, 'draw_distribution.R')
    distribution_file = os.path.join(out_dir, 'pairwise_alignment_result.png')
    title = 'Simulation {0} Generations, {1} HGTs on {2} Generations'\
            .format(str(format(generations_number, '.1e')), str(hgt_number),
                    str(format(hgt_generation, '.1e')))
    devnull = open(os.devnull, 'w')
    try:
        subprocess.call(['Rscript', r_script, final_result_file, distribution_file, title],
                        stdout=devnull, stderr=devnull)
    except OSError:
        sys.exit(1)


def all_process(item, seq_length, seq_number, generation_number,
                shape, scale, out_dir, hgt_generation, hgt_number):
    if item == 0:
        category_file = os.path.join(out_dir, 'all_gene_category.txt')
        with open(category_file, 'w') as f1:
            for i in range(seq_number):
                result_line = '{0}\tVGT\n'.format(str(i))
                f1.write(result_line)
        run_in_batch(seq_length, seq_number, generation_number, shape, scale, out_dir)
    else:
        simulate_hgt_event(generation_number, hgt_generation, hgt_number, seq_length,
                           seq_number, shape, scale, out_dir)
    pairwise_seq_align(out_dir, seq_number, item)
    draw_distribution(out_dir, generation_number, hgt_generation, hgt_number)


if __name__ == '__main__':
    for lina in range(10):
        the_item = 0
        the_shape = 2.3
        the_scale = 1.5
        the_generations_number = int(5 * 10 ** 5)
        the_seq_length = 951
        the_seq_number = 5000
        the_hgt_generation = int(5 * 10 ** 3)
        the_hgt_number = 200
        my_path = os.getcwd()
        the_out_dir = os.path.join(my_path, 'sim_seqs_{0}_{1}_{2}'.format(str(the_generations_number),
                                                                          str(the_item), str(lina + 1)))
        if not os.path.exists(the_out_dir):
            os.makedirs(the_out_dir)
        else:
            shutil.rmtree(the_out_dir)
            os.makedirs(the_out_dir)
        print('Start~')
        all_process(the_item, the_seq_length, the_seq_number, the_generations_number, the_shape, the_scale,
                    the_out_dir, the_hgt_generation, the_hgt_number)
        print('End~')
