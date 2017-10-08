#!/usr/bin/python
# -*- coding: UTF-8 -*-
# Introduction: Simulation, using with rHGT_sim.R,
#               draw_distribution.R, pyani and needle in EMBOSS.
#               Required Packages: Biopython, Numpy
# Created by galaxy on 2017/6/21 20:30

import os
import re
import sys
import random
import shutil
import datetime
import subprocess
import numpy as np
from multiprocessing import Pool, cpu_count
from collections import defaultdict, OrderedDict
from Bio.Emboss.Applications import NeedleCommandline
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def random_without_same(mi=0, ma=1, num=1):
    # Generate the ids of HGT genes
    temp = list(range(mi, ma))
    random.shuffle(temp)
    random_list = temp[0:num]
    return random_list


def random_aa_seq(length=0):
    # Generate the amino acid sequence at a given length
    all_aa = 'ACDEFGHIKLMNPQRSTVWY'
    result_seq = ''
    for i in range(length):
        result_seq += random.choice(all_aa)
    return result_seq


def generate_coding_seqs(length=0, number=0):
    # Randomly generate the protein coding sequences
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


def generate_mutation_rate(shape=0.0, scale=0.0, number=1):
    # Randomly generate the mutation rate
    mutation_rate_dict = OrderedDict()
    s = np.random.gamma(shape, scale, number)
    for i in range(number):
        mutation_rate_dict[i] = s[i] * (10 ** -7) * 1000
    return mutation_rate_dict


def each_generation(last_seq_dict=None, mutation_rate_dict=None):
    # Simulate the mutation in each generation.
    # If the simulated mutation rate less than the given rate, the mutation will happen.
    if last_seq_dict is None:
        last_seq_dict = {}
    if mutation_rate_dict is None:
        mutation_rate_dict = {}
    dna = 'ATCG'
    new_seq_dict = defaultdict()
    for id, seq_list in last_seq_dict.items():
        each_mutation_rate = mutation_rate_dict[id]
        random_index = random.random()
        if random_index < each_mutation_rate:
            position = random.randint(3, len(seq_list) - 4)
            seq_list[position] = random.choice(dna)
        new_seq_dict[id] = seq_list
    return new_seq_dict


def simulate_generations(first_seq_dict=None, generation_number=0,
                         mutation_rate_dict=None, out_dir='', id=0, item=0):
    # Simulate genomes evolution
    if mutation_rate_dict is None:
        mutation_rate_dict = {}
    if first_seq_dict is None:
        first_seq_dict = {}
    new_seq_dict = each_generation(first_seq_dict, mutation_rate_dict)
    for i in range(generation_number):
        new_seq_dict = each_generation(new_seq_dict, mutation_rate_dict)
    if item == 0:
        out_file = os.path.join(out_dir, 'Seq_{0}.fasta'.format(str(id)))
    else:
        out_file = os.path.join(out_dir, 'Seq_{0}_HGT.fasta'.format(str(id)))
    concatenated_seqs = ''
    with open(out_file, 'w') as f1:
        for seq_id in range(len(new_seq_dict)):
            seq_line = ''.join(str(e) for e in new_seq_dict[seq_id])
            each_line = '>{0}\n{1}\n'.format(str(seq_id), seq_line)
            concatenated_seqs += seq_line
            f1.write(each_line)
    genome_dir = os.path.join(out_dir, 'genome')
    if not os.path.exists(genome_dir):
        os.makedirs(genome_dir)
    genome_file = os.path.join(genome_dir, 'genome_{0}.fasta').format(str(id))
    genome_record = SeqRecord(Seq(concatenated_seqs), id=str(id), description='')
    SeqIO.write(sequences=genome_record, handle=genome_file, format='fasta')
    mutation_rate_file = os.path.join(out_dir, 'mt_rate_{0}.txt'.format(str(id)))
    with open(mutation_rate_file, 'w') as f3:
        for m, n in mutation_rate_dict.items():
            result_line = '{0}\t{1}\n'.format(str(m), str(n))
            f3.write(result_line)


def run_in_batch(seq_length=0, seq_number=0, generation_number=0,
                 shape=1.0, scale=1.0, out_dir=''):
    # run all functions
    origin_seq_dict = generate_coding_seqs(seq_length, seq_number)
    mutation_rate_dict = generate_mutation_rate(shape, scale, seq_number)
    for i in [0, 1]:
        simulate_generations(origin_seq_dict, generation_number,
                             mutation_rate_dict, out_dir, i, 0)


def run_hgt(seq_dict=None, hgt_generation=0, out_dir='', id=0):
    # Do the rest generations after the HGT event
    if seq_dict is None:
        seq_dict = {}
    mutation_rate_file = os.path.join(out_dir, 'mt_rate_{0}.txt'.format(str(id)))
    mutation_rate_dict = OrderedDict()
    with open(mutation_rate_file, 'r') as f1:
        for each_line in f1.readlines():
            a_list = each_line.strip().split('\t')
            mutation_rate_dict[int(a_list[0])] = float(a_list[1])
    simulate_generations(seq_dict, hgt_generation, mutation_rate_dict,
                         out_dir, id, item=1)


def simulate_hgt_event(generation_number=0, hgt_generation=0,
                       hgt_number=0, seq_length=0, seq_number=0,
                       shape=0.0, scale=0.0, out_dir=''):
    # Simulate HGT event, all in one event
    first_run_generations = generation_number - hgt_generation
    # Do generate
    run_in_batch(seq_length, seq_number, first_run_generations,
                 shape, scale, out_dir)
    # Read all sequences at the intermediate state
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
    # Generate HGT genes
    hgt_seq_id_list = random_without_same(0, seq_number, hgt_number)
    all_seq_id_file = os.path.join(out_dir, 'all_gene_category.txt')
    with open(all_seq_id_file, 'w') as f1:
        for each_id in range(seq_number):
            if each_id in hgt_seq_id_list:
                result_line = '{0}\tHGT\n'.format(str(each_id))
            else:
                result_line = '{0}\tVGT\n'.format(str(each_id))
            f1.write(result_line)
    # Gene Transfer
    for j in hgt_seq_id_list:
        dict_2[j] = dict_1[j]
    for k in [0, 1]:
        if k == 0:
            run_hgt(dict_1, hgt_generation, out_dir, k)
        else:
            run_hgt(dict_2, hgt_generation, out_dir, k)


def fetch_ani(out_dir='', item=0):
    # Calculate the ANI value of two generated genomes
    genome_dir = os.path.join(out_dir, 'genome')
    ani_out_dir = os.path.join(out_dir, 'ani')
    if item < 2:
        ani_cmd = 'average_nucleotide_identity.py -i {0} -o {1} ' \
                  '-m ANIm --workers 1 --force'.format(genome_dir,
                                                       ani_out_dir)
        devnull = open(os.devnull, 'w')
        try:
            subprocess.call(str(ani_cmd), shell=True,
                            stdout=devnull, stderr=devnull)
        except OSError:
            print('pyani calling failed, please check.')
            sys.exit(1)
    ani_file = os.path.join(ani_out_dir,
                            'ANIm_percentage_identity.tab')
    with open(ani_file, 'r') as f1:
        all_lines = f1.readlines()
        result_line = all_lines[1]
        ani = result_line.strip().split('\t')[-1]
    result_ani = round(float(ani) * 100, 2)
    return result_ani


def each_seq_align(each_id=0, record_list=list(), pair_aln_dir=''):
    # prepare pairwise sequence alignment files
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
        subprocess.call(str(needle_cline), shell=True,
                        stdout=devnull, stderr=devnull)
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
                gene_alignment_result = '{0}\t{1}\n'.format(str(each_id),
                                                            str(similarity))
    os.remove(result_file)
    with open(result_file, 'w') as f2:
        f2.write(gene_alignment_result)


def pairwise_seq_align(out_dir='', seq_number=0, item=0):
    # Call Needle program to do alignment in batch
    all_seq_dict = OrderedDict()
    if item == 0:
        item_line = ''
    else:
        item_line = '_HGT'
    seq_0 = os.path.join(out_dir, 'Seq_0{0}.fasta'.format(item_line))
    seq_1 = os.path.join(out_dir, 'Seq_1{0}.fasta'.format(item_line))
    for each_seq in [seq_0, seq_1]:
        for record in SeqIO.parse(each_seq, 'fasta'):
            final_record = SeqRecord(record.seq, record.id, description='')
            all_seq_dict.setdefault(str(record.id), []).append(final_record)
    pair_aln_dir = os.path.join(out_dir, 'pairwise_alignment')
    if not os.path.exists(pair_aln_dir):
        os.makedirs(pair_aln_dir)
    for each_id, record_list in all_seq_dict.items():
        each_seq_align(each_id, record_list, pair_aln_dir)
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


def draw_distribution(out_dir='', generations_number=0,
                      hgt_generation=0, hgt_number=0, ani=0.0):
    # Call a R script (draw_distribution.R) to draw sequence similarity distribution
    final_result_file = os.path.join(out_dir, 'pairwise_alignment_result.txt')
    r_script = os.path.join(my_path, 'draw_distribution.R')
    distribution_file = os.path.join(out_dir, 'pairwise_alignment_result.png')
    if not os.path.exists(distribution_file):
        title = 'Total: {0} G | {1} HGT: {2} G | ANI: {3}%'\
                .format(str(format(generations_number, '.1e')), str(hgt_number),
                        str(format(hgt_generation, '.1e')), str(ani))
        devnull = open(os.devnull, 'w')
        try:
            subprocess.call(['Rscript', r_script, final_result_file,
                             distribution_file, title],
                            stdout=devnull, stderr=devnull)
        except OSError:
            print('draw_distribution.R calling failed, please check.')
            sys.exit(1)


def recent_hgt_detection(cutoff=list(), out_dir=''):
    # Call rHGT_sim.R script to use EM model to detect HGT number and compare it with the real HGT number
    r_script = os.path.join(my_path, 'rHGT_sim.R')
    hgt_result_list = []
    resource_similarity_file = os.path.join(out_dir, 'pairwise_alignment_result.txt')
    real_hgt_num = 0
    with open(resource_similarity_file, 'r') as f:
        for each_line in f.readlines()[1:]:
            a_list = each_line.strip().split('\t')
            if a_list[1] == 'HGT':
                if cutoff[0] <= float(a_list[2]) < cutoff[1]:
                    real_hgt_num += 1
    for i in cutoff:
        tmp_dir = os.path.join(out_dir, 'tmp_{0}'.format(str(i)))
        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)
        tmp_similarity_file = os.path.join(tmp_dir, 'similarity_{0}.txt'.format(str(i)))
        tmp_result_file = os.path.join(tmp_dir, 'tmp_hgt_result_{0}.txt'.format(str(i)))
        shutil.copy(resource_similarity_file, tmp_similarity_file)
        devnull = open(os.devnull, 'w')
        try:
            subprocess.call(['Rscript', r_script, tmp_dir, tmp_result_file,
                             str(i)], stdout=devnull, stderr=devnull)
            with open(tmp_result_file, 'r') as f1:
                tmp_result_line = f1.readlines()[1]
                hgt_num = int(tmp_result_line.strip().split('\t')[1])
                hgt_result_list.append(hgt_num)
        except OSError:
            print('rHGT_sim.R calling failed, please check.')
            sys.exit(1)
        shutil.rmtree(tmp_dir)
    interval_hgt_num = abs(hgt_result_list[1] - hgt_result_list[0])
    hgt_result_file = os.path.join(out_dir, 'hgt_result.txt')
    hgt_line = ''
    with open(hgt_result_file, 'w') as f2:
        header = 'HGT({0})Num\tHGT({1})Num\tHGT({0}~{1})Num\tRealNum\n'.format(str(cutoff[0]),
                                                                               str(cutoff[1]))
        hgt_line += '{0}\t{1}\t{2}\t{3}\n'.format(str(hgt_result_list[0]),
                                                  str(hgt_result_list[1]),
                                                  str(interval_hgt_num),
                                                  str(real_hgt_num))
        f2.write(header + hgt_line)
    return hgt_line


def all_process(item=0, seq_length=0, seq_number=0, generation_number=0, cutoff=list(),
                shape=0.0, scale=0.0, out_dir='', hgt_generation=0, hgt_number=0):
    # Run all generations
    # item choice:
    # 0 no HGT happened
    # 1 Simulate HGT generations
    # 2 HGT number detection
    if item <= 1:
        if item == 0:
            category_file = os.path.join(out_dir, 'all_gene_category.txt')
            with open(category_file, 'w') as f1:
                for i in range(seq_number):
                    result_line = '{0}\tVGT\n'.format(str(i))
                    f1.write(result_line)
            run_in_batch(seq_length, seq_number, generation_number,
                         shape, scale, out_dir)
        else:
            simulate_hgt_event(generation_number, hgt_generation, hgt_number,
                               seq_length, seq_number, shape, scale, out_dir)
        ani = fetch_ani(out_dir, item)
        pairwise_seq_align(out_dir, seq_number, item)
        draw_distribution(out_dir, generation_number,
                          hgt_generation, hgt_number, ani)
        output = 'nothing here'
    else:
        output = recent_hgt_detection(cutoff, out_dir)
        ani = fetch_ani(out_dir, item)
    return output, ani


def each_run(item=0, seq_length=0, seq_number=0, generation_number=0,
             cutoff=list(), all_out_dir='', shape=0.0, scale=0.0, li=0,
             collections_dir='', hgt_generation=0, hgt_number=0):
    # Run a mode
    if item > 0:
        out_dir = os.path.join(all_out_dir,
                               'sim_seqs_{0}_1_{1}'.format(str(generation_number / (10 ** 5)),
                                                           str(li + 1)))
    else:
        out_dir = os.path.join(all_out_dir,
                               'sim_seqs_{0}_0_{1}'.format(str(generation_number / (10 ** 5)),
                                                           str(li + 1)))
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    result_line, ani = all_process(item, seq_length, seq_number, generation_number, cutoff,
                                   shape, scale, out_dir, hgt_generation, hgt_number)
    result_picture = os.path.join(out_dir, 'pairwise_alignment_result.png')
    collect_picture = os.path.join(collections_dir, '{0}.png'.format(str(li + 1)))
    if not os.path.exists(collect_picture):
        shutil.copy(result_picture, collect_picture)
    out_line = '{0}\t{1}'.format(str(ani), result_line)
    return out_line


def auto_run(item_list=list(), seq_length=0, seq_number=0, generation_number=0,
             cutoff=list(), shape=0.0, scale=0.0, all_out_dir='',
             hgt_generation=0, hgt_number=0, run_number=0):
    # auto run all processes
    if not os.path.exists(all_out_dir):
        os.makedirs(all_out_dir)
    collections_dir = os.path.join(all_out_dir, 'hgt_pictures_collection')
    for item in item_list:
        if item < 2:
            if not os.path.exists(collections_dir):
                os.makedirs(collections_dir)
            else:
                shutil.rmtree(collections_dir)
                os.makedirs(collections_dir)
            p = Pool(processes=cpu_count())
            for li in range(run_number):
                p.apply_async(each_run, args=(item, seq_length, seq_number, generation_number,
                                              cutoff, all_out_dir, shape, scale, li,
                                              collections_dir, hgt_generation, hgt_number,))
            p.close()
            p.join()
        else:
            a = float(generation_number / (10 ** 5))
            all_hgt_results_file = os.path.join(all_out_dir,
                                                'all_hgt_results_{0}g_{1}_{2}.txt'.format(str(a),
                                                                                          str(the_cutoff[0]),
                                                                                          str(the_cutoff[1])))
            all_hgt_lines = ''
            p = Pool(processes=cpu_count())
            results = []
            for li in range(run_number):
                res = p.apply_async(each_run, args=(item, seq_length, seq_number, generation_number,
                                                    cutoff, all_out_dir, shape, scale, li,
                                                    collections_dir, hgt_generation, hgt_number,))
                results.append(res)
            p.close()
            p.join()
            for i in range(len(results)):
                each_result_line = '{0}\t{1}'.format(str(i + 1), str(results[i].get()))
                all_hgt_lines += each_result_line
            with open(all_hgt_results_file, 'w') as file:
                header = 'RunID\tANI\tHGT({0})Num\tHGT({1})Num' \
                         '\tHGT({0}~{1})Num\tRealNum\n'.format(str(the_cutoff[0]),
                                                               str(the_cutoff[1]))
                file.write(header + all_hgt_lines)


if __name__ == '__main__':
    print('Start')
    start_time = datetime.datetime.now()
    my_path = os.getcwd()
    '''
    Running Mode: 0, 1, 2
    0：No HGT happened
    1：HGT simulation
    2：HGT number detection 
    '''
    # Parameters
    the_item_list = [2]
    the_shape = 2.3
    the_scale = 1.5
    the_seq_length = 951
    the_seq_number = 5000
    the_hgt_number = 250
    the_hgt_generation = 1.0
    the_run_number = 100
    # Do simulation
    for lina in [97.0, 97.5, 98.0, 98.5, 99.0]:
        the_cutoff = [lina, 100.5]
        for the_hgt_generation in [0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]:
            the_generation_number = 5.0
            the_out_dir = 'sim_{0}g_{1}g_{2}'.format(str(the_generation_number),
                                                     str(the_hgt_generation),
                                                     str(the_run_number))
            the_generation_number = int(the_generation_number * 10 ** 5)
            the_hgt_generation = int(the_hgt_generation * 10 ** 4)
            auto_run(the_item_list, the_seq_length, the_seq_number,
                     the_generation_number, the_cutoff, the_shape,
                     the_scale, the_out_dir, the_hgt_generation,
                     the_hgt_number, the_run_number)
    end_time = datetime.datetime.now()
    total_time = end_time - start_time
    minutes, seconds = divmod(total_time.seconds, 60)
    hours, minutes = divmod(minutes, 60)
    print("End: %02d:%02d:%02d" % (hours, minutes, seconds))
