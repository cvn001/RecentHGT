#!/usr/bin/python
# -*- coding: UTF-8 -*-
# Introduction: This script is used to prepare and run recent HGT detection automatically.
# Created by galaxy on 2016/10/19 11:19
#
# DEPENDENCIES
# =============================================================================
# Python (>=2.7):
# o Biopython (http://www.biopython.org)
# o Pyani (https://github.com/widdowquinn/pyani)
# R (>=3.0.1):
# o fitdistplus, ggplot2
# Software:
# o EMBOSS Needle (http://emboss.sourceforge.net/download/)
# o MUMmer (http://mummer.sourceforge.net/)
# =============================================================================
#
# USAGE
# =============================================================================
# calculate_ani.py [options]
#
# Options:
#   -h, --help            show this help message and exit
#   -o OUTDIRNAME, --outdir=OUTDIRNAME
#                         Output directory
#   -i INDIRNAME, --indir=INDIRNAME
#                         Input directory name
#   -v, --verbose         Give verbose output
#   -f, --force           Force file overwriting
#   --noclobber           Don't nuke existing files
#   -p, --part            Which part will be run? [0|1|2|3|4]
#   -t, --threads         How many threads will be used? [default all]
#   -g GFORMAT            Graphics output format(s) [pdf|png|jpg|svg]
#   -l, --logfile         Logfile location
# ==============================================================================
#
# The MIT License
#
# Copyright (c) 2017-2018 Northwest A&U University
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ===============================================================================


import os
import re
import sys
import time
import shutil
import subprocess
import logging.handlers
import traceback
import tarfile
import numpy as np
from argparse import ArgumentParser
from Bio.Emboss.Applications import NeedleCommandline
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool, cpu_count
from collections import defaultdict


def parse_cmdline():
    """
    Parse command-line arguments for script.
    :return: Input command-line arguments
    """
    parser = ArgumentParser(prog="main_process.py")
    parser.add_argument("-i", "--indir", dest="indirname", action="store", default=None,
                        help="Input directory name")
    parser.add_argument("-o", "--outdir", dest="outdirname", action="store", default=None,
                        help="Output directory")
    parser.add_argument("-t", "--threads", type=int, dest="threads", default=cpu_count(),
                        help="How many threads will be used? [default all]")
    parser.add_argument("-p", "--part", type=int, dest="part", default=0,
                        help="Which part will be run? [0|1|2|3|4|5]")
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", default=False,
                        help="Give verbose output")
    parser.add_argument("-l", "--logfile", dest="logfile", action="store", default=None,
                        help="Logfile location")
    parser.add_argument("-f", "--force", dest="force", action="store_true", default=False,
                        help="Force file overwriting")
    parser.add_argument("--noclobber", dest="noclobber", action="store_true", default=False,
                        help="Don't nuke existing files")
    parser.add_argument("-g", dest="gformat", action="store", default="pdf",
                        help="Graphics output format(s) [pdf|png|jpg|svg]")
    return parser.parse_args()


def last_exception():
    """ Returns last exception as a string, or use in logging.
    """
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type, exc_value, exc_traceback))


def make_outdir():
    """Make the output directory, if required.

    This is a little involved.  If the output directory already exists,
    we take the safe option by default, and stop with an error.  We can,
    however, choose to force the program to go on, in which case we can
    either clobber the existing directory, or not.  The options turn out
    as the following, if the directory exists:

    DEFAULT: stop and report the collision
    FORCE: continue, and remove the existing output directory
    NOCLOBBER+FORCE: continue, but do not remove the existing output
    """
    if os.path.exists(args.outdirname):
        if not args.force:
            logger.error("Output directory %s would overwrite existing " +
                         "files (exiting)", args.outdirname)
            sys.exit(1)
        else:
            logger.info("Removing directory %s and everything below it",
                        args.outdirname)
            if args.noclobber:
                logger.warning("NOCLOBBER: not actually deleting directory")
            else:
                shutil.rmtree(args.outdirname)
    logger.info("Creating directory %s", args.outdirname)
    try:
        os.makedirs(args.outdirname)  # We make the directory recursively
        # Depending on the choice of method, a subdirectory will be made for
        # alignment output files
    except OSError:
        # This gets thrown if the directory exists. If we've forced overwrite/
        # delete and we're not clobbering, we let things slide
        if args.noclobber and args.force:
            logger.info("NOCLOBBER+FORCE: not creating directory")
        else:
            logger.error(last_exception)
            sys.exit(1)


def load_strains_info(strain_file, form=1):
    """
    This function is used to load strains information.
    :param form: using different form to return useful information to the functions latter.
    :param strain_file: the path of the file contains the ids and names of strains.
    :return: a python dict
    """
    logger.info('Loading strain information, form {0}'.format(str(form)))
    strain_dict = defaultdict()
    try:
        with open(strain_file, 'r') as f:
            for a_line in f.readlines():
                a_list = a_line.strip().split('\t')
                strain_name = a_list[1].split(' ')[-1]
                strain_id = a_list[2]
                chr_id = a_list[3]
                if form == 1:
                    strain_dict[strain_id] = [strain_name, chr_id]
                elif form == 2:
                    strain_dict[strain_name] = [strain_id, chr_id]
    except IOError:
        logger.error("There is no file contains strain information or the file is locked, please check.")
        logger.error(last_exception())
        sys.exit(1)
    return strain_dict


def compress_delete_outdir(outdir):
    """Compress the contents of the passed directory to .tar.gz and delete."""
    # Compress output in .tar.gz file and remove raw output
    tar_fname = outdir + '.tar.gz'
    logger.info("\tCompressing output from %s to %s", outdir, tar_fname)
    with tarfile.open(tar_fname, "w:gz") as fh:
        fh.add(outdir)
    logger.info("\tRemoving output directory %s", outdir)
    shutil.rmtree(outdir)


def load_genome(genbank, strain_chr_id):
    pattern_1 = re.compile(r'LOCUS\s+(.*?)\s+')
    pattern_2 = re.compile(r'/db_xref="SEED:fig\|(.*)"')
    genome_dict = defaultdict()
    locus_value = 0
    with open(genbank, 'r') as f1:
        for each_line in f1.readlines():
            if 'LOCUS' in each_line:
                m = re.search(pattern_1, each_line.strip())
                locus = m.group(1)
                if locus == strain_chr_id:
                    locus_value = 0
                else:
                    locus_value = 1
            elif '/db_xref="SEED:' in each_line:
                n = re.search(pattern_2, each_line.strip())
                gene_id = n.group(1)
                genome_dict[gene_id] = locus_value  # key: *.peg.*, value: 0,1
    return genome_dict


def og_location(og_id, all_gene_dict, pair_gene_dir):
    if not os.path.exists(pair_gene_dir):
        logger.error("There is no directory contains gene file, please check.")
        logger.error(last_exception())
        sys.exit(1)
    tmp_gene_fasta = os.path.join(pair_gene_dir, og_id + '.fasta')
    re_pattern = re.compile(r'fig\|(\d+\.\d+\.peg\.\d+)\s.*')
    og_loc_value = 0  # 0: Both on Chromosome, 1: One on Plasmid, the other on Chromosome, 2: Both on Plasmid.
    for record in SeqIO.parse(tmp_gene_fasta, 'fasta'):
        m = re.search(re_pattern, record.description)
        gene_id = m.group(1)
        location_value = all_gene_dict[gene_id]
        og_loc_value += location_value
    return og_loc_value


def load_ani():
    ani_out_dir = os.path.join(args.outdirname, 'ANIm')
    if not os.path.exists(ani_out_dir):
        logger.info('IOError: try to open ANIm directory but failed, please check if it exists or renamed.')
        logger.error(last_exception())
        sys.exit(1)
    strain_ani_matrix_file = os.path.join(ani_out_dir, 'ANIm_percentage_identity.tab')
    matrix_strain_dict = defaultdict()
    strain_label_file = os.path.join(args.indirname, 'strain_info.txt')
    strain_dict = load_strains_info(strain_label_file, form=1)
    tmp_matrix_lines = ''
    with open(strain_ani_matrix_file, 'r') as f1:
        all_lines = f1.readlines()
        strain_list = all_lines[0].strip().split('\t')
        for strain_id in strain_list:
            strain_name = strain_dict[strain_id][0]
            strain_index = strain_list.index(strain_id)
            matrix_strain_dict[strain_name] = int(strain_index)  # key: strain name, e.g. CFN42, value: index, e.g. 0,1
        for each_line in all_lines[1:]:
            tmp_list = each_line.strip().split('\t')
            tmp_line = ''
            for i in tmp_list[1:]:
                ani = float('%.3f' % float(i)) * 100
                tmp_line += '{0}\t'.format(str(ani))
            tmp_matrix_lines += tmp_line.strip('\t') + '\n'
    tmp_matrix_file = os.path.join(args.outdirname, 'tmp.txt')
    with open(tmp_matrix_file, 'w') as f2:
        f2.write(tmp_matrix_lines)
    ani_matrix = np.loadtxt(tmp_matrix_file, delimiter='\t')
    os.remove(tmp_matrix_file)
    return ani_matrix, matrix_strain_dict


def each_needle_run(pair_gene_dir, tmp_gene_converted_dir, pair_gene_alignment_dir, og_id, strain_dict):
    """
    This function is used to call Needle program to do pairwise sequence alignment
    :param pair_gene_dir: each homologous gene directory
    :param tmp_gene_converted_dir: used to put some temporary files and will be deleted in the end
    :param pair_gene_alignment_dir: each orthologous gene pair-wised alignment directory
    :param og_id: each orthologous gene id
    :param strain_dict: inherit from load_strains_label function with strain information
    :return: the alignment result of each gene
    """
    if not os.path.exists(pair_gene_dir):
        logger.error("There is no directory contains gene file, please check.")
        logger.error(last_exception())
        sys.exit(1)
    tmp_gene_fasta = os.path.join(pair_gene_dir, og_id + '.fasta')
    converted_records = []
    re_pattern = re.compile(r'fig\|(\d+\.\d+)\.peg\.(\d+)\s(.*)')
    in_pattern = re.compile(r'Identity.*\((\d+\.\d+)%\)')
    annotation = ''
    og_list = []
    for record in SeqIO.parse(tmp_gene_fasta, 'fasta'):
        m = re.search(re_pattern, record.description)
        strain_id = m.group(1)
        gene_id = '{0}.peg.{1}'.format(strain_id, m.group(2))
        og_list.append(gene_id)
        annotation = m.group(3)
        record.id = strain_dict[strain_id][0]
        final_record = SeqRecord(record.seq, record.id, description='')
        converted_records.append(final_record)
    the_strain_fasta = os.path.join(tmp_gene_converted_dir, 'a.fasta')
    other_strain_fasta = os.path.join(tmp_gene_converted_dir, 'b.fasta')
    SeqIO.write(converted_records[0], the_strain_fasta, 'fasta')
    SeqIO.write(converted_records[1], other_strain_fasta, 'fasta')
    result_file = os.path.join(pair_gene_alignment_dir, "{0}.txt".format(og_id))
    needle_cline = NeedleCommandline()
    needle_cline.asequence = the_strain_fasta
    needle_cline.bsequence = other_strain_fasta
    needle_cline.gapopen = 10
    needle_cline.gapextend = 0.5
    needle_cline.outfile = result_file
    devnull = open(os.devnull, 'w')
    try:
        subprocess.call(str(needle_cline), shell=True, stdout=devnull, stderr=devnull)
    except OSError:
        logger.info('Try to call Needle program failed, please check if Needle has been installed successfully.')
        logger.error(last_exception())
        sys.exit(1)
    os.remove(the_strain_fasta)
    os.remove(other_strain_fasta)
    gene_alignment_result = ''
    with open(result_file, 'r') as f:
        for a_line in f.readlines():
            if 'Identity' in a_line:
                m = re.search(in_pattern, a_line.strip())
                identity = m.group(1)
                gene_alignment_result = '{0}\t[{1}|{2}]\t{3}\t{4}\n'.format(og_id, og_list[0],
                                                                            og_list[1], identity, annotation)
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
    pair_results = 'Orthologous\tIdentity\tAnnotation\n'
    for gene in gene_list:
        gene_alignment_result = each_needle_run(strain_pair_genes_dir, tmp_gene_converted_dir,
                                                strain_pair_result_dir, gene, strain_dict)
        pair_results += gene_alignment_result
    shutil.rmtree(tmp_gene_converted_dir)
    pair_result_collection_file = os.path.join(strain_results_dir, '{0}.txt'.format(strain_pair))
    with open(pair_result_collection_file, 'w') as f:
        f.write(pair_results)
    logger.info('{0} is over.'.format(strain_pair))


def first_part():
    """
    It is the first part. It is used to call pyani program to calculate ANI of each strain pair.
    :return: success message
    """
    logger.info('Part 1: Calculating average nucleotide identity (ANI) of each strain pair...')
    genbank_folder = os.path.join(args.indirname, 'genbank')
    fasta_folder = os.path.join(args.outdirname, 'fasta')
    if not os.path.exists(fasta_folder):
        os.makedirs(fasta_folder)
    i = 0
    for root, dirs, files in os.walk(genbank_folder):
        for each_file in files:
            strain_id = os.path.splitext(each_file)[0]
            gbk = os.path.join(genbank_folder, each_file)
            fasta = os.path.join(fasta_folder, strain_id + '.fasta')
            SeqIO.convert(gbk, 'genbank', fasta, 'fasta')
            i += 1
    logger.info('{0} genbank files have been converted to fasta.'.format(str(i)))
    ani_folder = os.path.join(args.outdirname, 'ANIm')
    subprocess.call(['average_nucleotide_identity.py', '-i', fasta_folder,
                     '-o', ani_folder, '-m', 'ANIm', '-g'])
    message = 'Average nucleotide identity analyses have been done.'
    return message


def second_part():
    """
    It is the second part. It is used to call Needle program to do pairwise sequence alignment.
    :return: success message
    """
    logger.info('Part 2: Processing pairwise sequence alignment...')
    strain_pair_dir = os.path.join(args.indirname, 'strain_pair_OG')
    output_dir = os.path.join(args.outdirname, 'strain_pair_OG_alignment')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    strain_label_file = os.path.join(args.indirname, 'strain_info.txt')
    p = Pool(args.threads)
    strain_dict = load_strains_info(strain_label_file, form=1)
    strain_results_collection_dir = os.path.join(args.outdirname, 'strain_pair_result')
    if not os.path.exists(strain_results_collection_dir):
        os.makedirs(strain_results_collection_dir)
    pairs = []
    for root, dirs, files in os.walk(strain_pair_dir):
        for each_dir in dirs:
            pairs.append(each_dir)
    logger.info('Running Needle sequence alignment in {0} threads, please wait...'.format(str(args.threads)))
    for each_pair in pairs:
        p.apply_async(each_strain_pair_run, args=(each_pair, strain_pair_dir, output_dir,
                                                  strain_dict, strain_results_collection_dir))
    p.close()
    p.join()
    message = 'All alignments have been finished.'
    logger.info("Compressing/deleting %s", output_dir)
    compress_delete_outdir(output_dir)
    return message


def third_part():
    """
    It is the third part. It is used to call R script to draw identity distribution pictures.
    :return: success message
    """
    logger.info('Part 3: Drawing alignment distribution pictures...')
    strain_result_dir = os.path.join(args.outdirname, 'strain_pair_result')
    all_result_dir = os.path.join(args.outdirname, 'strain_result')
    if not os.path.exists(all_result_dir):
        os.makedirs(all_result_dir)
    (ani_matrix, matrix_strain_dict) = load_ani()
    result_dict = defaultdict(list)
    for root, dirs, files in os.walk(strain_result_dir):
        for each_file in files:
            file_name = os.path.splitext(each_file)[0]
            file_path = os.path.join(strain_result_dir, each_file)
            pairs = file_name.split('_')
            if pairs[0] not in result_dict:
                result_dict[pairs[0]] = []
            pair_name = str(pairs[0]) + '~' + str(pairs[1])
            pair_ani = ani_matrix[matrix_strain_dict[pairs[0]]][matrix_strain_dict[pairs[1]]]
            if pair_ani < 95:
                with open(file_path) as f3:
                    for line in f3.readlines()[1:]:
                        result_line = '{0} ({1}%)\t{2}'.format(pair_name, str(pair_ani), line)
                        result_dict[pairs[0]].append(result_line)
    r_script = os.path.join(src_dir_name, 'draw_distribution.R')
    header_line = 'Pair\tOrthologous\tIdentity\tAnnotation\n'
    devnull = open(os.devnull, 'w')
    logger.info('Saving {0} pictures to {1} format.'.format(str(len(result_dict)), args.gformat))
    for each_strain, results in result_dict.items():
        strain_result_file = os.path.join(all_result_dir, each_strain + '.txt')
        with open(strain_result_file, 'w') as f4:
            f4.write(header_line)
            for each_result in results:
                f4.write(each_result)
        each_result_picture = os.path.join(all_result_dir, '{0}.{1}'.format(each_strain, args.gformat))
        try:
            subprocess.call(['Rscript', r_script, strain_result_file, each_result_picture],
                            stdout=devnull, stderr=devnull)
        except OSError:
            logger.info('Try to run {0} but failed, please check.'.format(r_script))
            logger.error(last_exception())
            sys.exit(1)
    message = 'All pictures have been saved in {0}'.format(all_result_dir)
    return message


def fourth_part():
    """
    It is the fourth part. It is used to call R script to infer the number of recent HGT genes.
    :return: success message
    """
    result_dir = args.outdirname
    logger.info('Part 4: Processing recent HGT detection...')
    strain_result_dir = os.path.join(result_dir, 'strain_pair_result')
    r_script = os.path.join(src_dir_name, 'rHGT_alpha.R')
    param_min = 50.0
    param_max = 98.5
    devnull = open(os.devnull, 'w')
    result_file = os.path.join(result_dir, 'recent_HGT_results.txt')
    try:
        subprocess.call(['Rscript', r_script, strain_result_dir, result_file,
                         str(param_min), str(param_max)], stdout=devnull, stderr=devnull)
    except OSError:
        logger.info('Try to run {0} but failed, please check.'.format(r_script))
        logger.error(last_exception())
        sys.exit(1)
    message = 'Recent HGT detections have been finished.'
    return message


def fifth_part():
    """
    It is the fifth part. It is used to call R script to draw the comparison between the number
    of rHGTs and specific location genes (chromosome and plasmid genes) to show the accuracy.
    :return: success message
    """
    logger.info('Part 5: Drawing the comparison pictures...')
    logger.info('Loading recent HGT results.')
    genome_dir = os.path.join(args.indirname, 'genbank')
    if not os.path.exists(genome_dir):
        logger.error('These is no directory contains genbank files, please check here.')
        sys.exit(1)
    strain_id_list = []
    strain_label_file = os.path.join(args.indirname, 'strain_info.txt')
    strain_dict = load_strains_info(strain_label_file, form=1)
    all_gene_dict = defaultdict()
    for root, dirs, files in os.walk(genome_dir):
        for each_file in files:
            fname = os.path.splitext(each_file)
            strain_id = fname[0]                # The fname[0] is a strain id, e.g. 379.140.
            strain_id_list.append(fname[0])
            genbank = os.path.join(genome_dir, each_file)
            strain_chr_id = strain_dict[strain_id][1]
            each_genome_dict = load_genome(genbank, strain_chr_id)
            all_gene_dict.update(each_genome_dict)
    result_dir = args.outdirname
    strain_result_dir = os.path.join(result_dir, 'strain_pair_result')
    hgt_result_file = os.path.join(args.outdirname, 'recent_HGT_results.txt')
    strain_pair_gene_dir = os.path.join(args.indirname, 'strain_pair_OG')
    detect_dict = defaultdict()
    with open(hgt_result_file, 'r') as f1:
        for each_line in f1.readlines()[1:]:
            a_list = each_line.strip().split('\t')
            detect_dict[a_list[0]] = a_list[1]  # key: strain_pair || Value: inferred HGT number
    pair_gene_location_dict = defaultdict()
    for root, dirs, files in os.walk(strain_result_dir):
        for each_file in files:
            fname = os.path.splitext(each_file)
            strain_pair = fname[0]
            pair_gene_dir = os.path.join(strain_pair_gene_dir, strain_pair)
            result_file = os.path.join(strain_result_dir, each_file)
            plasmid_gene = 0
            ambiguity_gene = 0
            chr_gene = 0
            with open(result_file, 'r') as f2:
                for each_line in f2.readlines()[1:]:
                    b_list = each_line.strip().split('\t')
                    identity = float(b_list[2])
                    if identity >= 98.5:
                        og_id = b_list[0]
                        og_loc_value = og_location(og_id, all_gene_dict, pair_gene_dir)
                        if og_loc_value == 0:
                            chr_gene += 1
                        elif og_loc_value == 1:
                            ambiguity_gene += 1
                        elif og_loc_value == 2:
                            plasmid_gene += 1
            pair_gene_location_dict[strain_pair] = '{0}|{1}|{2}'.format(str(chr_gene), str(plasmid_gene),
                                                                        str(ambiguity_gene))
    tmp_result_dict = defaultdict()
    with open(hgt_result_file, 'w') as f3:
        header_line = 'strain.pair\trHGT.number\tgene.number\n'
        result_lines = ''
        for each_pair, loc_result in pair_gene_location_dict.items():
            detect_num = detect_dict[each_pair]
            result_lines += '{0}\t{1}\t{2}\n'.format(each_pair, detect_dict[each_pair], loc_result)
            loc_list = loc_result.split('|')
            tmp_result_dict[each_pair] = [detect_num, loc_list[0], loc_list[1], loc_list[2]]
        f3.write(header_line + result_lines)
    tmp_order_list = ['rHGT', 'Chromosome', 'Plasmid', 'Ambiguity']
    logger.info('Drawing comparison pictures to show detection accuracy.')
    (ani_matrix, matrix_strain_dict) = load_ani()
    combine_result_dir = os.path.join(args.outdirname, 'combined_results')
    if not os.path.exists(combine_result_dir):
        os.makedirs(combine_result_dir)
    strain_list = []
    for each_strain_id in strain_dict.keys():
        each_strain_name = strain_dict[each_strain_id][0]
        strain_list.append(each_strain_name)
    r_script = os.path.join(src_dir_name, 'draw_rHGT_loc_genes.R')
    devnull = open(os.devnull, 'w')
    for query_strain in strain_list:
        query_strain_result_file = os.path.join(combine_result_dir, query_strain + '.txt')
        with open(query_strain_result_file, 'w') as f4:
            h = 'query.strain\ttype\tnumber\tani\n'
            l = ''
            for other_strain in strain_list:
                if other_strain != query_strain:
                    tmp_result_list = tmp_result_dict['{0}_{1}'.format(query_strain, other_strain)]
                    pair_ani = ani_matrix[matrix_strain_dict[query_strain]][matrix_strain_dict[other_strain]]
                    for i in tmp_result_list:
                        if i != '0':
                            index = tmp_result_list.index(i)
                            l += '{0} ({3})\t{1}\t{2}\t{3}\n'.format(other_strain, tmp_order_list[index], i, pair_ani)
            f4.write(h + l)
        comparison_pictures = os.path.join(combine_result_dir, '{0}.{1}'.format(query_strain, args.gformat))
        try:
            subprocess.call(['Rscript', r_script, query_strain_result_file, comparison_pictures],
                            stdout=devnull, stderr=devnull)
        except OSError:
            logger.info('Try to run {0} but failed, please check.'.format(r_script))
            logger.error(last_exception())
            sys.exit(1)
    message = 'All pictures have been saved in {0}'.format(combine_result_dir)
    return message


def auto_run():
    """
    This function is used to run all parts automatically.
    :return: success message
    """
    logger.info('Automatically run all processes...')
    message_1 = first_part()
    logger.info(message_1)
    message_2 = second_part()
    logger.info(message_2)
    message_3 = third_part()
    logger.info(message_3)
    message_4 = fourth_part()
    logger.info(message_4)
    message_5 = fifth_part()
    logger.info(message_5)
    done_message = 'All 5 parts have been done.'
    return done_message


def separate_run(part):
    """
    This function is used to choose one of four part to run by user.
    :param part: The part will be run.
    :return: success message
    """
    message = ''
    if part == 1:
        message = first_part()
    elif part == 2:
        message = second_part()
    elif part == 3:
        message = third_part()
    elif part == 4:
        message = fourth_part()
    elif part == 5:
        message = fifth_part()
    return message


if __name__ == '__main__':
    # Run as script
    # Parse command-line
    args = parse_cmdline()
    # Set up logging
    logger = logging.getLogger('main_process.py: %s' % time.asctime())
    t0 = time.time()
    src_dir_name = 'src'
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    # Was a logfile specified? If so, use it
    if args.logfile is not None:
        try:
            log_stream = open(args.logfile, 'w')
            err_handler_file = logging.StreamHandler(log_stream)
            err_handler_file.setFormatter(err_formatter)
            err_handler_file.setLevel(logging.INFO)
            logger.addHandler(err_handler_file)
        except IOError:
            logger.error("Could not open %s for logging", args.logfile)
            sys.exit(1)
    # Do we need verbosity?
    if args.verbose:
        err_handler.setLevel(logging.INFO)
    else:
        err_handler.setLevel(logging.WARNING)
    logger.addHandler(err_handler)
    # Report arguments, if verbose
    # logger.info("pyani version: %s", VERSION)
    logger.info(args)
    logger.info("command-line: %s", ' '.join(sys.argv))
    # Have we got an input and output directory? If not, exit.
    if args.indirname is None:
        logger.error("No input directory name (exiting)")
        sys.exit(1)
    logger.info("Input directory: %s", args.indirname)
    if args.outdirname is None:
        logger.error("No output directory name (exiting)")
        sys.exit(1)
    make_outdir()
    logger.info("Output directory: %s", args.outdirname)
    run_message = ''
    if args.part == 0:
        try:
            run_message = auto_run()
        except OSError:
            logger.info('Try to run all 4 parts automatically but failed, please check.')
    elif 1 <= args.part <= 5:
        try:
            run_message = separate_run(args.part)
        except OSError:
            logger.info('Try to run part {0} but failed, please check.'.format(args.part))
            logger.error(last_exception())
            sys.exit(1)
    else:
        logger.error('Part {0} is not valid, please choose one of [0|1|2|3|4|5].'.format(args.part))
        sys.exit(1)
    logger.info(run_message)
    # Report that we've finished
    logger.info("All jobs have been done: %s.", time.asctime())
    logger.info("Total time taken: %.2fs", (time.time() - t0))
