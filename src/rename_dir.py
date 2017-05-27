#!/usr/bin/python
# -*- coding: UTF-8 -*-
# Introduction: 
# Created by Xiangchen Li on 2017/4/22 1:39

import os
from collections import defaultdict

my_path = os.getcwd()
strain_info_file = os.path.join(my_path, 'strain_info.txt')
strain_dict = defaultdict()
with open(strain_info_file, 'r') as f1:
    for each_line in f1.readlines():
        a_list = each_line.strip().split('\t')
        strain_name = a_list[1].split(' ')[-1]
        strain_id = a_list[2]
        strain_dict[strain_id] = strain_name
pair_dir = os.path.join(my_path, 'strain_pair_OG')
for root, dirs, files in os.walk(pair_dir):
    for each_dir in dirs:
        dir_name_list = each_dir.split('_')
        a_strain = strain_dict[dir_name_list[0]]
        b_strain = strain_dict[dir_name_list[1]]
        new_dir = os.path.join(pair_dir, '{0}_{1}'.format(a_strain, b_strain))
        old_dir = os.path.join(pair_dir, each_dir)
        os.rename(old_dir, new_dir)
