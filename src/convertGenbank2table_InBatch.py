#!/usr/bin/python3
# -*- coding: UTF-8 -*-
# Introduction: 本程序用于
# Created by galaxy on 2016/9/8 10:50

import os

my_path = os.getcwd()
gbk_dir = os.path.join(my_path, 'gbk')
batch_lines = []
for root, dirs, files in os.walk(gbk_dir):
    for each_file in files:
        gbk_path = 'gbk/{0}'.format(each_file)
        each_cmd = 'python convertGenbank2table.py -g {0} -v 1'.format(gbk_path)
        batch_lines.append(each_cmd)
batch_file = os.path.join(my_path, 'convertGenbank2table.sh')
with open(batch_file, 'w') as f1:
    for each_batch_cmd in batch_lines:
        f1.write('{0}\n'.format(each_batch_cmd))
convert_cmd = 'sh convertGenbank2table.sh'
os.system(convert_cmd)


