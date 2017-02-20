#!/usr/bin/python
# -*- coding: UTF-8 -*-

import os.path
import sys

maxbit = sys.argv[1]
curDir = os.getcwd()
oldId = "all_I_2.0_c_{0}_m_maxbit_".format(maxbit)
newId = ""
for parent, dirnames, filenames in os.walk(curDir):  
    for filename in filenames:
        if filename.find(oldId) != -1:
            newName = filename.replace(oldId, newId)
            # print(filename, "---->", newName)
            os.rename(os.path.join(parent, filename), os.path.join(parent, newName))
# os.system("pause")
