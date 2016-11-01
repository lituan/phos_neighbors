#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
transform ent files to pdb files
change extention to '.pdb'
"""
import sys
import os
import re


def files_in_dir(directory):
    for root, dirs, files in os.walk(directory):
        for f in files:
            yield os.path.join(root, f)

for f in files_in_dir(sys.argv[-1]):
    f_dir, f_name = os.path.split(f)
    f_short_name, f_extention = os.path.splitext(f_name)
    if f_extention == '.ent':
        with open(f) as o_f:
            lines = o_f.readlines()
            new_name = f_short_name[3:] + '.pdb'
            new_name = os.path.join(f_dir, new_name)
            with open(new_name, 'w') as w_f:
                for line in lines:
                    print >> w_f, line
