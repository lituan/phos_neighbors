#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import Bio
from Bio.PDB import PDBList

def read_pdb_list(id_f):
    with open(id_f) as o_f:
        lines = o_f.readlines()
        lines = [line.split()[0] for line in lines if len(line.split()) > 0]
        return lines


def pdb_download(pdb_ids,dir_name):

    pdbl = PDBList()
    while len(pdb_ids) > 0:
        try:
            for i in pdb_ids:
                pdbl.retrieve_pdb_file(i,pdir=dir_name)

                f_name = os.path.join(dir_name,'pdb'+i+'.ent')
                f_name_new = os.path.join(dir_name,i+'.pdb')
                ent = open(f_name,'r')
                pdb = open(f_name_new,'w')
                pdb.writelines(ent.readlines())
                ent.close()
                pdb.close()
                os.remove(f_name)

                pdb_ids.pop(pdb_ids.index(i))
        except:
            continue



def main():
    pdb_ids = read_pdb_list(sys.argv[-1])

    f_path,f_name = os.path.split(sys.argv[-1])
    f_name,f_ext = os.path.splitext(f_name)
    dir_name = f_name+'_pdb'
    dir_name = os.path.join(f_path,f_name+'_pdb')
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)


    pdb_download(pdb_ids,dir_name)


if __name__ == "__main__":
    main()
