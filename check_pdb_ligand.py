#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
list ligands (not standard residues) in pdbfiles in a directory
"""
import os
import sys

from Bio.PDB import *

import lt


def get_ligand(pdb_f):
    lines = pdb_f.readlines()
    ligands = [l[7:10].strip() for l in lines if l.startswith('HET    ')]
    return ligands


def get_ligands(pdb_f):
    standard = ['His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Thr', 'Trp', 'Val',
                'Ala', 'Arg', 'Asp', 'Asn', 'Cys', 'Glu', 'Gln', 'Gly', 'Pro', 'Ser', 'Tyr']
    p = PDBParser()
    s = p.get_structure('', pdb_f)
    residues = [r.get_resname() for r in s.get_residues()]
    ligands = filter(lambda x: x not in standard, residues)
    return ligands


def check_ligands(pdb_f):
    wanted = ['TPO', 'SEP', 'PTR']
    ligands = get_ligands(pdb_f)
    ligands = filter(lambda x: x in wanted, ligands)
    return ligands

@lt.run_time
@lt.log
def main():
    with lt.open_file() as write_f:
        for f in lt.files_in_dir(sys.argv[-1]):
            pdb_id = os.path.split(f)[1].split('.')[0]
            ligands = check_ligands(f)
            if len(ligands) > 0:
                print >> write_f, pdb_id, '\t', ' '.join(ligands)

if __name__ == "__main__":
    main()

