#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
from collections import defaultdict

from Bio.PDB import NeighborSearch, PDBParser, Selection

import lt

POSITIVE_BOND_CUTOFF = 4.0


def pdb_neighbors(pdb_f, pdb_id):
    structure = PDBParser().get_structure(pdb_id, pdb_f)
    atom_list = Selection.unfold_entities(structure, 'A')
    ns = NeighborSearch(atom_list)
    center_res = [res for res in structure.get_residues() if res.get_resname() in ['PTR','SEP','TPO']]

    neighbors = []
    for res in center_res:
        if res.get_resname() == 'PTR':
            central_atoms = [atom for atom in res.child_list if atom.get_name() in ['O1P', 'O2P', 'O3P', 'OH']]
        elif res.get_resname() == 'SEP':
            central_atoms = [atom for atom in res.child_list if atom.get_name() in ['O1P', 'O2P', 'O3P', 'OG']]
        elif res.get_resname() == 'TPO':
            central_atoms = [atom for atom in res.child_list if atom.get_name() in ['O1P', 'O2P', 'O3P', 'OG1']]


        positive_atom_neighbors = [ns.search(a.get_coord(),POSITIVE_BOND_CUTOFF) for a in central_atoms]
        positive_atom_neighbors = [atom for atoms in positive_atom_neighbors for atom in atoms]
        positive_atom_neighbors = [atom for atom in positive_atom_neighbors if atom.get_name() in ['NE2','ND1','NZ','NE','NH2','NH1']]

        atom_neighbors = positive_atom_neighbors
        atom_neighbors = list(set(atom_neighbors))

        atom_neighbors = list(set(Selection.unfold_entities(atom_neighbors,'R')))
        atom_neighbors = [r for r in atom_neighbors if r != res]

        if len(atom_neighbors) > 0:
            res = res.get_resname()+'_'+str(res.get_id()[1])+'_'+res.get_parent().get_id()
            atom_neighbors = [n.get_resname()+'_'+str(n.get_id()[1])+'_'+n.get_parent().get_id()+'_'+'S' for n in atom_neighbors]
            neighbors.append((pdb_id,res,atom_neighbors))

    return neighbors

def main():
    res_neighbors = []
    for pdb_f in lt.files_in_dir(sys.argv[-1]):
        f_path,f_name = os.path.split(pdb_f)
        f_name,f_exten = os.path.splitext(f_name)
        if f_exten == '.pdb':
            pdb_id = f_name
            res_neighbors.extend(pdb_neighbors(pdb_f, pdb_id))

    with lt.open_file(file_name='salt_bridge_initial') as w_f:
        for pro, phos_res, neighbor_res in res_neighbors:
            print >> w_f, '{0:<8}{1:<15}{2}'.format(
                pro, phos_res, ',   '.join(neighbor_res))

    lt.pickle_dump(res_neighbors,'salt_bridge')

if __name__ == "__main__":
    main()
