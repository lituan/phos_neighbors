import os
import sys
from collections import defaultdict

from Bio.PDB import NeighborSearch, PDBParser, Selection

import lt

BOND_CUTOFF = 3.6
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

        atom_neighbors = [ns.search(a.get_coord(),BOND_CUTOFF) for a in central_atoms]
        atom_neighbors = [atom for atoms in atom_neighbors for atom in atoms]

        positive_atom_neighbors = [ns.search(a.get_coord(),POSITIVE_BOND_CUTOFF) for a in central_atoms]
        positive_atom_neighbors = [atom for atoms in positive_atom_neighbors for atom in atoms]
        positive_atom_neighbors = [atom for atom in positive_atom_neighbors if atom.get_name() in ['NE2','ND1','NZ','NE','NH2','NH1']]

        atom_neighbors.extend(positive_atom_neighbors)
        atom_neighbors = list(set(atom_neighbors))

        #filter self
        atom_neighbors = [atom for atom in atom_neighbors if not atom.get_parent() == res]

        # only consider those containing N or O
        atom_neighbors = [atom for atom in atom_neighbors if 'N' in atom.get_name() or 'O' in atom.get_name()]

        ## ignore water
        atom_neighbors = [atom for atom in atom_neighbors if not atom.get_parent().get_resname() == 'HOH']

        # filter main_chain O, they are not donor
        atom_neighbors = [atom for atom in atom_neighbors if not atom.get_name() == 'O']

        # filter O in N Q, they are not donor
        atom_neighbors = [atom for atom in atom_neighbors if not (atom.get_name() == 'OD1' and atom.get_parent().get_resname() == 'ASN')]
        atom_neighbors = [atom for atom in atom_neighbors if not (atom.get_name() == 'OE1' and atom.get_parent().get_resname() =='GLN')]

        # filter O in D E, they are not donor
        # atom_neighbors = [atom for atom in atom_neighbors if not (atom.get_name() == 'OD1' and atom.get_parent().get_resname() == 'ASP')]
        # atom_neighbors = [atom for atom in atom_neighbors if not (atom.get_name() == 'OD2' and atom.get_parent().get_resname() == 'ASP')]
        # atom_neighbors = [atom for atom in atom_neighbors if not (atom.get_name() == 'OE1' and atom.get_parent().get_resname() == 'GLU')]
        # atom_neighbors = [atom for atom in atom_neighbors if not (atom.get_name() == 'OE2' and atom.get_parent().get_resname() == 'GLU')]

        # ignore residues on the same chain of res using main-chain atom
        # atom_neighbors = [atom for atom in atom_neighbors if not (atom.get_name() == 'N' and atom.get_parent().get_parent() == res.get_parent())]

        ## filter non-standard residues
        STAND_RES = ['VAL', 'ILE', 'LEU', 'GLU', 'GLN', 'ASP', 'ASN', 'HIS', 'TRP', 'PHE', 'TYR', 'ARG',
                 'LYS', 'SER', 'THR', 'MET', 'ALA', 'GLY', 'PRO', 'CYS']
        for atom in atom_neighbors:
            if atom.get_parent().get_resname() not in STAND_RES:
                atom_neighbors = []

        ## filter same chain
        # for atom in atom_neighbors:
            # if atom.get_parent().get_parent() == res.get_parent():
                # atom_neighbors = []

        ## filter entry containing main_chain O of residues on different chain of res
        # for atom in atom_neighbors:
            # if atom.get_name() == 'N':
                # atom_neighbors = []

        atom_neighbors = list(set(Selection.unfold_entities(atom_neighbors,'R')))
        atom_neighbors = [r for r in atom_neighbors if r != res]

        if len(atom_neighbors) > 0:
            res = res.get_resname()+'_'+str(res.get_id()[1])+'_'+res.get_parent().get_id()
            atom_neighbors = [n.get_resname()+'_'+str(n.get_id()[1])+'_'+n.get_parent().get_id() for n in atom_neighbors]
            neighbors.append((pdb_id,res,atom_neighbors))

    return neighbors


def write_result(res_neighbors, file_suffix):
    pdbs = list(set([pdb_id for pdb_id,_,_ in res_neighbors]))
    with lt.open_file(file_name='pdb_id', inner_dir=file_suffix, file_suffix=file_suffix) as w_f:
        print >> w_f, 'num:', '\t', len(pdbs)
        for pdb in pdbs:
            print >> w_f, pdb

    res_neighbors = sorted(res_neighbors)
    with lt.open_file(file_name='res_neighbors', inner_dir=file_suffix, file_suffix=file_suffix) as w_f:
        for pdb_id, res, neighbors in res_neighbors:
            print >> w_f, '{0:<20}{1:<20}{2}'.format(
                pdb_id, res, '  '.join(neighbors))


"""
ignore water
filter non-standard residues
"""


@lt.run_time
def main():
    null = open(os.devnull, 'w')
    # sys.stderr = null
    """
    """

    res_neighbors = []
    for pdb_f in lt.files_in_dir(sys.argv[-1]):
        f_path,f_name = os.path.split(pdb_f)
        f_name,f_exten = os.path.splitext(f_name)
        if f_exten == '.pdb':
            pdb_id = f_name
            res_neighbors.extend(pdb_neighbors(pdb_f, pdb_id))

    write_result(res_neighbors, file_suffix='original')

    sta = [len(v) for k,v in res_neighbors.iteritems()]
    sta = lt.lis_sta(sta)
    with lt.open_file('phos_neighbors_num') as w_f:
        for k,v in sta.items():
            print >> '{0:<10}{1:<10}'.format(k,v)


if __name__ == "__main__":
    main()

