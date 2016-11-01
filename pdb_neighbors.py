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
        atom_neighbors = [atom for atom in atom_neighbors if not (atom.get_name() == 'OD1' and atom.get_parent().get_resname() == 'ASP')]
        atom_neighbors = [atom for atom in atom_neighbors if not (atom.get_name() == 'OD2' and atom.get_parent().get_resname() == 'ASP')]
        atom_neighbors = [atom for atom in atom_neighbors if not (atom.get_name() == 'OE1' and atom.get_parent().get_resname() == 'GLU')]
        atom_neighbors = [atom for atom in atom_neighbors if not (atom.get_name() == 'OE2' and atom.get_parent().get_resname() == 'GLU')]

        # ignore residues on the same chain of res using main-chain atom
        atom_neighbors = [atom for atom in atom_neighbors if not (atom.get_name() == 'N' and atom.get_parent().get_parent() == res.get_parent())]

        ## filter non-standard residues
        STAND_RES = ['VAL', 'ILE', 'LEU', 'GLU', 'GLN', 'ASP', 'ASN', 'HIS', 'TRP', 'PHE', 'TYR', 'ARG',
                 'LYS', 'SER', 'THR', 'MET', 'ALA', 'GLY', 'PRO', 'CYS']
        for atom in atom_neighbors:
            if atom.get_parent().get_resname() not in STAND_RES:
                atom_neighbors = []

        ## filter same chain
        for atom in atom_neighbors:
            if atom.get_parent().get_parent() == res.get_parent():
                atom_neighbors = []

        ## filter entry containing main_chain O of residues on different chain of res
        for atom in atom_neighbors:
            if atom.get_name() == 'N':
                atom_neighbors = []

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

    # do not distinguish between S and T
    # do not distinguish between G A V I L F P C M , because these residues all
    # use main chain to form hydrogen-bond
    HASH = {'VAL': 'G', 'ILE': 'G', 'LEU': 'G', 'GLU': 'G', 'GLN': 'Q', 'ASP': 'G', 'ASN': 'N', 'HIS': 'H', 'TRP': 'W', 'PHE': 'G',
            'TYR': 'Y', 'ARG': 'R', 'LYS': 'K', 'SER': 'S', 'THR': 'T', 'MET': 'G', 'ALA': 'G', 'GLY': 'G', 'PRO': 'G', 'CYS': 'G','HOH':'O'}
    neighbor_pdbs = {}
    for pdb_id, res, neighbors in res_neighbors:
        neighbors = [n.split('_')[0] for n in neighbors]
        neighbors = [HASH.get(n,'*') for n in neighbors]
        neighbors = ''.join(sorted(neighbors))
        if not neighbors in neighbor_pdbs.keys():
            neighbor_pdbs[neighbors] = [pdb_id]
        else:
            if not pdb_id in neighbor_pdbs[neighbors]:
                neighbor_pdbs[neighbors].append(pdb_id)
    neighbor_pdbs = [(len(pdbs), neighbors, pdbs)
                     for neighbors, pdbs in neighbor_pdbs.iteritems()]
    neighbor_pdbs = sorted(neighbor_pdbs, reverse=True)
    with lt.open_file(file_name='neighbor_pdbs', inner_dir=file_suffix, file_suffix=file_suffix) as w_f:
        for _, neighbors, pdbs in neighbor_pdbs:
            pdbs = [pdb for pdb in pdbs]
            print >> w_f, '{0:<20}{1}'.format(neighbors, ','.join(pdbs))

"""
first, delete empty
delete water
filter non-standard residues
filter same-chain residues
fifth, get four-neighbors pdb
"""


@lt.run_time
def main():
    null = open(os.devnull, 'w')
    # sys.stderr = null
    """
    write original center:neighbors
    write pdb_ids
    write neighbors:pdbs
    """

    res_neighbors = []
    for pdb_f in lt.files_in_dir(sys.argv[-1]):
        f_path,f_name = os.path.split(pdb_f)
        f_name,f_exten = os.path.splitext(f_name)
        if f_exten == '.pdb':
            pdb_id = f_name
            res_neighbors.extend(pdb_neighbors(pdb_f, pdb_id))

    write_result(res_neighbors, file_suffix='original')

    #neighbors = 1
    res_neighbors1 = [(pdb_id, res, neighbors) for pdb_id,
                       res, neighbors in res_neighbors if len(neighbors) == 1]
    write_result(res_neighbors1, file_suffix='1_neighbor')
    #neighbors = 2
    res_neighbors2 = [(pdb_id, res, neighbors) for pdb_id,
                       res, neighbors in res_neighbors if len(neighbors) == 2]
    write_result(res_neighbors2, file_suffix='2_neighbor')
    #neighbors = 3
    res_neighbors3 = [(pdb_id, res, neighbors) for pdb_id,
                       res, neighbors in res_neighbors if len(neighbors) == 3]
    write_result(res_neighbors3, file_suffix='3_neighbor')
    #neighbors = 4
    res_neighbors4 = [(pdb_id, res, neighbors) for pdb_id,
                       res, neighbors in res_neighbors if len(neighbors) == 4]
    write_result(res_neighbors4, file_suffix='4_neighbor')
    #neighbors = 5
    res_neighbors5 = [(pdb_id, res, neighbors) for pdb_id,
                       res, neighbors in res_neighbors if len(neighbors) == 5]
    write_result(res_neighbors5, file_suffix='5_neighbor')
    #neighbors = 6
    res_neighbors6 = [(pdb_id, res, neighbors) for pdb_id,
                       res, neighbors in res_neighbors if len(neighbors) == 6]
    write_result(res_neighbors6, file_suffix='6_neighbor')
    #neighbors > 6
    res_neighbors7 = [(pdb_id, res, neighbors) for pdb_id,
                       res, neighbors in res_neighbors if len(neighbors) > 6]
    write_result(res_neighbors7, file_suffix='7_neighbor')

    num1 = len(set([pdb_id for pdb_id, _, _ in res_neighbors1]))
    num2 = len(set([pdb_id for pdb_id, _, _ in res_neighbors2]))
    num3 = len(set([pdb_id for pdb_id, _, _ in res_neighbors3]))
    num4 = len(set([pdb_id for pdb_id, _, _ in res_neighbors4]))
    num5 = len(set([pdb_id for pdb_id, _, _ in res_neighbors5]))
    num6 = len(set([pdb_id for pdb_id, _, _ in res_neighbors6]))
    num7 = len(set([pdb_id for pdb_id, _, _ in res_neighbors7]))

    sta = [num1, num2, num3, num4, num5, num6, num7]
    with lt.open_file(file_suffix='sta') as w_f:
        print >> w_f, '{0:<20}{1}'.format('neighbors num', 'pdb_num')
        for i, num in enumerate(sta):
            print >> w_f, '{0:<20}{1}'.format(i + 1, num)
        print >> w_f,'{0:<20}{1}'.format('total',sum(sta))


main()

