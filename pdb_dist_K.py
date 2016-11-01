import os
import sys
from collections import defaultdict
import numpy as np

from Bio.PDB import NeighborSearch, PDBParser, Selection

import lt


def pdb_dist(pdb_f, pdb_id):
    structure = PDBParser().get_structure(pdb_id, pdb_f)
    atom_list = Selection.unfold_entities(structure, 'A')
    ns = NeighborSearch(atom_list)
    center = [a for a in structure.get_atoms() if a.get_parent().get_resname(
    ) in ['PTR', 'SEP', 'TPO'] and a.get_name() in ['O1P', 'O2P', 'O3P', 'OH', 'OG', 'OG1']]

    # neighbors = [a for a in structure.get_atoms() if a.get_parent().get_resname(
    # ) in ['ARG'] and a.get_name() in ['NE', 'NH2', 'NH1']]
    neighbors = [a for a in structure.get_atoms() if a.get_parent().get_resname(
    ) in ['LYS'] and a.get_name() in ['NZ']]

    def calc_dist(a, b):
        vector = a.coord - b.coord
        return np.sqrt(np.sum(vector * vector))

    dist = {}
    for c in center:
        for n in neighbors:
            value = calc_dist(c, n)
            key = str(c.get_parent()) + '_' + str(n.get_parent())
            if not key in dist.keys():
                dist[key] = [value]
            else:
                dist[key].append(value)

    for k, v in dist.items():
        dist[k] = min(v)

    dist = [v for k, v in dist.items()]

    return dist


@lt.run_time
def main():
    null = open(os.devnull, 'w')
    sys.stderr = null

    dist = []
    for pdb_f in lt.files_in_dir(sys.argv[-1]):
        pdb_id = pdb_f[-8:-4]
        p_dist = pdb_dist(pdb_f,pdb_id)
        dist += p_dist

    def myround(n):
        n = np.round(n, 1)
        a, b = str(n).split('.')
        n = a + '.' + b[0]
        return float(n)

    dist = sorted(dist)
    dist = [myround(d) for d in dist]
    dist = [d for d in dist if d < 10.0]

    sta = lt.lis_sta(dist)
    with lt.open_file(file_suffix='dist_sta') as w_f:
        for num, freq in sta:
            print >> w_f, '{0:<10}{1}'.format(num, freq)

main()
