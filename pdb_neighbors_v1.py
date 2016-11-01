import os
import sys
from collections import defaultdict

from Bio.PDB import NeighborSearch, PDBParser, Selection

import lt

BOND_CUTOFF = 3.5


def pdb_neighbors(pdb_f, pdb_id):
    structure = PDBParser().get_structure(pdb_id, pdb_f)
    atom_list = Selection.unfold_entities(structure, 'A')
    ns = NeighborSearch(atom_list)
    center = [(a, a.get_coord()) for a in structure.get_atoms() if a.get_parent().get_resname(
    ) in ['PTR', 'SEP', 'TPO'] and a.get_name() in ['O1P', 'O2P', 'O3P', 'OH', 'OG', 'OG1']]
    # if there are no phos-atomes, return
    if len(center) == 0:
        return ''
    neighbors = {}
    for a, c in center:
        # set neighbor distance cutoff
        neighbor_list = ns.search(c, BOND_CUTOFF)
        residue_list = list(set(Selection.unfold_entities(neighbor_list, 'R')))
        try:
            neighbors[Selection.unfold_entities(a, 'R')[0]] = [
                x for x in residue_list if not x == Selection.unfold_entities(a, 'R')[0]]
        except:
            continue
    # add residue id and chain id
    neighbors_full = dict([('_'.join([k.get_resname(), str(k.get_id()[1]), str(k.get_parent().get_id())]),
                            ['_'.join([vi.get_resname(), str(vi.get_id()[1]), str(vi.get_parent().get_id())]) for vi in v]) for k, v in neighbors.iteritems()])
    return neighbors_full


@lt.run_time
def main():
    phos_neighbors = {}
    good_pdb = []
    bad_pdb = []
    for pdb_f in lt.files_in_dir(sys.argv[-1]):
        f_path, f_name = os.path.split(pdb_f)
        f_id, f_ex = os.path.splitext(f_name)
        neighbors = pdb_neighbors(pdb_f, f_id)
        print neighbors
        if len(neighbors) > 0:
            for k, v in neighbors.iteritems():
                if len(v) > 0:
                    phos_neighbors[f_id] = neighbors
                    good_pdb.append(f_id)
                else:
                    bad_pdb.append(f_id)

    lt.write_list(good_pdb, 'good_pdb')
    lt.write_list(bad_pdb, 'bad_pdb')
    # delete res_id and clear empty
    phos_neighbors_v1 = {}
    for k, v in phos_neighbors.iteritems():
        phos_neighbors_v1[k] = {}
        for vk, vv in v.iteritems():
            phos_neighbors_v1[k][vk] = []
            vvn = []
            if len(vv) > 0:
                for vi in vv:
                    vvn.append(vi.split('_')[0])
                vvn = sorted(vvn)
                phos_neighbors_v1[k][vk] = vvn

    # delete non-standard residues
    stand_res = ['VAL', 'ILE', 'LEU', 'GLU', 'GLN', 'ASP', 'ASN', 'HIS', 'TRP', 'PHE', 'TYR', 'ARG',
                 'LYS', 'SER', 'THR', 'MET', 'ALA', 'GLY', 'PRO', 'CYS']
    phos_neighbors_v2 = {}
    for k, v in phos_neighbors_v1.iteritems():
        phos_neighbors_v2[k] = {}
        for vk, vv in v.iteritems():
            vvn = []
            for vi in vv:
                if vi in stand_res:
                    vvn.append(vi)
            if len(vvn) > 0:
                vvn = sorted(vvn)
                phos_neighbors_v2[k][vk] = vvn

    # delete uncomplete patterns
    phos_neighbors_v3 = {}
    for k, v in phos_neighbors_v2.iteritems():
        phos_neighbors_v3[k] = {}
        for vk, vv in v.iteritems():
            if len(vv) >= 3:
                phos_neighbors_v3[k][vk] = vv

    def write_phosneighbors(phos_neighbors, ofile):
        keys = phos_neighbors.keys()
        keys = sorted(keys)
        for k in keys:
            print >> ofile, k

        print >> ofile, '*' * 80
        for k in keys:
            v = phos_neighbors[k]
            print >> ofile, k
            for vk, vv in v.iteritems():
                print >> ofile, '    ', vk, '\t', ' '.join(vv)

        print >> ofile, '*' * 80
        phos_sta = {}
        for k, v in phos_neighbors.iteritems():
            for vk, vv in v.iteritems():
                if not tuple(vv) in phos_sta.keys():
                    phos_sta[tuple(vv)] = [k]
                else:
                    phos_sta[tuple(vv)].append(k)
        phos_sta_list = []
        for k, v in phos_sta.iteritems():
            phos_sta_list.append((len(v), k, v))
        phos_sta_list = sorted(phos_sta_list, reverse=True)
        for len_v, k, v in phos_sta_list:
            print >> ofile, '{0:150}{1:<8}{2}'.format(k, len_v, ','.join(v))

    ofile = lt.open_file('phos_neighbors_original')
    write_phosneighbors(phos_neighbors, ofile)
    ofile = lt.open_file('phos_neighbors_sorted')
    write_phosneighbors(phos_neighbors_v1, ofile)
    ofile = lt.open_file('phos_neighbors_clear')
    write_phosneighbors(phos_neighbors_v2, ofile)
    ofile = lt.open_file('phos_neighbors_clean')
    write_phosneighbors(phos_neighbors_v3, ofile)

main()
