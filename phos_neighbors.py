#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import lt
from collections import OrderedDict

def write_result(res_neighbors,filename):
    pro_ids = set([pro for pro,_,_ in res_neighbors])
    with lt.open_file(file_name='hbplus_salt_combine_'+filename+'_pdb') as w_f:
        for p in pro_ids:
            print >> w_f,p

    with lt.open_file(file_name='hbplus_salt_combine_'+filename) as w_f:
        for pro, phos_res, neighbor_res in res_neighbors:
            print >> w_f, '{0:<8}{1:<15}{2}'.format(
                pro, phos_res, ',   '.join(neighbor_res))

def write_sta_result(res_neighbors,filename):

    HASH = {'VAL': 'V', 'ILE': 'I', 'LEU': 'L', 'GLU': 'E', 'GLN': 'Q', 'ASP': 'D', 'ASN': 'N', 'HIS': 'H', 'TRP': 'W', 'PHE': 'F',
            'TYR': 'Y', 'ARG': 'R', 'LYS': 'K', 'SER': 'S', 'THR': 'T', 'MET': 'M', 'ALA': 'A', 'GLY': 'G', 'PRO': 'P', 'CYS': 'C','HOH':'O'}

    neighbor_pdbs = {}
    for pdb, res, neighbors in res_neighbors:
        neighbors = [n.split('_')[0] for n in neighbors]
        neighbors = [HASH.get(n,'*') for n in neighbors]
        neighbors = ''.join(sorted(neighbors))
        if not neighbors in neighbor_pdbs.keys():
            neighbor_pdbs[neighbors] = [(pdb,res)]
        else:
            if not (pdb,res)in neighbor_pdbs[neighbors]:
                neighbor_pdbs[neighbors].append((pdb,res))
    neighbor_pdbs = [(len(neighbors),neighbors, pdbs)
                     for neighbors, pdbs in neighbor_pdbs.iteritems()]
    neighbor_pdbs = sorted(neighbor_pdbs)

    with lt.open_file(file_name='hbplus_salt_combine_'+filename+'_sta') as w_f:
        for _,neighbors,pdbs in neighbor_pdbs:
            print >> w_f,'{0:<20}{1:<}'.format(neighbors,pdbs)


def main():

    res_neighbors = lt.pickle_load(sys.argv[-1])

    #delete water
    res_neighbors_dw = [(pdb,res,[n for n in neighbors if not 'HOH' in n]) for pdb,res,neighbors in res_neighbors]
    res_neighbors_dw = [(pdb,res,neighbors) for pdb,res,neighbors in res_neighbors_dw if len(neighbors) > 0]
    write_result(res_neighbors_dw,'1_delete_water')

    #filter entry containing hetero residues
    res_neighbors_fh = [(pdb,res,[n for n in neighbors if n.split('_')[-1] != 'H']) for pdb,res,neighbors in res_neighbors_dw]
    res_neighbors_fh = [(pdb,res,neighbors) for pdb,res,neighbors in res_neighbors_fh if len(neighbors) > 0]
    write_result(res_neighbors_fh,'2_dw_filter_hetero')

    #special case: phos_res neighboring residues involving main-chain interaction is not considered
    res_neighbors_pm = []
    for pdb,res,neighbors in res_neighbors_fh:
        res_id,res_chain = res.split('_')[1:3]
        new_neighbors = []
        for n in neighbors:
            words=n.split('_')
            if words[3] == 'M' and words[2] == res_chain and words[1] == res_id:
                pass
            else:
                new_neighbors.append(n)
        res_neighbors_pm.append((pdb,res,new_neighbors))
    write_result(res_neighbors_pm,'3_dw_fh_pm')

    #filter_entry containing only same-chain neighbors
    res_neighbors_fs = []
    for pdb,res,neighbors in res_neighbors_pm:
        res_chain = [res.split("_")[2]]
        neighbors_chain = [n.split("_")[2] for n in neighbors]
        if set(res_chain)== set(neighbors_chain):
            pass
        else:
            res_neighbors_fs.append((pdb,res,neighbors))
    write_result(res_neighbors_fs,'4_dw_fh_pm_fs')

    #change main-chain interaction residue as 'GLY'
    res_neighbors_cm = []
    for pdb,res,neighbors in res_neighbors_fs:
        new_neighbors = []
        for n in neighbors:
            words = n.split('_')
            if words[3] == 'M':
                words = ['GLY'] + words[1:3]
            else:
                words = words[0:3]
            new_neighbors.append('_'.join(words))
        res_neighbors_cm.append((pdb,res,new_neighbors))
    write_result(res_neighbors_cm,'4.0_dw_fh_fs_cm')
    write_sta_result(res_neighbors_cm,'4.0_dw_fh_pm_cm')


    #filter entry involves main-chain interaction
    res_neighbors_fm = []
    for pdb,res,neighbors in res_neighbors_fs:
        interaction_type = [n.split('_')[3] for n in neighbors]
        if 'M' in interaction_type:
            pass
        else:
            res_neighbors_fm.append((pdb,res,neighbors))
    write_result(res_neighbors_fm,'5_dw_fh_pm_fs_fm')
    write_sta_result(res_neighbors_fm,'5_dw_fh_pm_fs_fm')


    lt.pickle_dump(res_neighbors_dw,'hbplus_salt_combine_1_dw')
    lt.pickle_dump(res_neighbors_fh,'hbplus_salt_combine_2_dw_fh')
    lt.pickle_dump(res_neighbors_pm,'hbplus_salt_combine_3_dw_fh_pm')
    lt.pickle_dump(res_neighbors_fs,'hbplus_salt_combine_4_dw_fh_pm_fs')
    lt.pickle_dump(res_neighbors_cm,'hbplus_salt_combine_4.0_dw_fh_fs_cm')
    lt.pickle_dump(res_neighbors_fm,'hbplus_salt_combine_5_dw_fh_pm_fs_fm')


if __name__ == "__main__":
    main()
