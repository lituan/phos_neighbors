#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import lt
from collections import OrderedDict

def read_pickle():
    for f in lt.files_in_dir(sys.argv[-1]):
        if 'hbplus' in f:
            hbplus_neighbors = lt.pickle_load(f)
        if 'salt' in f:
            salt_neighbors = lt.pickle_load(f)

    hbplus_neighbors_dic = OrderedDict()
    for pdb_id,res,neighbors in hbplus_neighbors:
        hbplus_neighbors_dic[(pdb_id,res)] = neighbors
    salt_neighbors_dic = {}
    for pdb_id,res,neighbors in salt_neighbors:
        salt_neighbors_dic[(pdb_id,res)] = neighbors

    combine_neighbors = []
    for k,v in hbplus_neighbors_dic.iteritems():
        if k in salt_neighbors_dic.keys():
            v.extend(salt_neighbors_dic[k])
            v = set(v)
            combine_neighbors.append((k[0],k[1],v))
        else:
            combine_neighbors.append((k[0],k[1],v))

    return combine_neighbors




def main():
    res_neighbors = read_pickle()

    with lt.open_file(file_name='hbplus_salt_combine_initial') as w_f:
        for pro, phos_res, neighbor_res in res_neighbors:
            print >> w_f, '{0:<8}{1:<15}{2}'.format(
                pro, phos_res, ',   '.join(neighbor_res))

    lt.pickle_dump(res_neighbors, 'hbplus_salt_combine')


if __name__ == "__main__":
    main()
