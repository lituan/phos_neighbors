#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import lt
import operator

"""
find residues interacting with phos_res side-chain
"""

def read_hb2(hb2_f):
    with open(hb2_f) as o_f:
        pro = hb2_f[-8:-4]
        hb = {}
        hb_lines = {}
        lines = o_f.readlines()
        lines = [line for line in lines if len(line.split()) > 0]
        lines = [line.strip('\r\n') for line in lines]

        for line in lines:
            phos_res = ''
            neighbor_res = ''
            if line[20:23] == 'PTR' and line[34] in ['S', 's'] \
                    or line[20:23] == 'SEP' and line[34] in ['S', 's'] \
                    or line[20:23] == 'TPO' and line[34] in ['S', 's']:
                # format ILE_111_A
                phos_res = line[20:23] + '_' + \
                    str(int(line[15:19])) + '_' + line[14]

                neighbor_res = line[
                    6:9] + '_' + str(int(line[1:5])) + '_' + line[0] + '_' + line[33].capitalize()

            if line[6:9] == 'PTR' and line[33] in ['S', 's'] \
                    or line[6:9] == 'SEP' and line[33] in ['S', 's'] \
                    or line[6:9] == 'TPO' and line[33] in ['S', 's']:
                # format ILE_111_A
                phos_res = line[6:9] + '_' + \
                    str(int(line[1:5])) + '_' + line[0]

                neighbor_res = line[
                    20:23] + '_' + str(int(line[15:19])) + '_' + line[14] + '_' + line[34].capitalize()

            if phos_res != '' and neighbor_res != '':
                if not phos_res in hb.keys():
                    hb[phos_res] = [neighbor_res]
                    hb_lines[phos_res] = [line]
                else:
                    if not neighbor_res in hb[phos_res]:
                        hb[phos_res].append(neighbor_res)
                        hb_lines[phos_res].append(line)

        pro_hbs = []
        for phos_res, neighbor_res in hb.iteritems():
            hb_line = hb_lines[phos_res]
            pro_hbs.append((pro, phos_res, neighbor_res, hb_line))

    return pro_hbs


def write_initial(hbs):
    with lt.open_file(file_name='hbplus_initial') as w_f:
        for pro, phos_res, neighbor_res, _ in hbs:
            print >> w_f, '{0:<8}{1:<15}{2}'.format(
                pro, phos_res, ',   '.join(neighbor_res))
    with lt.open_file(file_name='hbplus_initial_proof') as w_f:
        for pro, phos_res, neighbor_res, lines in hbs:
            print >> w_f, '-' * 80
            print >> w_f, '{0:<8}{1:<15}{2}'.format(
                pro, phos_res, ',   '.join(neighbor_res))
            print >> w_f, '-' * 80
            for line in lines:
                print >> w_f, line


def main():
    hbs = []
    for hb2_f in lt.files_in_dir(sys.argv[-1]):
        if hb2_f[-4:] == '.hb2':
            hbs.extend(read_hb2(hb2_f))

    write_initial(hbs)
    hbs = [(pro, phos_res, neighbor_res)
           for pro, phos_res, neighbor_res, _ in hbs]

    lt.pickle_dump(hbs, 'hbplus')


if __name__ == "__main__":
    main()
