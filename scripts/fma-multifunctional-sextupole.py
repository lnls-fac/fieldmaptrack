#!/usr/bin/env python3

import sys
from lnls import dialog
import os

_default_dir = '/home/fac_files/data/sirius/si/magnet_modelling/si-s'

def read_analysis_file(fname):
    #print('reading file "' + fname + '" ... ')
    with open(fname, 'r') as fp:
        lines = fp.readlines()
    data = {}
    for line in lines:
        if line.startswith('n='):
            words = line.split()
            n = int(words[0][2:4])
            try:
                integ_normal_multipole = float(words[2])
            except:
                integ_normal_multipole = 0
            try:
                integ_skew_multipole = float(words[6])
            except:
                integ_skew_multipole = 0
            #print(n,integ_normal_multipole,integ_skew_multipole)
            data[n] = (integ_normal_multipole, integ_skew_multipole)
    return data

def get_data(text):
    ok, files = dialog.directories_dialog(path=_default_dir, name=text)
    if not ok: sys.exit()
    if len(files) > 1:
        print('please select only one file!')
        sys.exit(1)
    fname = files[0] + '/analysis.txt'
    if not os.path.exists(fname):
        print('selected directory does not contain ')
        sys.exit(1)

    return files[0], read_analysis_file(files[0] + '/analysis.txt')

def print_difference(fname, data0, data1):
    print(fname + ':')
    header = 'harm: normal0    skew0      | normal1    skew1      | normal_dif skew_dif  '
    print('-'*(1+len(fname)))
    print(header)
    print('-'*len(header))
    for h in data1.keys():
        n1, s1 = data1[h]
        try:
            n0, s0 = data0[h]
        except:
            n0, s0 = 0.0, 0.0
        print('n={0:02d}: {1:+.3e} {2:+.3e} | {3:+.3e} {4:+.3e} | {5:+.3e} {6:+.3e}'.format(h,n0,s0,n1,s1,n1-n0,s1-s0))

def run():
    fname, data_sx = get_data(text = 'Select folder containing fieldmap analysis of sextupolar-only excitation')
    while True:
        fname, data = get_data(text = 'Select folder containing fieldmap analysis of an additional excitation')
        print_difference(fname, data_sx, data)

if __name__ == "__main__":
    run()
