#!/usr/bin/env python-sirius

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
        if 'r0_for_relative_multipoles' in line:
            words = line.split()
            r0 = float(words[1])
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

            data[n] = (integ_normal_multipole, integ_skew_multipole)
    return data, r0

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

def print_difference(fname, data0, data1, order, mtype, r0):

    t = {'normal':0,'skew':1}[mtype]

    r0 = r0/1000.0
    mmult = data1[order][t]*(r0**order)
    m = {}
    dm = {}
    for h in data1.keys():
        dm[h] = (data1[h][0] - data0[h][0], data1[h][1] - data0[h][1]) 
        m[h] = ((data1[h][0] - data0[h][0])*(r0**h)/mmult, (data1[h][1] - data0[h][1])*(r0**h)/mmult)


    print(fname + ':')
    header = 'harm: normal0     skew0       | normal1     skew1       | dnormal     dskew       | rel_normal_dif_r0 rel_skew_dif_r0'
    print('-'*(1+len(fname)))
    print(header)
    print('-'*len(header))
    for h in data1.keys():
        n1, s1 = data1[h]
        n2, s2 = m[h]
        n3, s3 = dm[h]
        try:
            n0, s0 = data0[h]
        except:
            n0, s0 = 0.0, 0.0
        print('n={0:02d}: {1:+.4e} {2:+.4e} | {3:+.4e} {4:+.4e} | {5:+.4e} {6:+.4e} | {7:+.4e} {8:+.4e}'.format(h,n0,s0,n1,s1,n3,s3,n2,s2))

def run():
    fname, data = get_data(text = 'Select folder containing fieldmap analysis of sextupolar-only excitation')
    data_sx, r0_0 = data
    fname, data = get_data(text = 'Select folder containing fieldmap analysis of an additional excitation')
    data_m, r0_1 = data
    if r0_0 != r0_1: raise Exception('differete values for r0!')
    ok, order = dialog.input_dialog('please define main multipole order (0:dipole|1:quadrupole|...)',def_answer='0',name='main multipole order')
    if not ok: return
    ok, mtype = dialog.input_dialog('please define type of main multipole (normal|skew)',def_answer='normal',name='main multipole type')
    if not ok: return
    print_difference(fname, data_sx, data_m, int(order[0]), mtype[0], r0_0)

if __name__ == "__main__":

    run()
