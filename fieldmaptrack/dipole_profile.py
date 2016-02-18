import sys
import numpy as np
from lnls import dialog
import os
import math
import mathphys
import matplotlib.pyplot as plt

_default_dir = '/home/fac_files/data/sirius/bo/magnet_modelling/b'

def read_analysis_file(fname):
    #print('reading file "' + fname + '" ... ')
    with open(fname, 'r') as fp:
        lines = fp.readlines()
    data = []
    store_line = False
    for line in lines:
        if 'magnet_label' in line:
            words = line.strip().split()
            magnet_label = words[1]
        if line.startswith('beam_energy:'):
            energy = 1e9*float((line.split(':')[1]).strip().split(' ')[0])
        if 'len[m]' in line:
            store_line = True
            continue
        if store_line:
            d = [word.strip() for word in line.strip().split(',')]
            data.append([float(d[i]) for i in range(len(d)-1)])
    return magnet_label, energy, np.array(data)

def read_multipoles(pathdir):
    fname = pathdir + '/multipoles.txt'
    with open(fname, 'r') as fp:
        lines = fp.readlines()
    spos, multipoles = [], []
    for line in lines:
        if '#' in line: continue
        words = line.strip().split()
        m = [float(word) for word in words]
        multipoles.append(m)
    return multipoles

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
    magnet_label, energy, data = read_analysis_file(fname)
    return files[0], magnet_label, energy, data, read_multipoles(files[0])

def print_wiki_table(magnet_label, energy, data, multipoles):

    # --- table ---
    brho = mathphys.beam_optics.beam_rigidity(energy = energy)[0]
    lens, angle, quad, sext = np.array(data[:,0]), np.array(data[:,1]), data[:,3], data[:,4]
    field = -brho * (math.pi/180.0) * (angle / lens)
    print('|-')
    for i in range(len(field)):
        print('| {0:02d} || {1:.4f} || {2:.3f} || {3:+.3f} || {4:+.3f} || {5:+.3f}'.format(i+1, lens[i], angle[i], field[i], quad[i], sext[i]))
        print('|-')

    # --- segmented model plot ---
    lens = lens * 1000
    s,f = [0], [0]
    for i in range(len(lens)):
        s0 = s[-1]
        s.append(s0), f.append(field[i])
        s.append(s0+lens[i]), f.append(field[i])
        s.append(s0+lens[i]), f.append(0)
    plt.fill(s,f, 'lightblue')
    plt.fill([0, s[-1]], [0, 0], 'lightblue')

    # --- dipolar profile ---
    m = np.array(multipoles)
    spos, field = m[:,0], m[:,1]
    plt.plot(spos, field, 'r', linewidth=1.5)

    plt.xlabel('pos [mm]'), plt.ylabel('field [T]')
    plt.grid('on')
    plt.title('half segmented model of ' + magnet_label)
    plt.savefig('segmented_model_' + magnet_label + ".svg")
    #plt.show()


def run():
    pathdir, magnet_label, energy, data, multipoles = get_data(text = 'Select folder containing fieldmap analysis of dipole')
    print_wiki_table(magnet_label, energy, data, multipoles)

