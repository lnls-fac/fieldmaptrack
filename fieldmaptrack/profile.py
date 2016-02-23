import sys
import numpy as np
from lnls import dialog
import os
import math
import mathphys
import matplotlib.pyplot as plt

_default_dir = '/home/fac_files/data/sirius/'

def read_analysis_file(fname):
    #print('reading file "' + fname + '" ... ')
    with open(fname, 'r') as fp:
        lines = fp.readlines()
    model = []
    store_line = False
    for line in lines:
        if 'magnet_label' in line:
            words = line.strip().split()
            magnet_label = words[1]
        if 'main_monomial' in line:
            words = line.strip().replace(',','').split()
            magnet_type = {'0':['dipole','normal'], '1':['quadrupole','normal'], '2':['sextupole','normal']}[words[3]]
            if 'True' in words[4]: magnet_type[1] = 'skew'
        if line.startswith('beam_energy:'):
            energy = 1e9*float((line.split(':')[1]).strip().split(' ')[0])
        if 'len[m]' in line:
            harmonics = [int(word) for word in line.strip().replace('PolyB(n=','').replace(')','').split()[2:]]
            store_line = True
            continue
        if store_line:
            d = [word.strip() for word in line.strip().split(',')]
            model.append([float(d[i]) for i in range(len(d)-1)])
    return magnet_label, magnet_type, harmonics, energy, np.array(model)

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
    magnet_label, magnet_type, harmonics, energy, model = read_analysis_file(fname)
    return files[0], magnet_label, magnet_type, harmonics, energy, model, read_multipoles(files[0])

def print_wiki_table(magnet_label, magnet_type, harmonics, energy, data, multipoles):

    if magnet_type[1] is not 'normal': raise Exception('skew field not yet implemented')

    m = np.array(multipoles)

    if magnet_type[0] == 'dipole':
        lens, angle, quad, sext = np.array(data[:,0]), np.array(data[:,1]), data[:,3], data[:,4]
        brho = mathphys.beam_optics.beam_rigidity(energy = energy)[0]
        field = -brho * (math.pi/180.0) * (angle / lens)
        # --- table ---
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
        # --- runge-kutta profile ---
        spos, field = m[:,0], m[:,1]
        plt.plot(spos, field, 'r', linewidth=1.5)
        plt.xlim([0,1.1*s[-1]])
        plt.xlabel('pos [mm]'), plt.ylabel('field [T]')
        plt.grid('on')
        plt.title('half segmented model of ' + magnet_label)
        plt.savefig('segmented_model_' + magnet_label + ".svg")
        plt.show()
    elif magnet_type[0] == 'quadrupole':
        quad_idx = harmonics.index(1)
        quad = data[:,quad_idx+2]
        try:
            sext_idx = harmonics.index(2)
            sext = data[:,sext_idx+2]
        except:
            sext = 0 * quad
        lens, angle = np.array(data[:,0]), np.array(data[:,1])
        brho = mathphys.beam_optics.beam_rigidity(energy = energy)[0]
        field = -brho * (math.pi/180.0) * (angle / lens)
        # --- table ---
        print('|-')
        for i in range(len(field)):
            print('| {0:02d} || {1:.4f} || {2:.3f} || {3:+.3e} || {4:+.3e} || {5:+.3e}'.format(i+1, lens[i], angle[i], field[i], quad[i], sext[i]))
            print('|-')
        # --- segmented model plot ---
        lens = lens * 1000
        field = quad
        s,f = [0], [0]
        for i in range(len(lens)):
            s0 = s[-1]
            s.append(s0), f.append(field[i])
            s.append(s0+lens[i]), f.append(field[i])
            s.append(s0+lens[i]), f.append(0)
        plt.fill(s,f, 'lightblue')
        plt.fill([0, s[-1]], [0, 0], 'lightblue')
        # --- runge-kutta profile ---
        spos, field = m[:,0], -m[:,quad_idx+1]/brho
        plt.plot(spos, field, 'r', linewidth=1.5)
        plt.xlim([0,1.1*s[-1]])
        plt.xlabel('pos [mm]'), plt.ylabel('quadrupole strength [1/m²]')
        plt.grid('on')
        plt.title('half segmented model of ' + magnet_label)
        plt.savefig('segmented_model_' + magnet_label + ".svg")
        plt.show()


def run():
    pathdir, magnet_label, magnet_type, harmonics, energy, data, multipoles = get_data(text = 'Select folder containing fieldmap analysis of dipole')
    print_wiki_table(magnet_label, magnet_type, harmonics, energy, data, multipoles)
