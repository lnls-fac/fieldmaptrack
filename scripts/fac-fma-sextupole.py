#!/usr/bin/env python-sirius

import pickle
import os
import sys

def read_arguments():
    if len(sys.argv) != 3:
        print('NAME')
        print('       fma_sextupole.py [FOLDER1] [FOLDER2]')
        print('')
        print('DESCRIPTION')
        print('       Performs analysis of multifunction sextupoles.')
        print('       This scripts reads the analysis of a sextupole fieldmap with a particular function')
        print('       [FOLDER1] and compares it with that of the pure sextupole [FOLDER2].')
        print('')
        sys.exit(-1)
    else:
        folder1 = sys.argv[1]
        folder2 = sys.argv[2]
    return folder1, folder2

def read_data(folder1, folder2):
    files = os.listdir(folder1)
    for file in files:
        if '.pkl' in file:
            pickle_fname = os.path.join(folder1, file)
            with open(pickle_fname, 'rb') as pickle_file:
                config1 = pickle.load(pickle_file)
            break
    files = os.listdir(folder2)
    for file in files:
        if '.pkl' in file:
            pickle_fname = os.path.join(folder2, file)
            with open(pickle_fname, 'rb') as pickle_file:
                config2 = pickle.load(pickle_file)
            break
    return config1, config2


'''reads input arguments'''
folder1, folder2 = read_arguments()

'''reads analysis data'''
config1, config2 = read_data(folder1, folder2)


'''defines aux structures'''
normal_monomials1 = config1.multipoles_normal_field_fitting_monomials
skew_monomials1   = config1.multipoles_skew_field_fitting_monomials
normal_monomials2 = config2.multipoles_normal_field_fitting_monomials
skew_monomials2   = config2.multipoles_skew_field_fitting_monomials
r0 = config1.multipoles_r0 / 1000.0

'''calcs main integrated multipole'''
main_monomial = config1.multipoles.main_monomial
main_monomial_is_skew = config1.multipoles.main_monomial_is_skew
if main_monomial_is_skew:
    idx = skew_monomials1.index(main_monomial)
    main_multipole = config1.multipoles.skew_multipoles_integral[idx] * (r0**main_monomial)
else:
    idx = normal_monomials1.index(main_monomial)
    main_multipole = config1.multipoles.normal_multipoles_integral[idx] * (r0**main_monomial)


'''print header '''
print()
print('analysis of "' + folder1 + '" compared to "' + folder2 + '"')
print()
print('{0:<35s}: {1} mm'.format('multipoles calculated at r0', 1000*r0))
print('{0:<35s}: {1}'.format('main_monomial', main_monomial))
print('{0:<35s}: {1}'.format('is_skew?', main_monomial_is_skew))
print('{0:<35s}: {1:+.4e} {2}'.format('main integrated multipole', main_multipole, 'T.m**('+str(1-main_monomial)+')'))


'''prints relative normal multipoles'''
print()
for i in range(len(normal_monomials1)):
    monomial = normal_monomials1[i]
    m1 = config1.multipoles.normal_multipoles_integral[i]
    idx = normal_monomials2.index(monomial)
    m2 = config2.multipoles.normal_multipoles_integral[idx]
    multipole = (m1 - m2) * (r0**monomial)
    relative_multipole = multipole / main_multipole
    print('{0:<35s}: {1:+.4e}'.format('normal relative multipole n = {0:2d}'.format(monomial), relative_multipole))

'''prints relative skew multipoles'''
print()
for i in range(len(skew_monomials1)):
    monomial = skew_monomials1[i]
    m1 = config1.multipoles.skew_multipoles_integral[i]
    idx = skew_monomials2.index(monomial)
    m2 = config2.multipoles.skew_multipoles_integral[idx]
    multipole = (m1 - m2) * (r0**monomial)
    relative_multipole = multipole / main_multipole
    print('{0:<35s}: {1:+.4e}'.format('skew relative multipole n = {0:2d}'.format(monomial), relative_multipole))
