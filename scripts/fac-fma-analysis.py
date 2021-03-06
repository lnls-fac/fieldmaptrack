#!/usr/bin/env python-sirius

import os
import sys


def help():

    print('NAME')
    print('       fac-fma-analys.py - drive routine used to perform analysis of field maps from 3D magnet models.')
    print()
    print('SYNOPSIS')
    print('       fac-fma-analysis.py [COMMAND]')
    print()
    print('DESCRIPTION')
    print('       Accepts commands (single argument) that specify what type of action to do.')
    print('       Valid commands are (uselly invoked in order):')
    print()
    print('       help')
    print('              prints this text')
    print()
    print('       clean')
    print('              cleans all output files in current directory')
    print()
    print('       edit')
    print('              open gedit and edit input files')
    print()
    print('       run')
    print('              run the complete analysis')
    print()
    print('       rawfield')
    print('              does rawfield analysis')
    print()
    print('       trajectory')
    print('              does trajectory analysis (after rawfield)')
    print()
    print('       multipoles')
    print('              does multipoles analysis (after trajectory)')
    print()
    print('       model')
    print('              does model analysis (after multipoles)')
    print()
    print('       summary')
    print('              prints summary of results and plots analysis results')
    print()
    print('       profile')
    print('              prints and plots segmented models')
    print()
    print('       multifunctional-sextupole')
    print('              prints data on multifunctional sextupoles')
    print()



def profile():

    os.system('fac-fma-profile.py')

def edit():

    os.system('atom -n rawfield.in trajectory.in multipoles.in model.in &')

def clean():

    os.system('rm -rf *.out *.txt *.pdf *.pkl *.svg *~')

def run():

    print('1/6 - cleaning directory...'), os.system('fac-fma-analysis.py clean')
    print('2/6 - rawfield analysis...'), os.system('fac-fma-rawfield.py > rawfield.out')
    print('3/6 - trajectory calculation...'), os.system('fac-fma-trajectory.py > trajectory.out')
    print('4/6 - multipoles calculation...'), os.system('fac-fma-multipoles.py > multipoles.out')
    print('5/6 - model creation...'), os.system('fac-fma-model.py > model.out')
    print('6/6 - analysis summary and visualization...'), summary()

def summary():


    os.system('rm -rf analysis.txt')

    ''' reads and prints summary '''
    try:
        with open('rawfield.out', 'r') as fp:
            content = fp.read()
        print(content)
        os.system('cat rawfield.out >> analysis.txt')
    except:
        pass
    try:
        with open('trajectory.out', 'r') as fp:
            content = fp.read()
        print(content)
        os.system('cat trajectory.out >> analysis.txt')
    except:
        pass
    try:
        with open('multipoles.out', 'r') as fp:
            content = fp.read()
        print(content)
        os.system('cat multipoles.out >> analysis.txt')
    except:
        pass
    try:
        with open('model.out', 'r') as fp:
            content = fp.read()
        print(content)
        os.system('cat model.out >> analysis.txt')
    except:
        pass


    ''' unites all pdf figs '''
    try:
        os.system('ls *fig*.pdf 1> /dev/null 2> /dev/null && rm -rf analysis.pdf')
        os.system('ls *fig*.pdf 1> /dev/null 2> /dev/null && pdfunite *.pdf analysis.pdf')
        os.system('rm -rf *fig*.pdf')
    except:
        pass


    ''' shows summary figure '''
    # try:
    #     os.system('okular analysis.pdf 2> /dev/null')
    # except:
    #     pass


def rawfield():
    os.system('fac-fma-rawfield.py > rawfield.out')

def trajectory():
    os.system('fac-fma-trajectory.py > trajectory.out')

def multipoles():
    os.system('fac-fma-multipoles.py > multipoles.out')

def model():
    os.system('fac-fma-model.py > model.out')

def multifunctional_sextupole():
    os.system('fac-fma-multifunctional-sextupole.py')

if len(sys.argv) != 2:
    help()
if sys.argv[1] == 'help':
    help()
if sys.argv[1] == 'clean':
    clean()
if sys.argv[1] == 'run':
    run()
if sys.argv[1] == 'edit':
    edit()
if sys.argv[1] == 'rawfield':
    rawfield()
if sys.argv[1] == 'trajectory':
    trajectory()
if sys.argv[1] == 'multipoles':
    multipoles()
if sys.argv[1] == 'model':
    model()
if sys.argv[1] == 'summary':
    summary()
if sys.argv[1] == 'profile':
    profile()
if sys.argv[1] == 'multifunctional-sextupole':
    multifunctional_sextupole()
