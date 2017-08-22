#!/usr/bin/env python-sirius

import pickle
import fieldmaptrack.common_analysis

''' reads input file '''
import argparse
parser = argparse.ArgumentParser(description='Fieldmap analysis of raw field data')
parser.add_argument('-i', '--input-file', action='store', default = 'rawfield.in')
args = parser.parse_args()

''' analysis raw field data '''
config = fieldmaptrack.common_analysis.Config(args.input_file)
analysis = fieldmaptrack.common_analysis.get_analysis_symbol(config.magnet_type)
config = analysis.raw_fieldmap_analysis(config)

''' saves analysis into pickle file '''
fname = config.config_label + '.pkl'
with open(fname, 'wb') as output:
    pickle.dump(config, output, pickle.HIGHEST_PROTOCOL)
