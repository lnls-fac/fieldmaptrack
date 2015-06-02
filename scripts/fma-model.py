#!/usr/bin/env python3 

import pickle
import fieldmaptrack.common_analysis

''' reads input file '''
import argparse
parser = argparse.ArgumentParser(description='Fieldmap analysis of multipoles')
parser.add_argument('-i', '--input-file', action='store', default = 'model.in')
args = parser.parse_args()
config_input = fieldmaptrack.common_analysis.Config(args.input_file)

''' imports multipoles analysis '''
pickle_fname = config_input.config_label + '.pkl'
with open(pickle_fname, 'rb') as input:
    config = pickle.load(input)
for p in config_input.__dict__:
    setattr(config, p, getattr(config_input, p))

''' does model analysis '''
analysis = fieldmaptrack.common_analysis.get_analysis_symbol(config.magnet_type)
config = analysis.model_analysis(config)

''' saves analysis into pickle file '''
import pickle
with open(pickle_fname, 'wb') as output:
    pickle.dump(config, output, pickle.HIGHEST_PROTOCOL)
