#!/usr/bin/env python3 

import fieldmaptrack as fmaptrack
import fieldmaptrack.quadrupole_analysis as analysis
import numpy as np

config = analysis.Config()

config.config_label          = 'qc_model2'   # identification label
config.fmap_filename         = '/home/fac_files/data/sirius/si/magnet_modelling/qc/fieldmaps/2015-01-27 Quadrupolo_Anel_QC_Modelo 2_-12_12mm_-500_500mm.txt' # parameter
config.beam_energy           = 3.0          # [GeV] electron beam energy 
config.model_hardedge_length = 140          # [mm]  model hard-edge length of the magnet 
config.model_nominal_angle   = 0.0          # [deg] model nominal deflection angle of the magnet
config.traj_rk_s_step        = 1.0          # [mm]
config.traj_init_rx          = 0.0          # [mm]  init rx at center of dipole 
config.traj_center_sagitta_flag = False     # centers trajectory sagitta in good field region of the magnet
config.traj_force_midplane_flag = True      # forces trajectory on midplane (setting ry = py = 0)
config.traj_is_reference_traj   = False     # Rescale field so that nominal deflection is reached. Multipoles are calculated around this ref_traj

config.multipoles_main_monomials     = [1]
config.multipoles_perpendicular_grid = np.linspace(-11.5,11.5,41)  # grid of points on perpendicular line to ref trajectory [mm]
config.multipoles_normal_field_fitting_monomials  = (1,5,9,13,17) # monomials to include in the polynomial fit of multipoles 
config.multipoles_r0                 = 11.7                    # [mm] horizontal position at which polynomial fields are calculated relative to the principal multipole

config.model_segmentation            = [140.0/4.0, 140.0/4.0, ]
config.model_segmentation            = [140.0/2.0, ]

if __name__ == "__main__":
    
    print('QUADRUPOLE ANALYSIS')
    print('===================')
         
    print('{0:<35s} {1}'.format('label:', config.config_label))
    
    config = analysis.raw_fieldmap_analysis(config)
    config = analysis.trajectory_analysis(config)
    config = analysis.multipoles_analysis(config)
    config = analysis.model_analysis(config)
