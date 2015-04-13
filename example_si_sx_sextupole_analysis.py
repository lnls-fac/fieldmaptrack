#!/usr/bin/env python3 

import fieldmaptrack as fmaptrack
import fieldmaptrack.sextupole_analysis as analysis
import numpy as np

config = analysis.Config()

config.config_label          = 'sx_model1'   # identification label
config.fmap_filename         = '/home/fac_files/data/sirius/si/magnet_modelling/sx/fieldmaps/2015-02-03 Sextupolo_Anel_S_Modelo 1_-12_12mm_-500_500mm.txt'
config.beam_energy           = 3.0          # [GeV] electron beam energy 
config.model_hardedge_length = 150          # [mm]  model hard-edge length of the magnet 
config.model_nominal_angle   = 0.0          # [deg] model nominal deflection angle of the magnet
config.traj_rk_s_step        = 1.0          # [mm]
config.traj_init_rx          = 0.0          # [mm]  init rx at center of dipole 
config.traj_center_sagitta_flag = False     # centers trajectory sagitta in good field region of the magnet
config.traj_force_midplane_flag = True      # forces trajectory on midplane (setting ry = py = 0)
config.traj_is_reference_traj   = False     # Rescale field so that nominal deflection is reached. Multipoles are calculated around this ref_traj

config.multipoles_main_monomials     = [2,]
#config.multipoles_perpendicular_grid = np.linspace(-11.5,11.5,41)  # grid of points on perpendicular line to ref trajectory [mm]
config.multipoles_perpendicular_grid = np.linspace(-12,12,49)  # grid of points on perpendicular line to ref trajectory [mm]

config.multipoles_normal_field_fitting_monomials  = (0,2,4,6,8,10,14) # monomials to include in the polynomial fit of multipoles 
config.multipoles_r0                 = 11.7                    # [mm] horizontal position at which polynomial fields are calculated relative to the principal multipole

config.model_segmentation            = [150.0/2.0,] 

if __name__ == "__main__":
    
    print('SEXTUPOLE ANALYSIS')
    print('==================')
         
    print('{0:<35s} {1}'.format('label:', config.config_label))
    
    config = analysis.raw_fieldmap_analysis(config)
    config = analysis.trajectory_analysis(config)
    config = analysis.multipoles_analysis(config)
    config = analysis.model_analysis(config)
