#!/usr/bin/env python3 

import fieldmaptrack as fmaptrack
import fieldmaptrack.corrector_analysis as analysis
import numpy as np

config = analysis.Config()

config.config_label          = 'hcm_model2'   # identification label
config.fmap_filename         = '/home/fac_files/data/sirius/bo/magnet_modelling/cm/fieldmaps/2014-11-06 Corretora_Horizontal_Modelo 2_-20_20mm_-300_300mm.txt' # parameter
config.beam_energy           = 3.0          # [GeV] electron beam energy 
config.model_hardedge_length = 100          # [mm]  model hard-edge length of the magnet 
config.model_nominal_angle   = 0.0          # [deg] model nominal deflection angle of the magnet
config.traj_rk_nrpts         = 501          # [mm]  step in s for the 4th-order RK integration
config.traj_init_rx          = 0.0          # [mm]  init rx at center of dipole 
config.traj_center_sagitta_flag = False     # centers trajectory sagitta in good field region of the magnet
config.traj_force_midplane_flag = True      # forces trajectory on midplane (setting ry = py = 0)
config.traj_is_reference_traj   = False     # Rescale field so that nominal deflection is reached. Multipoles are calculated around this ref_traj

config.multipoles_main_monomials     = [0]
config.multipoles_perpendicular_grid = np.linspace(-18,18,57)  # grid of points on perpendicular line to ref trajectory [mm]
config.multipoles_normal_field_fitting_monomials  = (0,1,2,3,4,5,6,7,8,9,10)      # monomials to include in the polynomial fit of multipoles 
config.multipoles_r0                 = 17.5                    # [mm] horizontal position at which polynomial fields are calculated relative to the principal multipole

if __name__ == "__main__":
    
    print('SEXTUPOLE ANALYSIS')
    print('==================')
         
    print('{0:<35s} {1}'.format('label:', config.config_label))
    
    config = analysis.raw_fieldmap_analysis(config)
    config = analysis.trajectory_analysis(config)
    config = analysis.multipoles_analysis(config)
