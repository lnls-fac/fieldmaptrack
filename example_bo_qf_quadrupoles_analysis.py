#!/usr/bin/env python3

import fieldmaptrack as fmaptrack
import fieldmaptrack.quadrupole_analysis as analysis
import numpy as np

config = analysis.Config()

config.config_label          = 'qf_model4'  # identification label
config.fmap_filename         = '/home/fac_files/data/sirius/bo/magnet_modelling/qf/fieldmaps/2013-08-19 Quadrupolo_Booster_QF_Modelo 4_-32_32mm_-450_450mm.txt' # parameter
config.beam_energy           = 3.0          # [GeV] electron beam energy
config.model_hardedge_length = 200          # [mm]  model hard-edge length of the magnet
config.model_nominal_angle   = 0.0          # [deg] model nominal deflection angle of the magnet
config.traj_rk_nrpts         = 501          # [mm]  step in s for the 4th-order RK integration
config.traj_init_rx          = 0.0          # [mm]  init rx at center of dipole
config.traj_center_sagitta_flag = False     # centers trajectory sagitta in good field region of the magnet
config.traj_force_midplane_flag = True      # forces trajectory on midplane (setting ry = py = 0)
config.traj_is_reference_traj   = False     # Rescale field so that nominal deflection is reached. Multipoles are calculated around this ref_traj

config.multipoles_main_monomials     = [1]
config.multipoles_perpendicular_grid = np.linspace(-18,18,37)  # grid of points on perpendicular line to ref trajectory [mm]
config.multipoles_normal_field_fitting_monomials  = (1,3,5,7,9,11) # monomials to include in the polynomial fit of multipoles
config.multipoles_r0                 = 17.5                    # [mm] horizontal position at which polynomial fields are calculated relative to the principal multipole

if __name__ == "__main__":

    print('QUADRUPOLE ANALYSIS')
    print('===================')

    print('{0:<35s} {1}'.format('label:', config.config_label))

    config = analysis.raw_fieldmap_analysis(config)
    config = analysis.trajectory_analysis(config)
    config = analysis.multipoles_analysis(config)
