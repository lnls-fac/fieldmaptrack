#!/usr/bin/env python3

import fieldmaptrack as fmaptrack
import fieldmaptrack.dipole_analysis as analysis
import numpy as np

config = analysis.Config()

config.config_label          = 'b_model6'   # identification label
config.fmap_filename         = '/home/fac_files/data/sirius/bo/magnet_modelling/b/fieldmaps/2014-09-18_Dipolo_Booster_BD_Modelo_6_-80_35mm_-1000_1000mm.txt' # parameter
config.beam_energy           = 3.0          # [GeV] electron beam energy
config.model_hardedge_length = 1152         # [mm]  model hard-edge length of the magnet
config.model_nominal_angle   = 7.2          # [deg] model nominal deflection angle of the magnet
config.traj_rk_s_step        = 1.0          # [mm]  step in s for the 4th-order RK integration
config.traj_init_rx          = 9.045        # [mm]  init rx at center of dipole
config.traj_center_sagitta_flag = False     # centers trajectory sagitta in good field region of the magnet
config.traj_force_midplane_flag = True      # forces trajectory on midplane (setting ry = py = 0)
config.traj_is_reference_traj   = False     # Rescale field so that nominal deflection is reached. Multipoles are calculated around this ref_traj

config.multipoles_main_monomials     = [0,1,2]
config.multipoles_perpendicular_grid = np.linspace(-16,16,65)            # grid of points on perpendicular line to ref trajectory [mm]
config.multipoles_normal_field_fitting_monomials  = (0,1,2,3,4,5,6)      # monomials to include in the polynomial fit of multipoles
config.multipoles_r0                 = 17.5                              # [mm] horizontal position at which polynomial fields are calculated relative to the principal multipole

config.model_segmentation = (196, 192, 158, 34, 30, 158, 1)
#config.model_segmentation = (196+192+158+34+30+158+1,)

#if __name__ == "__main__":

print('DIPOLE ANALYSIS')
print('===============')

print('{0:<35s} {1}'.format('label:', config.config_label))

config = analysis.raw_fieldmap_analysis(config)
config = analysis.trajectory_analysis(config)
config = analysis.multipoles_analysis(config)
config = analysis.model_analysis(config)
