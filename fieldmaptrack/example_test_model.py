#!/usr/bin/env python3 

import fieldmaptrack
import fieldmaptrack.beam as beam
import fieldmaptrack.track as track
import fieldmaptrack.dipole_analysis as analysis
import fieldmaptrack.multipoles as multipoles
import math

import numpy as np

beam_energy      =  3      # [GeV]
hard_edge_length =  1000   # [mm]
peak_field       = -0.6    # [T]
quad_K           = -0.78   # [1/m^2]
min_rz           =  600    # [mm]
s_step           =  1.0    # [mm]
init_rx          =  0.0    # [mm]
init_ry          =  0.0    # [mm]
init_rz          =  0.0    # [mm]
 
perp_grid        = np.linspace(-5,5,11)
fit_monomials    = [0,1,2,3,4,5] 
 
gradient   = - quad_K * beam.Beam.calc_brho(energy = beam_energy)[0] / 1000.0  # [T/mm]

def dipole(rx, ry, rz):
    if rz >= -hard_edge_length/2.0 and rz <= hard_edge_length/2.0:    
        return (0,peak_field,0)
    else:
        return (0,0,0)

def dipole_with_gradient(rx, ry, rz):
    bx = gradient * ry
    by = peak_field + gradient * rx
    if rz >= -hard_edge_length/2.0 and rz <= hard_edge_length/2.0:    
        return (0,by,0)
    else:
        return (0,0,0)


ebeam    = beam.Beam(energy = beam_energy)
fieldmap = fieldmaptrack.FieldMap(field_function = dipole_with_gradient, rotation = 0.0)
ref_traj = track.Trajectory(beam = ebeam, fieldmap = fieldmap)
new_traj = track.Trajectory(beam = ebeam, fieldmap = fieldmap)

ref_traj.calc_trajectory(min_rz  = min_rz,  s_step  = s_step, 
                         init_rz = init_rz, init_rx = init_rx, init_ry = init_ry)
new_traj.calc_trajectory(min_rz  = min_rz,  s_step  = s_step, 
                         init_rz = init_rz, init_rx = init_rx + 10, init_ry = init_ry)

import matplotlib.pyplot as plt
plt.plot(ref_traj.rz, ref_traj.rx)
plt.plot(new_traj.rz, new_traj.rx)
plt.show()

sf = track.SerretFrenetCoordSystem(ref_traj, point_idx = 500);
idx, alpha, dx = sf.find_intersection(new_traj)
new_traj.px
    

multipoles = multipoles.Multipoles(trajectory = trajectory, perpendicular_grid=perp_grid, normal_field_fitting_monomials=fit_monomials)
multipoles.calc_multipoles(is_ref_trajectory_flag = False)

import matplotlib.pyplot as plt
plt.plot(trajectory.s, multipoles.normal_multipoles[1])
#plt.plot(trajectory.s, trajectory.by)
plt.show()

print('ok')

# 
# 
# config = analysis.Config()
# 
# # raw-field analysis
# config.config_label          = 'bc_model2_controlgap_50mm_modelo_mecanico_inclinado_450urad' # identification label
# config.fmap_filename         = '/home/fac_files/data/sirius/si/magnet_modelling/bc/fieldmaps/2014-10-07_Dipolo_Anel_BC_Modelo2_gap_lateral_50mm_modelo_mecanico_-50_50mm_-2000_2000mm.txt'
# config.beam_energy           = 3.0      # [GeV] electron beam energy 
# config.model_hardedge_length = 828.08   # [mm]  model hard-edge length of the magnet 
# config.model_nominal_angle   =  4.10351 # 2.76654  # [deg] model nominal deflection angle of the magnet
# config.traj_rk_s_step        = 0.1      # [mm]  step in s for the 4th-order RK integration
# config.traj_center_sagitta_flag = False # centers trajectory sagitta in good field region of the magnet
# config.traj_force_midplane_flag = True  # forces trajectory on midplane (setting ry = py = 0)
# config.traj_init_rx             = 0.0   # initial rx
# config.traj_load_filename       = None
# 
# config.multipoles_main_monomials     = [0,1]
# config.multipoles_perpendicular_grid = np.linspace(-10,10,41) # grid of points on perpendicular line to ref trajectory [mm]
# config.multipoles_normal_field_fitting_monomials  = (0,1,2,3,4,5,6,7,8,9,10)    # monomials to include in the polynomial fit of multipoles
# config.multipoles_r0                 = 11.7 # [mm] horizontal position at which polynomial fields are calculated relative to the principal multipole
# 
# config.model_segmentation = 6 * [828.08/6]
# #config.model_segmentation = (106.01, 106.01, 53.005, 53.005, 53.005, 53.005) #828.08,)
# 
# if __name__ == "__main__":
#     
#     print('DIPOLE ANALYSIS')
#     print('===============')
#          
#     print('{0:<35s} {1}'.format('label:', config.config_label))
#     
#     config = analysis.raw_fieldmap_analysis(config)
#     config = analysis.trajectory_analysis(config)
#     #config = analysis.multipoles_analysis(config)
#     #config = analysis.model_analysis(config)
#     
#     #print()
#     #print()
#     #import math
#     #brho = 10.00692271077752
#     #theta_x = (180.0/math.pi) * math.atan(config.traj.px[-1]/config.traj.pz[-1])
#     #print('deflection [deg]: ' + str(theta_x * 2))
#     #print('K [1/m^2]       : ' + str(-config.multipoles.normal_multipoles_integral[1]*2/0.82808/brho))