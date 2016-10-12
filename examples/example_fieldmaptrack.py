#!/usr/bin/env python3

import fieldmaptrack
import numpy
import matplotlib.pyplot as plt

# --- load fieldmap from file
fname = '2014-09-18_Dipolo_Booster_BD_Modelo_6_-80_35mm_-1000_1000mm.txt'
fmap = fieldmaptrack.FieldMap(fname = fname)

# --- print a summary of the fieldmap data
print('Fieldmap summary')
print('================')
print('{0:<25s}: {1:s}'.format('fieldmap label',fmap.fieldmap_label))
print('{0:<25s}: {1:s}'.format('timestamp',fmap.timestamp))
print('{0:<25s}: {1:s}'.format('current[A]',fmap.current))
print('{0:<25s}: {1:.1f}'.format('magnet length[mm]',fmap.length))
print('{0:<25s}: [{1:+.1f},{2:+.1f}], {3:d} points (step {4:.2f} mm)'.format('minx[mm]|maxx[mm]|nrptsx',fmap.rx_min,fmap.rx_max,fmap.rx_nrpts, fmap.rx_step))
print('{0:<25s}: [{1:+.1f},{2:+.1f}], {3:d} points (step {4:.2f} mm)'.format('minz[mm]|maxz[mm]|nrptsz',fmap.rz_min,fmap.rz_max,fmap.rz_nrpts, fmap.rx_step))
print('{0:<25s}: {1:+.1f}|{2:+.1f}'.format('min_bx|max_bx [T]', numpy.min(fmap.bx[0]), numpy.max(fmap.bx[0])))
print('{0:<25s}: {1:+.1f}|{2:+.1f}'.format('min_by|max_by [T]', numpy.min(fmap.by[0]), numpy.max(fmap.by[0])))
print('{0:<25s}: {1:+.1f}|{2:+.1f}'.format('min_bz|max_bz [T]', numpy.min(fmap.bz[0]), numpy.max(fmap.bz[0])))

# --- field interpolation (transverse field profile @ z = 0)
x = numpy.linspace(fmap.rx_min,fmap.rx_max,fmap.rx_nrpts)
points = numpy.array([(vx,0,0) for vx in x])
fields = fmap.interpolate_set(points.T)
by = fields[1,:]
plt.plot(x,by)
plt.xlabel('posx [mm]'), plt.ylabel('by [T]')
plt.show()

# --- field interpolation (longitudinal field profile @ x = 0)
z = numpy.linspace(fmap.rz_min,fmap.rz_max,fmap.rz_nrpts)
points = numpy.array([(0,0,vz) for vz in z])
fields = fmap.interpolate_set(points.T)
by = fields[1,:]
plt.plot(z,by)
plt.xlabel('posz [mm]'), plt.ylabel('by [T]')
plt.show()

# --- define electron beam and trajectory object
beam_energy = 3 # [GeV]
ebeam = fieldmaptrack.Beam(energy = beam_energy)
traj = fieldmaptrack.Trajectory(ebeam, fmap)

# --- do runge-kutta traj calculation
traj_init_rx = 9.045
traj_rk_s_step = 0.1
init_ry, init_rz = 0,0
min_rz = fmap.rz_max
traj.calc_trajectory(init_rx = traj_init_rx, initry = init_ry, init_rz = init_rz,
                     min_rz = min_rz, s_step = traj_rk_s_step)

# --- plot calculated trajectory
plt.plot(traj.rz, traj.rx)
plt.xlabel('rz[mm]'), plt.ylabel('rx[mm]')
plt.show()

# --- calculate multipoles around trajectory
perpendicular_grid = numpy.linspace(-12,12,65)
multipoles_normal_field_fitting_monomials = (0,1,2,3,4,5,6) # 0 - dipoles, 1 - quadrupole, etc
multipoles_skew_field_fitting_monomials = ()

multipoles = fieldmaptrack.Multipoles(trajectory = traj,
                                      perpendicular_grid = perpendicular_grid,
                                      normal_field_fitting_monomials = multipoles_normal_field_fitting_monomials,
                                      skew_field_fitting_monomials = multipoles_skew_field_fitting_monomials
                                     )
multipoles.calc_multipoles()

# --- plot  multipoles
normal_dipole = multipoles.normal_multipoles[0,:]
normal_quadrupole = multipoles.normal_multipoles[1,:]
normal_sextupole = multipoles.normal_multipoles[2,:]

plt.plot(traj.s, normal_dipole)
plt.xlabel('s[mm]'), plt.ylabel('normal dipole [T]')
plt.show();

plt.plot(traj.s, normal_quadrupole)
plt.xlabel('s[mm]'), plt.ylabel('normal quadrupole [T/m]')
plt.show();

plt.plot(traj.s, normal_sextupole)
plt.xlabel('s[mm]'), plt.ylabel('normal sextupole [T/m^2]')
plt.show(); 
