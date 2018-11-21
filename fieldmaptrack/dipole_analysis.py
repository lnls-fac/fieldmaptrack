import matplotlib.pyplot as plt
import numpy as np
import fieldmaptrack
import math
from   fieldmaptrack.common_analysis import *
import fieldmaptrack.track as track
import sys

class Config:

    def __init__(self):

        self.config_label = None
        self.config_timestamp = None
        self.fmap_filename = None
        self.model_hardedge_length = None
        self.traj_rk_s_step = None
        self.traj_rk_length = None
        self.traj_rk_nrpts = None
        self.traj_rk_min_rz = None
        self.traj_save = True
        self.traj_load_filename = None
        self.traj_init_rx = None
        self.traj_center_sagitta_flag = True
        self.traj_force_midplane_flag = True
        self.traj_is_reference_traj = False
        self.multipoles_r0 = None
        self.model_segmentation = None
        self.model_multipoles_integral = None

def calc_reference_trajectory_good_field_region(config):
    # calcs reference trajectory
    # ==========================
    # if it is the case, iterates the calculation of the trajectory
    # until its sagitta is centered in the good field region. This is
    # accomplished by changing the initial radial position of the trajectory
    # at the longitudinal center of the dipole.
    half_dipole_length = config.fmap.length / 2.0
    config.traj = fieldmaptrack.Trajectory(beam=config.beam, fieldmap=config.fmap)
    if config.traj_init_rx is not None:
        init_rx = config.traj_init_rx
    else:
        init_rx = 0.0
    if hasattr(config, 'traj_init_rz'):
        init_rz = config.traj_init_rz
    else:
        config.traj_init_rz = init_rz = 0.0
    if hasattr(config, 'traj_init_px'):
        init_px = config.traj_init_px * (math.pi/180.0)
    else:
        config.traj_init_px = init_px = 0.0
    init_ry = 0.0
    init_py = 0.0
    init_pz = math.sqrt(1.0 - init_px**2 - init_py**2)
    if config.traj_rk_s_step > 0.0:
        rk_min_rz = max(config.fmap.rz)
    else:
        rk_min_rz = abs(min(config.fmap.rz))
    # rk_min_rz = config.fmap.rz[-1]
    while True:
        config.traj.calc_trajectory(init_rx=init_rx, init_ry=init_ry, init_rz=init_rz,
                                    init_px=init_px, init_py=init_py, init_pz=init_pz,
                                    s_step         = config.traj_rk_s_step,
                                    s_length       = config.traj_rk_length,
                                    s_nrpts        = config.traj_rk_nrpts,
                                    min_rz         = rk_min_rz,
                                    force_midplane = config.traj_force_midplane_flag)
        config.traj_sagitta = config.traj.calc_sagitta(half_dipole_length)

        if not config.traj_center_sagitta_flag:
            break
        new_init_rx = config.traj_sagitta / 2.0
        change_init_rx = new_init_rx - init_rx
        if abs(change_init_rx) < 0.001:
            break
        else:
            init_rx = new_init_rx

    return config

def calc_reference_trajectory(config):

    nominal_deflection = -abs(config.model_nominal_angle/2.0)
    deflection = 0.0

    # search an upper energy for which deflection is lower than nominal
    energy1 = 1.005 * config.beam_energy
    while True:
        config.beam = fieldmaptrack.Beam(energy = energy1)
        config = calc_reference_trajectory_good_field_region(config)
        deflection1 = (180.0/math.pi)*math.atan(config.traj.px[-1]/config.traj.pz[-1])
        if deflection1 > nominal_deflection:
            break
        energy1 *= 1.005
    # search a lower energy for which deflection is higher than nominal
    energy2 = 0.995 * config.beam_energy
    while True:
        config.beam = fieldmaptrack.Beam(energy = energy2)
        config = calc_reference_trajectory_good_field_region(config)
        deflection2 = (180.0/math.pi)*math.atan(config.traj.px[-1]/config.traj.pz[-1])
        if deflection2 < nominal_deflection:
            break
        energy2 *= 0.995

    # does bisection
    iter = 0
    while True:
        energy = 0.5 * (energy1 + energy2)
        config.beam = fieldmaptrack.Beam(energy = energy)
        config = calc_reference_trajectory_good_field_region(config)
        deflection = (180.0/math.pi)*math.atan(config.traj.px[-1]/config.traj.pz[-1])
        if deflection > nominal_deflection:
            energy1 = energy
        else:
            energy2 = energy
        error = abs(100*(deflection - nominal_deflection)/nominal_deflection)
        if error < 0.001:
            break
        iter += 1
        if (iter > 20):
            raise Exception('max number of iterations while bisecting for correct deflection angle')

    energy = 0.5 * (energy1 + energy2)
    config.beam = fieldmaptrack.Beam(energy = energy)
    config.traj.bx = (config.beam_energy / energy) * config.traj.bx
    config.traj.by = (config.beam_energy / energy) * config.traj.by
    config.traj.bz = (config.beam_energy / energy) * config.traj.bz
    return config

def trajectory_analysis(config):

    if config.traj_load_filename is not None:
        # loads trajectory from file
        half_dipole_length = config.fmap.length / 2.0
        config.beam = fieldmaptrack.Beam(energy = config.beam_energy)
        config.traj = fieldmaptrack.Trajectory(beam=config.beam, fieldmap=config.fmap)
        config.traj.load(config.traj_load_filename)
        config.traj_sagitta = config.traj.calc_sagitta(half_dipole_length)
        config.traj_init_rz = config.traj.rz[0]
    else:
        if config.traj_is_reference_traj:
            # rescale field so that nominal deflection is reached
            config = calc_reference_trajectory(config)
        else:
            # calcs trajectory centered in good-field-region
            config.beam = fieldmaptrack.Beam(energy = config.beam_energy)
            config = calc_reference_trajectory_good_field_region(config)


    # prints basic information on the reference trajectory
    # ====================================================
    print('--- trajectory (rz > {0} mm) ---'.format(config.traj_init_rz))
    print(config.traj)
    print('{0:<35s} {1} mm'.format('sagitta:', config.traj_sagitta))

    if not config.interactive_mode:
        # saves trajectory in file
        config.traj.save(filename='trajectory.txt')
        # saves field on trajectory in file
        config.traj.save_field(filename='field_on_trajectory.txt')

    return config

def multipoles_analysis(config):

    # calcs multipoles around reference trajectory
    # ============================================
    config.multipoles = fieldmaptrack.Multipoles(trajectory=config.traj,
                                         perpendicular_grid=config.multipoles_perpendicular_grid,
                                         normal_field_fitting_monomials=config.multipoles_normal_field_fitting_monomials,
                                         skew_field_fitting_monomials=config.multipoles_skew_field_fitting_monomials)
    config.multipoles.calc_multipoles(is_ref_trajectory_flag = False)
    config.multipoles.calc_multipoles_integrals()
    config.multipoles.calc_multipoles_integrals_relative(config.multipoles.normal_multipoles_integral, main_monomial = 0, r0 = config.multipoles_r0, is_skew = False)

    # calcs effective length

    # main_monomial = config.normalization_monomial
    # monomials = config.multipoles.normal_field_fitting_monomials
    # idx_n = monomials.index(main_monomial)
    # idx_z = list(config.traj.s).index(0.0)
    # main_multipole_center = config.multipoles.normal_multipoles[idx_n,idx_z]
    # config.multipoles.effective_length = config.multipoles.normal_multipoles_integral[idx_n] / main_multipole_center

    main_monomial = config.normalization_monomial
    monomials = config.multipoles.normal_field_fitting_monomials
    idx_n = monomials.index(main_monomial)

    if hasattr(config, 'hardedge_half_region'):
        sel = config.traj.s < config.hardedge_half_region
        s = config.traj.s[sel]
        field = config.multipoles.normal_multipoles[idx_n,sel]
        integrated_field = np.trapz(field, s)
        hardedge_field = integrated_field / config.hardedge_half_region
        config.multipoles.effective_length = config.multipoles.normal_multipoles_integral[idx_n] / hardedge_field
    else:
        idx_z = list(config.traj.s).index(0.0)
        main_multipole_center = config.multipoles.normal_multipoles[idx_n,idx_z]
        config.multipoles.effective_length = config.multipoles.normal_multipoles_integral[idx_n] / main_multipole_center

    # saves multipoles to file
    if not config.interactive_mode: config.multipoles.save('multipoles.txt')

    # prints basic information on multipoles
    # ======================================
    print('--- multipoles on reference trajectory (rz > 0) ---')
    print(config.multipoles)

    if not config.interactive_mode:
        # plots normal multipoles
        config = plot_normal_multipoles(config)

        # plots skew multipoles
        config = plot_skew_multipoles(config)

        # plots residual normal field
        #config = plot_residual_field_in_curvilinear_system(config)
        config = plot_residual_normal_field(config)
        # plots residual skew field
        config = plot_residual_skew_field(config)

    return config
