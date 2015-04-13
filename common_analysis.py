import matplotlib.pyplot as plt
import numpy as np
import fieldmaptrack
import math
from fieldmaptrack.track import Trajectory
import fieldmaptrack.common_analysis

class Config:
    def __init__(self, fname = None):
        """loads parameters from a file"""
        if fname is not None:
            with open(fname, 'r') as fp:
                content = fp.read()
            lines = content.split('\n')
            for line in lines:
                line = line.strip()
                if len(line) == 0: continue
                if line[0] is not '#':
                    attribute, _, value = line.partition(' ')
                    value = value.strip()
                    strcode = 'self.{0} = {1}'.format(attribute, value)
                    exec(strcode)

def get_analysis_symbol(magnet_type):
    """returns analysis module according to magnet type"""
    if magnet_type == 'dipole':
        import fieldmaptrack.dipole_analysis as dipole_analysis
        return dipole_analysis
    if magnet_type == 'quadrupole':
        import fieldmaptrack.quadrupole_analysis as quadrupole_analysis
        return quadrupole_analysis
    if magnet_type == 'sextupole':
        import fieldmaptrack.sextupole_analysis as sextupole_analysis
        return sextupole_analysis
    if magnet_type == 'corrector':
        import fieldmaptrack.corrector_analysis as corrector_analysis
        return corrector_analysis

def raw_fieldmap_analysis(config):
    """does raw fieldmap analysis"""

    # loads fieldmap from file
    config.fmap = fieldmaptrack.FieldMap(config.fmap_filename)

    # plots basic data, longitudinal profile at (x,y) = (0,0)
    try:
        config.config_fig_number += 1
    except:
        config.config_fig_number = 1
    x,y = config.fmap.rz, config.fmap.by[config.fmap.ry_zero][config.fmap.rx_zero,:]
    plt.plot(x,y)
    plt.grid(True)
    plt.xlabel('rz [mm]'), plt.ylabel('by [T]')
    plt.title(config.config_label + '\n' + 'Longitudinal profile of vertical field')
    plt.savefig(config.config_label + '_fig{0:02d}_'.format(config.config_fig_number) + 'by-vs-rz.pdf')
    plt.clf()

    # plots basic data, transversal profile at (y,z) = (0,0)
    try:
        config.config_fig_number += 1
    except:
        config.config_fig_number = 1
    x, y = config.fmap.rx, config.fmap.by[config.fmap.ry_zero][:,config.fmap.rz_zero]
    plt.plot(x,y)
    plt.grid(True)
    plt.xlabel('rx [mm]'), plt.ylabel('by [T]')
    plt.title(config.config_label + '\n' + 'Transverse profile of vertical field')
    plt.savefig(config.config_label + '_fig{0:02d}_'.format(config.config_fig_number) +  'by-vs-rx.pdf')
    plt.clf()

    # calculates missing integrals
    if config.fmap_extrapolation_flag:
        config.fmap.field_extrapolation_analysis(threshold_field_fraction = config.fmap_extrapolation_threshold_field_fraction,
                                        polyfit_exponents = config.fmap_extrapolation_exponents)

    # prints basic raw information on the fieldmap
    print('--- fieldmap ---')
    print(config.fmap)

    return config

def calc_reference_trajectory(config):

    config.traj = Trajectory(beam=config.beam, fieldmap=config.fmap)
    max_z  = config.fmap.rz_max
    if config.traj_rk_nrpts is None:
        config.traj_rk_nrpts = int(math.ceil(max_z * 1.0 / config.traj_rk_s_step))
    else:
        rk_nrpts = config.traj_rk_nrpts
        config.traj_rk_s_step = max_z / (config.traj_rk_nrpts - 1.0)

    config.traj.s_step = config.traj_rk_s_step
    #config.traj.s  = np.array([config.traj.s_step * i for i in range(config.traj_rk_nrpts)])
    config.traj.rx = np.array([0 for i in range(config.traj_rk_nrpts)])
    config.traj.ry = np.array([0 for i in range(config.traj_rk_nrpts)])
    config.traj.rz = np.array([config.traj.s_step * i for i in range(config.traj_rk_nrpts)])
    config.traj.px = np.array([0 for i in range(config.traj_rk_nrpts)])
    config.traj.py = np.array([0 for i in range(config.traj_rk_nrpts)])
    config.traj.pz = np.array([1.0 for i in range(config.traj_rk_nrpts)])
    config.traj.s  = np.array(config.traj.rz)

    # rever!!! XRR
    #config.traj.bx = np.array([0 for i in range(config.traj_rk_nrpts)])
    #config.traj.by = np.array([0 for i in range(config.traj_rk_nrpts)])
    #config.traj.bz = np.array([0 for i in range(config.traj_rk_nrpts)])
    config.traj.bx = [0 for i in range(config.traj_rk_nrpts)]
    config.traj.by = [0 for i in range(config.traj_rk_nrpts)]
    config.traj.bz = [0 for i in range(config.traj_rk_nrpts)]

    return config

def trajectory_analysis(config):

    if config.traj_load_filename is not None:
        # loads trajectory from file
        config.beam = fieldmaptrack.Beam(energy = config.beam_energy)
        config.traj = fieldmaptrack.Trajectory(beam=config.beam, fieldmap=config.fmap)
        config.traj.load(config.traj_load_filename)
    else:
        # calcs trajectory centered in good-field-region
        config.beam = fieldmaptrack.Beam(energy = config.beam_energy)
        config = calc_reference_trajectory(config)


    # prints basic information on the reference trajectory
    # ====================================================
    print('--- trajectory (rz > 0) ---')
    print(config.traj)

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
    #main_monomial = {'corrector':0, 'dipole':0, 'quadrupole':1, 'sextupole':2}[config.magnet_type]
    main_monomial = config.normalization_monomial
    normalization_is_skew = config.normalization_is_skew
    if normalization_is_skew:
        config.multipoles.calc_multipoles_integrals_relative(config.multipoles.skew_multipoles_integral, main_monomial = main_monomial, r0 = config.multipoles_r0, is_skew = True)
    else:
        config.multipoles.calc_multipoles_integrals_relative(config.multipoles.normal_multipoles_integral, main_monomial = main_monomial, r0 = config.multipoles_r0, is_skew = False)

    # saves multipoles to file
    config.multipoles.save('multipoles.txt')

    # prints basic information on multipoles
    # ======================================
    print('--- multipoles on reference trajectory (rz > 0) ---')
    print(config.multipoles)

    # plots normal multipoles
    config = plot_normal_multipoles(config)
    # plots skew multipoles
    config = plot_skew_multipoles(config)

    # plots residual normal field
    config = plot_residual_normal_field(config)
    # plots residual skew field
    config = plot_residual_skew_field(config)

    return config

def plot_residual_field_in_curvilinear_system(config):

    """ NOT TO BE USED! """

    main_monomials = config.multipoles_main_monomials

    r0 = config.multipoles_r0/1000.0
    x = np.linspace(0,1.0 * r0,20)

    # by field reconstructed from fitted polynomials
    dby, by0 = 0*x, 0*x
    n   = config.multipoles.normal_field_fitting_monomials
    m   = config.multipoles.normal_multipoles_integral

    for i in range(len(n)):
        if n[i] not in main_monomials:
            dby += m[i] * (x ** n[i]) # [T.m]
        else:
            #dby += m[i] * (x ** n[i]) # [T.m]
            by0 += m[i] * (x ** n[i]) # [T.m]


    # by field calculated from fieldmap interpolation
    by = 0*x
    traj = Trajectory(beam=config.beam, fieldmap=config.fmap)
    for i in range(len(x)):
        traj.calc_trajectory(init_rx = 9.045 + 1000*x[i],
                             s_step = config.fmap.rz_step/2.0,
                             min_rz = config.fmap.rz[-1],
                             force_midplane = True)
        by[i] = np.trapz(y = traj.by, x = traj.s/1.0e3) # [T.m]
    dby2 = by - 1*by0

    try:
        config.config_fig_number += 1
    except:
        config.config_fig_number = 1

    kick1 = 1e6 * dby / config.beam.brho
    kick2 = 1e6 * dby2 / config.beam.brho
    plt.plot(x*1000,kick1)
    plt.plot(x*1000,kick2)
    #plt.gca().get_yaxis().get_major_formatter().set_powerlimits((0, 0))
    plt.xlabel('x [mm]')
    plt.ylabel('residual integrated field [urad]')
    plt.grid('on')
    #plt.plot(x*1000,dby_r02)
    plt.show()

    return config

def plot_normal_multipoles(config):

    for i in range(len(config.multipoles_normal_field_fitting_monomials)):
        n = config.multipoles_normal_field_fitting_monomials[i]
        try:
             config.config_fig_number += 1
        except:
            config.config_fig_number = 1
        x,y = config.traj.rz, config.multipoles.normal_multipoles[i,:]
        ylabel, title, fname = config.multipoles.get_multipole_labels('normal',n)
        plt.plot(x,y)
        plt.grid(True)
        plt.xlabel('s [mm]'), plt.ylabel(ylabel)
        plt.title(title + ' along trajectory' + '\n' + '(' + config.config_label + ')')
        plt.gca().yaxis.get_major_formatter().set_powerlimits((-3,3))
        plt.savefig(config.config_label + '_fig{0:02d}_'.format(config.config_fig_number) + fname + '.pdf')
        plt.clf()
    return config

def plot_skew_multipoles(config):

    for i in range(len(config.multipoles_skew_field_fitting_monomials)):
        n = config.multipoles_skew_field_fitting_monomials[i]
        try:
             config.config_fig_number += 1
        except:
            config.config_fig_number = 1
        x,y = config.traj.rz, config.multipoles.skew_multipoles[i,:]
        ylabel, title, fname = config.multipoles.get_multipole_labels('skew',n)
        plt.plot(x,y)
        plt.grid(True)
        plt.xlabel('s [mm]'), plt.ylabel(ylabel)
        plt.title(title + ' along trajectory' + '\n' + '(' + config.config_label + ')')
        plt.gca().yaxis.get_major_formatter().set_powerlimits((-3,3))
        plt.savefig(config.config_label + '_fig{0:02d}_'.format(config.config_fig_number) + fname + '.pdf')
        plt.clf()
    return config

def plot_residual_normal_field(config):

    main_monomials = config.normal_multipoles_main_monomials
    if not main_monomials:
        return config

    r0 = config.multipoles_r0/1000.0
    x = np.linspace(0,1.0*r0, 60)

    # by field reconstructed from fitted polynomials
    dby, by_nominal = 0*x, 0*x
    n   = config.multipoles.normal_field_fitting_monomials
    m   = config.multipoles.normal_multipoles_integral
    for i in range(len(n)):
        if n[i] not in main_monomials:
            dby += m[i] * (x ** n[i]) # [T.m]
        else:
            by_nominal += m[i] * (x ** n[i]) # [T.m]

    # by field integration
    s = config.traj.s
    by = np.zeros((len(s),len(x)))
    for i in range(len(s)):
        sf = fieldmaptrack.SerretFrenetCoordSystem(config.traj, i)
        points = sf.get_transverse_line(1000*x)
        field = config.traj.fieldmap.interpolate_set(points)
        by[i,:] = field[1,:]
    intby = 0*x
    for i in range(len(x)):
        intby[i] = np.trapz(y = by[:,i], x = s/1.0e3) # [T.m]
    dby2 = intby - 1*by_nominal


    # by field integration (rawfield)
    iz = list(config.fmap.rz).index(0.0)
    ix = list(config.fmap.rx).index(0.0)
    rz  = config.fmap.rz[iz:]/1000.0
    by  = config.fmap.by[0][ix:,iz:]
    rx  = config.fmap.rx[ix:]/1000.0
    intby = [0] * len(rx)
    for i in range(len(rx)):
        intby[i] = np.trapz(y = by[i,:], x = rz) # [T.m]
    by_nominal_rawfieldgrid = 0 * rx
    for i in range(len(n)):
        if n[i] in main_monomials:
            by_nominal_rawfieldgrid += m[i] * (rx ** n[i]) # [T.m]
    dby3 = np.array(intby) - 1*by_nominal_rawfieldgrid


    # plots curves
    try:
        config.config_fig_number += 1
    except:
        config.config_fig_number = 1
    kick1 = 1e6 * dby / config.beam.brho
    kick2 = 1e6 * dby2 / config.beam.brho
    kick3 = 1e6 * dby3 / config.beam.brho
    plt.plot(x*1000,kick1, label='sum of integrated multipoles')
    plt.plot(x*1000,kick2, label='integrated interpolated field on parallel lines')
    plt.plot(rx*1000,kick3, label='integrated raw field on parallel lines')
    #plt.xlim(1000*x[0],1000*x[-1])
    #plt.gca().get_yaxis().get_major_formatter().set_powerlimits((0, 0))
    plt.xlabel('x [mm]')
    plt.ylabel('residual integrated normal field [urad]')
    plt.grid('on')
    plt.legend(loc='upper left')
    plt.xlim(0,1000*r0)
    plt.savefig(config.config_label + '_fig{0:02d}_'.format(config.config_fig_number) + 'residual_normal_field.pdf')
    #plt.show()
    plt.clf()

    return config

def plot_residual_skew_field(config):

    main_monomials = config.skew_multipoles_main_monomials
    if not main_monomials:
        return config

    r0 = config.multipoles_r0/1000.0
    x = np.linspace(0,1.0*r0, 60)

    # bx field reconstructed from fitted polynomials
    dbx, bx_nominal = 0*x, 0*x
    n   = config.multipoles.skew_field_fitting_monomials
    m   = config.multipoles.skew_multipoles_integral
    for i in range(len(n)):
        if n[i] not in main_monomials:
            dbx += m[i] * (x ** n[i]) # [T.m]
        else:
            bx_nominal += m[i] * (x ** n[i]) # [T.m]

    # by field integration
    s = config.traj.s
    bx = np.zeros((len(s),len(x)))
    for i in range(len(s)):
        sf = fieldmaptrack.SerretFrenetCoordSystem(config.traj, i)
        points = sf.get_transverse_line(1000*x)
        field = config.traj.fieldmap.interpolate_set(points)
        bx[i,:] = field[0,:]
    intbx = 0*x
    for i in range(len(x)):
        intbx[i] = np.trapz(y = bx[:,i], x = s/1.0e3) # [T.m]
    dbx2 = intbx - 1*bx_nominal


    # by field integration (rawfield)
    iz = list(config.fmap.rz).index(0.0)
    ix = list(config.fmap.rx).index(0.0)
    rz  = config.fmap.rz[iz:]/1000.0
    bx  = config.fmap.bx[0][ix:,iz:]
    rx  = config.fmap.rx[ix:]/1000.0
    intbx = [0] * len(rx)
    for i in range(len(rx)):
        intbx[i] = np.trapz(y = bx[i,:], x = rz) # [T.m]
    bx_nominal_rawfieldgrid = 0 * rx
    for i in range(len(n)):
        if n[i] in main_monomials:
            bx_nominal_rawfieldgrid += m[i] * (rx ** n[i]) # [T.m]
    dbx3 = np.array(intbx) - 1*bx_nominal_rawfieldgrid


    # plots curves
    try:
        config.config_fig_number += 1
    except:
        config.config_fig_number = 1
    kick1 = 1e6 * dbx / config.beam.brho
    kick2 = 1e6 * dbx2 / config.beam.brho
    kick3 = 1e6 * dbx3 / config.beam.brho
    plt.plot(x*1000,kick1, label='sum of integrated multipoles')
    plt.plot(x*1000,kick2, label='integrated interpolated field on parallel lines')
    plt.plot(rx*1000,kick3, label='integrated raw field on parallel lines')
    #plt.xlim(1000*x[0],1000*x[-1])
    #plt.gca().get_yaxis().get_major_formatter().set_powerlimits((0, 0))
    plt.xlabel('x [mm]')
    plt.ylabel('residual integrated skew field [urad]')
    plt.grid('on')
    plt.legend(loc='upper left')
    plt.xlim(0,1000*r0)
    plt.savefig(config.config_label + '_fig{0:02d}_'.format(config.config_fig_number) + 'residual_skew_field.pdf')
    #plt.show()
    plt.clf()

    return config

def create_AT_model(config, segmentation):

    s     = np.copy(config.traj.s)
    p     = np.copy(config.multipoles.normal_multipoles)
    nr_monomials = len(p[:,0])
    nr_segments  = len(segmentation)

    seg_border     = np.append([0.0], np.cumsum(segmentation))
    seg_multipoles = np.zeros((nr_monomials, nr_segments))
    model_multipoles_integral = np.zeros((nr_segments, nr_monomials))
    for i in range(nr_segments-1):
        s1, s2 = seg_border[i], seg_border[i+1]
        # interpolates multipolos on borders of segment
        m1, m2 = np.zeros((nr_monomials, 1)), np.zeros((nr_monomials, 1))
        for j in range(nr_monomials):
            m1[j,:], m2[j,:] = np.interp(x = s1, xp=s, fp=p[j,:]), np.interp(x = s2, xp=s, fp=p[j,:])
        # calcs integral of multipoles within segment
        sel = (s >= s1) & (s <= s2)
        seg_s = np.append(np.append([s1], s[sel]), [s2])
        seg_m = np.append(np.append(m1, p[:,sel], axis=1), m2, axis=1)
        model_multipoles_integral[i,:] = np.trapz(x = seg_s/1000, y = seg_m)
    # last segment accumulates the remainder of the integrated multipoles
    s2  = seg_border[-2]
    m2 = np.zeros((nr_monomials, 1))
    for j in range(nr_monomials):
        m2[j,:] = np.interp(x = s2, xp=s, fp=p[j,:])
    sel = (s >= s2)
    seg_s = np.append([s2], s[sel])
    seg_m = np.append(m2, p[:,sel], axis=1)
    model_multipoles_integral[-1,:] = np.trapz(x = seg_s/1000, y = seg_m)

    config.model_multipoles_integral = model_multipoles_integral
    return config

def create_AT_model_orig(config, segmentation):

    s     = np.copy(config.traj.s)
    p     = np.copy(config.multipoles.normal_multipoles)
    nr_monomials = len(p[:,0])
    # interpolates multipoles on segmentation points
    s_seg = np.cumsum(segmentation)
    p_seg = np.zeros((nr_monomials, len(s_seg)))
    for i in range(nr_monomials):
        p_seg[i,:] = np.interp(x=s_seg, xp=s, fp=p[i,:])
    # adds interpolated points to data and sorts
    s = np.append(s, s_seg)
    p = np.append(p, p_seg, axis = 1)
    new_order = np.argsort(s)
    s = s[new_order]
    p = p[:,new_order]

    config.model_segmentation = segmentation
    config.model_multipoles_integral = np.zeros((len(s_seg), len(p[:,0])))
    left_s, right_s = 0, s_seg[0]
    for i in range(len(s_seg)-1):
        sel = (s >= left_s) & (s <= right_s)
        config.model_multipoles_integral[i,:] = np.trapz(x = s[sel]/1000, y = p[:,sel])
        left_s, right_s = right_s, s_seg[i+1]
    sel = (s >= s_seg[-1])
    config.model_multipoles_integral[-1,:] = np.trapz(x = s[sel]/1000, y = p[:,sel])

    return config

def model_analysis(config):

    # creates AT model
    config = create_AT_model(config, config.model_segmentation)

    # adds discrepancy of deflection angle as error in polynomb[0]
    l = np.array(config.model_segmentation) / 1000.0
    m = config.model_multipoles_integral.transpose() / (-config.beam.brho)
    mi = np.sum(m, axis=1)
    if 0 not in config.multipoles.normal_field_fitting_monomials:
        fmap_deflection  = 0.0
        error_deflection = 0.0
    else:
        fmap_deflection    = mi[0]
        nominal_deflection = abs(config.model_nominal_angle/2)*(math.pi/180.0)
        error_deflection = -(nominal_deflection - fmap_deflection)

    # prints info on model
    nr_monomials = len(config.multipoles.normal_field_fitting_monomials)

    monomials = []
    strapp = '{0:^6s} {1:^14s} '
    for i in range(nr_monomials):
        strapp += '{'+'{0}'.format(2+i)+':^14s} '
        monomials.append('PolynomB(n='+'{0:d}'.format(config.multipoles.normal_field_fitting_monomials[i])+')')

    if fmap_deflection != 0.0:
        m[0,:] *= nominal_deflection / fmap_deflection

    print('--- model polynom_b (rz > 0). units: [mm] for length, [rad] for angle and [m],[T] for polynom_b ---')
    print(strapp.format('len[mm]', 'angle[rad]', *monomials))

    fstr = '{0:<7.4f} {1:<+14.06e} '
    for i in range(m.shape[0]):
        fstr += '{'+str(i+2)+':<+14.6e} '
    for i in range(len(l)):
        if 0 not in config.multipoles.normal_field_fitting_monomials:
            val = [l[i]] + [0.0] + list(m[:,i] / l[i])
        else:
            val = [l[i]] + [m[0,i]] + list(m[:,i] / l[i])
            if i == len(l)-1:
                val[2] = error_deflection / l[i]
            else:
                val[2] = 0.0
        print(fstr.format(*val))

    return config

def model_analysis_orig(config):

    # creates AT model
    config = create_AT_model_orig(config, config.model_segmentation)

    # adds discrepancy of deflection angle as error in polynomb[0]
    l = np.array(config.model_segmentation) / 1000.0
    m = config.model_multipoles_integral.transpose() / (-config.beam.brho)
    mi = np.sum(m, axis=1)
    fmap_deflection    = mi[0]
    nominal_deflection = abs(config.model_nominal_angle/2)*(math.pi/180.0)
    error_polynomb = -(nominal_deflection - fmap_deflection)

    # prints info on model
    nr_monomials = len(config.multipoles.normal_field_fitting_monomials)

    monomials = []
    strapp = '{0:^6s} {1:^14s} '
    for i in range(nr_monomials):
        strapp += '{'+'{0}'.format(2+i)+':^14s} '
        monomials.append('PolynomB(n='+'{0:d}'.format(config.multipoles.normal_field_fitting_monomials[i])+')')

    m[0,:] *= nominal_deflection / fmap_deflection

    print('--- model polynom_b (rz > 0). units: [mm] for length, [rad] for angle and [m],[T] for polynom_b ---')
    print(strapp.format('len[mm]', 'angle[rad]', *monomials))

    fstr = '{0:<7.4f} {1:<+14.06e} '
    for i in range(m.shape[0]):
        fstr += '{'+str(i+2)+':<+14.6e} '
    for i in range(len(l)):
        val = [l[i]] + [m[0,i]] + list(m[:,i] / l[i])
        if i == len(l)-1:
            val[2] = error_polynomb / l[i]
        else:
            val[2] = 0.0
        print(fstr.format(*val))
    val = [sum(l)] + [sum(m[0,:])] + [0.0] + list(mi[1:])
    print('--- integrated polynom_b (rz > 0). units: [mm] for length, [rad] for angle and [m],[T] for polynom_b ---')
    #val[1] = nominal_deflection
    val[2] = error_polynomb / val[0]
    print(fstr.format(*val))

    return config
