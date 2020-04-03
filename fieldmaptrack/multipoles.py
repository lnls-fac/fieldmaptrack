"""."""
import numpy as _np

import mathphys as _mp

from .track import SerretFrenetCoordSystem as _SerretFrenetCoordSystem
class Multipoles:

    def __init__(self,
                 trajectory = None,
                 perpendicular_grid = None,
                 normal_field_fitting_monomials = None,
                 skew_field_fitting_monomials = None):
        self.trajectory = trajectory
        self.perpendicular_grid = perpendicular_grid
        if normal_field_fitting_monomials is None:
            self.normal_field_fitting_monomials = []
        else:
            self.normal_field_fitting_monomials = normal_field_fitting_monomials
        if skew_field_fitting_monomials is None:
            self.skew_field_fitting_monomials = []
        else:
            self.skew_field_fitting_monomials = skew_field_fitting_monomials

    def get_multipole_labels(self, type, n):

        if n == 0:
            title  = type.title() + ' dipolar field'
            ylabel = title + ' [T]'
        elif n == 1:
            title  = type.title() + ' quadrupolar field'
            ylabel = title + ' [T/m]'
        elif n == 2:
            title  = type.title() + ' sextupolar field'
            ylabel = title + ' [T/m$^\mathrm{2}$]'
        elif n == 3:
            title  = type.title() + ' octupolar field'
            ylabel = title + ' [T/m$^\mathrm{3}$]'
        elif n == 4:
            title  = type.title() + ' decapolar field'
            ylabel = title + ' [T/m$^\mathrm{4}$]'
        elif n == 5:
            title  = type.title() + ' duodecapolar field'
            ylabel = title + ' [T/m$^\mathrm{5}$]'
        else:
            title  = type.title() + ' 2*({0}+1)-polar field'.format(n)
            pot = '{0}'.format(n)
            ylabel = title + ' [T/m$^\mathrm{'+pot+'}$]'

        fname = type + '_' + 'b{0}rx{0}'.format(n)
        return ylabel, title, fname

    def calc_multipoles(self, is_ref_trajectory_flag = False):
        """ calculates multipoles ([T] and [m] units) around given trajectory

            inputs

            - is_ref_trajectory: True or False

            if trajectory is a reference trajectory then, for the dipolar
            term, the algorithm fits only the difference between the fieldmap
            and the field on the ref. trajectory (dipolar error fit)

            if not a reference trajectory the algorithm does not fit the dipolar
            term at all. It is set explicitly according to the field on the trajectory"""


        # checks if x == 0 is in the perpendicular grid. gets its index
        grid = list(self.perpendicular_grid)
        try:
            grid_zero = grid.index(0)
        except ValueError:
            print('perpendicular grid needs to contain origin point')
            raise ValueError

        s = self.trajectory.s
        grid_meter = _np.array(grid) * _mp.units.mm_2_meter
        normal_field_monomials = list(self.normal_field_fitting_monomials)
        skew_field_monomials   = list(self.skew_field_fitting_monomials)
        self.normal_multipoles = _np.zeros((len(normal_field_monomials), len(s)))
        self.skew_multipoles   = _np.zeros((len(skew_field_monomials), len(s)))

        if is_ref_trajectory_flag:
            reference_field = _np.zeros((3,len(s)))
            reference_field[0,:] = self.trajectory.bx
            reference_field[1,:] = self.trajectory.by
            reference_field[2,:] = self.trajectory.bz
        else:
            #monomials.remove(0)
            pass

        self.max_fit_error_normal = (0,0)
        self.max_fit_error_skew   = (0,0)
        for i in range(len(s)):
            #print(str(i) + '/' + str(len(s)))
            sf = fieldmaptrack.SerretFrenetCoordSystem(self.trajectory, i)
            points = sf.get_transverse_line(grid)
            fieldmap_field = self.trajectory.fieldmap.interpolate_set(points)
            if is_ref_trajectory_flag:
                # trajectory is a reference trajectory
                field = fieldmap_field - _np.tile(reference_field[:,i].reshape((3,1)), (1, len(grid)))
                self.max_fit_error_skew = max_error if max_error[0] > self.max_fit_error_skew[0] else self.max_fit_error_skew
                self.normal_multipoles[:,i], max_error = mathphys.functions.polyfit(grid_meter, field[1,:], normal_field_monomials)
                self.max_fit_error_normal = max_error if max_error[0] > self.max_fit_error_normal[0] else self.max_fit_error_normal
            else:
                # trajectory is not a reference trajectory
                field = fieldmap_field
                self.skew_multipoles[:,i], max_error = mathphys.functions.polyfit(grid_meter, field[0,:], skew_field_monomials, algorithm='*lstsq')
                self.max_fit_error_skew = max_error if max_error[0] > self.max_fit_error_skew[0] else self.max_fit_error_skew
                self.normal_multipoles[:,i], max_error = mathphys.functions.polyfit(grid_meter, field[1,:], normal_field_monomials, algorithm='*lstsq')
                self.max_fit_error_normal = max_error if max_error[0] > self.max_fit_error_normal[0] else self.max_fit_error_normal

    def calc_multipoles_integrals(self):
        normal_field_monomials = self.normal_field_fitting_monomials
        self.normal_multipoles_integral = _np.zeros(self.normal_multipoles.shape[0])
        x = self.trajectory.s * _mp.units.mm_2_meter
        # if RK integration was done ds < 0, we need to invert integration sign
        sign = +1.0 if self.trajectory.s[-1] > self.trajectory.s[0] else -1.0
        for i in range(len(normal_field_monomials)):
            yb = self.normal_multipoles[i,:]
            self.normal_multipoles_integral[i] = sign * _np.trapz(y = yb, x = x)
        skew_field_monomials = self.skew_field_fitting_monomials
        self.skew_multipoles_integral = _np.zeros(self.skew_multipoles.shape[0])
        x = self.trajectory.s * _mp.units.mm_2_meter
        for i in range(len(skew_field_monomials)):
            ya = self.skew_multipoles[i,:]
            self.skew_multipoles_integral[i] = sign * _np.trapz(y = ya, x = x)

    def calc_multipoles_integrals_relative(self, main_polynom, main_monomial, r0, is_skew = False):

        self.r0 = r0
        self.main_monomial = main_monomial
        self.main_monomial_is_skew = is_skew
        r0 = self.r0 * _mp.units.mm_2_meter
        if is_skew:
            main_idx = list(self.skew_field_fitting_monomials).index(main_monomial)
            main_multipole = main_polynom[main_idx] * r0 ** main_monomial
        else:
            main_idx = list(self.normal_field_fitting_monomials).index(main_monomial)
            main_multipole = main_polynom[main_idx] * r0 ** main_monomial
        self.skew_multipoles_integral_relative = _np.zeros(self.skew_multipoles_integral.shape)
        self.normal_multipoles_integral_relative = _np.zeros(self.normal_multipoles_integral.shape)
        for i in range(len(self.normal_field_fitting_monomials)):
            n = self.normal_field_fitting_monomials[i]
            self.normal_multipoles_integral_relative[i] = self.normal_multipoles_integral[i] * (r0 ** n) / main_multipole
        for i in range(len(self.skew_field_fitting_monomials)):
            n = self.skew_field_fitting_monomials[i]
            self.skew_multipoles_integral_relative[i]   = self.skew_multipoles_integral[i]   * (r0 ** n) / main_multipole

    def __str__(self):

        nrpts = len(self.perpendicular_grid)
        grid_min = min(self.perpendicular_grid)
        grid_max = max(self.perpendicular_grid)
        normal_field_monomials = self.normal_field_fitting_monomials
        skew_field_monomials = self.skew_field_fitting_monomials

        all_monomials = sorted(set(list(normal_field_monomials) + list(skew_field_monomials)))
        r = ''
        try:
            r +=   '{0:<35s} {1} mm'.format('effective_length:', 1000*self.effective_length)
        except:
            pass
        r += '\n{0:<35s} {1}'.format('perpendicular_grid:', '{0} points in [{1:+f},{2:+f}] mm'.format(nrpts, grid_min, grid_max))
        r += '\n{0:<35s} {1:.3f}/{2:.3f} G/G'.format('max_fitting_error_normal', 1e4*self.max_fit_error_normal[0], 1e4*abs(self.max_fit_error_normal[1]))
        r += '\n{0:<35s} {1:.3f}/{2:.3f} G/G'.format('max_fitting_error_skew', 1e4*self.max_fit_error_skew[0], 1e4*abs(self.max_fit_error_skew[1]))
        r += '\n{0:<35s} {1} mm'.format('r0_for_relative_multipoles', self.r0)
        r += '\n{0:<35s} {1}'.format('main_monomial', 'n = {0}, skew:{1}'.format(self.main_monomial, self.main_monomial_is_skew))
        #r += '\n{0:<35s} {1:^13s} {2:^13s} {5:^13s} | {3:^13s} {4:^13s} {6:^13s}'.format('                   ', 'MaxAbs_Nn', 'Integ_Nn', 'MaxAbs_Sn', 'Integ_Sn', 'Nn/N0(@r0)', 'Sn/S0(@r0)')
        r += '\n{0:<35s} {1:^13s} {2:^13s} {5:^13s} | {3:^13s} {4:^13s} {6:^13s}'.format('                   ', 'Nn(s=0)', 'Integ_Nn', 'Sn(s=0)', 'Integ_Sn', 'Nn/N0(@r0)', 'Sn/S0(@r0)')
        r += '\n{0:<35s} {1:^13s} {2:^13s} {5:^13s} | {3:^13s} {4:^13s} {6:^13s}'.format('<multipole_order n>', '[T/m^n]', '[T.m/m^n]', '[T/m^n]', '[T.m/m^n]', '[]', '[]')
        for i in range(len(all_monomials)):
            n = all_monomials[i]
            try:
                idx = normal_field_monomials.index(all_monomials[i])
                #max_poly_b   = '{0:^13.4e}'.format(max(_np.abs(self.normal_multipoles[idx,:])))
                max_poly_b   = '{0:^13.4e}'.format(_np.abs(self.normal_multipoles[idx,0]))
                integ_poly_b = '{0:^+13.4e}'.format(self.normal_multipoles_integral[idx])
                integ_poly_b_relative = '{0:^+13.4e}'.format(self.normal_multipoles_integral_relative[idx])
            except ValueError:
                max_poly_b, integ_poly_b, integ_poly_b_relative = '---','---','---'
            try:
                idx = skew_field_monomials.index(all_monomials[i])
                #max_poly_a   = '{0:^13.4e}'.format(max(_np.abs(self.skew_multipoles[idx,:])))
                max_poly_a   = '{0:^13.4e}'.format(_np.abs(self.skew_multipoles[idx,0]))
                integ_poly_a = '{0:^+13.4e}'.format(self.skew_multipoles_integral[idx])
                integ_poly_a_relative = '{0:^+13.4e}'.format(self.skew_multipoles_integral_relative[idx])
            except ValueError:
                max_poly_a, integ_poly_a, integ_poly_a_relative = '---','---','---'
            r += '\n{0:<35s} {1:^13s} {2:^13s} {5:^13s} | {3:^13s} {4:^13s} {6:^13s}'.format('n={0:02d}:'.format(n), max_poly_b, integ_poly_b, max_poly_a, integ_poly_a, integ_poly_b_relative, integ_poly_a_relative)

        return r

    def save(self, filename, is_skew=False):

        with open(filename, 'w') as fp:
            fp.write('# multipoles\n')
            fp.write('# s[mm] ')
            if is_skew:
                for j in range(len(self.skew_field_fitting_monomials)):
                    fp.write('polynom_a[n={0}][T/m^n] '.format(self.skew_field_fitting_monomials[j]))
            else:
                for j in range(len(self.normal_field_fitting_monomials)):
                    fp.write('polynom_b[n={0}][T/m^n] '.format(self.normal_field_fitting_monomials[j]))
            fp.write('\n')
            traj = self.trajectory
            for i in range(len(traj.s)):
                fp.write('{0:+.16e} '.format(traj.s[i]))
                if is_skew:
                    for j in range(len(self.skew_field_fitting_monomials)):
                        fp.write('{0:+.16e} '.format(self.skew_multipoles[j, i]))
                else:
                    for j in range(len(self.normal_field_fitting_monomials)):
                        fp.write('{0:+.16e} '.format(self.normal_multipoles[j, i]))
                fp.write('\n')
