"""Track module."""

import math
import numpy as np
import fieldmaptrack.fieldmap as fieldmap


class TrackException(Exception):
    """TrackException class."""

    pass


class SerretFrenetCoordSystem:
    """SerretFrenetCoordSystem class."""

    def __init__(self, trajectory, point_idx=0):
        """."""
        t, i = trajectory, point_idx  # syntactic-sugars
        self.s = t.s[i]                # s position
        self.p = np.array((t.rx[i], t.ry[i], t.rz[i]))   # (rx,ry,rz) position
        self.t = np.array((t.px[i], t.py[i], t.pz[i]))   # tangential versor
        self.t /= math.sqrt(np.sum(self.t**2))           # renormalization
        self.n = np.array((t.pz[i], t.py[i], -t.px[i]))  # normal versor
        tx, ty, tz = self.t  # syntactic-sugars
        nx, ny, nz = self.n  # syntactic-sugars
        # skew versor k = t x n:
        self.k = np.array((ty*nz-tz*ny, tz*nx-tx*nz, tx*ny-ty*nx))

    def get_transverse_line(self, grid):
        """."""
        points = np.zeros((3, len(grid)))
        for i in range(len(grid)):
            points[:, i] = self.p + grid[i] * self.n
        return points

    def find_intersection(self, trajectory):
        """."""
        s = trajectory.s
        rx = trajectory.rx
        rz = trajectory.rz

        hit_i, hit_alpha, hit_dx = 0, None, None
        for i in range(len(s)-1):
            M11 = rx[i+1] - rx[i]
            M12 = -self.n[0]
            M21 = rz[i+1] - rz[i]
            M22 = -self.n[2]
            B1 = self.p[0] - rx[i]
            B2 = self.p[2] - rz[i]
            detM = M11 * M22 - M12 * M21
            invM = [[M22/detM, -M12/detM], [-M21/detM, M22/detM]]
            alpha = invM[0][0] * B1 + invM[0][1] * B2
            dx = invM[1][0] * B1 + invM[1][1] * B2
            if hit_alpha is None or abs(alpha) < abs(hit_alpha):
                hit_i, hit_alpha, hit_dx = i, alpha, dx
        if hit_i == 0 and hit_alpha < 0.0:
            # left extrapolation
            return (hit_i, hit_alpha, hit_dx)
        elif hit_i == len(s)-1 and hit_alpha > 1.0:
            # right extrapolation
            return (hit_i, hit_alpha, hit_dx)
        else:
            # interpolation
            return (hit_i, hit_alpha, hit_dx)


class Trajectory:
    """Trajectory class."""

    def __init__(self,
                 beam=None,
                 fieldmap=None,
                 not_raise_range_exceptions=False):
        """."""
        self.beam = beam
        self.fieldmap = fieldmap
        self.not_raise_range_exceptions = not_raise_range_exceptions

    def calc_trajectory(self, **kwargs):
        """."""
        # return self.calc_trajectory_rk1(**kwargs)
        return self.calc_trajectory_rk4(**kwargs)

    def calc_reference_point(self):
        """Calc reference point.

        defined as intersection of two asymptoctic
        straight section lines.
        """
        x1, z1 = self.rx[-2], self.rz[-2]
        x2, z2 = self.rx[-1], self.rz[-1]
        a = (x2-x1)/(z2-z1)
        b = 0.5*((x1+x2) - a * (z1+z2))
        return (b, 0.0)

    def calc_force(self, alpha, p):
        """."""
        rx, ry, rz, px, py, pz = p
        # calcs magnetic field on current position
        try:
            bx, by, bz = self.fieldmap.interpolate(rx, ry, rz)
        except (fieldmap.OutOfRangeRx,
                fieldmap.OutOfRangeRxMin,
                fieldmap.OutOfRangeRy):
            rstr = 'extrapolation at ' + str((rx, ry, rz))
            if self.not_raise_range_exceptions:
                print(rstr)
                bx, by, bz = 0.0, 0.0, 0.0
            else:
                raise TrackException(rstr)
        except fieldmap.OutOfRangeRz:
            bx, by, bz = 0.0, 0.0, 0.0

        # calcs derivatives of the eqs. of motion
        derivs = np.zeros(p.shape)
        derivs[0] = px                              # drx/ds
        derivs[1] = py                              # dry/ds
        derivs[2] = pz                              # drz/ds
        derivs[3] = - alpha * (py * bz - pz * by)   # dpx/ds
        derivs[4] = - alpha * (pz * bx - px * bz)   # dpy/ds
        derivs[5] = - alpha * (px * by - py * bx)   # dpz/ds

        return (derivs, bx, by, bz)

    def calc_trajectory_rk4(self,
                            init_rx=0.0, init_ry=0.0, init_rz=0.0,
                            init_px=0.0, init_py=0.0, init_pz=1.0,
                            s_length=None,
                            s_nrpts=None,
                            s_step=None,
                            min_rz=None,
                            force_midplane=False,
                            **kwargs):
        """Calculate trajectory of charged particle in an magnetic field.

            Algorithm              1st-order Runge-Kutta with constant step
            Independent variable   s: trajectory arc-length
            Dependent variables    rx,ry,rz: rectangular coordinates [mm]
                                   px,py,pz: p0-normalized momenta [no unit]
                                     (px**2+py**2+pz**2 == 1)
                                     horiz. defl. angle is atan(px/pz) [rad]
                                     verti. defl. angle is atan(py/pz) [rad]

            Input      force_midplane: whether to force trajectory to midplane
                       (usefull, for exmaple, for dipoles)
                       init_rx,init_ry,init_rz: initial rectangular position
                       of the particle.

                       init_px,init_py,init_pz: initial normalized momenta
                       of the particle. (px,py,pz)=(0.0,0.0,1.0), for example,
                       represents an initial velocity along the longitudinal
                       direction z.

                       s_length,s_nrpts,s_step,min_rz: parameters that control
                       range of trajectory integrations. Any of subset of these
                       parameters may be used, so long as they are consistent.
                       s_nrpts:  number of points in trajectory
                       s_step:   fixed step size in longitudinal coordinates
                       s_length: trajectory length
                       min_rz:   integrates trajectory until while rz < min_rz
                                 (s_step has be to be set as well)

        Output     s: vector with longitudinal position of points on trajectory
                   rx,ry,rz: vector with rectangular coordinates of points on
                             trajectory
                   px,py,pz: vector with normalized momenta of points on
                             trajectory
                   bx,by,bz: vector with magnetic field at the points on
                             trajectory
        """
        # processes input parameters
        if min_rz is not None:
            if s_step is None:
                rstr = ('s_step has to be set for this opion of trajectory '
                        'integration')
                raise TrackException(rstr)
            if s_length is not None or s_nrpts is not None:
                rstr = ('Neither s_length nor s_nrpts parameters should be '
                        'set for this option of trajectory integration')
                raise TrackException(rstr)
            s_length = 0.0
            s_nrpts = 0
        else:
            if s_step is None:
                if s_length is None or s_nrpts is None:
                    rstr = ('Both s_length and s_nrpts need to be set for '
                            'this option of trajectory integration')
                    raise TrackException(rstr)
                s_step = s_length / (s_nrpts - 1.0)
            else:
                if s_length is None:
                    if s_nrpts is None:
                        rstr = ('s_nrpts has to be set when s_step is given '
                                'and s_length is not set')
                        raise TrackException(rstr)
                    else:
                        s_length = s_step * (s_nrpts + 1.0)
                else:
                    s_nrpts = 1 + int(s_length / s_step)

        if 'kicks' in kwargs:
            kicks_rz = kwargs['kicks'][0]
            kicks_hkicks = kwargs['kicks'][1]
            kicks_vkicks = kwargs['kicks'][2]
            kicks_idx = 0
        else:
            kicks_idx = -1  # no kicks

        # inits auxiliary data structures
        self.s_step = s_step
        self.force_midplane = force_midplane
        self.init_rx, self.init_ry, self.init_rz = init_rx, init_ry, init_rz
        self.init_px, self.init_py, self.init_pz = init_px, init_py, init_pz
        self.s = []
        self.rx, self.ry, self.rz = [], [], []
        self.px, self.py, self.pz = [], [], []
        self.bx, self.by, self.bz = [], [], []
        alpha = 1.0/(1000.0*self.beam.brho)/self.beam.beta

        # RK integratior proper
        s, h, i = 0.0, self.s_step, 0
        rx, ry, rz = init_rx, init_ry, init_rz
        px, py, pz = init_px, init_py, init_pz
        while True:

            # insert kicks if available
            if kicks_idx >= 0:
                if kicks_idx < len(kicks_rz) and rz >= kicks_rz[kicks_idx]:
                    px += kicks_hkicks[kicks_idx]
                    py += kicks_vkicks[kicks_idx]
                    kicks_idx += 1

            # forces midplane, if the case
            if self.force_midplane:
                ry, py = 0.0, 0.0

            # calcs derivatives of the eqs. of motion
            p0 = np.array((rx, ry, rz, px, py, pz))
            if self.force_midplane:
                p0[1] = p0[4] = 0.0
            k1, bx, by, bz = self.calc_force(alpha, p0)
            p1 = p0 + (h/2.0) * k1
            if self.force_midplane:
                p1[1] = p1[4] = 0.0
            k2, _, _, _ = self.calc_force(alpha, p1)
            p2 = p0 + (h/2.0) * k2
            if self.force_midplane:
                p2[1] = p2[4] = 0.0
            k3, _, _, _ = self.calc_force(alpha, p2)
            p3 = p0 + h * k3
            if self.force_midplane:
                p3[1] = p3[4] = 0.0
            k4, _, _, _ = self.calc_force(alpha, p3)

            # stores current position
            self.s.append(s)
            self.rx.append(rx), self.ry.append(ry), self.rz.append(rz)
            self.px.append(px), self.py.append(py), self.pz.append(pz)
            try:
                self.bx.append(bx[0]), self.by.append(by[0]), \
                    self.bz.append(bz[0])
            except:
                self.bx.append(bx), self.by.append(by), self.bz.append(bz)

            # propagates to next point
            rx, ry, rz, px, py, pz = p0 + (h/6.0) * (k1 + 2*k2 + 2*k3 + k4)
            s += self.s_step
            i += 1

            # tests if end of integration is reached
            if min_rz is not None:
                if self.s_step > 0.0:
                    if rz > min_rz:
                        break
                else:
                    if rz < min_rz:
                        break
                # if abs(rz) > abs(min_rz):
                #     print(rz, min_rz)
                #     break
            else:
                # integration is being done until <s_nrpts is reached>
                if i == s_nrpts:
                    break

        # converts python native lists into numpy arrays
        # (appending in native lists is faster than in nympy arrays)
        self.s = np.array(self.s)
        self.rx, self.ry, self.rz = \
            np.array(self.rx), np.array(self.ry), np.array(self.rz)
        self.px, self.py, self.pz = \
            np.array(self.px), np.array(self.py), np.array(self.pz)
        self.error_estimate = math.sqrt(self.px[-1]**2 + self.py[-1]**2 +
                                        self.pz[-1]**2) - 1.0

    def __str__(self):
        """."""
        bx, by, bz = \
            [abs(x) for x in self.bx], [abs(x) for x in self.by], \
            [abs(x) for x in self.bz]
        max_bx, max_by, max_bz = max(bx), max(by), max(bz)

        s_max_bx, rx_max_bx, ry_max_bx, rz_max_bx = \
            self.s[bx.index(max_bx)], self.rx[bx.index(max_bx)], \
            self.ry[bx.index(max_bx)], self.rz[bx.index(max_bx)]
        s_max_by, rx_max_by, ry_max_by, rz_max_by = \
            self.s[by.index(max_by)], self.rx[by.index(max_by)], \
            self.ry[by.index(max_by)], self.rz[by.index(max_by)]
        s_max_bz, rx_max_bz, ry_max_bz, rz_max_bz = \
            self.s[bz.index(max_bz)], self.rx[bz.index(max_bz)], \
            self.ry[bz.index(max_bz)], self.rz[bz.index(max_bz)]
        theta_x = math.atan(self.px[-1]/self.pz[-1])
        theta_y = math.atan(self.py[-1]/self.pz[-1])
        ref_point = self.calc_reference_point()
        dtheta_x = theta_x - math.atan(self.px[0]/self.pz[0])
        dtheta_y = theta_y - math.atan(self.py[0]/self.pz[0])

        r = ''
        r += '{0:<35s} {1:.6e} GeV'.format('beam_energy:', self.beam.energy)
        r += '\n{0:<35s} {1:+.4e} deg.'.format('horizontal_deflection_angle:',
                                               dtheta_x * (180.0/math.pi))
        r += '\n{0:<35s} {1:+.4e} deg.'.format('vertical_deflection_angle:',
                                               dtheta_y * (180.0/math.pi))
        r += '\n{0:<35s} {1:+.4e} deg.'.format('final_horizontal_angle:',
                                               theta_x * (180.0/math.pi))
        r += '\n{0:<35s} {1:+.4e} deg.'.format('final_vertical_angle:',
                                               theta_y * (180.0/math.pi))
        r += '\n{0:<35s} {1} mm'.format('trajectory_length:',
                                        self.s[-1]-self.s[0])
        r += '\n{0:<35s} {1}'.format('trajectory_nrpts:', len(self.s))
        r += '\n{0:<35s} {1} mm'.format('trajectory_s_step:', self.s_step)

        r += '\n{0:<35s} {1:+f} Tesla at (s,rx,ry,rz) = ({2},{3},{4},{5}) mm'.format('max_abs_bx@trajectory:', self.bx[bx.index(max_bx)], s_max_bx, rx_max_bx, ry_max_bx, rz_max_bx)
        r += '\n{0:<35s} {1:+f} Tesla at (s,rx,ry,rz) = ({2},{3},{4},{5}) mm'.format('max_abs_by@trajectory:', self.by[by.index(max_by)], s_max_by, rx_max_by, ry_max_bx, rz_max_by)
        r += '\n{0:<35s} {1:+f} Tesla at (s,rx,ry,rz) = ({2},{3},{4},{5}) mm'.format('max_abs_bz@trajectory:', self.bz[bz.index(max_bz)], s_max_bz, rx_max_bz, ry_max_bz, rz_max_bz)
        r += '\n{0:<35s} {1:+f} mm'.format('rx position of reference point:', ref_point[0])
        r += '\n{0:<35s} {1:+f} mm'.format('initial rx position of trajectory:', self.rx[0])
        return r

    def save(self, filename):

        with open(filename, 'w') as fp:
            fp.write('# trajectory\n')
            fp.write('# s[mm] rx[mm] ry[mm] rz[mm] px[rad] py[rad] pz[rad]\n')
            for i in range(len(self.s)):
                fp.write('{0:.16e} {1:+.16e} {2:+.16e} {3:+.16e} {4:+.16e} {5:+.16e} {6:+.16e}\n'.format(self.s[i],self.rx[i],self.ry[i],self.rz[i],self.px[i],self.py[i],self.pz[i]))

    def load(self, filename):
        lines = [line.strip() for line in open(filename)]
        s,rx,ry,rz,px,py,pz,bx,by,bz = [],[],[],[],[],[],[],[],[],[]
        for line in lines[2:]:
            words = line.split()
            s.append(float(words[0]))
            rx.append(float(words[1])), ry.append(float(words[2])), rz.append(float(words[3]))
            px.append(float(words[4])), py.append(float(words[5])), pz.append(float(words[6]))
        self.s = np.array(s)
        self.rx, self.ry, self.rz = np.array(rx), np.array(ry), np.array(rz)
        self.px, self.py, self.pz = np.array(px), np.array(py), np.array(pz)
        self.bx, self.by, self.bz = np.array(rx), np.array(ry), np.array(rz)
        # calcs field on reference trajectory
        for i in range(len(s)):
            self.bx[i], self.by[i], self.bz[i] = self.fieldmap.interpolate(self.rx[i], self.ry[i], self.rz[i])
        self.s_step = self.s[1] - self.s[0]

    def save_field(self, filename):

        with open(filename, 'w') as fp:
            fp.write('# field on trajectory\n')
            fp.write('# s[mm] rx[mm] ry[mm] rz[mm] bx[T] by[T] bz[T]\n')
            for i in range(len(self.s)):
                fp.write('{0:.16e} {1:+.16e} {2:+.16e} {3:+.16e} {4:+.16e} {5:+.16e} {6:+.16e}\n'.format(self.s[i], self.rx[i], self.ry[i], self.rz[i], self.bx[i], self.by[i], self.bz[i]))

    def calc_sagitta(self, half_dipole_length):

        rx = self.rx
        rz = self.rz

        if abs(rz[-1]) < half_dipole_length:
            raise TrackException('trajectory path does not exit dipole')

        i = 0
        while (abs(rz[i]) < half_dipole_length):
            i += 1
        sagitta = rx[0] - rx[i]
        return sagitta

    def find_intersection_point(self, csys):

        nx, nz = csys.n[0], csys.n[2]
        px, pz = csys.p[0], csys.p[2]
        r = (None, None, None)
        for i in range(len(self.s)-1):
            p1x, p1z = self.rx[i], self.rz[i]
            p2x, p2z = self.rx[i+1], self.rz[i+1]
            M11 = p2x - p1x
            M12 = nx
            M21 = p2z - p1z
            M22 = nz
            b1  = p1x - px
            b2  = p1z - pz
            det = M11*M22 - M12*M21
            iM11, iM12 = M22/det, -M12/det
            iM21, iM22 = -M21/det, M11/det
            ds = iM11 * b1 + iM12 * b2
            dx = iM21 * b1 + iM22 * b2
            if 0.0 <= ds <= 1.0:
                r = (i, ds, dx)
                break
        return r
