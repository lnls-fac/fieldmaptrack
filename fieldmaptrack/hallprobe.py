"""Hallsensor analysis module.

In this module classes are implemented that can be used to run interactive
fieldmap analysis and also to generate input files for the fac-fma-* scripts.
"""

import os as _os
import math as _math
import numpy as _np   # necessary for eval functions
import time as _time
import datetime as _datetime
import fieldmaptrack as _fmap
import matplotlib.pyplot as _plt
from copy import deepcopy as _dcopy
import mathphys as _mp


_path_base = '/home/fac_files/lnls-ima/'


defaults = {
    # default fieldmaptrack parameters.
    # These parameters are converted in FMapAnalysis object attributes
    # in the class constructor.

    'si-dipoles-b2': {
        # paths
        '_path_base': _path_base,
        '_path_repo': 'si-dipoles-b2/',
        '_path_model': 'model-08/',
        '_path_measf': 'measurement/magnetic/hallprobe/',
        '_path_dataset': '',
        # rawdata
        'magnet_type': 'dipole',
        'config_label': 'B2-analysis',
        # trajectory
        'model_nominal_angle': 4.0964,  # [degree]
        'model_nominal_KL': -0.95879438,  # from AT model
        # 'model_nominal_KL': -0.96018,  # from wiki spec
        'traj_init_rx': 8.153,  # [mm]
        'traj_init_px': 0.0,   # [deg]
        'traj_rk_s_step': 0.1,  # [mm]
        'model_nominal_refrx': 19.428,  # [mm]
        # multipoles
        'beam_energy': 3.0,
        'multipoles_normal_field_fitting_monomials': (0, 1, 2, 3, 4, 5, 6),
        'multipoles_skew_field_fitting_monomials': (),
        'multipoles_perpendicular_grid': '_np.linspace(-12,12,65)',
        'multipoles_r0': 12,  # [mm]
        'normalization_monomial': 0,
        'normalization_is_skew': False,
        'normal_multipoles_main_monomials': (0, 1),
        'skew_multipoles_main_monomials': (),
        'model_segmentation': (125, 55, 10, 5, 5, 5, 5, 10, 10, 175, 175,
                               20, 10, 15, 20, 30, 32, 32.5), },

    'si-dipoles-b1': {
        # paths
        '_path_base': _path_base,
        '_path_repo': 'si-dipoles-b1/',
        '_path_model': 'model-09/',
        '_path_measf': 'measurement/magnetic/hallprobe/',
        '_path_dataset': '',
        # rawdata
        'magnet_type': 'dipole',
        'config_label': 'B1-analysis',
        # trajectory
        'model_nominal_angle': 2.7553,  # [degree]
        # 'model_nominal_KL': -0.64400437024,  # from AT model (0.32% menor)
        'model_nominal_KL': -0.6460718,  # from AT model - dipole model 09
        'traj_init_rx': 8.527,  # [mm]
        'traj_init_px': 0.0,   # [deg]
        'traj_rk_s_step': 0.05,  # [mm]
        'model_nominal_refrx': 13.689,  # [mm]
        # multipoles
        'beam_energy': 3.0,
        'multipoles_normal_field_fitting_monomials': (0, 1, 2, 3, 4, 5, 6),
        'multipoles_skew_field_fitting_monomials': (),
        'multipoles_perpendicular_grid': '_np.linspace(-10,10,41)',
        'multipoles_r0': 12,  # [mm]
        'normalization_monomial': 0,
        'normalization_is_skew': False,
        'normal_multipoles_main_monomials': (0, 1),
        'skew_multipoles_main_monomials': (),
        'model_segmentation': (2, 3, 5, 5, 5, 10, 40, 150, 100,
                               50, 34, 16, 40, 40, 50), },
     }


class TemplateText:
    """."""

    rawfield = (
        '#  ==============================================',
        '# | fac-fma-rawfield.py input_file               |',
        '# | Accelerator Physics LNLS                     |',
        '# |                                              |',
        '# | Date: TIMESTAMP                |',
        "# | generated with 'fieldmaptrack.hallsensor.py' |",
        '#  ==============================================',
        '',
        '# --- Summary ---',
        '#',
        '# This is the input file for fac-fma-rawfield.py script',
        '# this script reads a fieldmap from a 3D magnet model, stores it',
        '# for latter analysis and prints and plots basic information on the',
        '# field map. It is used to quickly inspect the fieldmap',
        '',
        '',
        '# --- Input parameters ---',
        '',
        '# Each analysis has an identity label used for naming output files',
        '',
        "config_label                  'CONFIG_LABEL'",
        '',
        '',
        '# The next parameter specifies the type of magnet to be analysed.',
        '# each type may have its own particular algorithms to be applied',
        '',
        "magnet_type                   'MAGNET_TYPE'",
        '',
        '',
        '# the full name of the file that contains the field map',
        '',
        "fmap_filename                 'FMAP_FILENAME'",
        '',
        '# Runge-Kutta algorithm that is used for the integration of the eqs.',
        '# of motion needs to know what to do when trajectory reaches the',
        '# fieldmap bounds. It will either extrapolate the field along the',
        '# longitudinal (z) direction or consider it to have vanished. This',
        '# is controlled with the parameter below. Bear in mind that the ',
        '# calculation of extrapolation coefficients is very time-consuming',
        '# currently. As for the transverse directions (x and y), the RK ',
        '# algorithm will generate exceptions.',
        '',
        'fmap_extrapolation_flag  False',)

    trajectory = (
        '#  ==============================================',
        '# | fac-fma-trajectory.py input_file             |',
        '# | Accelerator Physics LNLS                     |',
        '# |                                              |',
        '# | Date: TIMESTAMP                |',
        "# | generated with 'fieldmaptrack.hallsensor.py' |",
        '#  ==============================================',
        '',
        '# --- Summary ---',
        '#',
        '# This is the input file for trajectory calculation based on a given',
        '# fieldmap which is performed with the script',
        "# 'fac-fma-trajectory.py'. A controllable fixed-size Runge-Kutta ",
        '# algorithm is used to integrate the equations of motion of a single',
        '# electron in the presence of the magnetic field as defined in the ',
        '# fieldmap.',
        '#',
        '# The implemented equations of motion are not approximated. Provided',
        '# a sufficiently fine RK step is chosen, this scripts may be used to',
        '# accurately obtain the trajectory of the electron with arbitrary ',
        '# energy.',
        '',
        '# Runge-Kutta algorithm used for the integration of the eqs. of ',
        '# motion needs to know what to do when trajectory reaches the ',
        '# fieldmap bounds. It will either extrapolate the field along the ',
        '# longitudinal (z) direction or consider it to have vanished.',
        '# As for the transverse directions (x and y), the RK algorithm will',
        '# generate exceptions.',
        '',
        '',
        '# --- Input parameters ---',
        '',
        '# Each analysis has an identity label used for naming output files',
        '',
        "config_label                  'CONFIG_LABEL'",
        '',
        '',
        '# beam energy',
        '',
        'beam_energy                   BEAM_ENERGY     # [GeV]',
        '',
        '',
        '# A trajectory can also be read from file. This is useful when the ',
        '# fieldmap of 3D models with errors are being analysed. In this ',
        '# case we want to use as reference trajectory a trajectory that was ',
        '# calculated from the 3D model without errors and saved to file.',
        "# If parameter 'traj_load_filename' is set to 'None' then a new",
        '# reference trajectory with be calculated with RK on the given',
        '# fieldmap.',
        '',
        'traj_load_filename            None',
        '',
        '',
        "# If parameter 'traj_is_reference_traj' is set to True the algorithm",
        '# will rescale the fieldmap so that the total trajectory deflection ',
        '# will exactly match the nominal deflection',
        '',
        'traj_is_reference_traj        False',
        ('model_nominal_angle           MODEL_NOMINAL_ANGLE  '
         ' # [deg] nominal deflection angle of the magnet.'),
        '',
        '# There is the option to restrain the trajectory to the midplane ',
        '# (y = 0 mm) of the magnet',
        '',
        'traj_force_midplane_flag      True',
        '',
        '',
        '# There is the option to serach for a initial rx position that will',
        '# result in a trajectory that is centered in the good-field region ',
        '# of the magnet (around rx == 0)',
        '',
        'traj_center_sagitta_flag      False',
        '',
        '',
        '# The RK algorithm always integrates the trajectory from the center',
        '# of the magnet (z = s = 0 mm). The limits of the RK integration may',
        '# be specified in various ways:',
        "#   If only 'traj_rk_s_step' is given then the algorithm will ",
        '#   integrate until the z coordinate of the particle reaches the ',
        '#   fieldmap bound.',
        '',
        'traj_init_rx                  TRAJ_INIT_RX  # [mm]',
        'traj_init_px                  TRAJ_INIT_PX  # [deg]',
        'traj_rk_s_step                TRAJ_RK_S_STEP  # [mm]',
        'traj_rk_length                None  # [mm]',
        'traj_rk_nrpts                 None',
        '',
        '# Whether to save trajectory to an ASCII file',
        '',
        'traj_save                     True',)

    multipoles = (
        '#  ==============================================',
        '# | fac-fma-multipoles.py input_file             |',
        '# | Accelerator Physics LNLS                     |',
        '# |                                              |',
        '# | Date: TIMESTAMP                |',
        "# | generated with 'fieldmaptrack.hallsensor.py' |",
        '#  ==============================================',
        '',
        '# --- Summary ---',
        '#',
        "# This is the input file for the 'fac-fma-multipoles.py' script",
        '# this script calculates the multipoles around the reference ',
        '# trajectory.',
        '',
        '',
        '# --- Input parameters ---',
        '',
        '# Each analysis has an identity label used for naming output files',
        '',
        "config_label                  'CONFIG_LABEL'",
        '',
        '',
        '# the multipoles (m1,m2,...) to be calculated are defined by a list',
        '# of position x exponents (n1,n2,...):',
        '#   By = m1 * x^n1 + m2 * x^n2 + ...',
        '',
        ('multipoles_normal_field_fitting_monomials        '
         'MULTIPOLES_NORMAL_FIELD_FITTING_MONOMIALS   # monomials to be '
         'included in the polynomial fit of multipoles'),
        ('multipoles_skew_field_fitting_monomials          '
         'MULTIPOLES_SKEW_FIELD_FITTING_MONOMIALS'),
        '',
        '# Grid of perpendicular points around each point of the reference '
        '# trajectory for the polynomial fit of By and Bx',
        '',
        ('multipoles_perpendicular_grid           '
         'MULTIPOLES_PERPENDICULAR_GRID  # grid of points on perpendicular '
         'line to ref trajectory [mm]'),
        '',
        '# After multipole coeffs are calculated, their normalized strengths',
        '# at perp. position r0 are calculated (as defined in tracy)',
        '',
        ('multipoles_r0                     MULTIPOLES_R0  # [mm] horizontal'
         ' position at which polynomial fields are calculated relative to the'
         ' principal multipole'),
        'normalization_monomial            NORMALIZATION_MONOMIAL',
        'normalization_is_skew             NORMALIZATION_IS_SKEW',
        '',
        '# Integrated residual field (converted to kick angle) calculated ',
        '# from fitted multipoles and from integrated fieldmap are compared.',
        '# The parameter below lists the monomials which are supposed to ',
        '# define the main field. The rest makes up for the residual field',
        '',
        ('normal_multipoles_main_monomials        '
         'NORMAL_MULTIPOLES_MAIN_MONOMIALS'),
        ('skew_multipoles_main_monomials          '
         'SKEW_MULTIPOLES_MAIN_MONOMIALS'), )

    model = (
        '#  ==============================================',
        '# | fac-fma-model.py input_file                  |',
        '# | Accelerator Physics LNLS                     |',
        '# |                                              |',
        '# | Date: TIMESTAMP                |',
        "# | generated with 'fieldmaptrack.hallsensor.py' |",
        '#  ==============================================',
        '',
        '# --- Summary ---',
        '#',
        '# This script integrates fitted multipoles at each segment of the ',
        '# hard-edge model',
        '',
        '# --- Input parameters ---',
        '',
        '# each analysis has an identity label used for naming output files',

        "config_label                      'CONFIG_LABEL'",
        '',
        '',
        '# list with lengths of model segments',
        '',
        'model_segmentation            MODEL_SEGMENTATION', )


class FMapAnalysis(_fmap.common_analysis.Config):
    """Double-sided FieldMap analysis class."""

    def __init__(self, magnet=None, fmap_fname=None):
        """."""
        super().__init__(fname=None)
        self._magnet = magnet
        self._set_defaults()
        self.fmap_filename = self._get_fname(magnet, fmap_fname)
        # -- default rawfield parameters
        self.fmap_extrapolation_flag = False
        self.interactive_mode = True
        # -- default trajectory parameters
        self.traj_rk_length = None
        self.traj_rk_nrpts = None
        self.traj_center_sagitta_flag = False
        self.traj_force_midplane_flag = True
        self.traj_is_reference_traj = False
        self.traj_load_filename = None

    @property
    def energy(self):
        """."""
        return self.beam_energy

    @energy.setter
    def energy(self, value):
        """."""
        self._set_energy(value)

    @property
    def deflection_angle(self):
        """."""
        return self._get_deflection_angle()

    @property
    def reference_point(self):
        """."""
        return self._get_reference_point()

    def analysis_rawfield(self):
        """."""
        analysis = _fmap.common_analysis.get_analysis_symbol(self.magnet_type)
        analysis.raw_fieldmap_analysis(self)

    def analysis_trajectory(self):
        """."""
        analysis = _fmap.common_analysis.get_analysis_symbol(self.magnet_type)
        analysis.trajectory_analysis(self)

    def analysis_multipoles(self):
        """."""
        analysis = _fmap.common_analysis.get_analysis_symbol(self.magnet_type)
        analysis.multipoles_analysis(self)

    def search_energy(self, print_flag=False, angle_tol=0.001):
        """."""
        def calc_dangle():
            angle = abs(self.deflection_angle)
            dangle = 1000 * (_np.pi/180.0) * \
                (angle - self.model_nominal_angle)
            return dangle, angle

        _fmt_search_energy = \
            '{:<8s} {:+010.4f}A n={:02d}: ' + \
            '{:10.8f} GeV - {:07.4f}|{:07.4f} dang: {:+09.4f} mrad'

        if not hasattr(self, 'fmap'):
            self.analysis_rawfield()
        f = 1.0
        dangle = float('+inf')
        iter = 1
        while abs(dangle) > angle_tol:
            self.energy = self.energy * f
            # print('h:', self.beam_energy)
            # if hasattr(self, '_configN'):
            #     print('h:', self._configN.beam_energy)
            self.analysis_trajectory()
            dangle, angle = calc_dangle()
            if print_flag:
                print(_fmt_search_energy.format(
                    self._magnet, float(self.fmap.current), iter,
                    self.energy, angle, self.model_nominal_angle, dangle))
            f = angle / self.model_nominal_angle
            iter += 1
        return iter

    def search_x0(self, print_flag=False, angle_tol=0.001):
        """."""
        def calc_dangle():
            angle = abs(self.deflection_angle)
            dangle = 1000 * (_np.pi/180.0) * \
                (angle - self.model_nominal_angle)
            return dangle, angle

        _fmt_search_x0 = \
            '{:<8s} {:+010.4f}A n={:02d}: ' + \
            '{:6.3f} mm - {:07.4f}|{:07.4f} dang: {:+09.4f} mrad'

        if not hasattr(self, 'fmap'):
            self.analysis_rawfield()
        dangle = float('+inf')
        iter = 1
        x, a = [], []
        while abs(dangle) > angle_tol:
            self.analysis_trajectory()
            dangle, angle = calc_dangle()
            x.append(self.traj_init_rx)
            a.append(angle)
            if len(x) == 1:
                self.traj_init_rx += 0.1
            else:
                z = _np.polyfit(a, x, len(x)-1)
                p = _np.poly1d(z)
                x0 = p(self.model_nominal_angle)
                self.traj_init_rx = x0

            if print_flag:
                print(_fmt_search_x0.format(
                    self._magnet, float(self.fmap.current), iter,
                    self.traj_init_rx, angle, self.model_nominal_angle, dangle))
            iter += 1
        return iter

    def create_input_rawfield(self, path, force=False, path_subdir=None):
        """."""
        text = TemplateText.rawfield
        text = self._process_text(list(text))
        FMapAnalysis._create_file(
            text, path, '/rawfield.in', force, path_subdir)
        # return text

    def create_input_trajectory(self, path, force=False, path_subdir=None):
        """."""
        text = TemplateText.trajectory
        text = self._process_text(list(text))
        FMapAnalysis._create_file(
            text, path, '/trajectory.in', force, path_subdir)
        # return text

    def create_input_multipoles(self, path, force=False, path_subdir=None):
        """."""
        text = TemplateText.multipoles
        text = self._process_text(list(text))
        FMapAnalysis._create_file(
            text, path, '/multipoles.in', force, path_subdir)
        # return text

    def create_input_model(self, path, force=False, path_subdir=None):
        """."""
        text = TemplateText.model
        text = self._process_text(list(text))
        FMapAnalysis._create_file(
            text, path, '/model.in', force, path_subdir)
        # return text

    def create_input(self, path, force=False, path_subdir=None):
        """."""
        self.create_input_rawfield(
            path, force=force, path_subdir=path_subdir)
        self.create_input_trajectory(
            path, force=force, path_subdir=path_subdir)
        self.create_input_multipoles(
            path, force=force, path_subdir=path_subdir)
        self.create_input_model(
            path, force=force, path_subdir=path_subdir)

    @staticmethod
    def load_output_trajectory(path, path_subdir=None):
        """."""
        data = dict()
        if path_subdir is None:
            fname = path + 'trajectory.out'
        else:
            fname = path + path_subdir + 'trajectory.out'
        with open(fname, 'r') as f:
            lines = f.readlines()
        for line in lines:
            if 'horizontal_deflection_angle' in line:
                data['horizontal_deflection_angle'] = float(line.split()[1])
            elif 'rx position of reference' in line:
                data['rx position of reference'] = float(line.split()[5])
            elif 'initial rx position of trajectory' in line:
                data['initial rx position of trajectory'] = \
                    float(line.split()[5])
            elif 'beam_energy' in line:
                data['beam_energy'] = float(line.split()[1])
        return data

    @staticmethod
    def load_output_multipoles(path, path_subdir=None):
        """."""
        if path_subdir is None:
            fname = path + 'multipoles.out'
        else:
            fname = path + path_subdir + 'multipoles.out'
        with open(fname, 'r') as f:
            lines = f.readlines()
        harms = []
        mpoles_normal = []
        mpoles_skew = []
        for line in lines:
            if 'n=' in line:
                words = line.split()
                h = int(words[0].replace('n=', '').replace(':', ''))
                mn = float(words[2].replace('---', 'NaN'))
                ms = float(words[6].replace('---', 'NaN'))
                harms.append(h)
                mpoles_normal.append(mn)
                mpoles_skew.append(ms)
        return mpoles_normal, mpoles_skew, harms

    @staticmethod
    def load_output_model(path, path_subdir=None):
        """."""
        if path_subdir is None:
            fname = path + 'model.out'
        else:
            fname = path + path_subdir + 'model.out'
        with open(fname, 'r') as f:
            lines = f.readlines()
        lens = []
        angs = []
        for line in lines:
            if 'angle' not in line:
                words = line.split()
                le = float(words[0])
                an = float(words[2])
                lens.append(le)
                angs.append(an)
        return lens, angs

    def get_defaults(self):
        """Return default attributes."""
        if 'B2-' in self._magnet:
            return defaults['si-dipoles-b2']
        elif 'B1-' in self._magnet:
            return defaults['si-dipoles-b1']
        else:
            raise ValueError('defaults not defined for magnet!!!')

    def _process_text(self, text):
        defs = self.get_defaults()
        keyvals = {k.upper(): str(defs[k]) for k in defs}
        for i in range(len(text)):
            # replace with current attribute values
            text[i] = text[i].replace('TIMESTAMP', self._timestamp())
            text[i] = text[i].replace('FMAP_FILENAME', self.fmap_filename)
            text[i] = text[i].replace('BEAM_ENERGY', str(self.energy))
            text[i] = text[i].replace(
                'TRAJ_RK_S_STEP', str(self.traj_rk_s_step))
            text[i] = text[i].replace(
                'TRAJ_INIT_PX', str(self.traj_init_px))
            # replace with default attribute values
            for k in keyvals:
                if k in text[i]:
                    text[i] = text[i].replace(k, keyvals[k])
        return text

    def _get_deflection_angle(self):
        return self.calc_deflection_angle(self.traj)

    def _set_energy(self, value):
        self.beam_energy = value

    def _set_defaults(self):
        defs = self.get_defaults()
        for k, v in defs.items():
            if isinstance(v, str):
                try:
                    expr = 'setattr(self, "' + k + '", ' + v + ')'
                    # print('e1:', expr)
                    eval(expr)
                except (NameError, SyntaxError, TypeError):
                    expr = 'setattr(self, "' + k + '", "' + v + '")'
                    # print('e2:', expr)
                    eval(expr)
            else:
                expr = 'setattr(self, "' + k + '", ' + str(v) + ')'
                # print('e3:', expr)
                eval(expr)

    def _get_fname(self, magnet, fmap_fname):
        self._path_magnet = self._magnet + '/'
        fname = \
            self._path_base + \
            self._path_repo + \
            self._path_model + \
            self._path_measf + \
            self._path_magnet + \
            self._path_dataset + \
            fmap_fname
        return fname

    # def _get_reference_point(self):
    #     return self.traj.calc_reference_point()

    @staticmethod
    def calc_asymptotic_line_coeffs(traj):
        """."""
        a1, b1 = _np.polyfit(traj.rz[-5:],
                             traj.rx[-5:], 1)
        return a1, b1

    def _get_reference_point(self):
        a1, b1 = FMapAnalysis.calc_asymptotic_line_coeffs(self.traj)
        return (b1, 0.0)

    @staticmethod
    def _timestamp(now=None):
        """Get formatted timestamp ."""
        if now is None:
            now = _time.time()
        st = _datetime.datetime.fromtimestamp(now).strftime(
            '%Y-%m-%d-%H:%M:%S')
        st = st + '.{0:03d}'.format(int(1000*(now-int(now))))
        return st

    @staticmethod
    def _create_file(text, path, fname, force, path_subdir):
        if path_subdir is not None:
            path += path_subdir
        if not _os.path.exists(path):
            _os.makedirs(path)
        else:
            if not force:
                raise _os.FileExistsError()
        with open(path + fname, 'w') as fp:
            for line in text:
                fp.write(line + '\n')

    @staticmethod
    def calc_deflection_angle(traj):
        """."""
        px1, pz1 = traj.px[0], traj.pz[0]
        px2, pz2 = traj.px[-1], traj.pz[-1]
        t1 = _math.atan(px1/pz1)
        t2 = _math.atan(px2/pz2)
        return (180/_math.pi) * (t2 - t1)


class DoubleFMapAnalysis(FMapAnalysis):
    """Double-sided FieldMap analysis class.

    To be used for tracking for Z > 0, as well as Z < 0 regions.
    """

    def __init__(self, **kwargs):
        """."""
        super().__init__(**kwargs)

    @property
    def configN(self):
        """."""
        return self._configN

    def analysis_trajectory(self):
        """."""
        super().analysis_trajectory()
        analysis = _fmap.common_analysis.get_analysis_symbol(self.magnet_type)
        # if not hasattr(self, '_configN'):
        #     self._configN = _dcopy(self)
        #     self._configN.fmap = self.fmap  # to save space
        #     self._configN.traj_rk_s_step = -self.traj_rk_s_step
        self._configN = _dcopy(self)
        self._configN.fmap = self.fmap  # to save space
        self._configN.traj_rk_s_step = -self.traj_rk_s_step
        analysis.trajectory_analysis(self._configN)

    def analysis_multipoles(self):
        """."""
        super().analysis_multipoles()
        analysis = _fmap.common_analysis.get_analysis_symbol(self.magnet_type)
        self._configN = analysis.multipoles_analysis(self._configN)

    def _get_deflection_angle(self):
        a1 = self.calc_deflection_angle(self.traj)
        a2 = self.calc_deflection_angle(self._configN.traj)
        return a1 - a2

    def _set_energy(self, value):
        self.beam_energy = value
        if hasattr(self, '_configN'):
            self._configN.beam_energy = value

    def _get_reference_point(self):
        a1, b1 = FMapAnalysis.calc_asymptotic_line_coeffs(self.traj)
        a2, b2 = FMapAnalysis.calc_asymptotic_line_coeffs(self.configN.traj)
        z = - (b2 - b1)/(a2 - a1)
        x = ((a1*z + b1) + (a2*z + b2))/2.0
        return (x, z)


_fmap_files = None


def get_fmap_files_B2(magnet=None, current=None):
    """."""
    global _fmap_files

    if _fmap_files is None:
        _s = '_Hall_X=-64_28mm_Z=-1200_1200mm_Imc='
        _fmap_files = (
            '2018-11-10_B2-001' + _s + '381.7A_ID=474.dat',
            '2018-11-10_B2-001' + _s + '401.8A_ID=475.dat',
            '2018-11-10_B2-001' + _s + '421.9A_ID=476.dat',
            '2018-11-17_B2-002' + _s + '381.7A_ID=558.dat',
            '2018-11-17_B2-002' + _s + '401.8A_ID=559.dat',
            '2018-11-17_B2-002' + _s + '421.9A_ID=560.dat',
            '2018-11-13_B2-003' + _s + '381.7A_ID=510.dat',
            '2018-11-13_B2-003' + _s + '401.8A_ID=511.dat',
            '2018-11-13_B2-003' + _s + '421.9A_ID=512.dat',
            '2018-11-10_B2-004' + _s + '381.7A_ID=478.dat',
            '2018-11-10_B2-004' + _s + '401.8A_ID=479.dat',
            '2018-11-10_B2-004' + _s + '421.9A_ID=480.dat',
            '2018-11-01_B2-005' + _s + '381.7A_ID=383.dat',
            '2018-11-01_B2-005' + _s + '401.8A_ID=384.dat',
            '2018-11-01_B2-005' + _s + '421.9A_ID=385.dat',
            '2018-11-08_B2-006' + _s + '381.7A_ID=448.dat',
            '2018-11-08_B2-006' + _s + '401.8A_ID=449.dat',
            '2018-11-08_B2-006' + _s + '421.9A_ID=450.dat',
            '2018-11-16_B2-007' + _s + '381.7A_ID=536.dat',
            '2018-11-16_B2-007' + _s + '401.8A_ID=537.dat',
            '2018-11-16_B2-007' + _s + '421.9A_ID=538.dat',
            '2018-11-08_B2-008' + _s + '381.7A_ID=444.dat',
            '2018-11-08_B2-008' + _s + '401.8A_ID=443.dat',
            '2018-11-08_B2-008' + _s + '421.9A_ID=445.dat',
            '2018-11-23_B2-009' + _s + '381.7A_ID=607.dat',
            '2018-11-23_B2-009' + _s + '401.8A_ID=608.dat',
            '2018-11-23_B2-009' + _s + '421.9A_ID=609.dat',
            '2018-11-17_B2-010' + _s + '381.7A_ID=548.dat',
            '2018-11-17_B2-010' + _s + '401.8A_ID=549.dat',
            '2018-11-17_B2-010' + _s + '421.9A_ID=550.dat',
            '2018-11-11_B2-011' + _s + '381.7A_ID=498.dat',
            '2018-11-11_B2-011' + _s + '401.8A_ID=499.dat',
            '2018-11-11_B2-011' + _s + '421.9A_ID=500.dat',
            '2018-11-07_B2-013' + _s + '381.7A_ID=433.dat',
            '2018-11-07_B2-013' + _s + '401.8A_ID=434.dat',
            '2018-11-07_B2-013' + _s + '421.9A_ID=435.dat',
            '2018-11-30_B2-014' + _s + '381.7A_ID=641.dat',
            '2018-11-30_B2-014' + _s + '401.8A_ID=642.dat',
            '2018-11-30_B2-014' + _s + '421.9A_ID=643.dat',
            '2018-11-17_B2-015' + _s + '381.7A_ID=544.dat',
            '2018-11-17_B2-015' + _s + '401.8A_ID=545.dat',
            '2018-11-17_B2-015' + _s + '421.9A_ID=546.dat',
            '2018-11-21_B2-016' + _s + '381.7A_ID=576.dat',
            '2018-11-21_B2-016' + _s + '401.8A_ID=577.dat',
            '2018-11-21_B2-016' + _s + '421.9A_ID=578.dat',
            '2018-11-24_B2-017' + _s + '381.7A_ID=613.dat',
            '2018-11-24_B2-017' + _s + '401.8A_ID=614.dat',
            '2018-11-24_B2-017' + _s + '421.9A_ID=615.dat',
            '2018-11-24_B2-018' + _s + '381.7A_ID=617.dat',
            '2018-11-24_B2-018' + _s + '401.8A_ID=618.dat',
            '2018-11-24_B2-018' + _s + '421.9A_ID=619.dat',
            '2018-11-05_B2-019' + _s + '381.7A_ID=408.dat',
            '2018-11-05_B2-019' + _s + '401.8A_ID=409.dat',
            '2018-11-05_B2-019' + _s + '421.9A_ID=410.dat',
            '2018-11-30_B2-021' + _s + '381.7A_ID=649.dat',
            '2018-11-30_B2-021' + _s + '401.8A_ID=650.dat',
            '2018-11-30_B2-021' + _s + '421.9A_ID=651.dat',
            '2018-11-30_B2-022' + _s + '381.7A_ID=680.dat',
            '2018-11-30_B2-022' + _s + '401.8A_ID=681.dat',
            '2018-11-30_B2-022' + _s + '421.9A_ID=683.dat',
            '2018-11-26_B2-023' + _s + '381.7A_ID=627.dat',
            '2018-11-26_B2-023' + _s + '401.8A_ID=628.dat',
            '2018-11-26_B2-023' + _s + '421.9A_ID=629.dat',
            '2018-11-03_B2-024' + _s + '381.7A_ID=390.dat',
            '2018-11-03_B2-024' + _s + '401.8A_ID=391.dat',
            '2018-11-03_B2-024' + _s + '421.9A_ID=392.dat',
            '2018-11-06_B2-025' + _s + '381.7A_ID=423.dat',
            '2018-11-06_B2-025' + _s + '401.8A_ID=424.dat',
            '2018-11-06_B2-025' + _s + '421.9A_ID=425.dat',
            '2018-11-10_B2-026' + _s + '381.7A_ID=484.dat',
            '2018-11-10_B2-026' + _s + '401.8A_ID=485.dat',
            '2018-11-10_B2-026' + _s + '421.9A_ID=486.dat',
            '2018-11-23_B2-027' + _s + '381.7A_ID=603.dat',
            '2018-11-23_B2-027' + _s + '401.8A_ID=604.dat',
            '2018-11-23_B2-027' + _s + '421.9A_ID=605.dat',
            '2018-11-21_B2-028' + _s + '381.7A_ID=563.dat',
            '2018-11-21_B2-028' + _s + '401.8A_ID=564.dat',
            '2018-11-21_B2-028' + _s + '421.9A_ID=565.dat',
            '2018-11-30_B2-029' + _s + '381.7A_ID=669.dat',
            '2018-11-30_B2-029' + _s + '401.8A_ID=670.dat',
            '2018-11-30_B2-029' + _s + '421.9A_ID=671.dat',
            '2018-11-04_B2-030' + _s + '381.7A_ID=399.dat',
            '2018-11-04_B2-030' + _s + '401.8A_ID=400.dat',
            '2018-11-04_B2-030' + _s + '421.9A_ID=401.dat',
            '2018-10-30_B2-031' + _s + '381.7A_ID=366.dat',
            '2018-10-30_B2-031' + _s + '401.8A_ID=367.dat',
            '2018-10-30_B2-031' + _s + '421.9A_ID=368.dat',
            '2018-11-22_B2-032' + _s + '381.7A_ID=583.dat',
            '2018-11-22_B2-032' + _s + '401.8A_ID=584.dat',
            '2018-11-22_B2-032' + _s + '421.9A_ID=585.dat',
            '2018-11-22_B2-033' + _s + '381.7A_ID=589.dat',
            '2018-11-22_B2-033' + _s + '401.8A_ID=590.dat',
            '2018-11-22_B2-033' + _s + '421.9A_ID=591.dat',
            '2018-11-24_B2-034' + _s + '381.7A_ID=621.dat',
            '2018-11-24_B2-034' + _s + '401.8A_ID=622.dat',
            '2018-11-24_B2-034' + _s + '421.9A_ID=623.dat',
            '2018-11-30_B2-036' + _s + '381.7A_ID=657.dat',
            '2018-11-30_B2-036' + _s + '401.8A_ID=658.dat',
            '2018-11-30_B2-036' + _s + '421.9A_ID=659.dat',
            '2018-11-09_B2-037' + _s + '381.7A_ID=459.dat',
            '2018-11-09_B2-037' + _s + '401.8A_ID=460.dat',
            '2018-11-09_B2-037' + _s + '421.9A_ID=461.dat',
            '2018-11-30_B2-038' + _s + '381.7A_ID=634.dat',
            '2018-11-30_B2-038' + _s + '401.8A_ID=635.dat',
            '2018-11-30_B2-038' + _s + '421.9A_ID=636.dat',
            '2018-11-30_B2-040' + _s + '381.7A_ID=653.dat',
            '2018-11-30_B2-040' + _s + '401.8A_ID=654.dat',
            '2018-11-30_B2-040' + _s + '421.9A_ID=655.dat',
            '2018-11-30_B2-042' + _s + '381.7A_ID=663.dat',
            '2018-11-30_B2-042' + _s + '401.8A_ID=664.dat',
            '2018-11-30_B2-042' + _s + '421.9A_ID=665.dat',
            '2018-11-14_B2-043' + _s + '381.7A_ID=518.dat',
            '2018-11-14_B2-043' + _s + '401.8A_ID=519.dat',
            '2018-11-14_B2-043' + _s + '421.9A_ID=520.dat',
            '2018-11-05_B2-044' + _s + '381.7A_ID=404.dat',
            '2018-11-05_B2-044' + _s + '401.8A_ID=405.dat',
            '2018-11-05_B2-044' + _s + '421.9A_ID=406.dat',
            '2018-11-09_B2-045' + _s + '381.7A_ID=464.dat',
            '2018-11-09_B2-045' + _s + '401.8A_ID=465.dat',
            '2018-11-09_B2-045' + _s + '421.9A_ID=466.dat',
            '2018-11-21_B2-046' + _s + '381.7A_ID=570.dat',
            '2018-11-21_B2-046' + _s + '401.8A_ID=571.dat',
            '2018-11-21_B2-046' + _s + '421.9A_ID=572.dat',
            )

    if magnet is None and current is None:
        return _fmap_files

    files = []
    for f in _fmap_files:
        if magnet is None or magnet in f:
            if current is None or current in f:
                files.append(f)
    return tuple(files)


def get_magnets_B2():
    """."""
    fnames = get_fmap_files_B2()
    magnets = []
    for f in fnames:
        magnet = f[11:11+6]
        if magnet not in magnets:
            magnets.append(magnet)
    return sorted(magnets)


def get_currents_B2():
    """."""
    fnames = get_fmap_files_B2()
    currents = []
    for f in fnames:
        current = f[53:53+6]
        if current not in currents:
            currents.append(current)
    return sorted(currents)


def get_fmap_files_B1(magnet=None, current=None):
    """."""
    global _fmap_files

    if _fmap_files is None:
        _s = '_Hall_X=-32_32mm_Z=-800_800mm_Imc='
        _fmap_files = (
             '2018-12-07_B1-002' + _s + '381.7A_ID=799.dat',
             '2018-12-07_B1-002' + _s + '401.8A_ID=800.dat',
             '2018-12-07_B1-002' + _s + '403.6A_ID=801.dat',
             '2018-12-07_B1-002' + _s + '421.9A_ID=802.dat',
             '2018-12-15_B1-003' + _s + '381.7A_ID=981.dat',
             '2018-12-15_B1-003' + _s + '401.8A_ID=982.dat',
             '2018-12-15_B1-003' + _s + '403.6A_ID=983.dat',
             '2018-12-15_B1-003' + _s + '421.9A_ID=984.dat',
             '2018-12-17_B1-004' + _s + '381.7A_ID=1021.dat',
             '2018-12-17_B1-004' + _s + '401.8A_ID=1022.dat',
             '2018-12-17_B1-004' + _s + '403.6A_ID=1023.dat',
             '2018-12-17_B1-004' + _s + '421.9A_ID=1024.dat',
             '2018-12-13_B1-005' + _s + '381.7A_ID=912.dat',
             '2018-12-13_B1-005' + _s + '401.8A_ID=913.dat',
             '2018-12-13_B1-005' + _s + '403.6A_ID=914.dat',
             '2018-12-13_B1-005' + _s + '421.9A_ID=915.dat',
             '2018-12-15_B1-006' + _s + '381.7A_ID=970.dat',
             '2018-12-15_B1-006' + _s + '401.8A_ID=971.dat',
             '2018-12-15_B1-006' + _s + '403.6A_ID=972.dat',
             '2018-12-15_B1-006' + _s + '421.9A_ID=973.dat',
             '2018-12-10_B1-009' + _s + '381.7A_ID=827.dat',
             '2018-12-10_B1-009' + _s + '401.8A_ID=828.dat',
             '2018-12-10_B1-009' + _s + '403.6A_ID=829.dat',
             '2018-12-10_B1-009' + _s + '421.9A_ID=830.dat',
             '2018-12-07_B1-010' + _s + '381.7A_ID=794.dat',
             '2018-12-07_B1-010' + _s + '401.8A_ID=795.dat',
             '2018-12-07_B1-010' + _s + '403.6A_ID=796.dat',
             '2018-12-07_B1-010' + _s + '421.9A_ID=797.dat',
             '2018-12-05_B1-011' + _s + '381.7A_ID=756.dat',
             '2018-12-05_B1-011' + _s + '401.8A_ID=757.dat',
             '2018-12-05_B1-011' + _s + '403.6A_ID=758.dat',
             '2018-12-05_B1-011' + _s + '421.9A_ID=759.dat',
             '2018-12-17_B1-012' + _s + '381.7A_ID=1007.dat',
             '2018-12-17_B1-012' + _s + '401.8A_ID=1008.dat',
             '2018-12-17_B1-012' + _s + '403.6A_ID=1009.dat',
             '2018-12-17_B1-012' + _s + '421.9A_ID=1010.dat',
             '2018-12-10_B1-013' + _s + '381.7A_ID=840.dat',
             '2018-12-10_B1-013' + _s + '401.8A_ID=841.dat',
             '2018-12-10_B1-013' + _s + '403.6A_ID=842.dat',
             '2018-12-10_B1-013' + _s + '421.9A_ID=843.dat',
             '2018-12-08_B1-014' + _s + '381.7A_ID=815.dat',
             '2018-12-08_B1-014' + _s + '401.8A_ID=816.dat',
             '2018-12-08_B1-014' + _s + '403.6A_ID=817.dat',
             '2018-12-08_B1-014' + _s + '421.9A_ID=818.dat',
             '2018-12-16_B1-015' + _s + '381.7A_ID=1001.dat',
             '2018-12-16_B1-015' + _s + '401.8A_ID=1002.dat',
             '2018-12-16_B1-015' + _s + '403.6A_ID=1003.dat',
             '2018-12-16_B1-015' + _s + '421.9A_ID=1004.dat',
             '2018-12-11_B1-016' + _s + '381.7A_ID=856.dat',
             '2018-12-11_B1-016' + _s + '401.8A_ID=857.dat',
             '2018-12-11_B1-016' + _s + '403.6A_ID=858.dat',
             '2018-12-11_B1-016' + _s + '421.9A_ID=859.dat',
             '2018-12-03_B1-017' + _s + '381.7A_ID=712.dat',
             '2018-12-03_B1-017' + _s + '401.8A_ID=713.dat',
             '2018-12-03_B1-017' + _s + '403.6A_ID=714.dat',
             '2018-12-03_B1-017' + _s + '421.9A_ID=715.dat',
             '2018-12-15_B1-018' + _s + '381.7A_ID=963.dat',
             '2018-12-15_B1-018' + _s + '401.8A_ID=964.dat',
             '2018-12-15_B1-018' + _s + '403.6A_ID=965.dat',
             '2018-12-15_B1-018' + _s + '421.9A_ID=966.dat',
             '2018-12-10_B1-019' + _s + '381.7A_ID=822.dat',
             '2018-12-10_B1-019' + _s + '401.8A_ID=823.dat',
             '2018-12-10_B1-019' + _s + '403.6A_ID=824.dat',
             '2018-12-10_B1-019' + _s + '421.9A_ID=825.dat',
             '2018-12-14_B1-020' + _s + '381.7A_ID=947.dat',
             '2018-12-14_B1-020' + _s + '401.8A_ID=948.dat',
             '2018-12-14_B1-020' + _s + '403.6A_ID=949.dat',
             '2018-12-14_B1-020' + _s + '421.9A_ID=950.dat',
             '2018-12-12_B1-021' + _s + '381.7A_ID=899.dat',
             '2018-12-12_B1-021' + _s + '401.8A_ID=900.dat',
             '2018-12-12_B1-021' + _s + '403.6A_ID=901.dat',
             '2018-12-12_B1-021' + _s + '421.9A_ID=902.dat',
             '2018-12-04_B1-022' + _s + '381.7A_ID=720.dat',
             '2018-12-04_B1-022' + _s + '401.8A_ID=721.dat',
             '2018-12-04_B1-022' + _s + '403.6A_ID=722.dat',
             '2018-12-04_B1-022' + _s + '421.9A_ID=723.dat',
             '2018-12-12_B1-023' + _s + '381.7A_ID=890.dat',
             '2018-12-12_B1-023' + _s + '401.8A_ID=891.dat',
             '2018-12-12_B1-023' + _s + '403.6A_ID=892.dat',
             '2018-12-12_B1-023' + _s + '421.9A_ID=893.dat',
             '2018-12-16_B1-024' + _s + '381.7A_ID=988.dat',
             '2018-12-16_B1-024' + _s + '401.8A_ID=989.dat',
             '2018-12-16_B1-024' + _s + '403.6A_ID=990.dat',
             '2018-12-16_B1-024' + _s + '421.9A_ID=991.dat',
             '2018-12-15_B1-025' + _s + '381.7A_ID=958.dat',
             '2018-12-15_B1-025' + _s + '401.8A_ID=959.dat',
             '2018-12-15_B1-025' + _s + '403.6A_ID=960.dat',
             '2018-12-15_B1-025' + _s + '421.9A_ID=961.dat',
             '2018-12-05_B1-026' + _s + '381.7A_ID=747.dat',
             '2018-12-05_B1-026' + _s + '401.8A_ID=748.dat',
             '2018-12-05_B1-026' + _s + '403.6A_ID=749.dat',
             '2018-12-05_B1-026' + _s + '421.9A_ID=750.dat',
             '2018-12-15_B1-027' + _s + '381.7A_ID=976.dat',
             '2018-12-15_B1-027' + _s + '401.8A_ID=977.dat',
             '2018-12-15_B1-027' + _s + '403.6A_ID=978.dat',
             '2018-12-15_B1-027' + _s + '421.9A_ID=979.dat',
             '2018-12-13_B1-028' + _s + '381.7A_ID=922.dat',
             '2018-12-13_B1-028' + _s + '401.8A_ID=923.dat',
             '2018-12-13_B1-028' + _s + '403.6A_ID=924.dat',
             '2018-12-13_B1-028' + _s + '421.9A_ID=925.dat',
             '2018-12-16_B1-029' + _s + '381.7A_ID=996.dat',
             '2018-12-16_B1-029' + _s + '401.8A_ID=997.dat',
             '2018-12-16_B1-029' + _s + '403.6A_ID=998.dat',
             '2018-12-16_B1-029' + _s + '421.9A_ID=999.dat',
             '2018-12-10_B1-030' + _s + '381.7A_ID=833.dat',
             '2018-12-10_B1-030' + _s + '401.8A_ID=834.dat',
             '2018-12-10_B1-030' + _s + '403.6A_ID=835.dat',
             '2018-12-10_B1-030' + _s + '421.9A_ID=836.dat',
             '2018-12-11_B1-031' + _s + '381.7A_ID=848.dat',
             '2018-12-11_B1-031' + _s + '401.8A_ID=849.dat',
             '2018-12-11_B1-031' + _s + '403.6A_ID=850.dat',
             '2018-12-11_B1-031' + _s + '421.9A_ID=851.dat',
             '2018-12-06_B1-032' + _s + '381.7A_ID=783.dat',
             '2018-12-06_B1-032' + _s + '401.8A_ID=784.dat',
             '2018-12-06_B1-032' + _s + '403.6A_ID=785.dat',
             '2018-12-06_B1-032' + _s + '421.9A_ID=786.dat',
             '2018-12-17_B1-033' + _s + '381.7A_ID=1016.dat',
             '2018-12-17_B1-033' + _s + '401.8A_ID=1017.dat',
             '2018-12-17_B1-033' + _s + '403.6A_ID=1018.dat',
             '2018-12-17_B1-033' + _s + '421.9A_ID=1019.dat',
             '2018-12-07_B1-034' + _s + '381.7A_ID=804.dat',
             '2018-12-07_B1-034' + _s + '401.8A_ID=805.dat',
             '2018-12-07_B1-034' + _s + '403.6A_ID=806.dat',
             '2018-12-07_B1-034' + _s + '421.9A_ID=807.dat',
             '2018-12-13_B1-035' + _s + '381.7A_ID=930.dat',
             '2018-12-13_B1-035' + _s + '401.8A_ID=931.dat',
             '2018-12-13_B1-035' + _s + '403.6A_ID=932.dat',
             '2018-12-13_B1-035' + _s + '421.9A_ID=933.dat',
             '2018-12-13_B1-036' + _s + '381.7A_ID=917.dat',
             '2018-12-13_B1-036' + _s + '401.8A_ID=918.dat',
             '2018-12-13_B1-036' + _s + '403.6A_ID=919.dat',
             '2018-12-13_B1-036' + _s + '421.9A_ID=920.dat',
             '2018-12-12_B1-037' + _s + '381.7A_ID=883.dat',
             '2018-12-12_B1-037' + _s + '401.8A_ID=884.dat',
             '2018-12-12_B1-037' + _s + '403.6A_ID=885.dat',
             '2018-12-12_B1-037' + _s + '421.9A_ID=886.dat',
             '2018-12-08_B1-038' + _s + '381.7A_ID=809.dat',
             '2018-12-08_B1-038' + _s + '401.8A_ID=810.dat',
             '2018-12-08_B1-038' + _s + '403.6A_ID=811.dat',
             '2018-12-08_B1-038' + _s + '421.9A_ID=812.dat',
             '2018-12-06_B1-039' + _s + '381.7A_ID=776.dat',
             '2018-12-06_B1-039' + _s + '401.8A_ID=777.dat',
             '2018-12-06_B1-039' + _s + '403.6A_ID=778.dat',
             '2018-12-06_B1-039' + _s + '421.9A_ID=779.dat',
             '2018-12-05_B1-040' + _s + '381.7A_ID=767.dat',
             '2018-12-05_B1-040' + _s + '401.8A_ID=768.dat',
             '2018-12-05_B1-040' + _s + '403.6A_ID=769.dat',
             '2018-12-05_B1-040' + _s + '421.9A_ID=770.dat',
             '2018-12-12_B1-041' + _s + '381.7A_ID=905.dat',
             '2018-12-12_B1-041' + _s + '401.8A_ID=906.dat',
             '2018-12-12_B1-041' + _s + '403.6A_ID=907.dat',
             '2018-12-12_B1-041' + _s + '421.9A_ID=908.dat',
             '2018-12-05_B1-042' + _s + '381.7A_ID=741.dat',
             '2018-12-05_B1-042' + _s + '401.8A_ID=742.dat',
             '2018-12-05_B1-042' + _s + '403.6A_ID=743.dat',
             '2018-12-05_B1-042' + _s + '421.9A_ID=744.dat',
             '2018-12-14_B1-043' + _s + '381.7A_ID=939.dat',
             '2018-12-14_B1-043' + _s + '401.8A_ID=940.dat',
             '2018-12-14_B1-043' + _s + '403.6A_ID=941.dat',
             '2018-12-14_B1-043' + _s + '421.9A_ID=942.dat',
             '2018-12-04_B1-046' + _s + '381.7A_ID=735.dat',
             '2018-12-04_B1-046' + _s + '401.8A_ID=736.dat',
             '2018-12-04_B1-046' + _s + '403.6A_ID=737.dat',
             '2018-12-04_B1-046' + _s + '421.9A_ID=738.dat',
    )

    if magnet is None and current is None:
        return _fmap_files

    files = []
    for f in _fmap_files:
        if magnet is None or magnet in f:
            if current is None or current in f:
                files.append(f)
    return tuple(files)


def get_magnets_B1():
    """."""
    fnames = get_fmap_files_B1()
    magnets = []
    for f in fnames:
        magnet = f[11:11+6]
        if magnet not in magnets:
            magnets.append(magnet)
    return sorted(magnets)


def get_currents_B1():
    """."""
    fnames = get_fmap_files_B1()
    currents = []
    for f in fnames:
        current = f[51:51+6]
        if current not in currents:
            currents.append(current)
    return sorted(currents)


def load_analysis_result(folder, dipole_type, plots=None):
    """."""
    if dipole_type == 'B1':
        b2d = defaults['si-dipoles-b1']
        magnets = get_magnets_B1()
        currents = get_currents_B1()
    else:
        b2d = defaults['si-dipoles-b2']
        magnets = get_magnets_B2()
        currents = get_currents_B2()

    path_base = b2d['_path_base'] + b2d['_path_repo'] + b2d['_path_model'] + \
                'analysis/hallprobe/production/' + folder

    data = {
        'angle': dict(),
        'angleP': dict(),
        'angleN': dict(),
        'dangle': dict(),
        'energy': dict(),
        'rx0': dict(),
        'refrxP': dict(),
        'refrxN': dict(),
        'dip_normal': dict(),
        'dip_skew': dict(),
        'quad_normal': dict(),
        'quad_skew': dict(),
        'sext_normal': dict(),
        'sext_skew': dict(),
        'dip_normalP': dict(),
        'dip_skewP': dict(),
        'quad_normalP': dict(),
        'quad_skewP': dict(),
        'sext_normalP': dict(),
        'sext_skewP': dict(),
        'dip_normalN': dict(),
        'dip_skewN': dict(),
        'quad_normalN': dict(),
        'quad_skewN': dict(),
        'sext_normalN': dict(),
        'sext_skewN': dict(),
        }
    convd = {
        'initial rx position of trajectory': 'rx0',
        'beam_energy': 'energy'}
    for current in currents:
        for p in data:
            data[p][current] = []
        for magnet in magnets:
            path = path_base + magnet + '/' + current.replace('.', 'p') + '/'
            if dipole_type == 'B1':
                fname = get_fmap_files_B1(magnet, current)[0]
            else:
                fname = get_fmap_files_B2(magnet, current)[0]
            f = DoubleFMapAnalysis(magnet=magnet, fmap_fname=fname)
            dataP = f.load_output_trajectory(path, 'z-positive/')
            dataN = f.load_output_trajectory(path, 'z-negative/')

            angle = \
                dataP['horizontal_deflection_angle'] - \
                dataN['horizontal_deflection_angle']
            data['angle'][current].append(angle)
            data['angleP'][current].append(
                dataP['horizontal_deflection_angle'])
            data['angleN'][current].append(
                dataN['horizontal_deflection_angle'])
            # v = dataP['rx position of reference'] - f.model_nominal_refrx
            v = dataP['rx position of reference']
            data['refrxP'][current].append(v)
            # v = dataN['rx position of reference'] - f.model_nominal_refrx
            v = dataN['rx position of reference']
            data['refrxN'][current].append(v)

            for k in convd:
                dP = dataP[k]
                dN = dataN[k]
            if dP == dN:
                data[convd[k]][current].append(dP)
            else:
                raise ValueError('{} != {}'.format(dP, dN))

            mpoles_normalP, mpoles_skewP, harms = \
                f.load_output_multipoles(path, 'z-positive/')
            mpoles_normalN, mpoles_skewN, harms = \
                f.load_output_multipoles(path, 'z-negative/')
            mpoles_normal = _np.array(mpoles_normalP) + _np.array(mpoles_normalN)
            mpoles_skew = _np.array(mpoles_skewP) + _np.array(mpoles_skewN)
            data['dip_normal'][current].append(mpoles_normal[0])
            data['dip_skew'][current].append(mpoles_skew[0])
            data['quad_normal'][current].append(mpoles_normal[1])
            data['quad_skew'][current].append(mpoles_skew[1])
            data['sext_normal'][current].append(mpoles_normal[2])
            data['sext_skew'][current].append(mpoles_skew[2])

            data['dip_normalP'][current].append(mpoles_normalP[0])
            data['dip_skewP'][current].append(mpoles_skewP[0])
            data['quad_normalP'][current].append(mpoles_normalP[1])
            data['quad_skewP'][current].append(mpoles_skewP[1])
            data['sext_normalP'][current].append(mpoles_normalP[2])
            data['sext_skewP'][current].append(mpoles_skewP[2])

            data['dip_normalN'][current].append(mpoles_normalN[0])
            data['dip_skewN'][current].append(mpoles_skewN[0])
            data['quad_normalN'][current].append(mpoles_normalN[1])
            data['quad_skewN'][current].append(mpoles_skewN[1])
            data['sext_normalN'][current].append(mpoles_normalN[2])
            data['sext_skewN'][current].append(mpoles_skewN[2])

    for current in currents:
        a = _np.array(data['angle'][current])
        am = _np.mean(a)
        # dm = 100*(a - am)/am
        dm = 100*(a - (-f.model_nominal_angle))/((-f.model_nominal_angle))
        data['dangle'][current] = dm

    if not plots:
        return magnets, currents, data

    # --- energy ---
    if not plots or 'energy' in plots:
        for current in currents:
            print(current, _np.mean(data['energy'][current]))
            _plt.plot(data['energy'][current], 'o', label=current)
        _plt.legend()
        _plt.xlabel('Magnet index')
        _plt.ylabel('Energy [GeV]')
        _plt.show()

    # --- angle ---
    if not plots or 'angle' in plots:
        for current in currents:
            _plt.plot(data['angle'][current], 'o', label=current)
            n = len(data['angle'][current])
        _plt.plot([0, n-1], [-f.model_nominal_angle, ]*2, '--k', label='spec')
        _plt.legend()
        _plt.xlabel('Magnet index')
        _plt.ylabel('Angle [deg]')
        _plt.show()

    # --- angle error ---
    if not plots or 'dangle' in plots:
        dat = []
        for current in currents:
            v = data['dangle'][current]
            _plt.plot(v, 'o', label=current)
            n = len(v)
            print('current:{}, avg_angle_error:{:+f} %'.format(current,
                                                           _np.mean(v)))
            dat = dat + list(v)
        _plt.plot([0, n-1], [-0.05, -0.05], '--k', label='spec')
        _plt.plot([0, n-1], [+0.05, +0.05], '--k')
        _plt.legend()
        _plt.xlabel('Magnet index')
        _plt.ylabel('Angle Deviation from Average for each Current[%]')
        _plt.show()
    print('avg_angle_error: {:+f} %'.format(_np.mean(dat)))

    # --- quadrupolar error ---
    if not plots or 'quad' in plots:
        nominal_KL = b2d['model_nominal_KL']
        nominal_GL = dict()
        dat = []
        for current in currents:
            energy = _np.mean(data['energy'][current])
            brho, *_ = _mp.beam_optics.beam_rigidity(energy=energy)
            nominal_GL = -nominal_KL * brho
            GL = data['quad_normal'][current]
            v = 100*(GL - nominal_GL)/nominal_GL
            dat = dat + list(v)
            print('current:{}, avg_quad_error:{:+f} %'.format(current,
                                                           _np.mean(v)))
            _plt.plot(v, label=current)
            n = len(data['quad_normal'][current])
        print('avg_quad_error: {:+f} %'.format(_np.mean(dat)))
        _plt.plot([0, n-1], [-0.1, -0.1], '--k', label='spec')
        _plt.plot([0, n-1], [+0.1, +0.1], '--k')
        _plt.xlabel('Magnet index')
        _plt.ylabel('Quadrupolar Error [%]')
        _plt.grid()
        _plt.legend()
        _plt.show()

    # --- reference point rx variation ---
    if not plots or 'refrx' in plots:
        s1 = ['ob', 'or', 'og', 'oy']
        s2 = ['sb', 'sr', 'sg', 'sy']
        datP = []
        datN = []
        for i in range(len(currents)):
            v = 1e0*_np.array(data['refrxP'][currents[i]])
            datP = datP + list(v)
            _plt.plot(v, s1[i], label=currents[i] + ' P')
            v = 1e0*_np.array(data['refrxN'][currents[i]])
            datN = datN + list(v)
            _plt.plot(v, s2[i], label=currents[i] + ' N')
        print('average ref rx (pos): ', _np.mean(datP))
        print('average ref rx (neg): ', _np.mean(datN))
        print('average ref rx      : ', _np.mean(datN+datP))
        _plt.legend(fontsize=20)
        _plt.grid()
        _plt.xlabel('Magnet Index', fontsize=20)
        _plt.ylabel('Reference rx point [mm]', fontsize=20)
        _plt.show()


def search_for_deflection_angle(dipole_type):
    """."""
    if dipole_type == 'B1':
        init_energies = [2.840677756097561, 2.987685682926829, 3.0008746829268294, 3.134536902439024]
        magnets = get_magnets_B1()
        currents = get_currents_B1()
        init_rx = 8.285 + 0.313
    else:
        init_energies = [2.8426315121951222 , 2.990131219512195, 3.0131593524884295, 3.137303658536585]
        magnets = get_magnets_B2()
        currents = get_currents_B2()
        init_rx = 7.920 + 0.245 - 0.012

    fstr = 'magnet:{}, current:{} => nr_iter:{:02d}, energy:{:8.6f} GeV'
    for i in range(len(currents)):
        curr = currents[i]
        for magnet in magnets:
            if dipole_type == 'B1':
                files = get_fmap_files_B1(magnet, curr)
            else:
                files = get_fmap_files_B2(magnet, curr)
            fa = DoubleFMapAnalysis(magnet=magnet, fmap_fname=files[0])
            fa.traj_init_rx = init_rx
            fa.energy = init_energies[i]
            n = fa.search_energy(True)
            print(fstr.format(magnet, curr, n, fa.energy))


def search_for_deflection_angle_vary_x0(c2e, dipole_type):
    """."""
    if dipole_type == 'B1':
        magnets = get_magnets_B1()
        currents = get_currents_B1()
        init_rx = 8.285 + 0.313 - 0.012
    else:
        magnets = get_magnets_B2()
        currents = get_currents_B2()
        init_rx = 7.920 + 0.245 - 0.012

    fstr = 'magnet:{}, current:{} => nr_iter:{:02d}, energy:{:8.6f} GeV, x0:{:6.4f} mm'
    for i in range(len(currents)):
        curr = currents[i]
        for magnet in magnets:
            if dipole_type == 'B1':
                files = get_fmap_files_B1(magnet, curr)
            else:
                files = get_fmap_files_B2(magnet, curr)
            fa = DoubleFMapAnalysis(magnet=magnet, fmap_fname=files[0])
            fa.traj_init_rx = init_rx
            fa.energy = c2e[curr]
            n = fa.search_x0(True)
            print(fstr.format(magnet, curr, n, fa.energy, fa.traj_init_rx))


def load_search_deflection_angle_file(fname=None):
    """."""
    if fname is None:
        fname = 'search-energies.txt'
    with open(fname) as fp:
        text = fp.readlines()
    data = dict()
    for line in text:
        words = line.split(' ')
        magnet = words[0].split(':')[1].replace(',','')
        current = words[1].split(':')[1]
        energy = float(words[4].split(':')[1].split(' ')[0])
        if current not in data:
            data[current] = {magnet: energy}
        elif magnet not in data[current]:
            data[current][magnet] = energy
        else:
            raise ValueError()

    # --- average energy ---
    for current in data:
        avg = _np.mean(_np.array([data[current][m] for m in data[current]]))
        print('{} : {}'.format(current, avg))

    currents = tuple(data.keys())
    magnets = tuple(data[currents[0]].keys())
    energies = tuple(tuple(data[c][m] for m in magnets) for c in currents)
    return currents, magnets, energies


def plot_results_search_deflection_angle(fname):
    """."""
    currents, magnets, energies = load_search_deflection_angle_file(fname=fname)
    fst = 'current:{}  energy_avg:{:8.6f} Gev  energy_std:{:5.3f} %'
    for i in range(len(energies)):
        em = _np.mean(energies[i])
        es = _np.std(energies[i])
        ev = 100*(energies[i] - em)/em
        _plt.plot(ev, 'o', label=currents[i])
        print(fst.format(currents[i], em, 100*es/em))
    _plt.plot([0, len(ev)-1], [+0.05, ]*2, '--k', label='spec')
    _plt.plot([0, len(ev)-1], [-0.05, ]*2, '--k',)
    _plt.xlabel('Magnet Index')
    _plt.ylabel('Nominal deflection energy spread [%]')
    _plt.grid()
    _plt.legend()
    _plt.show()


def generate_inputs(c2e, x0, dipole_type='B1'):
    """."""
    folder = 'x0-' + x0 + 'mm/'
    init_rx = float(x0.replace('p','.'))
    if dipole_type == 'B1':
        b2d = defaults['si-dipoles-b1']
    else:
        b2d = defaults['si-dipoles-b2']
    path_base = b2d['_path_base'] + b2d['_path_repo'] + b2d['_path_model'] + \
                'analysis/hallprobe/production/' + folder
    if dipole_type == 'B1':
        magnets = get_magnets_B1()
        currents = get_currents_B1()
    else:
        magnets = get_magnets_B2()
        currents = get_currents_B2()
    for magnet in magnets:
        print('creating input files for magnet {}'.format(magnet))
        for curr in currents:
            path = path_base + magnet + '/' + curr.replace('.', 'p') + '/'
            if dipole_type == 'B1':
                fname = get_fmap_files_B1(magnet, curr)[0]
            else:
                fname = get_fmap_files_B2(magnet, curr)[0]
            f = DoubleFMapAnalysis(magnet=magnet, fmap_fname=fname)
            default_s_step = f.get_defaults()['traj_rk_s_step']
            f.energy = c2e[curr]
            f.traj_init_rx = init_rx
            # positive
            f.traj_rk_s_step = +abs(default_s_step)
            f.create_input(
                path=path, force=True, path_subdir='z-positive/')
            # negative
            f.traj_rk_s_step = -abs(default_s_step)
            f.create_input(
                path=path, force=True, path_subdir='z-negative/')


def save_readme_files(c2e, folder, dipole_type):
    """."""
    header = (
        'Sirius '+dipole_type+' Dipoles Integrated Principal Multipoles',
        '=================================================',
        '',
        'As calculated in SIDE-half Runge-Kutta trajectory,',
        ('defined by measured fieldmap with magnet excitated with current '
         'of CURRENT,'),
        'corresponding to nominal particle energy of ENERGY GeV.',
        '',
        ('  Dipole   |  Angle [°]   |  Dint [T.m]  |   Gint [T]   |  Sint '
         '[T/m]  |'),
        ('           |              |              |              '
         '|              |'))
    sfmt = ('|{0:^10s}| {1:^+12.5f} | {2:^+12.5f} | {3:^+12.5f} '
            '| {4:^+12.5f} |\n')

    magnets, currents, data = load_analysis_result(folder, dipole_type, False)

    for current in currents:
        for side in ('Zpositive', 'Znegative'):
            with open('README-' + current + '-' + side + '.md', 'w') as fp:
                # header
                for line in header:
                    line = line.replace('CURRENT', current)
                    line = line.replace('SIDE', side)
                    line = line.replace('ENERGY', str(c2e[current]))
                    fp.write(line + '\n')
                # data
                for i in range(len(magnets)):
                    magnet = magnets[i]
                    angle = data['angleP'][current][i] if \
                        side == 'Zpositive' else data['angleN'][current][i]
                    dip = data['dip_normalP'][current][i] if \
                        side == 'Zpositive' else \
                        data['dip_normalN'][current][i]
                    quad = data['quad_normalP'][current][i] if \
                        side == 'Zpositive' else \
                        data['quad_normalN'][current][i]
                    sext = data['sext_normalP'][current][i] if \
                        side == 'Zpositive' else \
                        data['sext_normalN'][current][i]
                    fp.write(sfmt.format(magnet, angle, dip, quad, sext))


def calc_average_angles(folder, dipole_type):
    """."""
    if dipole_type == 'B1':
        b2d = defaults['si-dipoles-b1']
        magnets = get_magnets_B1()
        currents = get_currents_B1()
    else:
        b2d = defaults['si-dipoles-b2']
        magnets = get_magnets_B2()
        currents = get_currents_B2()

    # def_angle = b2d['model_nominal_angle']
    path_base = b2d['_path_base'] + b2d['_path_repo'] + b2d['_path_model'] + \
                'analysis/hallprobe/production/' + folder

    data = dict()
    aP, aN = [], []
    for current in currents:
        for p in data:
            data[p][current] = []
        for magnet in magnets:
            path = path_base + magnet + '/' + current.replace('.', 'p') + '/'
            if dipole_type == 'B1':
                fname = get_fmap_files_B1(magnet, current)[0]
            else:
                fname = get_fmap_files_B2(magnet, current)[0]
            f = DoubleFMapAnalysis(magnet=magnet, fmap_fname=fname)
            lens, angsP = f.load_output_model(path, 'z-positive/')
            lens, angsN = f.load_output_model(path, 'z-negative/')
            aP.append(angsP)
            aN.append(angsN)
    aP = _np.array(aP)
    aN = _np.array(aN)

    # aritmetic average
    aP_avg = _np.mean(aP, 0)
    aN_avg = _np.mean(aN, 0)
    a_avg = 0.5*(aP_avg + aN_avg)
    # # harmonic average (curvature radius)
    # aP_avg = 1.0/_np.mean(1.0/aP, 0)
    # aN_avg = 1.0/_np.mean(1.0/aN, 0)
    # a_avg = 1.0/((1.0/aP_avg + 1.0/aN_avg)/2.0)

    # f = def_angle / a_avg
    # print('default angle: {:.16f} deg'.format(def_angle))
    print('Idx Len[m]       Angle[deg]')
    for i in range(len(a_avg)):
        print('{:02d}: {:.5f} {:.16f}'.format(i+1, lens[i], a_avg[i]))
        # print('{:.16f}, '.format(a_avg[i]), end='')
    # print('')
    # print(lens)
    return 1000*_np.array(lens), a_avg


def load_traj(fname):
    """."""
    with open(fname, 'r') as f:
        lines = f.readlines()
    rx, rz, s = [], [], []
    for line in lines:
        if '#' in line:
            continue
        words = line.split()
        s.append(float(words[0]))
        rx.append(float(words[1]))
        rz.append(float(words[3]))
    return rz, rx, s


def plot_trajectories(folder, dipole_type):
    """."""
    if dipole_type == 'B1':
        b2d = defaults['si-dipoles-b1']
        magnets = get_magnets_B1()
        currents = get_currents_B1()
        trajfname = 'trajectory-b1-pos.in'
    else:
        b2d = defaults['si-dipoles-b2']
        magnets = get_magnets_B2()
        currents = get_currents_B2()
        trajfname = 'trajectory-b2-pos.in'

    # def_angle = b2d['model_nominal_angle']
    path_base = b2d['_path_base'] + b2d['_path_repo'] + b2d['_path_model'] + \
                'analysis/hallprobe/production/' + folder

    c2c = {
        '381p7A' : 'b',
        '401p8A' : 'r',
        '403p6A' : 'g',
        '421p9A' : 'y',
    }

    # plot all dipole RK trajs
    data = dict()
    legs = dict()
    for current in currents:
        for p in data:
            data[p][current] = []
        for magnet in magnets:
            path = path_base + magnet + '/' + current.replace('.', 'p') + '/z-positive/trajectory.txt'
            rz, rx, _ = load_traj(path)
            for k in c2c:
                if k in path:
                    if k not in legs:
                        _plt.plot(rz, rx, c2c[k], label=k)
                        legs[k] = None
                    else:
                        _plt.plot(rz, rx, c2c[k])

    # plot average traj
    rz, rx, _ = load_traj('./production/'+folder+trajfname)
    _plt.plot(rz, rx, 'ok', label='avg')

    _plt.legend()
    _plt.show()


def calc_average_rk_traj(folder, dipole_type, plt=None):

    if dipole_type == 'B1':
        b2d = defaults['si-dipoles-b1']
        magnets = get_magnets_B1()
        currents = get_currents_B1()
        trajfname = 'trajectory-b1-pos.in'
    else:
        b2d = defaults['si-dipoles-b2']
        magnets = get_magnets_B2()
        currents = get_currents_B2()
        trajfname = 'trajectory-b2-pos.in'

    # def_angle = b2d['model_nominal_angle']
    path_base = b2d['_path_base'] + b2d['_path_repo'] + b2d['_path_model'] + \
                'analysis/hallprobe/production/' + folder

    data = dict()
    rz_avg, rz_avg = None, None
    for current in currents:
        for p in data:
            data[p][current] = []
        for magnet in magnets:
            path = path_base + magnet + '/' + current.replace('.', 'p') + '/z-positive/trajectory.txt'
            rz, rx, s = load_traj(path)
            if plt is not None:
                plt.plot(rz, rx, 'g')
            if rz_avg is None:
                rz_avg = _np.array(rz)
                rx_avg = _np.array(rx)
                n = 1
            else:
                rz_avg += rz
                rx_avg += rx
                n += 1
    rz_avg /= n
    rx_avg /= n
    return s, rx_avg, rz_avg


def gen_trajectory(x0, le, an, s_step, s_max):
    """."""
    # print(le)
    # print(an)
    pi = _math.pi
    an = (pi/180) * _np.array(an)
    rho = _np.array(le)/_np.array(an)
    s = [0.0]
    p = [_np.array([0,x0])]
    v = [_np.array([1,0])]
    C, S = _math.cos(pi/2), _math.sin(pi/2)
    r90 = _np.array([[C, +S], [-S, C]])
    for i in range(len(an)):
        u = _np.dot(r90, v[-1])
        # print(u)
        n = -u * rho[i]
        o = p[-1] - n;
        nps = _math.ceil(le[i]/s_step)
        # print(nps)
        av = _np.linspace(0, an[i], nps)
        v0 = v[-1]
        s0 = s[-1]
        for j in range(1, len(av)):
            C, S = _math.cos(av[j]), _math.sin(av[j])
            m = _np.array([[C, S], [-S, C]])
            nr = _np.dot(m, n)
            p.append(o + nr)
            v.append(_np.dot(m, v0))
            # s.append(s0 + av[j] * rho[i])
            s.append(s[-1] + s_step)
    while s[-1] < s_max:
        s.append(s[-1] + s_step)
        v.append(v[-1])
        p.append(p[-1] + v[-1]*s_step)
    rz = [va[0] for va in p]
    rx = [va[1] for va in p]
    pz = [va[0] for va in v]
    px = [va[1] for va in v]
    return s, rx, rz, px, pz


def save_trajectory(folder, dipole_type, s, rx, rz, px, pz):
    """."""
    fmts = '{:+.14e} {:+.14e} {:+.14e} {:+.14e} {:+.14e} {:+.14e} {:+.14e} \n'
    if dipole_type == 'B1':
        fnameP = 'production/' + folder + 'trajectory-b1-pos.in'
        fnameN = 'production/' + folder + 'trajectory-b1-neg.in'
    elif dipole_type == 'B2':
        fnameP = 'production/' + folder + 'trajectory-b2-pos.in'
        fnameN = 'production/' + folder + 'trajectory-b2-neg.in'

    with open(fnameP, 'w') as f:
        f.write('# trajectory\n')
        f.write('# s[mm] rx[mm] ry[mm] rz[mm] px py pz\n')
        for i in range(len(s)):
            f.write(fmts.format(s[i], rx[i], 0.0, rz[i], px[i], 0.0, pz[i]))

    with open(fnameN, 'w') as f:
        f.write('# trajectory\n')
        f.write('# s[mm] rx[mm] ry[mm] rz[mm] px py pz\n')
        for i in range(len(s)):
            f.write(fmts.format(-s[i], rx[i], 0.0, -rz[i], -px[i], 0.0, pz[i]))


def save_reference_trajectory(dipole_type):
    """."""
    if dipole_type == 'B2':
        x0_folder = 'x0-8p153mm/'
        trajfname = 'trajectory-b2-pos.in'
    elif dipole_type == 'B1':
        x0_folder = 'x0-8p527mm/'
        trajfname = 'trajectory-b1-pos.in'

    s_rk, rx_rk, rz_rk = calc_average_rk_traj(x0_folder, dipole_type, _plt)
    le, an = calc_average_angles(x0_folder, dipole_type)
    s, rx, rz, px, pz = gen_trajectory(rx_rk[0], le, an, (s_rk[1]-s_rk[0]), s_rk[-1])
    save_trajectory(x0_folder, dipole_type, s, rx, rz, px, pz)


def plot_reference_trajectory(dipole_type):
    """."""
    if dipole_type == 'B2':
        x0_folder = 'x0-8p153mm/'
        trajfname = 'trajectory-b2-pos.in'
    elif dipole_type == 'B1':
        x0_folder = 'x0-8p527mm/'
        trajfname = 'trajectory-b1-pos.in'
    s_rk, rx_rk, rz_rk = calc_average_rk_traj(x0_folder, dipole_type, _plt)
    le, an = calc_average_angles(x0_folder, dipole_type)
    s, rx, rz, px, pz = gen_trajectory(rx_rk[0], le, an, (s_rk[1]-s_rk[0]), s_rk[-1])
    z_avg, x_avg, _ = load_traj('./production/' + x0_folder + trajfname)


    _plt.plot(rz_rk, rx_rk, 'b', label='RK (Avg)')
    _plt.plot(rz, rx, 'r', label='SegModel (Avg)')
    _plt.xlabel('Z [mm]')
    _plt.ylabel('X [mm]')
    _plt.legend()
    _plt.show()
