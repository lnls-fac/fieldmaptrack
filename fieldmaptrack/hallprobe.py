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
from copy import deepcopy as _dcopy


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
        'traj_init_rx': 7.92,  # [mm]
        'traj_init_px': 0.0,   # [deg]
        'traj_rk_s_step': 0.1,  # [mm]
        'model_nominal_refrx': 19.239210,  # [mm]
        # multipoles
        'beam_energy': 3.0,
        'multipoles_normal_field_fitting_monomials': (0, 1, 2, 3, 4, 5, 6),
        'multipoles_skew_field_fitting_monomials': (),
        'multipoles_perpendicular_grid': '_np.linspace(-12,12,65)',
        'multipoles_r0': 17.5,  # [mm]
        'normalization_monomial': 0,
        'normalization_is_skew': False,
        'normal_multipoles_main_monomials': (0, 1, 2),
        'skew_multipoles_main_monomials': (),
        'model_segmentation': (196, 192, 182, 10, 10, 13, 17, 20, 30, 50), },

    'si-dipoles-b1': {}, }


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
            '{:10.8f} GeV - {:07.4f}|{:07.4f} deg: {:+09.4f} mrad'

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

    def get_defaults(self):
        """Return default attributes."""
        if 'B2-' in self._magnet:
            return defaults['si-dipoles-b2']

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
                raise FileExistsError()
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
