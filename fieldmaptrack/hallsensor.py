#!/usr/bin/env python-sirius
"""Run analysis."""

import math as _math
import fieldmaptrack as _fmap
from copy import deepcopy as _dcopy


defaults = {
    'si-dipoles-b2': {
        # paths
        '_path_repo': 'si-dipoles-b2/',
        '_path_model': 'model-08/',
        '_path_measf': 'measurement/magnetic/hallprobe/',
        '_path_dataset': '',
        # rawdata
        'magnet_type': 'dipole',
        # trajectory
        'model_nominal_angle': 4.0964,
        'traj_init_rx': 7.92,  # [mm]
        'traj_rk_s_step': 0.1,  # [mm]
         }}


class FMapAnalysis(_fmap.common_analysis.Config):
    """."""

    _path_base = '/home/fac_files/lnls-ima/'

    def __init__(self, magnet=None, fmap_fname=None):
        """."""
        super().__init__(fname=None)
        self._set_defaults(magnet)
        self.fmap_filename = self._get_fname(magnet, fmap_fname)
        # -- default rawfield parameters
        self.config_label = magnet
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

    def analysis_rawfield(self):
        """."""
        analysis = _fmap.common_analysis.get_analysis_symbol(self.magnet_type)
        analysis.raw_fieldmap_analysis(self)

    def analysis_trajectory(self):
        """."""
        analysis = _fmap.common_analysis.get_analysis_symbol(self.magnet_type)
        analysis.trajectory_analysis(self)

    def search_energy(self, print_flag=False, angle_tol=0.001):
        """."""
        def calc_dangle():
            angle = abs(self.deflection_angle)
            dangle = 1000 * (_math.pi/180.0) * \
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
                    self.config_label, float(self.fmap.current), iter,
                    self.energy, angle, self.model_nominal_angle, dangle))
            f = angle / self.model_nominal_angle
            iter += 1
        return iter

    @property
    def deflection_angle(self):
        """."""
        return self._get_deflection_angle()

    @staticmethod
    def _calc_deflection_angle(traj):
        px1, pz1 = traj.px[0], traj.pz[0]
        px2, pz2 = traj.px[-1], traj.pz[-1]
        t1 = _math.atan(px1/pz1)
        t2 = _math.atan(px2/pz2)
        return (180/_math.pi) * (t2 - t1)

    def _get_deflection_angle(self):
        return self._calc_deflection_angle(self.traj)

    def _set_energy(self, value):
        self.beam_energy = value

    def _set_defaults(self, magnet):
        if 'B2-' in magnet:
            defs = defaults['si-dipoles-b2']
        for k, v in defs.items():
            setattr(self, k, v)

    def _get_fname(self, magnet, fmap_fname):
        self._path_magnet = magnet + '/'
        fname = \
            self._path_base + \
            self._path_repo + \
            self._path_model + \
            self._path_measf + \
            self._path_magnet + \
            self._path_dataset + \
            fmap_fname
        return fname


class DoubleFMapAnalysis(FMapAnalysis):
    """."""

    def __init__(self, **kwargs):
        """."""
        super().__init__(**kwargs)

    def analysis_trajectory(self):
        """."""
        super().analysis_trajectory()
        analysis = _fmap.common_analysis.get_analysis_symbol(self.magnet_type)
        if not hasattr(self, '_configN'):
            self._configN = _dcopy(self)
            self._configN.fmap = self.fmap  # to save space
            self._configN.traj_rk_s_step *= -1.0
        self._configN = analysis.trajectory_analysis(self._configN)

    def _get_deflection_angle(self):
        a1 = self._calc_deflection_angle(self.traj)
        a2 = self._calc_deflection_angle(self._configN.traj)
        return a1 - a2

    def _set_energy(self, value):
        self.beam_energy = value
        if hasattr(self, '_configN'):
            self._configN.beam_energy = value
