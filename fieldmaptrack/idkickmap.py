"""IDKickmap module."""


import numpy as _np

from . import common_analysis as _common_analysis
from .beam import Beam as _Beam
from .fieldmap import FieldMap as _FieldMap
from .track import Trajectory as _Trajectory


class IDKickMap:
    """."""

    def __init__(self, fname=None):
        """."""
        self.id_length = None  # [m]
        self.posx = None  # [m]
        self.posy = None  # [m]
        self.kickx = None  # [T².m²]
        self.kicky = None  # [T².m²]
        self.fposx = None  # [m]
        self.fposy = None  # [m]
        self.fmap_config = None
        
        # load
        if fname:
            self.load(fname)
    
    def load(self, fname):
        """."""
        self.fname = fname
        self.id_length, \
        self.posx, self.posy, \
        self.kickx, self.kicky, \
        self.fposx, self.fposy = IDKickMap._load(self.fname)

    def fmap_calc_kickmap(self, fmap_fname, posx, posy, beam_energy=3.0, rk_s_step=0.2):
        """."""
        self.posx = _np.array(posx)
        self.posy = _np.array(posy)
        self.fmap_config = \
            IDKickMap._create_fmap_config(fmap_fname, beam_energy=beam_energy, rk_s_step=rk_s_step)
    
        brho = self.fmap_config.beam.brho
        self.id_length = (self.fmap_config.fmap.rz[-1] - self.fmap_config.fmap.rz[0])/1e3
        self.kickx = _np.full((len(self.posy), len(self.posx)), _np.inf)
        self.kicky = _np.full((len(self.posy), len(self.posx)), _np.inf)
        self.fposx = _np.full((len(self.posy), len(self.posx)), _np.inf)
        self.fposy = _np.full((len(self.posy), len(self.posx)), _np.inf)
        for i, ryi in enumerate(self.posy):
            for j, rxi in enumerate(self.posx):
                self.fmap_config = IDKickMap._calc_trajectory(self.fmap_config, init_rx=1e3*rxi, init_ry=1e3*ryi)
                pxf, pyf = self.fmap_config.traj.px[-1], self.fmap_config.traj.py[-1]
                rxf, ryf = self.fmap_config.traj.rx[-1], self.fmap_config.traj.ry[-1]
                print('rx = {:.01f} mm, ry = {:.01f} : px = {:.01f} urad, py = {:.01f} urad'.format(1e3*rxi, 1e3*ryi, 1e6*pxf, 1e6*pyf))
                self.kickx[i, j] = pxf * brho**2
                self.kicky[i, j] = pyf * brho**2
                self.fposx[i, j] = rxf / 1e3
                self.fposy[i, j] = ryf / 1e3

    @staticmethod
    def fmap_calc_trajectory(fmap_fname, init_rx, init_ry, init_px=0, init_py=0, beam_energy=3.0, rk_s_step=0.2):
        """."""
        fmap_config = \
            IDKickMap._create_fmap_config(fmap_fname, beam_energy=beam_energy, rk_s_step=rk_s_step)
        fmap_config = IDKickMap._calc_trajectory(fmap_config, init_rx, init_ry, init_px, init_py)
        return fmap_config
        
    def __str__(self):
        """."""

        rst = ''
        # header
        rst += '# Author: FAC script'
        rst += '\n# '
        rst += '\n# Total Length of Longitudinal Interval [m]'
        rst += '\n{}'.format(self.id_length)
        rst += '\n# Number of Horizontal Points'
        rst += '\n{}'.format(len(self.posx))
        rst += '\n# Number of Vertical Points'
        rst += '\n{}'.format(len(self.posy))
        
        rst += '\n# Total Horizontal 2nd Order Kick [T2m2]'
        rst += '\nSTART'
        # first line
        rst += '\n{:11s} '.format('')
        for rxi in self.posx:
            rst += '{:+011.5f} '.format(rxi)
        # table
        for i, ryi in enumerate(self.posy):
            rst += '\n{:+011.5f} '.format(ryi)
            for j, rxi in enumerate(self.posx):
                rst += '{:+11.4e} '.format(self.kickx[i, j])

        rst += '\n# Total Vertical 2nd Order Kick [T2m2]'
        rst += '\nSTART'
        # first line
        rst += '\n{:11s} '.format('')
        for rxi in self.posx:
            rst += '{:+011.5f} '.format(rxi)
        # table
        for i, ryi in enumerate(self.posy):
            rst += '\n{:+011.5f} '.format(ryi)
            for j, rxi in enumerate(self.posx):
                rst += '{:+11.4e} '.format(self.kicky[i, j])

        rst += '\n# Horizontal Final Position [m]'
        rst += '\nSTART'
        # first line
        rst += '\n{:11s} '.format('')
        for rxi in self.posx:
            rst += '{:+011.5f} '.format(rxi)
        # table
        for i, ryi in enumerate(self.posy):
            rst += '\n{:+011.5f} '.format(ryi)
            for j, rxi in enumerate(self.posx):
                rst += '{:+11.4e} '.format(self.fposx[i, j])

        rst += '\n# Vertical Final Position [m]'
        rst += '\nSTART'
        # first line
        rst += '\n{:11s} '.format('')
        for rxi in self.posx:
            rst += '{:+011.5f} '.format(rxi)
        # table
        for i, ryi in enumerate(self.posy):
            rst += '\n{:+011.5f} '.format(ryi)
            for j, rxi in enumerate(self.posx):
                rst += '{:+11.4e} '.format(self.fposy[i, j])

        return rst

    @staticmethod
    def _load(fname):
        """."""
        with open(fname) as fp:
            lines = fp.readlines()

        tables = []
        params = []
        for line in lines:
            line = line.strip()
            if line.startswith('START'):
                pass
            elif line.startswith('#'):
                pass
            else:
                data = [float(val) for val in line.split()]
                if len(data) == 1:
                    params.append(data[0])
                elif len(data) == int(params[1]):
                    posx = _np.array(data)
                else:
                    # print(data)
                    # return
                    tables.append(data)

        id_length = params[0]
        nrpts_y = int(params[2])
        tables = _np.array(tables)
        posy = tables[:nrpts_y, 0]
        tables = tables[:, 1:]
        
        kickx = tables[0*nrpts_y:1*nrpts_y, :]
        kicky = tables[1*nrpts_y:2*nrpts_y, :]
        fposx = tables[2*nrpts_y:3*nrpts_y, :]
        fposy = tables[3*nrpts_y:4*nrpts_y, :]

        return id_length, posx, posy, kickx, kicky, fposx, fposy

    @staticmethod
    def _create_fmap_config(fmap_fname, beam_energy, rk_s_step):
        config = _common_analysis.Config()
        config.config_label = 'id-3gev'
        config.magnet_type = 'insertion-device'
        config.fmap_filename = fmap_fname
        config.fmap_extrapolation_flag = False
        config.not_raise_range_exceptions = False
        config.interactive_mode = False
        config.not_raise_range_exceptions = True

        transforms = dict()
        config.fmap = _FieldMap(
            config.fmap_filename,
            transforms=transforms,
            not_raise_range_exceptions=config.not_raise_range_exceptions)

        config.traj_load_filename = None
        config.traj_is_reference_traj = True
        config.traj_init_rz = min(config.fmap.rz)
    
        config.traj_rk_s_step = rk_s_step
        config.traj_rk_length = None
        config.traj_rk_nrpts = None
        config.traj_force_midplane_flag = False
        config.interactive_mode = True

        config.beam_energy = beam_energy
    
        config.beam = _Beam(energy=config.beam_energy)
        config.traj = _Trajectory(
            beam=config.beam,
            fieldmap=config.fmap,
            not_raise_range_exceptions=config.not_raise_range_exceptions)
        config.traj_init_rz = min(config.fmap.rz)

        return config

    @staticmethod
    def _calc_trajectory(config,
        init_rx=0, init_ry=0,
        init_px=0, init_py=0,
        ):
        """."""
        config.traj_init_rx = init_rx
        config.traj_init_ry = init_ry
        config.traj_init_px = init_px
        config.traj_init_py = init_py
    
        analysis = _common_analysis.get_analysis_symbol(config.magnet_type)
        config = analysis.trajectory_analysis(config)

        return config
