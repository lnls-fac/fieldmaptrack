import os as _os
from fieldmaptrack.beam import Beam
from fieldmaptrack.fieldmap import FieldMap
from fieldmaptrack.fieldmap import FieldMapSet
from fieldmaptrack.track import SerretFrenetCoordSystem
from fieldmaptrack.track import Trajectory
from fieldmaptrack.track import TrackException
from fieldmaptrack.multipoles import Multipoles

from . import common_analysis

with open(_os.path.join(__path__[0], 'VERSION'), 'r') as _f:
    __version__ = _f.read().strip()
