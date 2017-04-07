import math
import numpy as np
from scipy import interpolate

class OutOfRange(Exception):
    pass
class OutOfRangeRx(OutOfRange):
    pass
class OutOfRangeRy(OutOfRange):
    pass
class OutOfRangeRz(OutOfRange):
    pass
class OutOfRangeRxMax(OutOfRangeRx):
    pass
class OutOfRangeRxMin(OutOfRangeRx):
    pass
class OutOfRangeRyMax(OutOfRangeRy):
    pass
class OutOfRangeRyMin(OutOfRangeRy):
    pass
class OutOfRangeRzMax(OutOfRangeRz):
    pass
class OutOfRangeRzMin(OutOfRangeRz):
    pass
class IrregularFieldMap(Exception):
    pass

min_rx = 0
max_rx = 0

class FieldMapSet:

    def __init__(self):
        self.fieldmaps = []

    def add(self, fieldmap):
        self.fieldmaps.append(fieldmap)

    def interpolate(self, rx_global, ry_global, rz_global):
        bx, by, bz = 0.0, 0.0, 0.0
        for fm in self.fieldmaps:

            try:
                tbx,tby,tbz = fm.interpolate(rx_global,ry_global,rz_global)
            except (fieldmap.OutOfRangeRx, fieldmap.OutOfRangeRxMin, fieldmap.OutOfRangeRy):
                raise TrackException('extrapolation at ' + str((rx_global,ry_global,rz_global)))
            except fieldmap.OutOfRangeRz:
                tbx,tby,tbz = 0.0,0.0,0.0

            #tbx, tby, tbz = fm.interpolate(rx_global, ry_global, rz_global)

            bx += tbx
            by += tby
            bz += tbz
        return (bx,by,bz)

    def interpolate_set(self, points):

        field = np.zeros(points.shape)
        for i in range(points.shape[1]):
            field[:,i] = self.interpolate(*points[:,i])
        return field

class FieldMap:

    def __init__(self,
                 fname = None,
                 field_function = None,
                 polyfit_exponents = None,
                 rotation = 0.0,
                 translation = (0,0),
                 rescale_factor = 1.0):

        self.filename = fname
        self.fieldmap_label = None
        self.field_function = field_function
        self.rotation = rotation
        self.translation = translation
        self.rescale_factor = rescale_factor

        self.rx,       self.ry,       self.rz        = None, None, None # unique and sorted 1D numpy arrays with values of corresponding coordinates
        self.rx_nrpts, self.ry_nrpts, self.rz_nrpts  = None, None, None # number of distinct values of each coordinate in the fieldmap
        self.rx_min,   self.ry_min,   self.rz_min    = None, None, None # minimum values of coordinates
        self.rx_max,   self.ry_max,   self.rz_max    = None, None, None # maximum values of coordinates
        self.rx_step,  self.step,     self.rz_step   = None, None, None # spacing between consecutive coordinate values

        if self.field_function is None and self.filename is None:
            raise IrregularFieldMap('source of field not defined!')
        if self.field_function is not None and self.filename is not None:
            raise IrregularFieldMap('two sources of field defined!')

        # bx, by and bz field components of the fieldmap
        # ----------------------------------------------
        # these are 3D fieldmaps stored as lists whose elements are 2D fieldmaps.
        # (one for each y plane. usually only one at the midplane y = 0)
        # 2D fieldmaps are two dimensional numpy arrays: the first index (row) runs through different x values
        # wheres the second index (column) runs through differency longitudinal z coordinate
        #
        # Ex:    If there is only the midplave 2D fieldmap, then
        #
        #        by_on_axis = self.by[0][self.rx_zero,:]
        #
        #        is a 1D numpy array containg the longitudinal vertical profile of the field at x == 0.
        #
        self.bx, self.by, self.bz = None, None, None

        # In order to estimate what field integrals are missing from the fielmap (outside its region)
        # (1/z)^n, n>=2 polynomial interpolations are performed in the vicinities of both upstream and downstream
        # limits of the fieldmap. The interpolated polynomial is used to extrapolate the asymptotic decaying of the
        # field and to thus obtain numerical estimates of the field integrals outside the fieldmap.
        # these estimates are stored in the lists below.
        #

        # reads fieldmap from givem filename
        if self.filename is not None:
            self.read_fieldmap(self.filename)

    def read_fieldmap(self, fname):

        with open(fname, 'r') as fp:
            content = fp.read()

        ''' finds index of data section start '''
        idx = content.find('Z[mm]')
        idx = content.find('\n', idx+1)
        idx = content.find('\n', idx+1)

        ''' data section '''
        raw_data = np.fromstring(content[idx+1:], dtype=float, sep=' ')
        data = raw_data.view()
        data.shape = (-1,6)
        # position data
        self.rx = np.unique(data[:,0])
        self.ry = np.unique(data[:,1])
        self.rz = np.unique(data[:,2])
        self.rx_min, self.rx_max, self.rx_nrpts = min(self.rx), max(self.rx), len(self.rx)
        self.ry_min, self.ry_max, self.ry_nrpts = min(self.ry), max(self.ry), len(self.ry)
        self.rz_min, self.rz_max, self.rz_nrpts = min(self.rz), max(self.rz), len(self.rz)
        if self.rx_nrpts * self.ry_nrpts * self.rz_nrpts != data.shape[0]:
            raise IrregularFieldMap('not a rectangular grid')
        self.rx_step = (self.rx_max - self.rx_min) / (self.rx_nrpts - 1.0) if self.rx_nrpts > 1 else 0.0
        self.ry_step = (self.ry_max - self.ry_min) / (self.ry_nrpts - 1.0) if self.ry_nrpts > 1 else 0.0
        self.rz_step = (self.rz_max - self.rz_min) / (self.rz_nrpts - 1.0) if self.rz_nrpts > 1 else 0.0
        self.rx_zero = np.where(self.rx == 0)[0][0]
        self.ry_zero = np.where(self.ry == 0)[0][0]
        try:
            self.rz_zero = np.where(self.rz == 0)[0][0]
        except:
            print('data points do not contain z=0 !')
            self.rz_zero = 0
            for i in range(len(self.rz)):
                if abs(self.rz[i]) < abs(self.rz[self.rz_zero]):
                    self.rz_zero = i

        # field data
        self.bx, self.by, self.bz = data[:,3].view(), data[:,4].view(), data[:,5].view()
        self.bx.shape, self.by.shape, self.bz.shape = (-1,self.rx_nrpts), (-1,self.rx_nrpts), (-1,self.rx_nrpts)

        # rescale field
        if self.rescale_factor != 1.0:
            self.bx *= self.rescale_factor
            self.by *= self.rescale_factor
            self.bz *= self.rescale_factor

        # lookup tables for field interpolation
        kind = 'cubic'
        self.bxf = interpolate.interp2d(self.rx, self.rz, self.bx, kind=kind)
        self.byf = interpolate.interp2d(self.rx, self.rz, self.by, kind=kind)
        self.bzf = interpolate.interp2d(self.rx, self.rz, self.bz, kind=kind)

        self.bx = [np.transpose(self.bx)]
        self.by = [np.transpose(self.by)]
        self.bz = [np.transpose(self.bz)]

        ''' header section '''
        lines = content[:idx].split('\n')
        for line in lines:

            # empty line or comment
            if not line or (line[0] == '#'):
                continue
            words = line.split()
            if not words:
                continue

            cmd = words[0].lower()
            if cmd == 'fieldmap_name:' or cmd == 'nome_do_mapa:':
                self.fieldmap_label = ' '.join(words[1:])
                continue
            if cmd == 'timestamp:' or cmd == 'data_hora:':
                self.timestamp = ' '.join(words[1:])
                continue
            if cmd == 'filename:' or cmd == 'nome_do_arquivo:':
                self.filename = ' '.join(words[1:])
                continue
            if cmd == 'nr_magnets:' or cmd == 'numero_de_imas:':
                self.nr_magnets = int(words[1])
                continue
            if cmd == 'magnet_name:' or cmd == 'nome_do_ima:':
                self.magnet_label = ' '.join(words[1:])
                continue
            if cmd == 'rotation[deg]:':
                self.rotation = float(words[1]) * (math.pi/180.0)
                continue
            if cmd == 'gap[mm]:':
                try:
                    self.gap = float(words[1]) #[mm]
                except ValueError:
                    self.gap = None
                continue
            if cmd == 'control_gap[mm]:' or cmd == 'gap_controle[mm]:':
                try:
                    self.control_gap = float(words[1]) #[mm]
                except:
                    self.control_gap = None
                continue
            if cmd == 'magnet_length[mm]:' or cmd == 'comprimento[mm]:':
                self.length = float(words[1]) #[mm]
                continue
            if cmd == 'current_main[a]:' or cmd == 'corrente[a]:':
                try:
                    self.current = words[1]#[A]
                except ValueError:
                    self.current = None
                continue
            if cmd == 'current_qs[a]:':
                try:
                    if float(words[1]) != 0:
                        self.current = words[1]#[A]
                    self.current_qs = words[1]
                except ValueError:
                    pass
                continue
            if cmd == 'current_ch[a]:':
                try:
                    if float(words[1]) != 0:
                        self.current = words[1]#[A]
                    self.current_ch = words[1]
                except ValueError:
                    pass
                continue
            if cmd == 'current_cv[a]:':
                try:
                    if float(words[1]) != 0:
                        self.current = words[1]#[A]
                    self.current_cv = words[1]
                except ValueError:
                    pass
                continue
            if cmd == 'ni_main[a.esp]:':
                try:
                    self.ni = float(words[1]) #[]
                except:
                    self.ni = None
                continue
            if cmd == 'ni_trim[a.esp]:':
                try:
                    self.ni_trim = float(words[1]) #[]
                except:
                    self.ni_trim = None
                continue
            if cmd == 'current_trim[a]:':
                try:
                    self.current_trim = words[1]#[A]
                except ValueError:
                    self.current_trim = None
                continue
            if cmd == 'corrente_aux[a]:':
                try:
                    self.current_aux = words[1]#[A]
                except ValueError:
                    self.current_aux = None
                continue
            if cmd == 'ni_aux[a.esp]:':
                try:
                    self.ni_trim = float(words[1]) #[]
                except:
                    self.ni_trim = None
                continue

    def interpolate_set(self, points):

        field = np.zeros(points.shape)
        for i in range(points.shape[1]):
            field[:,i] = self.interpolate(*points[:,i])
        return field

    def interpolate(self, rx_global, ry_global, rz_global):

        #global min_rx, max_rx

        # converts to local coordinate system
        C, S = math.cos(self.rotation), math.sin(self.rotation)
        rx =  C * (rx_global - self.translation[0]) + S * (rz_global - self.translation[1])
        ry =  ry_global
        rz = -S * (rx_global - self.translation[0]) + C * (rz_global - self.translation[1])


        # if rx < min_rx:
        #     min_rx = rx
        #     print(min_rx)
        #
        # if rx > max_rx:
        #     max_rx = rx
        #     print(max_rx)

        if rx < self.rx_min:
            #print('temporary in fieldmap.py!')
            raise OutOfRangeRxMin('rx = {0:f} < rx_min = {1:f} [mm]'.format(rx, self.rx_min))
            #return (0,0,0)

        if rx > self.rx_max:
            #print('temporary in fieldmap.py!')
            raise OutOfRangeRxMax('rx = {0:f} > rx_max = {1:f} [mm]'.format(rx, self.rx_max))
            #return (0,0,0)


        field = (self.bxf(rx, rz), self.byf(rx, rz), self.bzf(rx, rz))

        # converts field back to global coordinates
        bx =  C * field[0] - S * field[2]
        bz =  S * field[0] + C * field[2]
        field = (bx, field[1], bz)

        return field

    def __str__(self):
        r = ''
        r += '{0:<35s} {1}'.format('timestamp:', self.timestamp)
        r += '\n{0:<35s} {1}'.format('filename:', self.filename)
        r += '\n{0:<35s} {1}'.format('magnet_label:', self.magnet_label)
        r += '\n{0:<35s} {1} mm'.format('magnet_length:', self.length)
        r += '\n{0:<35s} {1} A'.format('main_coil_current:', self.current)
        try:
            r += '\n{0:<35s} {1} A'.format('main_coil_NI:', self.ni)
        except:
            pass
        try:
            r += '\n{0:<35s} {1} A'.format('trim_coil_current:', self.current_trim)
        except:
            pass
        try:
            r += '\n{0:<35s} {1} A'.format('trim_coil_NI:', self.ni_trim)
        except:
            pass
        try:
            r += '\n{0:<35s} {1} mm'.format('magnetic_gap:', self.gap)
        except:
            pass
        try:
            r += '\n{0:<35s} {1} mm'.format('control_gap:', self.control_gap)
        except:
            pass

        if self.ry_nrpts == 1:
            r += '\n{0:<35s} {3} point in [{1},{2}] mm (step of {4:f} mm)'.format('ry_grid:', self.ry_min, self.ry_max, self.ry_nrpts, self.ry_step)
        else:
            r += '\n{0:<35s} {3} points in [{1},{2}] mm (step of {4:f} mm)'.format('ry_grid:', self.ry_min, self.ry_max, self.ry_nrpts, self.ry_step)
        if self.rx_nrpts == 1:
            r += '\n{0:<35s} {3} point in [{1},{2}] mm (step of {4:f} mm)'.format('rx_grid:', self.rx_min, self.rx_max, self.rx_nrpts, self.rx_step)
        else:
            r += '\n{0:<35s} {3} points in [{1},{2}] mm (step of {4:f} mm)'.format('rx_grid:', self.rx_min, self.rx_max, self.rx_nrpts, self.rx_step)
        if self.rz_nrpts == 1:
            r += '\n{0:<35s} {3} point in [{1},{2}] mm (step of {4:f} mm)'.format('rz_grid:', self.rz_min, self.rz_max, self.rz_nrpts, self.rz_step)
        else:
            r += '\n{0:<35s} {3} points in [{1},{2}] mm (step of {4:f} mm)'.format('rz_grid:', self.rz_min, self.rz_max, self.rz_nrpts, self.rz_step)
        r += '\n{0:<35s} (min:{1:+8.5f} max:{2:+8.5f}) (min:{3:+8.5f} max:{4:+8.5f}) Tesla'.format('by@(all)(axis):',
        np.amin(self.by[self.ry_zero]), np.amax(self.by[self.ry_zero]), min(self.by[self.ry_zero][self.rx_zero]), max(self.by[self.ry_zero][self.rx_zero]))
        r += '\n{0:<35s} (min:{1:+8.5f} max:{2:+8.5f}) (min:{3:+8.5f} max:{4:+8.5f}) Tesla'.format('bx@(all)(axis):',
        np.amin(self.bx[self.ry_zero]), np.amax(self.bx[self.ry_zero]), min(self.bx[self.ry_zero][self.rx_zero]), max(self.bx[self.ry_zero][self.rx_zero]))
        r += '\n{0:<35s} (min:{1:+8.5f} max:{2:+8.5f}) (min:{3:+8.5f} max:{4:+8.5f}) Tesla'.format('bz@(all)(axis):',
        np.amin(self.bz[self.ry_zero]), np.amax(self.bz[self.ry_zero]), min(self.bz[self.ry_zero][self.rx_zero]), max(self.bz[self.ry_zero][self.rx_zero]))

        return r
