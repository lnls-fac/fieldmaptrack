import math
import mathphys.constants as consts
import mathphys.units as units
import mathphys.constants as consts

class Beam:

    def __init__(self, energy, current = 0):

        self.energy = energy
        self.current = current
        self.brho, self.velocity, self.beta, self.gamma = Beam.calc_brho(self.energy)

    @staticmethod
    def calc_brho(energy):
        electron_rest_energy_GeV = units.joule_2_eV(consts.electron_rest_energy) / 1e9
        gamma    = energy/electron_rest_energy_GeV
        beta     = math.sqrt(((gamma-1.0)/gamma)*((gamma+1.0)/gamma))
        velocity = consts.light_speed * beta
        brho     = beta * (energy * 1e9) / consts.light_speed
        return brho, velocity, beta, gamma

    def __str__(self):
        r = ''
        r += '{0:<10s} {1:f} GeV'.format('energy:', self.energy)
        r += '\n{0:<10s} {1:f}'.format('gamma:', self.gamma)
        r += '\n{0:<10s} 1 - {1:e}'.format('beta:', 1.0-self.beta)
        r += '\n{0:<10s} {1:.0f} - {2:f} m/s'.format('velocity:', consts.light_speed, consts.light_speed - self.velocity)
        r += '\n{0:<10s} {1:f} T.m'.format('brho:', self.brho)
        return r
