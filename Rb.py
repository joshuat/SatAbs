from math import pi
from Constants import u, c, h, eV_to_J

### Rubidium ###
class Rb():
    mass = 85.4678 * u   # kg (mass of rubidium)


class Transition():
    def __init__(self, E, linewidth):
        self.energy = E

        self.linewidth = linewidth

    def wavelength(self):
        return c/self.frequency()

    def frequency(self):
        return self.energy/h

    def Isat(self):
        return pi*h*c*self.linewidth/(3*self.wavelength()**3)

    def clone(self, freq_offset):
        return Transition( (self.frequency() + freq_offset)*h, self.linewidth )
        

# Rb85
class Rb85(Rb):
    mass = 84.911789732 * u # kg (from Steck)
    abundance = 0.7217

    D1 = Transition(1.559590695*eV_to_J,  # J - Energy
                    5.7500e6)           # Hz - linewidth
    D2 = Transition(1.589049139*eV_to_J,  # J - Energy
                    6.0666e6)           # Hz - linewidth


# Rb87
class Rb87(Rb):
    m = 87 * u
    abundance = 0.2783

    D1 = Transition(1.559591016*eV_to_J, # J - Energy
                    5.7500e6)            # Hz - linewidth
    D2 = Transition(1.589049462*eV_to_J, # J - Energy
                     6.0666e6)           # Hz - linewidth

    D2_F0 = D2.clone(-302.0738e6)         # MHz - freq offset
    D2_F1 = D2.clone(-229.8518e6)         # MHz - freq offset
    D2_F2 = D2.clone(-72.9112e6)          # MHz - freq offset
    D2_F3 = D2.clone(193.7407e6)          # MHz - freq offset


