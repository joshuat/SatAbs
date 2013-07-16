from math import pi

### Physical Constants ###
c = 299792458.       # m/s (Speed of light)
kB = 1.3806488e-23  # m^2 kg s^-2 K^-1 (Boltzmann's Constant)
u = 1.660538921e-27 # kg (Atomic mass unit/Dalton)
h = 6.62606957e-34  # J s (Plank constant)
hbar = h / (2*pi)

### Rubidium ###
class Rb():
    mass = 85.4678 * u   # kg (mass of rubidium)


class AtomicLine():
    def __init__(self, freq, linewidth):
        self.freq = freq
        self.linewidth = linewidth

    def wavelength(self):
        return c/self.freq

    def Isat(self):
        return pi*h*c*self.linewidth/(3*self.wavelength()**3)

# Rb85
class Rb85(Rb):
    mass = 84.911789732 * u # kg (from Steck)
    abundance = 0.7217

    D1 = AtomicLine(377.107385690e12,  # Hz - transition frequency (from Steck)
                    5.7500e6)           # Hz - linewidth
    D2 = AtomicLine(384.230406373e12,  # Hz - transition frequency (from Steck)
                    6.0666e6)           # Hz - linewidth


# Rb87
class Rb87(Rb):
    m = 87 * u
    abundance = 0.2783

    D1 = AtomicLine(377.107463380e12,  # Hz - transition frequency (from Steck)
                    5.7500e6)           # Hz - linewidth
    D2 = AtomicLine(384.230484468e12,  # Hz - transition frequency (from Steck)
                    6.0666e6)           # Hz - linewidth


