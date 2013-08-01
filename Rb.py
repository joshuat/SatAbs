from math import pi
from Constants import u, c, h, eV_to_J

### Rubidium ###
class Rb():
    mass = 85.4678 * u   # kg (mass of rubidium)


class Transition():
    def __init__(self, E, linewidth, strength = 0):
        self.energy = E

        self.linewidth = linewidth

        self.strength = strength

    def wavelength(self):
        return c/self.frequency()

    def frequency(self):
        return self.energy/h

    def Isat(self):
        return pi*h*c*self.linewidth/(3*self.wavelength()**3)

    def clone(self, freq_offset, strength):
        return Transition( (self.frequency() + freq_offset)*h, self.linewidth, strength = strength )
        

# Rb85
class Rb85(Rb):
    mass = 84.911789732 * u # kg (from Steck)
    abundance = 0.7217

    D1 = Transition(1.559590695*eV_to_J,  # J - Energy
                    5.7500e6)           # Hz - linewidth
    D2 = Transition(1.589049139*eV_to_J,  # J - Energy
                    6.0666e6)           # Hz - linewidth

    transitions = dict(
# From F = 2 to F = something
    D2_F2_F1 = D2.clone(1.7708439228e9-113.208e6, 3/10.0), # Hz - freq offset
    D2_F2_F2 = D2.clone(1.7708439228e9-83.835e6, 7/18.0), # Hz - freq offset
    D2_F2_F3 = D2.clone(1.7708439228e9-20.435e6, 14/45.0),  # Hz - freq offset
#    D2_F2_F4 = D2.clone(1.7708439228e9+100.205e6, 0), This transition is not allowed.

    # From F = 3 to F = something
#    D2_F2_F1 = D2.clone(-1.264888516e9-113.208e6, 0), This transition is not allowed.
    D2_F3_F2 = D2.clone(-1.264888516e9-83.835e6, 5/63.0),   # Hz - freq offset
    D2_F3_F3 = D2.clone(-1.264888516e9-20.435e6, 5/18.0),    # Hz - freq offset
    D2_F3_F4 = D2.clone(-1.264888516e9+100.205e6, 9/14.0))   # Hz - freq offset

# Rb87
class Rb87(Rb):
    m = 87 * u
    abundance = 0.2783

    D1 = Transition(1.559591016*eV_to_J, # J - Energy
                    5.7500e6)            # Hz - linewidth
    D2 = Transition(1.589049462*eV_to_J, # J - Energy
                     6.0666e6)           # Hz - linewidth

    transitions = dict(
    # From F = 1 to F = something
    D2_F1_F0 = D2.clone(4.271676631815181e9-302.0738e6, 1/6.0), # Hz - freq offset
    D2_F1_F1 = D2.clone(4.271676631815181e9-229.8518e6, 5/12.0), # Hz - freq offset
    D2_F1_F2 = D2.clone(4.271676631815181e9-72.9112e6, 5/12.0),  # Hz - freq offset
#    D2_F1_F3 = D2.clone(4.271676631815181e9+193.7407e6, 0), This transition is not allowed.

    # From F = 2 to F = something
#    D2_F2_F0 = D2.clone(-2.563005979089109e9-302.0738e6, 0), This transition is not allowed.
    D2_F2_F1 = D2.clone(-2.563005979089109e9-229.8518e6, 1/20.0),   # Hz - freq offset
    D2_F2_F2 = D2.clone(-2.563005979089109e9-72.9112e6, 1/4.0),    # Hz - freq offset
    D2_F2_F3 = D2.clone(-2.563005979089109e9+193.7407e6, 7/10.0))   # Hz - freq offset


