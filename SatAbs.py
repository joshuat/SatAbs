from math import sqrt, pi, exp
from Constants import kB, h, hbar, c
from numpy import arange, array, append
from Rb import Rb85, Rb87

import matplotlib.pyplot as plt

def crossSection(freq_laser, freq_resonance, linewidth, I, Is):
    detuning = freq_laser - freq_resonance
    sigma0 = hbar*freq_resonance*linewidth/(2*Is) # Low intensity, on resonance
    return sigma0/(1+4*(detuning/linewidth)**2+I/Is)

def boltzmannDistribution(v, T, m, n0=1):
  return n0 * (m / sqrt(2*pi*kB*T)) * exp(-m*v*v / (2*kB*T))

def gaussianIntensity(I0, r, waist):
  return I0 * exp(-2 * r**2 / waist**2)

def lorentzian(freq, freq0, linewidth):
  return 1 / (1 + 4*((freq-freq0)/linewidth)**2)

def beer(I0, x, D):
  return I0*exp(-D(x)) 

# Atoms
T = 330 # Kelvin

# Laser
P = 20e-3 # W
beam_waist = 0.01 # m -realistic albeit largish
laser_linewidth = 300e3 # Hz - FWHM, from MOGlabs website
I0 = 2*P/(pi*beam_waist**2) # W/m^2

df = 5e6
fs = arange(-600e6, 600e6, df)
dv = 2
vs = arange(-600, 600, dv)
I_on_Isat = 10
alpha0 = 2e6
i=0
#transitions = Rb85.transitions.values()
transitions2 = [[Rb85.transitions.get('D2_F3_F3')], [Rb85.transitions.get('D2_F3_F4')],[Rb85.transitions.get('D2_F3_F2')], [Rb85.transitions.get('D2_F3_F4'), Rb85.transitions.get('D2_F3_F3')]]
#transitions2 = [[Rb85.transitions.get('D2_F3_F4'), Rb85.transitions.get('D2_F3_F3'), Rb85.transitions.get('D2_F3_F2')]]
for transitions in transitions2:
    gs2 = array([])
    transition_indicator = array([])
    for f in fs:
        freq_L = Rb85.transitions.get('D2_F3_F4').frequency() + f
        gs = array([])
        for v in vs:
            p0_minus_p1 = 0

            lorentz = 1
            for tran in transitions:
                freq_A = tran.frequency()*(1-v/c)
                freq_A_inv = tran.frequency()*(1+v/c) # From the perspective of the counterpropagating beam.

                p0_minus_p1 += 1 / (1+2*lorentzian(freq_L, freq_A, tran.linewidth)*I_on_Isat) + \
                              1 / (1+2*lorentzian(freq_L, freq_A_inv, tran.linewidth)*I_on_Isat)
                p0_minus_p1 *= tran.strength

                lorentz *= lorentzian(freq_L, freq_A, tran.linewidth*sqrt(1+2*I_on_Isat))*dv
                                                       # power broadened ^
            p0_minus_p1 *= boltzmannDistribution(v, T, Rb87.mass) # The v gets squared in here so it's the same for both beams

            dk = h*freq_L*alpha0*p0_minus_p1*lorentz

            gs = append(gs, dk)

        gs2 = append(gs2, gs.sum())

        val = 0
        for tran in transitions:
            if abs(tran.frequency() - freq_L) < 1*df:
                val =  gs2[gs2.size-1]
        transition_indicator = append(transition_indicator, val)


    plt.plot(fs, gs2/gs2.max(), label=str(i))
    i+=1
    plt.plot(fs, transition_indicator/transition_indicator.max())

plt.legend()
plt.xlabel('Detuning (Hz)')
plt.ylabel('Absorption')
plt.show()
