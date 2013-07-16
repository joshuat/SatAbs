from math import sqrt, pi, exp
from Constants import kB, h, hbar, c, Rb87, Rb85
from numpy import arange, array, append

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
T = 300 # Kelvin

# Laser
P = 20e-3 # W
beam_waist = 0.01 # m -realistic albeit largish
laser_linewidth = 300e3 # Hz - FWHM, from MOGlabs website
I0 = 2*P/(pi*beam_waist**2) # W/m^2

fs = arange(-800e6, 800e6, 1e6)
dv = 1
vs = arange(-600, 600, dv)
I_on_Isat = 10
gs2 = array([])
alpha0 = 2e6
for f in fs:
    freq_L1 = Rb87.D2.freq - f
    freq_L2 = Rb87.D2.freq + f
    gs = array([])
    for v in vs:
        freq_A = Rb87.D2.freq*(1-v/c)
        p0_minus_p1 = 1 / (1+2*lorentzian(freq_L1, freq_A, Rb87.D2.linewidth)*I_on_Isat) + \
                      1 / (1+2*lorentzian(freq_L2, freq_A, Rb87.D2.linewidth)*I_on_Isat)
        p0_minus_p1 *= boltzmannDistribution(v, T, Rb87.mass)
        dk = h*freq_L1*alpha0*p0_minus_p1*lorentzian(freq_L1, freq_A, Rb87.D2.linewidth)*dv
        gs = append(gs, dk)

    gs2 = append(gs2, gs.sum())

plt.plot(fs, gs2)
plt.xlabel('Detuning (Hz)')
plt.ylabel('Absorption')
plt.show()
