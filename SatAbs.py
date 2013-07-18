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
T = 300 # Kelvin

# Laser
P = 20e-3 # W
beam_waist = 0.01 # m -realistic albeit largish
laser_linewidth = 300e3 # Hz - FWHM, from MOGlabs website
I0 = 2*P/(pi*beam_waist**2) # W/m^2

fs = arange(-800e6, 800e6, 5e6)
dv = 2
vs = arange(-600, 600, dv)
I_on_Isat = 10
alpha0 = 2e6
transitions = [Rb87.D2_F3, Rb87.D2_F2, Rb87.D2_F1, Rb87.D2_F0]


gs2 = array([])
for f in fs:
    freq_L = Rb87.D2.frequency() + f
    gs = array([])
    for v in vs:
        p0_minus_p1 = 0

        for tran in transitions:
            freq_A = tran.frequency()*(1-v/c)
            freq_A_inv = tran.frequency()*(1+v/c) # From the perspective of the counterpropagating beam.

            p0_minus_p1 += 1 / (1+2*lorentzian(freq_L, freq_A, tran.linewidth)*I_on_Isat) + \
                          1 / (1+2*lorentzian(freq_L, freq_A_inv, tran.linewidth)*I_on_Isat)

        p0_minus_p1 *= boltzmannDistribution(v, T, Rb87.mass) # The v gets squared in here so it's the same for both beams

        dk = 1
        for tran in transitions:
            freq_A = tran.frequency()*(1-v/c)

            dk *= h*freq_L*alpha0*p0_minus_p1*lorentzian(freq_L, freq_A, tran.linewidth*sqrt(1+2*I_on_Isat))*dv
                                                                           # power broadened ^

        gs = append(gs, dk)

    gs2 = append(gs2, gs.sum())

plt.plot(fs, gs2)

plt.xlabel('Detuning (Hz)')
plt.ylabel('Absorption')
plt.show()

#for f in [(Rb87.D2_F2.frequency() - Rb87.D2_F3.frequency())/2.]:
#    freq_L = Rb87.D2_F3.frequency() + f
#    gs = array([])
#    for v in vs:
#        p0_minus_p1 = 0

#        for tran in transitions:
#            freq_A = tran.frequency()*(1-v/c)
#            freq_A_inv = tran.frequency()*(1+v/c) # From the perspective of the counterpropagating beam.

#            p0_minus_p1 += 1 / (1+2*lorentzian(freq_L, freq_A, tran.linewidth)*I_on_Isat) + \
#                           1 / (1+2*lorentzian(freq_L, freq_A_inv, tran.linewidth)*I_on_Isat)

#        p0_minus_p1 *= boltzmannDistribution(v, T, Rb87.mass) # The v gets squared in here so it's the same for both beams
#        gs = append(gs, p0_minus_p1)

#    plt.plot(vs, gs)
#plt.xlabel('Velocity (m/s)')
#plt.ylabel('P0 - P1')
#plt.show()
