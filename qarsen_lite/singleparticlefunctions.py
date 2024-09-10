import numpy as np
from numpy.polynomial.hermite import hermval
import math
import matplotlib.pyplot as plt
from matplotlib import cm


########### 1D ###########
def gwp_1D(x, xc:float, pc:float, alpha:float, gamma:float):
    """Gaussian wavepacket in 1D."""

    return np.exp(np.imag(gamma)) * (2*np.real(alpha)/np.pi)**0.25 * np.exp(-alpha*(x-xc)**2 + 1j*pc*(x-xc) + 1j*gamma)



def harmonic_osc_1D(x, n, m, k, xc):
    """Generates the eigenfunction of the harmonic oscillator system.
    Arguments
    x: is a space coordinate.
    n: is the vibrational quantum number, for this case only v=0 is considered.
    m: is the (reduced) mass of the system.
    k: is the force constant of the harmonic potential.
    xc: is the equilibrium separation.
    """

    hermite_sum = np.zeros(n+1)
    hermite_sum[-1] = 1

    return 1/(2**n * math.factorial(n))**0.5 * (((m*k)**0.5)/np.pi)**0.25 * np.exp(-(x-xc)**2 * ((m*k)**0.5)/2) * hermval((m*k)**0.25 * (x-xc),hermite_sum)