import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as const
from astropy import units as u
from astropy.analytic_functions import blackbody_lambda, blackbody_nu

def compute_planck(freq, temp):
    """
    Input: Frequency [s^-1 numpy array], temp [K float]
    Output: Brightness [erg s^-1 cm^-2 Hz^-1 ster^-1 numpy array]
    All constants in cgs

    Expression for brightness taken from Rybicki &
    Lightman eq. 1.51
    """
    brightness = ((2*planck_c*freq**3.)/(light_speed**2))*\
                    1./(np.exp(planck_c*freq/(boltzmann_c*temp))-1)

    return brightness

def plot_curve(freq, flux):
    plt.plot(freq, flux, color='k', lw=2)
    #plt.gca().invert_xaxis()
    plt.xscale('log')
    #plt.yscale('log')
    plt.show()

def total_brightness_is(freq, flux):
    """
    Integrate Planck spectrum over frequency to get
    total brightness. Because of how the arrays are
    defined need to reverse them for the integration
    """
    return np.trapz(flux[::-1], freq[::-1])

def peak_ofcurve_is(freq, flux):


light_speed = 2.9979e10  # cm s^-1
planck_c = 6.626e-27     #
boltzmann_c = 1.38e-16   #
if __name__ == '__main__':
    temp = 10e5
    freq = np.linspace(10e20, 10e4, 10e6)
    flux = compute_planck(freq, temp)

    print total_brightness_is(freq, flux)
    plot_curve(freq, flux)
