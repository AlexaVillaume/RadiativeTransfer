import sys
import numpy as np
import matplotlib.pyplot as plt

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

def compute_planck_tderive(freq, temp):
    """
    Compute the curve of the derivative of the Planck curve
    with respect to temperature.

    Expression for this taken from Rybicki & Lightman eq. 1.55
    """
    exp1 = (2*(planck_c**2)*(freq**4))/((light_speed**2)*(boltzmann_c)*(temp**2))
    exp2 = np.exp(planck_c*freq/(boltzmann_c*temp))

    return exp1*(exp2/(exp2-1)**2)

def plot_curve(freq, flux, peak, temp):
    plt.plot(freq, flux, color='k', lw=2, label='Temp = {}'.format(temp))
    plt.axvline(peak['peak'], ls='--', label='x = {}'.format(peak['x']))
    #plt.gca().invert_xaxis()
    plt.xscale('log')
    #plt.yscale('log')
    plt.legend(frameon=False, loc='upper left')
    plt.show()

def total_brightness_is(freq, flux):
    """
    Integrate Planck spectrum over frequency to get
    total brightness. Because of how the arrays are
    defined need to reverse them for the integration
    """
    return np.trapz(flux[::-1], freq[::-1])

def peak_ofcurve_is(freq, flux, temp):
    peak = np.where(flux == max(flux))
    peak_freq = freq[peak]

    const = planck_c/(boltzmann_c * temp)
    return {'peak': peak_freq, 'x': peak_freq*const}

light_speed = 2.9979e10  # cm s^-1
planck_c = 6.626e-27     #
boltzmann_c = 1.38e-16   #
if __name__ == '__main__':
    temp = 10e5
    freq = np.linspace(10e20, 10e4, 10e6)
    flux = compute_planck(freq, temp)

    print total_brightness_is(freq, flux)
    peak = peak_ofcurve_is(freq, flux, temp)
    plot_curve(freq, flux, peak, temp)

    dBdT = compute_planck_tderive(freq, temp)
    i = np.where(~np.isnan(dBdT))
    peak = peak_ofcurve_is(freq[i], dBdT[i], temp)
    plot_curve(freq, dBdT, peak, temp)
