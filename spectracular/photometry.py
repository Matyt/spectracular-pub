import numpy as np
import scipy.integrate
import scipy.interpolate


def synthetic_photometry(spectrum, passband):
    spi = scipy.interpolate.interp1d(
        spectrum.wave, spectrum.flux, kind='cubic',
        axis=-1, copy=True, bounds_error=False, fill_value=(0, 0), assume_sorted=True)
    tci = passband._transmission_function
    tcd = passband._transmission_curve
    
    f = lambda w: tci(w) * spi(w)
    i = scipy.integrate.quad(f, tcd[0,0], tcd[-1,0])[0]
    #ti = scipy.integrate.quad(tci, tcd[0,0], tcd[tcd[:,0] <= spectrum.wave[-1].value,0][-1])[0]
    iti = i #iti = i / ti
    F0 = passband._zeropoint_AB.value
    mag = - 2.5 * np.log10(iti / F0)
    
    return (i, iti, mag)

def level_spectrum_to_photometry(spectrum, passband, mag):
    f = synthetic_photometry(spectrum, passband)
    ref_flux = passband._zeropoint_AB.value * np.power(10, mag * (-1/2.5))
    
    scale = ref_flux / f[1]
    return spectrum * scale
    