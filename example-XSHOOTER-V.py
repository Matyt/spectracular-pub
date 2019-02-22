import numpy as np
import astropy.units as u

from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation

import spectracular as us

U_FLUX_DENSITY = u.Unit('erg cm-2 s-1 AA-1')


plots = True

####### LOAD YOUR DATA #######

prefixpath = './example_data/'
rawdata = {
    'A': np.loadtxt(prefixpath + 'A_VIS.txt'),
    'B': np.loadtxt(prefixpath + 'B_VIS.txt'),
    'C': np.loadtxt(prefixpath + 'C_VIS.txt'),
}

sc = SkyCoord(ra=214.013691*u.deg, dec=-56.91315*u.deg)

EarthLocation._get_site_registry(force_download=True)
    # update the astropy's sitename list
    # and force it, because Paranal is not on the included "off-line list"
vlt = EarthLocation.of_site('Paranal Observatory')

mjds = {
    'A': Time(58282.9828701409, format='mjd'),
    'B': Time(58330.0597107003, format='mjd'),
    'C': Time(58355.9962764888, format='mjd'),
}

descriptions = {
    'A': 'A: highly magnified',
    'B': 'B: slightly magnified',
    'C': 'C: baseline',
}

labels = {
    'A': 'A magn (20180613)',
    'B': 'B near baseline (20180731)',
    'C': 'C baseline (20180826)',
}

colours = {
    'A': 'orange',
    'B': 'violet',
    'C': 'blue'
}

obs = {key:
    us.Observation.from_arrays(
        skycoord=sc, time=mjds[key], earthlocation=vlt,
        w=rawdata[key][:,0] * u.AA,
        f=rawdata[key][:,1] * U_FLUX_DENSITY,
        e=rawdata[key][:,2] * U_FLUX_DENSITY,
        name=key, description=descriptions[key],
        label=labels[key], colour=colours[key]
    )
    for key in ('A', 'B', 'C')
}

magsV = {
    'A': 13.6,
    'B': 15.85,
    'C': 16.25
}

####### SOME CLEANING AND CUTS #######

for s in obs.values():
    s.mask_wave(lambda w: w <= 5800 * u.AA)
    s.mask_flux(lambda f: f <= 0)

that_visible_spikes_mask = lambda w, f, e: (
    (
          ((6400 * u.AA <= w) & (w <= 6450 * u.AA))
        & (1.1e-14 * U_FLUX_DENSITY <= f)
    ) |
    (
          ((8950 * u.AA <= w) & (w <= 9050 * u.AA))
        & (1.6e-14 * U_FLUX_DENSITY <= f)
    ) |
    (
          ((9800 * u.AA <= w) & (w <= 9900 * u.AA))
        & (1.7e-14 * U_FLUX_DENSITY <= f)
    )
)
obs['A'].mask_generalised(that_visible_spikes_mask)

if plots:
    us.plotter.express_plot(obs.values(), title="Cut'n'cleaned spectra")
    
####### PHOTOMETRIC CALIBRATION #######

cal = {key:
    us.level_spectrum_to_photometry(obs[key], us.passbands['V'], magsV[key])
    for key in obs.keys()
}

if plots:
    us.plotter.express_plot(cal.values(), title='Photometrically calibrated spectra')

####### MASK TELLURICS AND SODIUM LINE #######

cor = {key: obj.get_copy() for key, obj in cal.items()}
for s in cor.values():
    for tell in us.known_tellurics.values():
        s.mask_wave(tell.get_selector())
        
    sodium = us.ranges_of_interest['Na D1+D2']
    s.mask_wave(sodium.get_selector())

if plots:
    to_plot = list(cor.values()) + list(us.known_tellurics.values()) + [us.ranges_of_interest['Na D1+D2']]
    us.plotter.express_plot(to_plot, title='After applying range masks')

####### BARYCENTRIC CORRECTION #######

for k, s in cor.items():
    v_cor = s.barycentric_correction(return_vcor=True)
    print('Barycentric correction for "{}": applying Doppler shift by {}'.format(k, v_cor))

####### RESAMPLING WAVELENGTHS TO MATCH BASELINE #######

wave = cor['C'].wave
res = {key: obj.interpolation(wave) for key, obj in cor.items()}

####### SPECTRA ARITHMETIC #######

### SIMPLE DIFFERENCE:
spec_diff = res['A'] - res['C']

if plots:
    us.plotter.express_plot([spec_diff], title=r'Simple difference ($A-C$)')
    
    us.plotter.express_plot([spec_diff], xlim=us.ranges_of_interest['CaT'].wavelength_range,
                            title=r'Simple difference ($A-C$) – Calcium Triplet range')

### ALGEBRAIC SOLUTION:
muA, muC, muAerr, muCerr = 12.0, 1.1, 0.1, 0.1
    # let us use some guess values
source, lens = us.decompose(res['A'], res['C'], muA, muC, muAerr, muCerr)

source.s_colour, lens.s_colour = 'green', 'red'
if plots:
    us.plotter.express_plot([source, lens], title=r'Algebraic solution (from $A$ and $C$)')
    
    us.plotter.express_plot([source, lens], xlim=us.ranges_of_interest['CaT'].wavelength_range,
                            title=r'Algebraic solution (from $A$ and $C$) – Calcium Triplet range')

####### SAVING #######

np.savetxt('simple_difference.txt', spec_diff.get_as_array())
np.savetxt('algebraic_source.txt', source.get_as_array())
np.savetxt('algebraic_lens.txt', lens.get_as_array())
