import numpy as np
import scipy
import astropy
import specutils
import pickle
import gc

from astropy.wcs import WCS
from astropy.io import fits
import astropy.units as u
import specutils.extinction
import astropy.constants as const
from astropy.table import QTable

from specutils.extinction import reddening
import scipy.optimize

from .config import coelho_path, processed_models_path

extinction_models = {'ccm89', 'od94', 'gcc09', 'f99', 'fm07'}

U_FLUX = u.Unit('erg AA-1 cm-2 s-1')


def fit_1Coelho(spectrum, CMM, teff, logg, mfeh, mafe, Gsigma=0., v_r=0.*u.km/u.s, binning=10, full_output=False):
    cm = CMM.get_Coelho_model(teff, logg, mfeh, mafe, Gsigma)
    if cm is None:
        print("Did not found Coelho model {:.0f} {:.1f} {:.1f} {:.1f} {:.3f}".format(teff, logg, mfeh, mafe, Gsigma))
        return None
    w, f, e = spectrum.wave, spectrum.flux, spectrum.err
    
    mmask = (w < CMM.wave[0]) | (CMM.wave[-1] < w) | np.isnan(f)
    w, f, e = w[~mmask], f[~mmask], e[~mmask]
    
    res = lambda Av, dist2red: (
        f - dist2red * cm.get_Spectrum(w, v_r=v_r, Av=Av, binning=10, interpolation_mode='extrapolate').flux
    )

    Av0 = 1
    ix = round(len(f) / 2)
    #print(w[ix])
    #print(f[ix], cm.get_Spectrum([w[ix].value] * u.AA, v_r=v_r, Av=Av0, binning=1, interpolation_mode='extrapolate').flux[0])
    dist2red0 = f[ix] / cm.get_Spectrum([w[ix].value] * u.AA, v_r=v_r, Av=Av0, binning=1, interpolation_mode='extrapolate').flux[0]
    sol = scipy.optimize.leastsq(
        func=lambda x: np.divide(res(x[0], x[1]), e),
        x0=[Av0, dist2red0.value]
    )
    x = sol[0]
    Av, dist2red = x[0], x[1] * dist2red0.unit
    
    chi2red = np.sum(np.divide(res(Av, dist2red), e)**2) / (len(w))
    resf = res(Av, dist2red)
    
    del cm
    del res
    
    if full_output:
        return Av, dist2red, chi2red, len(w), w, resf
    else:
        return Av, dist2red, chi2red, len(w)

def fit_Coelho_rows(spectrum, CMM, rows, Gsigma=0., v_r=0.*u.km/u.s, binning=10, save_interval=10, savefilename='__last_fit_Coelho_rows'):
    result_rows = []
    rowslen = len(rows)
    print('Iterating over {} models:'.format(rowslen))
    for ix, row in enumerate(rows):
        gc.collect()
        
        print(ix+1, end=' ')
        teff, logg, mfeh, mafe = row['Teff'], row['logg'], row['FeH'], row['alphaFe']
        Av, dist2red, chi2red, lenw = fit_1Coelho(spectrum, CMM, teff, logg, mfeh, mafe, Gsigma=0., v_r=v_r, binning=binning)
        result_rows.append((teff, logg, mfeh, mafe, chi2red, Av, dist2red, lenw))
        
        if ix % save_interval == 0:
            tab = QTable(rows=result_rows,
                   names=('Teff', 'logg', 'FeH', 'alphaFe', 'chi2', 'Av', 'dist2red', 'lenw')
                  ).write(savefilename + '.partial' + '.ecsv', format='ascii.ecsv', overwrite=True)
            del tab
    result_tab = QTable(rows=result_rows, names=('Teff', 'logg', 'FeH', 'alphaFe', 'chi2', 'Av', 'dist2red', 'lenw'))
    result_tab.write(savefilename + '.ecsv', format='ascii.ecsv', overwrite=True)
    
    print('\n' 'ended fitting!')
    return result_tab
