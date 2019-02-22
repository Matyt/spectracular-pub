import numpy as np
import scipy
import astropy
import specutils
import pickle

from astropy.wcs import WCS
from astropy.io import fits
import astropy.units as u
import specutils.extinction
import astropy.constants as const
from astropy.table import QTable

from specutils.extinction import reddening
import scipy.optimize

from .config import coelho_path, processed_models_path
from .spectrum import Spectrum

extinction_models = {'ccm89', 'od94', 'gcc09', 'f99', 'fm07'}

U_FLUX = u.Unit('erg AA-1 cm-2 s-1')


class CoelhoModel():
    wave = np.linspace(2500, 2500 + 0.02 * 325000, 325001) * u.AA
    
    def __init__(self, flux, teff, logg, mfeh, mafe, Gsigma=0., Av=0.):
        self.flux = flux
        self.teff = teff
        self.logg = logg
        self.mfeh = mfeh
        self.mafe = mafe
        self.Gsigma = Gsigma
        self.Av = Av
        
        self.iflux = lambda interpolation_mode, binning: scipy.interpolate.interp1d(
            self.wave, np.convolve(self.flux, np.ones((binning,))/binning, mode='same') * U_FLUX, kind='linear',
            axis=-1, copy=True, bounds_error=False, fill_value=interpolation_mode, assume_sorted=True)
    
    @property
    def err(self):
        return np.zeros_like(self.flux) * U_FLUX
    
    def interpolation(self, new_wave):
        w, f, e = new_wave, self.iflux(new_wave) * U_FLUX, np.zeros_like(new_wave.value) * U_FLUX
        tab = Spectrum(QTable([w, f, e], names=('WAVE', 'FLUX', 'ERR')))
        return tab
    
    def get_Spectrum(self, wave, v_r, Av, ext_model='fm07', binning=10, interpolation_mode=None):
        deDoppler_factor = np.sqrt((1 - v_r / const.c) / (1 + v_r / const.c))
        labwave = wave * deDoppler_factor
        labflux = self.iflux(interpolation_mode, binning)(labwave) * U_FLUX
        obsflux = np.divide(labflux, reddening(wave, a_v=Av, model=ext_model))
        e = np.zeros_like(obsflux.value) * U_FLUX
        tab = Spectrum(QTable([wave, obsflux, e], names=('WAVE', 'FLUX', 'ERR')),
                       description='DESCRIPTION',
                       name='NAME',
                       label=r'Coelho($T$={Teff:d}, $\log\ {{g}}$={logg:.1f}, [Fe/H]={mfeh:.1f}, [$\alpha$/Fe]={mafe:.1f}); $A_V$={Av:.1f}, $v_r$={vr:.2f} km/s'.format(Teff=round(self.teff), logg=self.logg, mfeh=self.mfeh, mafe=self.mafe, Av=Av, vr=v_r.value))
        return tab

import fnmatch
import os

class CoelhoModelManager():
    wave = np.linspace(2500, 2500 + 0.02 * 325000, 325001) * u.AA
    Coelho_orig_subd = 's_coelho14_highres/'
    
    def __init__(self):
        self.orig_main_path = coelho_path
        self.main_path = processed_models_path
        
        self.try_pickling()
        self.tab = self.tabularize_pickles('highres_0.000/')
    
    @staticmethod
    def try_pickling():
        print('Will now try to process Coelho original files, if they are not processed')
        print('This step might take a long time (5–10 minutes, probably)...')
        dump_orig_Coelho(coelho_path, CoelhoModelManager.Coelho_orig_subd, processed_models_path)
        print('...done!')
    
    @staticmethod
    def tabularize_pickles(subd, overwrite=False):
        main_path = processed_models_path
        if os.path.isfile(main_path + 'CoelhoQTable.votable') and not overwrite:
            tab = QTable.read(main_path + 'CoelhoQTable.votable')
            return tab
        
        data_rows = []
        for file in os.listdir(main_path + subd):
            with open(main_path+subd+file, 'rb') as f:
                p = pickle.load(f)

                teff = p['teff']
                logg = p['logg']
                mfeh = p['mfeh']
                mafe = p['mafe']
                
                data_rows.append((file, teff, logg, mfeh, mafe))
        
        colnames = ('filename', 'Teff', 'logg', 'FeH', 'alphaFe')
        #colnames = ('filename', 'Teff', 'logg', '[Fe/H]', '[alpha/Fe]')
        #colnames=('filename', r'$T_\text{eff}$', r'$\log{g}$', r'[Fe/H]', r'[$\alpha$/Fe]')
        tab = QTable(rows=data_rows, names=colnames)
        tab.meta['searched_path'] = main_path + subd
        tab.write(main_path + 'CoelhoQTable.votable', format='votable', overwrite=True)
        
        return tab

    @staticmethod
    def read_pickled(filepath, Gsigma):
        with open(filepath, 'rb') as f:
            p = pickle.load(f)
            
            flux = p['flux']
            teff = p['teff']
            logg = p['logg']
            mfeh = p['mfeh']
            mafe = p['mafe']
            Gsigma = Gsigma
            
            del p
            
            Cmodel = CoelhoModel(flux, teff, logg, mfeh, mafe, Gsigma=Gsigma)
            
            return Cmodel
    
    def get_Coelho_filename(self, teff, logg, mfeh, mafe):
        #+ 'plc'?
        signc = lambda x: 'p' if x >= 0 else 'm'
        fname = 't{teff:05d}_g{logg:+01.1f}_{mfeh_sign}{mfeh10a:02.0f}{mafe_sign}{mafe10a:02.0f}_hr'.format(
            teff=teff, logg=logg, mfeh_sign=signc(mfeh), mfeh10a=10*abs(mfeh), mafe_sign=signc(mafe), mafe10a=10*abs(mafe)
        )
        #print(fname)
        opath = self.main_path + self.Coelho_orig_subd
        for ofile in os.listdir(opath):
            if fnmatch.fnmatch(ofile, '{}*'.format(fname)):
                return ofile[:-5]
    
    def get_Coelho_model(self, teff, logg, mfeh, mafe, Gsigma):
        #tried_name = self.get_Coelho_filename(teff, logg, mfeh, mafe) + '.pickle'
        def du(arg, dunit):
            return arg if not u.Quantity(arg).unit.is_unity() else arg * dunit
        
        matchingmask = (
            (self.tab['Teff'] == du(teff, u.K).value) &
            (self.tab['logg'] == du(logg, u.dex).value) &
            (self.tab['FeH'] == du(mfeh, u.dex).value) &
            (self.tab['alphaFe'] == du(mafe, u.dex).value)
        )
        rows = self.tab[matchingmask]
        tried_name = None if 0 == len(rows) else rows[0]['filename']
        
        if tried_name is not None:
            filepath = self.main_path + 'highres_{Gsigma:.3f}/'.format(Gsigma=Gsigma) + tried_name
            Cmodel = self.read_pickled(filepath, Gsigma=Gsigma)
            return Cmodel
    
    def select_rows(self, teffs=None, loggs=None, mfehs=None, mafes=None):
        default_units = ((teffs, u.K), (loggs, u.dex), (mfehs, u.dex), (mafes, u.dex))
        for arg, dunit in default_units:
            if arg is None:
                continue
            if u.Quantity(arg).unit.is_unity():
                pass
                #arg *= dunit
        matchingmask = np.ones((len(self.tab),), dtype=bool)
        query_parameters = {'Teff': teffs, 'logg': loggs, 'FeH': mfehs, 'alphaFe': mafes}
        for colname, testset in query_parameters.items():
            if testset is not None:
                matchingmask &= np.isin(self.tab[colname], testset)
        selection = self.tab[matchingmask]
        
        return selection
    
import pickle
import os

def read_orig_Coelho(fpath, dpath):
    with fits.open(fpath) as f:
        head = f[0].header

        #wave = np.linspace(2500, 2500 + 0.02 * 325000, 325001)
        flux = f[0].data
        teff = float(head['TEFF'])
        logg = float(head['LOG_G'])
        mfeh = float(head['FEH'])
        mafe = float(head['AFE'])
        
        Cmodel = {
            #'wave': wave,
            'flux': flux,
            'teff': teff,
            'logg': logg,
            'mfeh': mfeh,
            'mafe': mafe
        }
        
        with open(dpath, 'wb') as d:
            pickle.dump(Cmodel, d)

def dump_orig_Coelho(orig_main_path, Coelho_orig_subd, main_path, overwrite=False):
    opath = orig_main_path + Coelho_orig_subd
    dpath = main_path + 'highres_0.000/'
    
    mkdir_p(dpath)
    
    for ofile in os.listdir(opath):
        dfull = os.path.join(dpath, ofile[:-4] + 'pickle')
        if os.path.isfile(dfull) and not overwrite:
            continue
        ofull = os.path.join(opath, ofile)
        with open(ofull) as o:
            read_orig_Coelho(ofull, dfull)

gaussian = lambda mu, sigma, x: np.exp(- (x - mu)**2 / (2 * sigma**2)) / np.sqrt(2*np.pi * sigma**2)

n3s = lambda sigma: round(3 * sigma / 0.02)
pxgauss3 = lambda sigma: gaussian(0, sigma, np.linspace(-n3s(sigma)*0.02, n3s(sigma)*0.02, 2*n3s(sigma) + 1))
pxgauss3 = lambda sigma: pxgauss3(sigma) / np.sum(pxgauss3(sigma))

def resample_pxprofile(wave, flux, pxprofile, npx):
    nf = []
    lenwave = len(wave)
    lenpxpr = len(pxprofile)
    for ix, w in enumerate(wave):
        m0, m1 = max(0, ix-npx), min(lenwave, ix+npx+1)
        
        wint, fint = wave[m0:m1], flux[m0:m1]
        pint = pxprofile[m0-(ix-npx):(m1-(ix+npx+1) if m1-(ix+npx+1) < 0 else lenpxpr)]
        
        f = np.multiply(fint, pint)
        nf.append(np.sum(f) / np.sum(pint))
        
    return np.array(nf)

def resample_sigma_Coelho(sigma, overwrite=False, silent=False):
    opath = Coelho_main_path + Coelho_0000_subd
    npath = Coelho_main_path + 'highres_{:.3f}/'.format(sigma)
    
    mkdir_p(npath)
    
    lenx = len(os.listdir(opath))
    for ix, ofile in enumerate(os.listdir(opath)):
        nfull = os.path.join(npath, ofile)
        if os.path.isfile(nfull) and not overwrite:
            continue
        if not silent:
            print('{}/{}: {}'.format(ix, lenx, nfull), end='')
        
        ofull = os.path.join(opath, ofile)
        wave, flux, teff, logg, mfeh, mafe = read_pickled_Coelho(ofull)
        flux = resample_pxprofile(wave, flux, pxgauss3(sigma), n3s(sigma))
        with open(nfull, 'wb') as f:
            Cmodel = {
                #'wave': wave,
                'flux': flux,
                'teff': teff,
                'logg': logg,
                'mfeh': mfeh,
                'mafe': mafe
            }
            pickle.dump(Cmodel, f)
        
        if not silent:
            print(' dumped')

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        pass
        #if exc.errno == errno.EEXIST and os.path.isdir(path):
        #    pass
        #else:
        #    raise
