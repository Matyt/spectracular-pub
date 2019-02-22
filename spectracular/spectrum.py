import numpy as np
import scipy

import astropy.constants as const
from astropy.table import QTable


class Spectrum(QTable):
    def __init__(self, table, description='DESCRIPTION', name='NAME', label='LABEL', colour=None):
        super().__init__(table)
        
        self.s_description = description
        self.s_name = name
        self.s_label = label
        self.s_colour = None
        self.w_mask = np.zeros(self['WAVE'].shape, dtype=bool)
        self.f_mask = np.zeros(self['FLUX'].shape, dtype=bool)
        self.e_mask = np.zeros(self['ERR'].shape, dtype=bool)
        self.g_mask = np.zeros(self['WAVE'].shape, dtype=bool)
    
    @staticmethod
    def from_arrays(w, f, e, **kwargs):
        a = np.argsort(w)
        return Spectrum(QTable([w[a], f[a], e[a]], names=('WAVE', 'FLUX', 'ERR')), **kwargs)
    
    @property
    def wave(self):
        return self['WAVE'][~self.mask]
    
    @property
    def flux(self):
        return self['FLUX'][~self.mask]
    
    @property
    def err(self):
        return self['ERR'][~self.mask]
    
    @property
    def color(self):
        return self.s_colour
    
    @property
    def mask(self):
        return self.w_mask | self.f_mask | self.e_mask | self.g_mask
    
    def mask_wave(self, mask_condition):
        self.w_mask = self.w_mask | mask_condition(self['WAVE'])
    
    def unmask_wave(self, unmask_condition):
        self.w_mask[unmask_condition(self['WAVE'])] = False
    
    def mask_flux(self, mask_condition):
        self.f_mask = self.f_mask | mask_condition(self['FLUX'])
    
    def unmask_flux(self, unmask_condition):
        self.f_mask[unmask_condition(self['FLUX'])] = False
    
    def mask_err(self, mask_condition):
        self.e_mask = self.e_mask | mask_condition(self['ERR'])
    
    def unmask_err(self, unmask_condition):
        self.e_mask[unmask_condition(self['ERR'])] = False
        
    def mask_generalised(self, mask_condition):
        self.g_mask = self.g_mask | mask_condition(self['WAVE'], self['FLUX'], self['ERR'])
    
    def unmask_generalised(self, unmask_condition):
        self.g_mask[unmask_condition(self['WAVE'], self['FLUX'], self['ERR'])] = False
    
    def doppler(self, v_r):
        doppler_factor = np.sqrt((1 + v_r / const.c) / (1 - v_r / const.c))
        self['WAVE'] *= doppler_factor
    
    def _interpolation(self, new_wave):
        funflux = scipy.interpolate.interp1d(
            self.wave, self.flux, kind='linear',
            axis=-1, copy=True, bounds_error=False, fill_value=None, assume_sorted=True)
        funerr = scipy.interpolate.interp1d(
            self.wave, self.err, kind='linear',
            axis=-1, copy=True, bounds_error=False, fill_value=None, assume_sorted=True)
        w, f, e = new_wave, funflux(new_wave) * self.flux.unit, funerr(new_wave) * self.err.unit
        
        return w, f, e
    
    def interpolation(self, new_wave):
        w, f, e = self._interpolation(new_wave)
        tab = Spectrum(QTable([w, f, e], names=('WAVE', 'FLUX', 'ERR')))
        return tab
    
    def _scale(self, factor):
        w, f, e = self.wave, self.flux * factor, self.err * abs(factor)
        return w, f, e
    
    def __mul__(self, factor):
        w, f, e = self._scale(factor)
        tab = Spectrum(QTable([w, f, e], names=('WAVE', 'FLUX', 'ERR')))
        return tab
    
    def __rmul__(self, factor):
        return self.__mul__(factor)
    
    def __add__(self, addend):
        if np.any(self.wave != addend.wave):
            raise ValueError('Wavelengths do not match (-> use interpolation)')
        
        def combine(f1, e1, f2, e2):
            f = f1 + f2
            e = np.sqrt(e1**2 + e2**2)
            return f, e
        
        w = self.wave
        f, e = combine(self.flux, self.err, addend.flux, addend.err)
        tab = Spectrum(QTable([w, f, e], names=('WAVE', 'FLUX', 'ERR')))
        return tab
    
    def __sub__(self, subtrahend):
        return self.__add__((-1) * subtrahend)
    
    def get_copy(self):
        return self * 1. # tricky
    
    def get_as_3arrays(self):
        w, f, e = self['WAVE'], self['FLUX'], self['ERR']
        return w, f, e
    
    def get_as_array(self):
        w, f, e = self.get_as_3arrays()
        arr = np.transpose(np.array([w, f, e]))
        return arr
    