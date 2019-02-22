import numpy as np
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.table import QTable

from .spectrum import Spectrum


class Observation(Spectrum):    
    _time = None
    _skycoord = None
    _earthlocation = None
    
    def __init__(self, skycoord, time, earthlocation, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
        assert type(skycoord) is SkyCoord
        assert type(time) is Time
        assert type(earthlocation) is EarthLocation
        
        self._time = time
        self._skycoord = skycoord
        self._earthlocation = earthlocation
    
    @staticmethod
    def from_arrays(skycoord, time, earthlocation, w, f, e, *args, **kwargs):
        a = np.argsort(w)
        return Observation(skycoord, time, earthlocation,
                           QTable([w[a], f[a], e[a]], names=('WAVE', 'FLUX', 'ERR')),
                           *args, **kwargs)
    
    def interpolation(self, new_wave):
        w, f, e = self._interpolation(new_wave)
        tab = Observation(
            self._skycoord, self._time, self._earthlocation,
            QTable([w, f, e], names=('WAVE', 'FLUX', 'ERR')),
            description='interpolated'
        )
        return tab
    
    def __mul__(self, factor):
        w, f, e = self._scale(factor)
        tab = Observation(
            self._skycoord, self._time, self._earthlocation,
            QTable([w, f, e], names=('WAVE', 'FLUX', 'ERR')),
            description='scaled'
        )
        return tab
    
    def barycentric_correction(self, return_vcor=False):
        v_cor = self._skycoord.radial_velocity_correction(
            obstime=self._time,
            location=self._earthlocation
        )
        self.doppler(v_cor)
        
        if return_vcor:
            return v_cor
