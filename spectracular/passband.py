import numpy as np
import scipy


class Passband(object):
    _name = None
    _description = None
    _label = None
    
    _transmission_curve = None
    _transmission_function = None
    
    _zeropoint_AB = None
    
    def __init__(self, transcurve, zeropoint_AB,
                 name=None, description=None, label=None):
        self._name = name
        self._description = description
        self._label = label
        
        self._transmission_curve = transcurve
        self._transmission_function = scipy.interpolate.interp1d(
            transcurve[:,0], transcurve[:,1],
            kind='cubic', axis=-1, copy=True,
            bounds_error=False,
            fill_value=(0, 0), # expected no observed transmission outside the domain
            assume_sorted=True)
        
        self._zeropoint_AB = zeropoint_AB
    
    @staticmethod
    def from_SVO_filter_profile(transcurve_file, zeropoint_AB, *args, **kwargs):
        transcurve = np.loadtxt(transcurve_file)
        return Passband(transcurve, zeropoint_AB, *args, **kwargs)
    