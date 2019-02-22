
class WavelengthRange(object):
    _name = None
    _description = None
    _label = None
    
    _wmin = None
    _wmax = None
    
    def __init__(self, wmin, wmax,
                 name=None, description=None, label=None):
        self._name = name
        self._description = description
        self._label = label
        
        self._wmin = wmin
        self._wmax = wmax
        
    @property
    def name(self):
        return self._name
    
    @property
    def description(self):
        return self._description
    
    @property
    def label(self):
        return self._label
    
    @property
    def wavelength_range(self):
        return self._wmin, self._wmax
    
    def get_selector(self):
        return lambda w: (self._wmin <= w) & (w <= self._wmax)
