
class SpectralLine(object):
    _name = None
    _description = None
    _label = None
    
    _wavelength = None
    
    def __init__(self, wavelength,
                 name=None, description=None, label=None):
        self._name = name
        self._description = description
        self._label = label
        
        self._wavelength = wavelength
        
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
    def wavelength(self):
        return self._wavelength
        