from .observation import Observation
from .spectrum import Spectrum
from .spectrum_loader import try_load_fits

from .spectral_line import SpectralLine
from .wavelength_range import WavelengthRange
from .tellurics import TelluricRange

from .photometry import level_spectrum_to_photometry

from .algebra import decompose
from . import plotter

from .coelho_models import CoelhoModelManager, resample_sigma_Coelho

from . import modelling

print('Spectracular core initialized')
print('Now processing config...')

from .config import passbands, spectral_lines, known_tellurics, ranges_of_interest

print('Spectracular config processed')

print('Spectracular ready to go!')
