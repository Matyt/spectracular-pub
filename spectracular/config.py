import astropy.units as u

from .passband import Passband
from .spectral_line import SpectralLine
from .wavelength_range import WavelengthRange
from .tellurics import TelluricRange

_U_FLUX = u.Unit('erg cm-2 s-1')


####### MODELS #######

    # point to a directory containing the s_coelho14_highres/ catalog
coelho_path = '/home/stud/mzielinski/scr/ws/coelho_models/'
    # point to a directory where spectracular will save resampled models
processed_models_path = coelho_path

####### PASSBANDS #######

__filterpath1 = '/home/stud/mzielinski/scr/ws/spectracular_dir/passband_data/'
_passbands = [
    Passband.from_SVO_filter_profile(
        transcurve_file=__filterpath1+'GAIA_GAIA0.G.dat',
        zeropoint_AB=3.173e-9*_U_FLUX,
        name='G',
        description='Gaia G',
        label='G'),
    Passband.from_SVO_filter_profile(
        transcurve_file=__filterpath1+'Misc_APASS.sdss_r.dat',
        zeropoint_AB=2.904e-9*_U_FLUX,
        name='r',
        description='SDSS r',
        label='r'),
    Passband.from_SVO_filter_profile(
        transcurve_file=__filterpath1+'Misc_APASS.sdss_i.dat',
        zeropoint_AB=1.967e-9*_U_FLUX,
        name='i',
        description='SDSS i',
        label='i'),
    Passband.from_SVO_filter_profile(
        transcurve_file=__filterpath1+'GCPD_Johnson.V.dat',
        zeropoint_AB=3.662e-9*_U_FLUX,
        name='V',
        description='Johnson V',
        label='G'),
]

####### SPECTRAL LINES #######
_spectral_lines = [
    SpectralLine(6562.8518 * u.AA,
                 name='H-alpha',
                 description='H-alpha line (6562.8518 AA)',
                 label=r'H$\alpha$'),
    SpectralLine(4861.2786 * u.AA,
                 name='H-beta',
                 description='H-beta line (4861.2786 AA)',
                 label=r'H$\beta$'),
    SpectralLine(4340.462 * u.AA,
                 name='H-gamma',
                 description='H-gamma line (4340.462 AA)',
                 label=r'H$\gamma$'),
    SpectralLine(8498.018 * u.AA,
                 name='CaT 1',
                 description='1st line of the infrared Ca II triplet (8498.018 AA)',
                 label=r'Ca II'),
    SpectralLine(8542.089 * u.AA,
                 name='CaT 2',
                 description='2nd line of the infrared Ca II triplet (8542.089 AA)',
                 label=r'Ca II'),
    SpectralLine(8662.140 * u.AA,
                 name='CaT 3',
                 description='3rd line of the infrared Ca II triplet (8662.140 AA)',
                 label=r'Ca II'), 
]

####### VARIOUS WAVELENGTH RANGES #######

_known_tellurics = [
    TelluricRange(6863 * u.AA, 6919 * u.AA,
                  name='O2 ~6900',
                  description='O2 lines (6863-6919)',
                  label=r'O$_2$'),
    TelluricRange(7160 * u.AA, 7300 * u.AA,
                  name='H2O ~7200',
                  description='H2O lines (7160-7300)',
                  label=r'H$_2$O'),
    TelluricRange(7588 * u.AA, 7700 * u.AA,
                  name='O2 ~7600',
                  description='H2O lines (7588-7700)',
                  label=r'H$_2$O'),
    TelluricRange(8120 * u.AA, 8350 * u.AA,
                  name='H2O ~8200',
                  description='H2O lines (8120-8350)',
                  label=r'H$_2$O'),
    TelluricRange(8900 * u.AA, 9900 * u.AA,
                  name='H2O ~9000',
                  description='H2O range (8900-9900)',
                  label=r'H$_2$O'),
]

_ranges_of_interest = [
    WavelengthRange(5882 * u.AA, 5905 * u.AA,
                    name='Na D1+D2',
                    description='Na Fraunhofer D1+D2 lines (5882-5905)',
                    label=r'Na'),
    WavelengthRange(8450 * u.AA, 8700 * u.AA,
                    name='CaT',
                    description='infrared Calcium II triplet lines (8498, 8542, 8662)',
                    label=r'CaT'),
]

####### PROCESS CONFIGS #######

# Conversions list -> dict
try:
    step = 'passbands'
    passbands = {obj._name: obj for obj in _passbands}
    
    step = 'spectral_lines'
    spectral_lines = {obj.name: obj for obj in _spectral_lines}
    
    step = 'known_tellurics'
    known_tellurics = {obj.name: obj for obj in _known_tellurics}
    
    step = 'ranges_of_interest'
    ranges_of_interest = {obj.name: obj for obj in _ranges_of_interest}
except Exception as e:
    print("Caught an exception during initial configuration of Spectracular")
    print("You should possibly check the {} config variable".format(step))
    raise    
