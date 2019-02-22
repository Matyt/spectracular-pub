import astropy.units as u
from astropy.table import QTable

from .spectrum import Spectrum


def try_load_fits(path):
    t = QTable.read(path, format='fits')
    q = QTable(data=[t[0][c] for c in t.colnames], names=t.colnames, meta=t.meta)
    q = q[['WAVE', 'FLUX', 'ERR']]
    q.replace_column('WAVE', q['WAVE'].to(u.AA))
    return Spectrum(q)
