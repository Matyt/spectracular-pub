import numpy as np

from astropy.table import QTable
from .spectrum import Spectrum


def decompose(specA, specB, muA, muB, muAerr, muBerr):

    if np.any(specA.wave != specB.wave):
        raise ValueError('Wavelengths do not match (-> use interpolation)')
    
    fA, fB = specA.flux, specB.flux
    eA, eB = specA.err, specB.err
    
    mudiff, mudifferr2 = muA - muB, muAerr**2 + muBerr**2
    
    fS = (fA - fB) / mudiff
    eS2 = (
        np.multiply((fA / mudiff)**2,
                    np.divide(eA, fA)**2 + (muAerr**2 + muBerr**2) / mudiff**2)
        +
        np.multiply((fB / mudiff)**2,
                    np.divide(eB, fB)**2 + (muAerr**2 + muBerr**2) / mudiff**2)
    )
    
    fL = (fA * (-muB) + fB * muA) / mudiff
    eL2 = (
        np.multiply((fA * (-muB / mudiff))**2,
                    np.divide(eA, fA)**2 + np.divide(muBerr, muB)**2 + (muAerr**2 + muBerr**2) / mudiff**2)
        +
        np.multiply((fB * (muA / mudiff))**2,
                    np.divide(eB, fB)**2 + np.divide(muAerr, muA)**2 + (muAerr**2 + muBerr**2) / mudiff**2)
    )

    w = specA.wave
    
    sour = Spectrum(QTable([w, fS, np.sqrt(eS2)], names=('WAVE', 'FLUX', 'ERR')),
                    description='algebraic source')
    lens = Spectrum(QTable([w, fL, np.sqrt(eL2)], names=('WAVE', 'FLUX', 'ERR')),
                    description='algebraic lens')
    return sour, lens
