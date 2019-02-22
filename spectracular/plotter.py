import numpy as np

import astropy.units as u

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

from .spectrum import Spectrum
from .spectral_line import SpectralLine
from .tellurics import TelluricRange
from .wavelength_range import WavelengthRange


def plot_any(obj, *args, **kwargs):
    if issubclass(type(obj), Spectrum):
        plot_spectrum(obj, *args, **kwargs)
    elif issubclass(type(obj), SpectralLine):
        plot_spectralline(obj, *args, **kwargs)
    elif issubclass(type(obj), TelluricRange):
        plot_telluricrange(obj, *args, **kwargs)
    elif issubclass(type(obj), WavelengthRange):
        plot_telluricrange(obj, facecolor='#ffd1b2', *args, **kwargs)
    else:
        raise TypeError('Did not found plotting instructions for type {}'.format(type(obj)))

def plot_spectrum(spec, ax=None,
                  s='-', erralpha=0.5, label=None, colour=None,
                  wavemask=None, scale=1, add_fluxshift=0, add_waveshift=0,
                  *args, **kwargs):
    
    if ax is None:
        ax = plt.gca()
    
    w, f, e = spec.wave + add_waveshift, spec.flux * scale + add_fluxshift, spec.err * scale
    
    if wavemask is not None:
        mask = wavemask(w)
        [w, f, e] = np.ma.masked_where([mask, mask, mask], [w, f, e])
    
    if label is None:
        label = spec.s_label
    
    if colour is None:
        colour = spec.s_colour
    
    ax.plot(w, f, s, color=colour, label=label, *args, **kwargs)
    ax.fill_between(w, f-e, f+e, facecolor=colour, alpha=erralpha)

def plot_spectralline(line,
                      ax=None, xlim=None, ylim=None,
                      x_unit=None, annotate=True,
                      c='k', ls=':', lw=0.5, *args, **kwargs):
    if ax is None:
        ax = plt.gca()
    if xlim is None:
        xlim = ax.get_xlim()
    if ylim is None:
        ylim = ax.get_ylim()
    if x_unit is None:
        x_unit = u.AA

    xmin, xmax = xlim
    ymin, ymax = ylim
    w = line.wavelength
    w = w.to(x_unit).value
        
    ax.axvline(w, c=c, ls=ls, lw=lw, *args, **kwargs)
    if annotate:
        ax.annotate(line._label,
                    xy=(w, ymin + (ymax-ymin) * 0.006),
                    xytext=(w + (xmax-xmin) * 0.0006, ymin + (ymax-ymin) * 0.006))

def plot_telluricrange(trange,
                       ax=None, xlim=None, ylim=None,
                       x_unit=None, annotate=True,
                       facecolor='#c0c0c0', alpha=0.3, *args, **kwargs):
    if ax is None:
        ax = plt.gca()
    if xlim is None:
        xlim = ax.get_xlim()
    if ylim is None:
        ylim = ax.get_ylim()
    if x_unit is None:
        x_unit = u.AA

    xmin, xmax = xlim
    ymin, ymax = ylim
    
    w1, w2 = trange._wmin, trange._wmax
    w1, w2 = w1.to(x_unit).value, w2.to(x_unit).value
    
    if (xmin < w1) and (w2 < xmax):
        rect = Rectangle((w1, ymin), w2-w1, ymax-ymin)
        if annotate:
            ax.annotate(trange._label,
                        xy=(w2, ymin + (ymax-ymin) * 0.96),
                        xytext=(w2 + (xmax-xmin) * 0.0006, ymin + (ymax-ymin) * 0.96))
            boxes = [rect]
        else:
            boxes = []
            
        pc = PatchCollection(boxes,
                             facecolor=facecolor, alpha=alpha,
                             edgecolor=None)
        ax.add_collection(pc)

def express_plot(things, title=None, xlim=None):
    plt.figure()
    ax = plt.gca()
    for obj in things:
        plot_any(obj, ax=ax)
    
    if xlim is not None:
        ax.set_xlim(xlim[0].value, xlim[1].value)
    
    if title is not None:
        plt.title(title)
    
    plt.show()
