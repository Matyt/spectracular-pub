import numpy as np
import astropy.units as u
from astropy.table import QTable
import matplotlib.pyplot as plt

import spectracular as us

U_FLUX_DENSITY = u.Unit('erg cm-2 s-1 AA-1')


####### LOAD A SPECTRUM TO BE MODELLED #######

rawdata = np.loadtxt('simple_difference.txt')

spec_diff = us.Spectrum.from_arrays(
    w=rawdata[:,0] * u.AA,
    f=rawdata[:,1] * U_FLUX_DENSITY,
    e=rawdata[:,2] * U_FLUX_DENSITY,
    name='spec_diff', description='simple difference (A – C)',
    label=r'$A-C$', colour='green'
    )

####### PREPARE COELHO MODELS MANAGER #######

cmm = us.CoelhoModelManager()

####### SHOW SOME INFO ABOUT THE MODELS #######

print('')
print('Effective temperatures in the set:')
print(np.unique(list(cmm.tab['Teff'])))
print('')
print('Metallicities in the set:')
print('\t[Fe/H]: ', np.unique(list(cmm.tab['FeH'])))
print('\t[alpha/Fe]: ', np.unique(list(cmm.tab['alphaFe'])))
print('')

####### CHOOSE A SUBSET OF MODELS #######

rows = cmm.select_rows(teffs=[4750, 5000, 5250], mfehs=[-0.8, -0.5, -0.3])

print('Chosen subset has {} rows'.format(len(rows)))

####### PERFORM SOME FITTING #######

v_r = -20.34 * u.km/u.s
fitres = us.modelling.fit_Coelho_rows(spec_diff, cmm, rows, v_r=v_r, binning=10, savefilename='fit_table')
    # results will be saved with 'ascii.ecsv' format

####### EGZAMINE FIT RESULTS #######

fittab = QTable.read('fit_table.ecsv')
fittab.sort('chi2')
print('Fit results table:')
print(fittab)

    # check the best (allegedly) parameters
teff, logg, mfeh, mafe, chi2, Av, scale, _ = list(fittab[0])

    # get the spectrum and residuals
cm = cmm.get_Coelho_model(teff, logg, mfeh, mafe, Gsigma=0.)
model = cm.get_Spectrum(cmm.wave, v_r=v_r, Av=Av, binning=10)

fit = us.modelling.fit_1Coelho(spec_diff, cmm, teff, logg, mfeh, mafe, v_r=v_r, full_output=True)
Av, dist2red, resw, resf = fit[0], fit[1], fit[4], fit[5]

####### A PLOT #######

spec_diff.s_colour = 'green'
model.s_colour = 'blue'

plt.figure()
ax = plt.gca()

us.plotter.plot_spectrum(model, ax, scale=dist2red, alpha=0.5)
us.plotter.plot_spectrum(spec_diff, ax, label=r'magnified source ($A-C$)', alpha=0.1)
ax.plot(resw, resf, label='residuals (obs – model)')
ax.axhline(0, ls=':', c='k')

plt.legend()
plt.xlabel(r'Wavelength [$\AA$]')
plt.ylabel(r'Flux density [erg / cm$^2$ / s / $\AA$]')
plt.show()

####### EGZAMINE SOME OTHER (PERHAPS BETTER) PARAMETER SET #######

teff, logg, mfeh, mafe, Av = 5000, 2.0, -0.5, 0.0, 2.5
cm = cmm.get_Coelho_model(teff, logg, mfeh, mafe, Gsigma=0.)
model = cm.get_Spectrum(cmm.wave, v_r=v_r, Av=Av, binning=10)

fit = us.modelling.fit_1Coelho(spec_diff, cmm, teff, logg, mfeh, mafe, v_r=v_r, full_output=True)
Av, dist2red, resw, resf = fit[0], fit[1], fit[4], fit[5]

spec_diff.s_colour = 'green'
model.s_colour = 'orange'

plt.figure()
ax = plt.gca()

us.plotter.plot_spectrum(model, ax, scale=dist2red, alpha=0.5)
us.plotter.plot_spectrum(spec_diff, ax, label=r'magnified source ($A-C$)', alpha=0.1)
ax.plot(resw, resf, label='residuals (obs – model)')
ax.axhline(0, ls=':', c='k')

x1, x2 = us.ranges_of_interest['CaT'].wavelength_range
plt.xlim(x1.value, x2.value)
plt.legend()
plt.xlabel(r'Wavelength [$\AA$]')
plt.ylabel(r'Flux density [erg / cm$^2$ / s / $\AA$]')
plt.show()
