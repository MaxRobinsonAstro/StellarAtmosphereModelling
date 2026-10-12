"""
Spectra for HIP_36966
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

h = 6.62607015e-27
c = 2.99792458e10  
k = 1.380649e-16

def planck_law(wavelength, temperature):
   
    exponent = (h * c) / (wavelength * 1e-7 * k * temperature)
    sr = 3.3 * 1e-27 * (2 * h * c**2) / ((wavelength*1e-7)**5 * (np.exp(exponent) - 1))
    return sr

files = ['HIP_36966_UVB.fits','HIP_36966_VIS.fits','HIP_36966_NIR.fits']
waves, fluxes, e = [], [], []

for f in files:
    with fits.open(f) as s:
        data = s[1].data
        waves.append(data['WAVE'][0])
        fluxes.append(data['FLUX'][0])
        e.append(data['ERR'][0])

plt.figure(figsize=(12,5))
colours = ['cyan','orange','tomato']
for i in range(3):
    plt.plot(waves[i],fluxes[i],linewidth = 0.1, alpha = 0.8, color = colours[i])
    plt.ylim(0,4e-12)
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'Flux ($erg\ cm^{-2} s^{-1} \AA^{-1}$)')
             
arg1 = np.argmin(np.abs(waves[1] - waves[0][-1]))
waves[1] = waves[1][arg1:]
fluxes[1] = fluxes[1][arg1:]
e[1] = e[1][arg1:]

arg2 = np.argmin(np.abs(waves[2] - waves[1][-1]))
waves[2] = waves[2][arg2:]
fluxes[2] = fluxes[2][arg2:]
e[2] = e[2][arg2:]

print(waves[1][arg1])
print(waves[2][arg2])
wave = np.concatenate(waves)
flux = np.concatenate(fluxes)
err = np.concatenate(e)

plt.style.use('dark_background')
plt.rcParams.update({
    'font.family' : 'serif',
    'font.serif' : 'cmr10',
    'mathtext.fontset' : 'cm',
    'axes.formatter.use_mathtext' : True,
    'axes.unicode_minus' : False,
    'font.size' : 12 })

plt.figure(figsize=(12,5))
for i in range(3):
    plt.plot(waves[i],fluxes[i],linewidth = 0.3, alpha = 0.8, color = colours[i])
    plt.ylim(0,4e-12)
plt.plot(wave,planck_law(wave, 7500), c = 'violet')
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'Flux ($erg\ cm^{-2} s^{-1} \AA^{-1}$)')
plt.xlim(250,600)

plt.figure(figsize=(12,5))
plt.plot(wave,flux,linewidth = 0.1, c = 'cyan', alpha = 0.8)
plt.ylim(0,4e-12)
plt.axvline(364.7, linewidth = 0.5, c = 'tomato')
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'Flux ($erg\ cm^{-2} s^{-1} \AA^{-1}$)')

plt.figure(figsize=(12,5))
plt.plot(wave,err, c = 'tomato')
plt.ylim(0,4e-12)
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'Flux Error ($erg\ cm^{-2} s^{-1} \AA^{-1}$)')
plt.show()