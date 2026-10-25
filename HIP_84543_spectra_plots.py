"""
Plots for HIP_84543
"""
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

def background():
    if input('Dark plot? y or n:') == 'y':
        plt.style.use('dark_background')
    else:
        plt.style.use('default')
    
background()  
plt.rcParams.update({
    'font.family' : 'serif',
    'font.serif' : 'cmr10',
    'mathtext.fontset' : 'cm',
    'axes.formatter.use_mathtext' : True,
    'axes.unicode_minus' : False,
    'font.size' : 12 })

input_file = Path(__file__).parent / 'HIP_84543.txt'
wave, flux, err = np.loadtxt(input_file, delimiter = '\t', unpack = True)

"""plt.figure(figsize=(12,5))
colours = ['cyan','orange','red']
for i in range(3):
    plt.plot(waves[i],fluxes[i],linewidth = 0.1, alpha = 0.8, color = colours[i])
    plt.ylim(0,0.2e-12)
#plt.plot(wave,planck_law(wave, 5500), c ='magenta')
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'Flux ($erg\ cm^{-2} s^{-1} \AA^{-1}$)')"""


plt.figure(figsize=(12,5))
plt.plot(wave,flux,linewidth = 0.1, c = 'cyan', alpha = 0.8)
#plt.plot(lam,f+w, c = 'violet', linewidth = 0.1, alpha = 0.8)
plt.ylim(0,0.5e-12)
plt.axvline(364.7, linewidth = 0.5, c = 'r')
plt.xlim(0,2500)
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'Flux ($erg\ cm^{-2} s^{-1} \AA^{-1}$)')

plt.figure(figsize=(12,5))
plt.plot(wave,err, c = 'red')
plt.ylim(0,0.5e-12)

plt.figure(figsize=(12,5))
#plt.plot(lam,f+w, c = 'cyan', linewidth = 0.1, alpha = 0.8)
plt.axvline(364.7, linewidth = 0.5, c = 'r')
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'Flux ($erg\ cm^{-2} s^{-1} \AA^{-1}$)')
plt.xlim(0,2500)

plt.figure(figsize=(12,5))
plt.plot(wave,flux,linewidth = 1, c = 'cyan', alpha = 0.8)
plt.ylim(0,.2e-12)
plt.axvline(364.7, linewidth = 0.5, c = 'r')
#plt.yscale('log')
#plt.xscale('log')
plt.xlim(910,1130)
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'Flux ($erg\ cm^{-2} s^{-1} \AA^{-1}$)')
plt.legend()

plt.show()

