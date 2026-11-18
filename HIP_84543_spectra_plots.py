"""
Plots for HIP_84543
"""
import os
os.environ['PATH'] += os.pathsep + r'C:\Users\Owner\AppData\Local\Programs\MiKTeX\miktex\bin\x64' 
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append(r'GitHub/StellarAtmosphereModelling')
from HIP_84543_spectra import wave, flux, flux_slitloss, o

def background():
    if input('Dark plot? y or n:') == 'y':
        plt.style.use('dark_background')
    else:
        plt.style.use('default')
    
background()  
plt.rcParams.update({
    'text.usetex' : True,
    'font.family' : 'serif',
    'font.serif' : ['Computer Modern Roman'],
    'axes.formatter.use_mathtext' : True,
    'axes.unicode_minus' : False,
    'font.size' : 12 })


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
plt.plot(wave,flux_slitloss,linewidth = 0.1, c = 'cyan', alpha = 0.8)
plt.ylim(0,1.5e-12)
plt.axvline(364.7, linewidth = 0.5, c = 'r')
plt.xlim(wave[o],2500)
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'Flux ($erg\ cm^{-2} s^{-1} \AA^{-1}$)')



plt.figure(figsize=(12,5))
plt.plot(wave,flux,linewidth = 0.1, c = 'cyan', alpha = 0.8)
#plt.ylim(0,.2e-12)
plt.ylim(1e-14,1e-12)
plt.axvline(364.7, linewidth = 0.5, c = 'r')
plt.yscale('log')
#plt.xscale('log')
plt.xlim(o,800)
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'Flux ($erg\ cm^{-2} s^{-1} \AA^{-1}$)')
plt.legend()

plt.show()

