"""
Spectra for HIP_84543
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

h = 6.62607015e-27  # in erg cgs
c = 2.99792458e10  
k = 1.380649e-16

def planck_law(wavelength, temperature):
   
    exponent = (h * c) / (wavelength * 1e-7 * k * temperature)
    sr = 1.66 * 1e-27 * (2 * h * c**2) / ((wavelength*1e-7)**5 * (np.exp(exponent) - 1))
    return sr

files = ['C:/Users/mrobinson/OneDrive - Maynooth University/Project/HIP_84543/HIP_84543_UVB.fits','C:/Users/mrobinson/OneDrive - Maynooth University/Project/HIP_84543/HIP_84543_VIS.fits','C:/Users/mrobinson/OneDrive - Maynooth University/Project/HIP_84543/HIP_84543_NIR.fits']
waves, fluxes, e = [], [], []

for f in files:
    with fits.open(f) as s:
        data = s[1].data
        waves.append(data['WAVE'][0])
        fluxes.append(data['FLUX'][0])
        e.append(data['ERR'][0])

plt.figure(figsize=(12,5))
for i in range(3):
    plt.plot(waves[i],fluxes[i],linewidth = 0.1, alpha = 0.8)
    plt.ylim(0,0.5e-12)
    plt.grid()

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

model = np.loadtxt("C:/Users/mrobinson/OneDrive - Maynooth University/Project/HIP_84543/mp_t6500_g+4.5_m+0.00_a+0.00_c+0.00_v0_r10000_resam.txt")
lam = np.loadtxt("C:/Users/mrobinson/OneDrive - Maynooth University/Project/HIP_84543/bosz2024_wave_r10000.txt")

plt.figure(figsize=(12,5))
for i in range(3):
    plt.plot(waves[i],fluxes[i],linewidth = 0.1, alpha = 0.8)
    plt.ylim(0,0.5e-12)
    plt.grid()
plt.plot(wave,planck_law(wave, 5500))

plt.figure(figsize=(12,5))
plt.plot(wave,flux,linewidth = 0.1, c = 'black', alpha = 0.8)
plt.ylim(0,0.5e-12)
plt.axvline(364.7, linewidth = 0.5, c = 'r')
plt.grid()

plt.figure(figsize=(12,5))
plt.plot(wave,err, c = 'red')
plt.ylim(0,0.5e-12)
plt.grid()

w = model[:,0]
f = model[:,1]

plt.figure(figsize=(12,5))
plt.plot(lam,f+w)
#plt.ylim(0,0.5e-12)
#plt.axvline(364.7, linewidth = 0.5, c = 'r')
plt.xscale('log')
plt.grid()

plt.show()
