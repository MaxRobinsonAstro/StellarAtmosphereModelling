"""
Spectra for HIP_84543
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy import special
import csv
from pathlib import Path

h = 6.62607015e-27  # in erg cgs
c = 2.99792458e10  
k = 1.380649e-16

def planck_law(wavelength, temperature):
    """Planck's Law gives the blackbody curve"""
    exponent = (h * c) / (wavelength * 1e-7 * k * temperature)
    sr = 6 * 1e-27 * (2 * h * c**2) / ((wavelength*1e-7)**5 * (np.exp(exponent) - 1))
    return sr

airmass = 1.5
sigma_dimm = 2.16

def dimm(w, f, sigma_dimm, airmass, width):
    """Calculates T (slit loss factor"""
    sigma_lam = sigma_dimm * (w/500)**(-0.2) * airmass**(0.6)
    T = special.erf(np.sqrt(np.log(2)) * width / sigma_lam)
    return f / T


files = ['OneDrive - Maynooth University/Project/HIP_84543/HIP_84543_UVB.fits',
         'OneDrive - Maynooth University/Project/HIP_84543/HIP_84543_VIS.fits',
         'OneDrive - Maynooth University/Project/HIP_84543/HIP_84543_NIR.fits'] #fits files from 3 Xshooter arms

waves, fluxes, e = [], [], [] #2d arrays for wavelength, flux and error

for f in files:
    with fits.open(f) as s:
        data = s[1].data #data is contained in first extension
        waves.append(data['WAVE'][0])
        fluxes.append(data['FLUX'][0])
        e.append(data['ERR'][0])

"""Plotting arms before concatenating"""
colours = ['cyan','green','red']
plt.figure(figsize=(12,5))
for i in range(3):
    plt.plot(waves[i],fluxes[i],linewidth = 0.1, alpha = 0.8, color = colours[i])
    plt.ylim(0,2e-12)
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'Flux ($erg\ cm^{-2} s^{-1} \AA^{-1}$)')

plt.savefig("GitHub/StellarAtmosphereModelling/figures/HIP_84543_separate_before_slit")

arg1 = np.argmin(np.abs(waves[1] - waves[0][-1])) #index of VIS arm with same wavelength as last wavelength of UVB arm
waves[1] = waves[1][arg1:] #all VIS arrays from arg1 onward
fluxes[1] = fluxes[1][arg1:]
e[1] = e[1][arg1:]

arg2 = np.argmin(np.abs(waves[2] - waves[1][-1])) #index of NIR arm with same wavelength as last wavelength of VIS arm
waves[2] = waves[2][arg2:]  #all NIR arrays from arg1 onward
fluxes[2] = fluxes[2][arg2:]
e[2] = e[2][arg2:]

o = 300
waves[0] = waves[0][o:]
fluxes[0] = fluxes[0][o:]
e[0] = e[0][o:]

wave = np.concatenate(waves) #concatenating each arm into one array
flux = np.concatenate(fluxes)
err = np.concatenate(e)


fluxes_slitloss = fluxes.copy()
fluxes_slitloss[0] = dimm(waves[0], fluxes_slitloss[0], sigma_dimm, airmass, 1) #applying DIMM seeing to the fluxes for different slit widths on each arm
fluxes_slitloss[1] = dimm(waves[1], fluxes_slitloss[1], sigma_dimm, airmass, 0.9)
fluxes_slitloss[2] = dimm(waves[2], fluxes_slitloss[2], sigma_dimm, airmass, 0.6)

flux_slitloss = np.concatenate(fluxes_slitloss)

diff1 = np.abs((fluxes[0][-1] - fluxes[1][0]) / fluxes[0][-1] * 100) #percentage difference between arms
diff2 = np.abs((fluxes[1][-1] - fluxes[2][0]) / fluxes[1][-1] * 100)
print(f'Before Slitloss UVB-VIS: {diff1:.2f}%, VIS-NIR: {diff2:.2f}%')

diff1 = np.abs((fluxes_slitloss[0][-1] - fluxes_slitloss[1][0]) / fluxes_slitloss[0][-1] * 100) #percentage difference between arms
diff2 = np.abs((fluxes_slitloss[1][-1] - fluxes_slitloss[2][0]) / fluxes_slitloss[1][-1] * 100)
print(f'After Slitloss UVB-VIS: {diff1:.2f}%, VIS-NIR: {diff2:.2f}%')





output_file = Path(__file__).parent / 'HIP_84543.txt'

with open(output_file , 'w') as f:
    writer = csv.writer(f, delimiter = '\t')
    writer.writerows(zip(wave,flux,err))

model = np.loadtxt("OneDrive - Maynooth University/Project/HIP_84543/mp_t6500_g+4.5_m+0.00_a+0.00_c+0.00_v0_r10000_resam.txt") * 2 * 1e-19 #bosz model
lam = np.loadtxt("OneDrive - Maynooth University/Project/HIP_84543/bosz2024_wave_r10000.txt") / 10 #converting to nm

plt.style.use('dark_background')
plt.rcParams.update({
    'text.usetex' : True,
    'font.family' : 'serif',
    'font.serif' : 'cmr10',
    'axes.formatter.use_mathtext' : True,
    'axes.unicode_minus' : False,
    'font.size' : 12 })

plt.figure(figsize=(12,5))
for i in range(3):
    plt.plot(waves[i],fluxes[i],linewidth = 0.1, alpha = 0.8, color = colours[i])
    plt.ylim(0,2e-12)
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'Flux ($erg\ cm^{-2} s^{-1} \AA^{-1}$)')

plt.figure(figsize=(12,5))
for i in range(3):
    plt.plot(waves[i],fluxes_slitloss[i],linewidth = 0.1, alpha = 0.8, color = colours[i])
plt.ylim(0,2e-12)
plt.plot(wave,planck_law(wave, 5500), c ='magenta')
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'Flux ($erg\ cm^{-2} s^{-1} \AA^{-1}$)')
#plt.yscale('log')

plt.savefig("GitHub/StellarAtmosphereModelling/figures/HIP_84543_separate_after_slit")

w = model[:,0]
f = model[:,1]

plt.figure(figsize=(12,5))
plt.plot(wave,flux,linewidth = 0.1, c = 'cyan', alpha = 0.8)
plt.plot(lam,f+w, c = 'violet', linewidth = 0.1, alpha = 0.8)
plt.ylim(0,2e-12)
plt.axvline(364.7, linewidth = 0.5, c = 'r')
plt.xlim(0,2500)
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'Flux ($erg\ cm^{-2} s^{-1} \AA^{-1}$)')

plt.figure(figsize=(12,5))
plt.plot(wave,err, c = 'red')
plt.ylim(0,0.5e-12)

plt.figure(figsize=(12,5))
plt.plot(lam,f+w, c = 'cyan', linewidth = 0.1, alpha = 0.8)
plt.axvline(364.7, linewidth = 0.5, c = 'r')
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'Flux ($erg\ cm^{-2} s^{-1} \AA^{-1}$)')
plt.xlim(0,2500)

plt.figure(figsize=(12,5))
plt.plot(wave,flux,linewidth = 1, c = 'cyan', alpha = 0.8)
plt.ylim(0,1e-12)
plt.axvline(364.7, linewidth = 0.5, c = 'r')
#plt.yscale('log')
#plt.xscale('log')
plt.xlim(910,1130)
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'Flux ($erg\ cm^{-2} s^{-1} \AA^{-1}$)')
plt.legend()

plt.savefig("GitHub/StellarAtmosphereModelling/figures/HIP_84543_VIStoNIR")

plt.show()