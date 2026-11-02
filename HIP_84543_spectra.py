"""
Spectra for HIP_84543
"""
#import pyphot
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
    sr = 8 * 1e-28 * (2 * h * c**2) / ((wavelength*1e-7)**5 * (np.exp(exponent) - 1))
    return sr

airmass = 1.5
sigma_dimm = 2.16

def dimm(wave, sigma_dimm, airmass, width):
    """Calculates T (slit loss factor"""
    sigma_lam = sigma_dimm * (wave/500)**(-0.2) * airmass**(0.6)
    return special.erf(np.sqrt(np.log(2)) * width / sigma_lam)


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
colours = ['cyan','orange','red']
plt.figure(figsize=(12,5))
for i in range(3):
    plt.plot(waves[i],fluxes[i],linewidth = 0.1, alpha = 0.8, color = colours[i])
    plt.ylim(0,0.5e-12)
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'Flux ($erg\ cm^{-2} s^{-1} \AA^{-1}$)')

arg1 = np.argmin(np.abs(waves[1] - waves[0][-1])) #index of VIS arm with same wavelength as last wavelength of UVB arm
waves[1] = waves[1][arg1:] #all VIS arrays from arg1 onward
fluxes[1] = fluxes[1][arg1:]
e[1] = e[1][arg1:]

arg2 = np.argmin(np.abs(waves[2] - waves[1][-1])) #index of NIR arm with same wavelength as last wavelength of VIS arm
waves[2] = waves[2][arg2:]  #all NIR arrays from arg1 onward
fluxes[2] = fluxes[2][arg2:]
e[2] = e[2][arg2:]

fluxes[0] = fluxes[0] / dimm(waves[0], sigma_dimm, airmass, 1) #applying DIMM seeing to the fluxes for different slit widths on each arm
fluxes[1] = fluxes[1] / dimm(waves[1], sigma_dimm, airmass, 0.9)
fluxes[2] = fluxes[2] / dimm(waves[2], sigma_dimm, airmass, 0.6)

diff1 = np.abs((fluxes[0][-1] - fluxes[1][0]) / fluxes[0][-1] * 100) #percentage difference between arms
diff2 = np.abs((fluxes[1][-1] - fluxes[2][0]) / fluxes[1][-1] * 100)
print(f'UVB-VIS: {diff1:.2f}%, VIS-NIR: {diff2:.2f}%')

print(waves[1][arg1])
print(waves[2][arg2])
wave = np.concatenate(waves) #concatenating each arm into one array
flux = np.concatenate(fluxes)
err = np.concatenate(e)

output_file = Path(__file__).parent / 'HIP_84543.txt'

with open(output_file , 'w') as f:
    writer = csv.writer(f, delimiter = '\t')
    writer.writerows(zip(wave,flux,err))

model = np.loadtxt("OneDrive - Maynooth University/Project/HIP_84543/mp_t6500_g+4.5_m+0.00_a+0.00_c+0.00_v0_r10000_resam.txt") * 7 * 1e-20 #bosz model
lam = np.loadtxt("OneDrive - Maynooth University/Project/HIP_84543/bosz2024_wave_r10000.txt") / 10 #converting to nm

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
    plt.plot(waves[i],fluxes[i],linewidth = 0.1, alpha = 0.8, color = colours[i])
    plt.ylim(0,1e-12)
plt.plot(wave,planck_law(wave, 5500), c ='magenta')
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'Flux ($erg\ cm^{-2} s^{-1} \AA^{-1}$)')

w = model[:,0]
f = model[:,1]

plt.figure(figsize=(12,5))
plt.plot(wave,flux,linewidth = 0.1, c = 'cyan', alpha = 0.8)
plt.plot(lam,f+w, c = 'violet', linewidth = 0.1, alpha = 0.8)
plt.ylim(0,1e-12)
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

plt.show()

"""lib = pyphot.get_library()

wave = np.asarray(wave, dtype=float) * 10
flux = np.asarray(flux, dtype=float)

f = lib['GROUND_JOHNSON_V']
sed = f.get_flux(wave,flux)
sed1 = f.Vega_zero_flux
mag = -2.5 * np.log10(f.get_flux(wave,flux) / f.Vega_zero_flux)
print(mag)"""
