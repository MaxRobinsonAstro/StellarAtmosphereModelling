"""
Spectra for HIP_84543
"""
import pyphot
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy import special

h = 6.62607015e-27  # in erg cgs
c = 2.99792458e10  
k = 1.380649e-16

def planck_law(wavelength, temperature):
   
    exponent = (h * c) / (wavelength * 1e-7 * k * temperature)
    sr = 8 * 1e-28 * (2 * h * c**2) / ((wavelength*1e-7)**5 * (np.exp(exponent) - 1))
    return sr

airmass = 1.194

sigma_dimm = 1.406
def dimm(wave, sigma_dimm, airmass, width):
    sigma_lam = sigma_dimm * (wave/500)**(-0.2) * airmass**(0.6)
    
    return special.erf(np.sqrt(np.log(2)) * width / sigma_lam)


files = ['OneDrive - Maynooth University/Project/HIP_84543/HIP_84543_UVB.fits',
         'OneDrive - Maynooth University/Project/HIP_84543/HIP_84543_VIS.fits',
         'OneDrive - Maynooth University/Project/HIP_84543/HIP_84543_NIR.fits']
waves, fluxes, e = [], [], []

for f in files:
    with fits.open(f) as s:
        data = s[1].data
        waves.append(data['WAVE'][0])
        fluxes.append(data['FLUX'][0])
        e.append(data['ERR'][0])

colours = ['cyan','orange','red']
plt.figure(figsize=(12,5))
for i in range(3):
    plt.plot(waves[i],fluxes[i],linewidth = 0.1, alpha = 0.8, color = colours[i])
    plt.ylim(0,0.5e-12)
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

fluxes[0] = fluxes[0] * (1 - dimm(waves[0], sigma_dimm, airmass, 1))
fluxes[1] = fluxes[1] * (1 - dimm(waves[1], sigma_dimm, airmass, 0.9))
fluxes[2] = fluxes[2] * (1 - dimm(waves[2], sigma_dimm, airmass, 0.6))

print(fluxes[1][-1],fluxes[2][0])
print((fluxes[1][-1] - fluxes[2][0]) / fluxes[1][-1] * 100)

print(waves[1][arg1])
print(waves[2][arg2])
wave = np.concatenate(waves)
flux = np.concatenate(fluxes)
err = np.concatenate(e)

model = np.loadtxt("OneDrive - Maynooth University/Project/HIP_84543/mp_t6500_g+4.5_m+0.00_a+0.00_c+0.00_v0_r10000_resam.txt") * 7 * 1e-20
lam = np.loadtxt("OneDrive - Maynooth University/Project/HIP_84543/bosz2024_wave_r10000.txt") / 10

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
    plt.ylim(0,0.2e-12)
plt.plot(wave,planck_law(wave, 5500), c ='magenta')
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'Flux ($erg\ cm^{-2} s^{-1} \AA^{-1}$)')

w = model[:,0]
f = model[:,1]

plt.figure(figsize=(12,5))
plt.plot(wave,flux,linewidth = 0.1, c = 'cyan', alpha = 0.8)
plt.plot(lam,f+w, c = 'violet', linewidth = 0.1, alpha = 0.8)
plt.ylim(0,0.5e-12)
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
plt.ylim(0,.2e-12)
plt.axvline(364.7, linewidth = 0.5, c = 'r')
#plt.yscale('log')
#plt.xscale('log')
plt.xlim(910,1130)
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'Flux ($erg\ cm^{-2} s^{-1} \AA^{-1}$)')
plt.legend()

plt.show()

lib = pyphot.get_library()

wave = np.asarray(wave, dtype=float) * 10
flux = np.asarray(flux, dtype=float)

f = lib['GROUND_JOHNSON_V']
sed = f.get_flux(wave,flux)
sed1 = f.Vega_zero_flux
mag = -2.5 * np.log10(f.get_flux(wave,flux) / f.Vega_zero_flux)
print(mag)
