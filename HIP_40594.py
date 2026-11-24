"""
Spectra for HIP_40594
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy import special
import sys
sys.path.append(r'GitHub/StellarAtmosphereModelling')
from Selection import df0
from scipy import stats

h = 6.62607015e-27  # in erg cgs
c = 2.99792458e10  
k = 1.380649e-16

def planck_law(wavelength, temperature, dist):
    """Planck's Law gives the blackbody curve"""
    exponent = (h * c) / (wavelength * 1e-7 * k * temperature)
    sr = (1 / (4 * np.pi * (dist * 3.0857e10)**2)) * (2 * h * c**2) / ((wavelength*1e-7)**5 * (np.exp(exponent) - 1))
    return sr



def dimm(w, f, sigma_dimm, airmass, width):
    """Calculates T (slit loss factor"""
    sigma_lam = sigma_dimm * (w/500)**(-0.2) * airmass**(0.6)
    T = special.erf(np.sqrt(np.log(2)) * width / sigma_lam)
    return f / T


files = ['OneDrive - Maynooth University/Project/HIP_40594/HIP_40594_UVB.fits',
         'OneDrive - Maynooth University/Project/HIP_40594/HIP_40594_VIS.fits',
         'OneDrive - Maynooth University/Project/HIP_40594/HIP_40594_NIR.fits'] #fits files from 3 Xshooter arms

waves, fluxes, e = [], [], [] #2d arrays for wavelength, flux and error

for f in files:
    with fits.open(f) as s:
        data = s[1].data #data is contained in first extension
        waves.append(data['WAVE'][0])
        fluxes.append(data['FLUX'][0])
        e.append(data['ERR'][0])
        
        date = s[0].header['DATE']
        date_obs = s[0].header['DATE-OBS']
        obj = s[0].header['OBJECT']
        start = s[0].header['MJD-OBS']
        end = s[0].header['MJD-END']
        airmass1 = s[0].header['HIERARCH ESO TEL AIRM START']
        airmass2 = s[0].header['HIERARCH ESO TEL AIRM END']
        seeing1 = s[0].header['HIERARCH ESO TEL AMBI FWHM START']
        seeing2 = s[0].header['HIERARCH ESO TEL AMBI FWHM END']
        
            
        print(f'date: {date}')
        print(f'date_obs: {date_obs}')
        print(f'obj: {obj}')
        print(f'start: {start}')
        print(f'end: {end}')
        print(f'airmass1: {airmass1}')
        print(f'airmass2: {airmass2}')
        print(f'seeing1: {seeing1}')
        print(f'seeing2: {seeing2}')

"""Plotting arms before concatenating"""

plt.style.use('default')
plt.rcParams.update({
    'text.usetex' : True,
    'font.family' : 'serif',
    'font.serif' : 'cmr10',
    'axes.formatter.use_mathtext' : True,
    'axes.unicode_minus' : False,
    'font.size' : 12 })


colours = ['blue','green','red']
plt.figure(figsize=(12,5))
for i in range(3):
    plt.plot(waves[i],fluxes[i],linewidth = 0.1, alpha = 0.8, color = colours[i])
plt.ylim(0,8e-12)
plt.xlim(200,2500)
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'Flux ($erg\ cm^{-2} s^{-1} \AA^{-1}$)')

plt.savefig("GitHub/StellarAtmosphereModelling/figures/HIP_40594_separate_before_slit")

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


airmass = (airmass2 + airmass1) / 2
sigma_dimm0 = (0.756 + 0.808) / 2
sigma_dimm1 = (0.756 + 0.808) / 2
sigma_dimm2 = (0.756 + 0.808 + 0.891) / 3


fluxes_slitloss = fluxes.copy()
fluxes_slitloss[0] = dimm(waves[0], fluxes_slitloss[0], sigma_dimm0, airmass, 1) #applying DIMM seeing to the fluxes for different slit widths on each arm
fluxes_slitloss[1] = dimm(waves[1], fluxes_slitloss[1], sigma_dimm1, airmass, 0.9)
fluxes_slitloss[2] = dimm(waves[2], fluxes_slitloss[2], sigma_dimm2, airmass, 0.6)

flux_slitloss = np.concatenate(fluxes_slitloss)

diff1 = np.abs((fluxes[0][-1] - fluxes[1][0]) / fluxes[0][-1] * 100) #percentage difference between arms
diff2 = np.abs((fluxes[1][-1] - fluxes[2][0]) / fluxes[1][-1] * 100)
print(f'Before Slitloss UVB-VIS: {diff1:.2f}%, VIS-NIR: {diff2:.2f}%')

diff1 = np.abs((fluxes_slitloss[0][-1] - fluxes_slitloss[1][0]) / fluxes_slitloss[0][-1] * 100) #percentage difference between arms
diff2 = np.abs((fluxes_slitloss[1][-1] - fluxes_slitloss[2][0]) / fluxes_slitloss[1][-1] * 100)
print(f'After Slitloss UVB-VIS: {diff1:.2f}%, VIS-NIR: {diff2:.2f}%')

star = df0.loc[df0['Name'] == 'HIP 40594'].iloc[0]

plt.figure(figsize=(12,5))
plt.plot(wave,err, c = 'red')
plt.ylim(0,0.5e-12)

B_mag = star['BTmag']
V_mag = star['VTmag']
R_mag = star['Rmag']

print('Catalogue:')
print(f'BTmag: {B_mag}')
print(f'VTmag: {V_mag}')
print(f'Rmag: {R_mag}')

filters = ['OneDrive - Maynooth University/Project/BESSEL_B.fit',
           'OneDrive - Maynooth University/Project/BESSEL_V.fit',
           'OneDrive - Maynooth University/Project/BESSEL_R.fit']
names = ['B', 'V', 'R']

transmit_f = []

for i in range(len(filters)):
    with fits.open(filters[i]) as s:
        transmit_f.append(s[0].data)
        wave_f = np.linspace(375,740,len(transmit_f[i]))

T_lam = transmit_f[0] + transmit_f[1] + transmit_f[2] 

flux_binned = np.zeros(len(wave_f))
min_phot = np.where(wave == 375)[0][0]
max_phot = np.where(wave == 740)[0][0]
flux_binned, wave_binned, binnumber = stats.binned_statistic(wave[min_phot:max_phot], flux_slitloss[min_phot:max_phot], bins = len(wave_f), statistic = 'mean')

for j in range(len(wave_binned) - 1):
    wave_binned[j] = (wave_binned[j] + wave_binned[j+1]) / 2
wave_binned = wave_binned[:-1]


magB = -2.5 * np.log10((np.sum(flux_binned * transmit_f[0])/np.sum(transmit_f[0])) / 632.5e-11)
magV = -2.5 * np.log10((np.sum(flux_binned * transmit_f[1])/np.sum(transmit_f[1])) / 363.1e-11)
magR = -2.5 * np.log10((np.sum(flux_binned * transmit_f[2])/np.sum(transmit_f[2])) / 217.7e-11)
  
print('Calculated:')
print(f'Bmag: {magB}')
print(f'Vmag: {magV}')
print(f'Rmag: {magR}')

S = 10**(-0.4 * (np.abs(B_mag - magB)))
flux_phot = flux_slitloss * S
flux_binned_phot = flux_binned * S
fluxes_phot = [[a * S for a in fluxes_slitloss[b]] for b in range(len(fluxes_slitloss))]

magB_scaled = -2.5 * np.log10((np.sum(flux_binned_phot * transmit_f[0])/np.sum(transmit_f[0])) / 632.5e-11)
magV_scaled = -2.5 * np.log10((np.sum(flux_binned_phot * transmit_f[1])/np.sum(transmit_f[1])) / 363.1e-11)
magR_scaled = -2.5 * np.log10((np.sum(flux_binned_phot * transmit_f[2])/np.sum(transmit_f[2])) / 217.7e-11)
print('Scaled by B:')
print(f'Bmag: {magB_scaled}')
print(f'Vmag: {magV_scaled}')
print(f'Rmag: {magR_scaled}')


fig, ax1 = plt.subplots(figsize=(10,5))
ax2 = ax1.twinx()
colours1 = ['red', 'blue', 'green']
for i in range(3):
    ax1.plot(waves[i],fluxes[i],linewidth = 0.2, alpha = 0.3, color = colours1[i])
    ax1.plot(waves[i],fluxes_slitloss[i],linewidth = 0.2, alpha = 0.6, color = colours1[i])
    ax1.plot(waves[i],fluxes_phot[i],linewidth = 0.2, alpha = 0.9, color = colours1[i])
    ax2.plot(wave_f, transmit_f[i], lw = 0.5, color = colours[i], label = names[i])
temp = star['Teff']
dist = star['Dist']
ax1.plot(wave,planck_law(wave, temp, dist), c ='magenta', lw = '0.4', label = 'Blackbody Curve')
ax1.set_ylim(0,1e-11)
ax2.set_ylim(0, 1.2)
ax1.set_xlim(o, 900)
ax1.set_xlabel(r'Wavelength (nm)')
ax1.set_ylabel(r'Flux ($erg\ cm^{-2} s^{-1} \AA^{-1}$)')
ax2.set_ylabel(r'Filter Transmission')
ax1.legend()
ax2.legend()

plt.savefig("GitHub/StellarAtmosphereModelling/figures/HIP_40594_corrections")

plt.figure(figsize=(10,5))
plt.plot(wave,flux,linewidth = 0.3, c = 'blue', alpha = 0.8, label = 'uncorrected')
plt.plot(wave,flux_slitloss,linewidth = 0.3, c = 'purple', alpha = 0.8, label = 'slit loss corrected')
plt.plot(wave,flux_phot,linewidth = 0.3, c = 'red', alpha = 0.8, label = 'slit loss + flux corrected')
plt.ylim(0,1e-11)
plt.axvline(364.7, linewidth = 0.5, c = 'r')
plt.xlim(o,600)
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'Flux ($erg\ cm^{-2} s^{-1} \AA^{-1}$)')
plt.legend()

plt.savefig("GitHub/StellarAtmosphereModelling/figures/HIP_40594__Balmer_jump_corrections")

plt.show()