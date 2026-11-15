import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append(r'GitHub/StellarAtmosphereModelling')
from HIP_84543_spectra import wave, flux
from Selection import df0
from astropy.io import fits
from scipy import stats


fig, ax1 = plt.subplots(figsize=(12,5))
ax1.plot(wave,flux,linewidth = 0.25, c = 'cyan', alpha = 0.8)

ax1.set_ylim(0,2e-12)
#ax1.axvline(364.7, linewidth = 0.5, c = 'r')
ax1.set_xlim(350,800)
ax1.set_xlabel(r'Wavelength (nm)')
ax1.set_ylabel(r'Flux ($erg\ cm^{-2} s^{-1} \AA^{-1}$)')

ax2 = ax1.twinx()

filters = ['OneDrive - Maynooth University/Project/BESSEL_B.fit',
           'OneDrive - Maynooth University/Project/BESSEL_V.fit',
           'OneDrive - Maynooth University/Project/BESSEL_R.fit']
colours = ['cyan','green','red']
names = ['B', 'V', 'R']

transmit_f = []

for i in range(len(filters)):
    with fits.open(filters[i]) as s:
        transmit_f.append(s[0].data)
        wave_f = np.linspace(375,740,len(transmit_f[i]))
        ax2.plot(wave_f, transmit_f[i], lw = 0.8, color = colours[i], label = names[i])

ax2.set_ylim(0, 1.2)
ax2.legend()

star = df0.loc[df0['Name'] == 'HIP 84543'].iloc[0]
B_mag = star['BTmag']
V_mag = star['VTmag']
R_mag = star['Rmag']

print('Catalogue:')
print(f'BTmag: {B_mag}')
print(f'VTmag: {V_mag}')
print(f'Rmag: {R_mag}')

T_lam = transmit_f[0] + transmit_f[1] + transmit_f[2] 

ax2.plot(wave_f, T_lam, lw = 0.8)

flux_binned = np.zeros(len(wave_f))
min_phot = np.where(wave == 375)[0][0]
max_phot = np.where(wave == 740)[0][0]
flux_binned, wave_binned, binnumber = stats.binned_statistic(wave[min_phot:max_phot], flux[min_phot:max_phot], bins = len(wave_f), statistic = 'mean')

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

plt.figure(figsize=(12,5))
plt.plot(wave,flux, c = 'red', lw = 0.1, alpha = 0.5)

plt.ylim(0,1.5e-12)
#ax1.axvline(364.7, linewidth = 0.5, c = 'r')
plt.xlim(0,2500)
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'Flux ($erg\ cm^{-2} s^{-1} \AA^{-1}$)')

S = 10**(-0.4 * (np.abs(B_mag - magB)))
flux_phot = flux * S
flux_binned_phot = flux_binned * S

magB_scaled = -2.5 * np.log10((np.sum(flux_binned_phot * transmit_f[0])/np.sum(transmit_f[0])) / 632.5e-11)
magV_scaled = -2.5 * np.log10((np.sum(flux_binned_phot * transmit_f[1])/np.sum(transmit_f[1])) / 363.1e-11)
magR_scaled = -2.5 * np.log10((np.sum(flux_binned_phot * transmit_f[2])/np.sum(transmit_f[2])) / 217.7e-11)
print('Scaled by B:')
print(f'Bmag: {magB_scaled}')
print(f'Vmag: {magV_scaled}')
print(f'Rmag: {magR_scaled}')

plt.figure(figsize=(12,5))
plt.plot(wave,flux_phot, c = 'red', lw = 0.1, alpha = 0.5)

plt.ylim(0,1.5e-12)
#ax1.axvline(364.7, linewidth = 0.5, c = 'r')
plt.xlim(0,2500)
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'Flux ($erg\ cm^{-2} s^{-1} \AA^{-1}$)')



plt.show()
