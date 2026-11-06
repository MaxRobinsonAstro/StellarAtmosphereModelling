import pyphot
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append(r'GitHub/StellarAtmosphereModelling')
from HIP_84543_spectra import wave, flux
from Selection import df0

lib = pyphot.get_library()

fig, ax1 = plt.subplots(figsize=(12,5))
ax1.plot(wave,flux,linewidth = 0.1, c = 'cyan', alpha = 0.8)

ax1.set_ylim(0,2e-12)
#ax1.axvline(364.7, linewidth = 0.5, c = 'r')
ax1.set_xlim(0,2500)
ax1.set_xlabel(r'Wavelength (nm)')
ax1.set_ylabel(r'Flux ($erg\ cm^{-2} s^{-1} \AA^{-1}$)')

ax2 = ax1.twinx()

filters = ['TYCHO_B_MvB',
           'TYCHO_V_MvB',
           '2MASS_Ks']

for name in filters:
    f = lib[name]
    wave_f = f.wavelength/10
    transmit_f = f.transmit / np.max(f.transmit)
    ax2.plot(wave_f, transmit_f, label = f.name, lw = 0.7)

ax2.set_ylim(0, 1.2)
ax2.legend()

star = df0.loc[df0['Name'] == 'HIP 84543'].iloc[0]
B_mag = star['BTmag']
V_mag = star['VTmag']
K_mag = star['Ksmag']

print(f'BTmag: {B_mag}')
print(f'VTmag: {V_mag}')
print(f'Ksmag: {K_mag}')

F_syn = []
for name in filters:
    f = lib[name]
    #f.info(show_zeropoints = True)
    fluxes = f.get_flux(wave * 10, flux)
    mag = -2.5 * np.log10(fluxes.value) - f.Vega_zero_mag
    F_syn.append(mag)
    print(f'{f.name}: {mag}')

F_syn = np.array(F_syn)
F_cat = np.array([B_mag, V_mag, K_mag])

S = 2 * np.sum(F_cat * F_syn) / np.sum(F_syn)**2
flux_S = flux * S


plt.figure(figsize=(12,5))
plt.plot(wave,flux_S,linewidth = 0.1, c = 'cyan', alpha = 0.8)

plt.ylim(0,2e-12)
#ax1.axvline(364.7, linewidth = 0.5, c = 'r')
plt.xlim(0,2500)
plt.xlabel(r'Wavelength (nm)')
plt.ylabel(r'Flux ($erg\ cm^{-2} s^{-1} \AA^{-1}$)')

for name in filters:
    f = lib[name]
    #f.info(show_zeropoints = True)
    fluxes = f.get_flux(wave * 10, flux_S)
    mag = -2.5 * np.log10(fluxes.value) - f.Vega_zero_mag
    print(f'{f.name}: {mag}')

plt.show()


