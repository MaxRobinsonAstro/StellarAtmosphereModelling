import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

files = ['OneDrive - Maynooth University/Project/HIP_84543/HIP_84543_UVB.fits',
         'OneDrive - Maynooth University/Project/HIP_84543/HIP_84543_VIS.fits',
         'OneDrive - Maynooth University/Project/HIP_84543/HIP_84543_NIR.fits'] #fits files from 3 Xshooter arms

s, data, waves, fluxes, e = [], [], [], [], [] #2d arrays for wavelength, flux and error

for f in files:
    s.append(fits.open(f))


for i in range(len(s)):
    data.append(s[i][1].data)

for j in range(3):
    date = s[j][0].header['DATE']
    date_obs = s[j][0].header['DATE-OBS']
    obj = s[j][0].header['OBJECT']
    start = s[j][0].header['MJD-OBS']
    end = s[j][0].header['MJD-END']
    airmass1 = s[j][0].header['HIERARCH ESO TEL AIRM START']
    airmass2 = s[j][0].header['HIERARCH ESO TEL AIRM END']
    seeing1 = s[j][0].header['HIERARCH ESO TEL AMBI FWHM START']
    seeing2 = s[j][0].header['HIERARCH ESO TEL AMBI FWHM END']
    
        
    print(date)
    print(date_obs)
    print(obj)
    print(start)
    print(end)
    print(airmass1)
    print(airmass2)
    print(seeing1)
    print(seeing2)

for i in range(len(data)):
    waves.append(data[i]['WAVE'][0])
    fluxes.append(data[i]['FLUX'][0])
    e.append(data[i]['ERR'][0])