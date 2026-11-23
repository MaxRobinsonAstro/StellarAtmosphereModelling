import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.nddata import CCDData
import ccdproc
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def frames(file_type, file_amount):
    files = []
    for i in range(1, file_amount + 1):
        files.append(f'OneDrive - Maynooth University\Project\AB Aurigae\{file_type}-{i:04d}.fits')
    return [CCDData.read(f, unit = 'adu') for f in files]

def sci_frames(file_type, file_amount):
    files = []
    for i in range(1, file_amount + 1):
        files.append(f'OneDrive - Maynooth University\Project\AB Aurigae\AB Aurigae-{i:04d}{file_type}.fits')
    return [CCDData.read(f, unit = 'adu') for f in files]

def flat_frames(file_type, file_amount):
    files = []
    for i in range(1, file_amount + 1):
        files.append(f'OneDrive - Maynooth University\Project\AB Aurigae\Flats_12_11-{i:04d}-{file_type}.fits')
    return [CCDData.read(f, unit = 'adu') for f in files]

plt.style.use('dark_background')
plt.rcParams.update({
    'text.usetex' : True,
    'font.family' : 'serif',
    'font.serif' : 'cmr10',
    'axes.formatter.use_mathtext' : True,
    'axes.unicode_minus' : False,
    'font.size' : 12 })

"""Bias files"""

MasterBias = ccdproc.combine(frames('Bias_12_11', 5), 'MasterBias.fits', 'median', overwrite_output = True)
plt.figure(figsize = (5,5))
plt.imshow(MasterBias.data, cmap = 'magma', origin = 'lower', norm = LogNorm())
plt.colorbar()

"""Dark files"""

MasterDark = ccdproc.combine(frames('Dark_12_11', 5), 'MasterDark.fits', 'median', overwrite_output = True)
plt.figure(figsize = (5,5))
dark_image = MasterDark.data / 150
MasterDark = CCDData(dark_image, unit = 'adu')
plt.imshow(dark_image, cmap = 'magma', origin = 'lower', norm = LogNorm())
plt.colorbar()

"""Ha_flat files"""

MasterB_flat = ccdproc.combine(flat_frames('Ha', 3), 'MasterR_flat.fits', 'median', overwrite_output = True)
plt.figure(figsize = (5,5))
plt.imshow(MasterB_flat.data, cmap = 'magma', origin = 'lower', norm = LogNorm())
plt.colorbar()

"""B science files"""
n = 3
R_sci = sci_frames('Ha', n)

BiasSubDark = ccdproc.subtract_bias(MasterDark, MasterBias)
BiasSubR_flat = ccdproc.subtract_bias(MasterB_flat, MasterBias)

BiasSubR_sci = []
for i in range(n):
    BiasSubR_sci.append(ccdproc.subtract_bias(R_sci[i], MasterBias))

DarkSubR_flat = ccdproc.subtract_dark(BiasSubR_flat, BiasSubDark, dark_exposure = (2.0 * u.second),
                                    data_exposure = (0.6 * u.second), scale = True)

DarkSubR_sci = []
for i in range(n):
    DarkSubR_sci.append(ccdproc.subtract_dark(BiasSubR_sci[i], BiasSubDark, dark_exposure = (2.0 * u.second),
                                    data_exposure = (2.0 * u.second), scale = True))

FinalR_sci = []
for i in range(n):
    FinalR_sci.append(ccdproc.flat_correct(DarkSubR_sci[i], DarkSubR_flat))

#FinalR_sci = np.asarray(FinalR_sci)
#FinalR_sci = CCDData(FinalR_sci, unit = 'adu')

MasterR_sci = ccdproc.combine(FinalR_sci, 'MasterR_sci.fits', 'median', overwrite_output = True)

#MasterR_sci = ccdproc.cosmicray_lacosmic(MasterR_sci, sigclip = 3.0)
#MasterR_sci.write('OneDrive - Maynooth University\Project\OHP_sample_data\MasterR_sci.fits', overwrite = True)

image = MasterR_sci.data
plt.figure(figsize = (5,5))
plt.hist(image.flat, 1000, log = True)

vmin, vmax = np.percentile(image, [1, 99.9])
vmin = 1
plt.figure(figsize = (5,5))
plt.imshow(image, cmap = 'magma', origin = 'lower', norm = LogNorm(vmin, vmax))
plt.colorbar()
plt.show()