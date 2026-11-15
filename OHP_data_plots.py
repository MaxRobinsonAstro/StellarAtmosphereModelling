from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

plt.style.use('dark_background')
plt.rcParams.update({
    'font.family' : 'serif',
    'font.serif' : 'cmr10',
    'mathtext.fontset' : 'cm',
    'axes.formatter.use_mathtext' : True,
    'axes.unicode_minus' : False,   
    'font.size' : 12 })


hdul = fits.open(r'OneDrive - Maynooth University\Project\OHP_sample_data\MasterR_sci.fits')
image = hdul[0].data

vmin, vmax = np.percentile(image, [1, 99])

plt.figure(figsize=(7,7))
plt.imshow(image, cmap = 'twilight_shifted', origin = 'lower', norm = LogNorm(vmin, vmax))
plt.colorbar(label = 'pixel value')
plt.title('Dumbbell Nebula R Filter')
plt.xlabel('x')
plt.ylabel('y')

plt.show()