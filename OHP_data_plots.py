from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

hdul = fits.open(r'OneDrive - Maynooth University\Project\OHP_sample_data\MasterR_sci.fits')
image = hdul[0].data

vmin, vmax = np.percentile(image, [1, 99])

plt.figure(figsize=(5,5))
plt.imshow(image, cmap = 'twilight_shifted', origin = 'lower', norm = LogNorm(vmin, vmax))
plt.colorbar(label = 'pixel value')
plt.title('Dumbbell Nebula R Filter')
plt.xlabel('x')
plt.ylabel('y')

plt.show()