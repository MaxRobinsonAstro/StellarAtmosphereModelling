import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.nddata import CCDData
import ccdproc
import matplotlib.pyplot as plt

bias_files = ['OneDrive - Maynooth University\Project\OHP_sample_data\Bias-001.fit',
              'OneDrive - Maynooth University\Project\OHP_sample_data\Bias-002.fit',
              'OneDrive - Maynooth University\Project\OHP_sample_data\Bias-003.fit',
              'OneDrive - Maynooth University\Project\OHP_sample_data\Bias-004.fit',
              'OneDrive - Maynooth University\Project\OHP_sample_data\Bias-005.fit',]

bias = []
i = 0
for b in bias_files:
    bias.append(fits.getdata(b))
    i += 1
    
bias = CCDData(bias, unit = 'adu')

MasterBias = ccdproc.combine(bias, 'MasterBias.fit', 'median', overwrite_output = True)

flat = fits.getdata('OneDrive - Maynooth University\Project\OHP_sample_data\\flat.fit')
MasterFlat = CCDData(flat, unit = 'adu')
#MasterFlat = ccdproc.combine(flat, 'MasterFlat.fit', 'median', overwrite_output = True)

dark = fits.getdata('OneDrive - Maynooth University\Project\OHP_sample_data\dark_600s.fit')
MasterDark = CCDData(dark, unit = 'adu')
#MasterDark = ccdproc.combine(dark, 'MasterDark.fit', 'median', overwrite_output = True)

sci = fits.getdata('OneDrive - Maynooth University\Project\OHP_sample_data\CAS_A_SII_900s.fit')
sci = CCDData(sci, unit = 'adu')

BiasSubDark = ccdproc.subtract_bias(MasterDark, MasterBias)
BiasSubFlat = ccdproc.subtract_bias(MasterFlat, MasterBias)
BiasSubSci = ccdproc.subtract_bias(sci, MasterBias)

DarkSubFlat = ccdproc.subtract_dark(BiasSubFlat, BiasSubDark, dark_exposure = (600.0 * u.second),
                                    data_exposure = (13.0 * u.second), scale = True)

DarkSubSci = ccdproc.subtract_dark(BiasSubSci, BiasSubDark, dark_exposure = (600.0 * u.second),
                                    data_exposure = (900.0 * u.second), scale = True)

FinalSci = ccdproc.flat_correct(DarkSubSci, DarkSubFlat)
FinalSci = np.asarray(FinalSci)
FinalSci = CCDData(FinalSci, unit = 'adu')

MasterScience = ccdproc.cosmicray_lacosmic(FinalSci, sigclip = 3.0)
MasterScience.write('OneDrive - Maynooth University\Project\OHP_sample_data\MasterScience.fit')

"""plt.figure(figsize = (5,5))
plt.imshow(, cmap = 'gray')
plt.colorbar()
plt.show()"""