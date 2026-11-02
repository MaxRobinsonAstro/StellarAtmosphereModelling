import pyphot
import numpy as np
from HIP_84543_spectra import wave, flux

lib = pyphot.get_library()
print(lib)
f = lib['GROUND_JOHNSON_V']

mag = -2.5 * np.log10(f.get_flux(wave,flux) / f.Vega_zero_flux)





