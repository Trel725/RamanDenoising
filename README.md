# RamanDenoising

This is a python translation of MATLAB code, available in supplementary information for article  
Barton, Sinead J., Tomas E. Ward, and Bryan M. Hennelly. "Algorithm for optimal denoising of Raman spectra." Analytical Methods 10.30 (2018): 3759-3769.  

The bottleneck is rewritten in Numba which gave about 88x speedup (should be comparable with original implementation).

**Warning**: since this code was translated, it looks painfully non-pythonic. Freakin' MATLAB!

# Usage

You can use that code like

```python
from mlesg import MLESG
# given x - intensities of Raman spectrum,
# wvn - wavenumbers of the spectrum
peaks = find_peaks(x, prominence=20, distance=10)[0]
filtered = MLESG(x, wvn, peaks)
plt.plot(wvn, filtered)

# if you have example of noise, you can slightly improve filtering
# by providing the noise to the function as MLESG(x, wvn, peaks, mu=noise)
```
