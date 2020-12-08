# RamanDenoising

This is a python translation of MATLAB code, available in supplementary information for article  
Barton, Sinead J., Tomas E. Ward, and Bryan M. Hennelly. "Algorithm for optimal denoising of Raman spectra." Analytical Methods 10.30 (2018): 3759-3769.  

The bottleneck is rewritten in Numba which gave about 88x speedup (should be comparable with original implementation).

**Warning**: since this code was translated, it looks painfully non-pythonic. *Freakin' MATLAB!*