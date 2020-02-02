MasSpOT - The Sinkhorn algorithm for spectra comparison
=======================================================

A tool for comparing MS spectra using Sinkhorn algorithm [1]. Basing application is clustering of MSI pictures.

# Installation

1. Obtain masserstein package: https://github.com/mciach/masserstein
2. Obtain [eigen library](eigen.tuxfamily.org) and build it under MasSpOT/eigen.
3. Run: python setup.py install.

# Usage
Examples/Clustering example.ipynb file contains an exemplary usage on MSI picture of [mouse cerebellum](https://www.ebi.ac.uk/metabolights/MTBLS487).

# References
[1] Chizat, L., Peyré, G., Schmitzer, B., & Vialard, F.-X. (2018). Scaling algorithms for unbalanced optimal transport problems. Mathematics of Computation, 87(314), 2563–2609.

[2] Bond NJ, Koulman A, Griffin JL, Hall Z. massPix: an R package for annotation and interpretation of mass spectrometry imaging data for lipidomics. Metabolomics. 2017;13(11) 128.
