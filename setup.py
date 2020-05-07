from setuptools import find_packages, setup
from distutils.core import setup, Extension

module1 = Extension('MasSpOTCppToPy',
                    sources = ['MasSpOT/optimal_transport.cpp'],
                    include_dirs = ['MasSpOT/eigen'])


setup(
    name='MasSpOT',
    packages=find_packages(),
    version='0.1.0',
    description='Comparing MS spectra using Sinkhorn algorithm.',
    author='Grzegorz Skoraczynski',
    license='MIT',
    install_requires=["cffi", "matplotlib", "scikit-learn"],
    python_requires='>=3.6, <4',
    ext_modules=[module1],
)
