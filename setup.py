#!/usr/bin/env python-sirius

import pathlib

from setuptools import setup


def get_abs_path(relative):
    """."""
    return str(pathlib.Path(__file__).parent / relative)


with open(get_abs_path("README.md"), "r") as _f:
    _long_description = _f.read().strip()

with open(get_abs_path("VERSION"), "r") as _f:
    __version__ = _f.read().strip()

with open(get_abs_path("requirements.txt"), "r") as _f:
    _requirements = _f.read().strip().split("\n")

setup(
    name='fieldmaptrack',
    version=__version__,
    author='lnls-fac',
    description='Fieldmap analysis utilities',
    long_description=_long_description,
    url='https://github.com/lnls-fac/fieldmaptrack',
    download_url='https://github.com/lnls-fac/fieldmaptrack',
    license='MIT License',
    classifiers=[
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering'
    ],
    packages=['fieldmaptrack'],
    install_requires=_requirements,
    package_data={'fieldmaptrack': ['VERSION']},
    scripts=[
        'scripts/fac-fma-analysis.py',
        'scripts/fac-fma-model.py',
        'scripts/fac-fma-multifunctional-sextupole.py',
        'scripts/fac-fma-multipoles.py',
        'scripts/fac-fma-profile.py',
        'scripts/fac-fma-rawfield.py',
        'scripts/fac-fma-sextupole.py',
        'scripts/fac-fma-trajectory.py',
    ],
    zip_safe=False,
)
