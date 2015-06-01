#!/usr/bin/env python3

from setuptools import setup

with open('VERSION','r') as _f:
    __version__ = _f.read().strip()

setup(
    name='fieldmaptrack',
    version=__version__,
    author='lnls-fac',
    description='Fieldmap analysis utilities',
    url='https://github.com/lnls-fac/fieldmaptrack',
    download_url='https://github.com/lnls-fac/fieldmaptrack',
    license='MIT License',
    classifiers=[
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering'
    ],
    packages=['fieldmaptrack'],
    package_data={'fieldmaptrack': ['VERSION']},
    zip_safe=False,
)
