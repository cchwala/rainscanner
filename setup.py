# -*- coding: utf-8 -*-


import os
from setuptools import setup


# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "rainscanner",
    version = "0.1.0",
    author = "Christian Chwala",
    author_email = "christian.chwala@kit.edu",
    description = "Python tools for processing SELEX rainscanner data",
    license = "BSD",
    keywords = "precipitation radar x-band wradlib rainscanner",
    url = "...",
    packages=['rainscanner'],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Atmospheric Science",
        "License :: OSI Approved :: BSD License",
        'Programming Language :: Python :: 2.7',
    ],
    # A list of all available classifiers can be found at 
    # https://pypi.python.org/pypi?%3Aaction=list_classifiers
    install_requires=[
        'wradlib',
        'numpy',
        'pandas',
        'xarray',
        'tqdm',
        'xmltodict'
    ],
)
