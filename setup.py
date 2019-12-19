#!/usr/bin/env python

import numpy as np
from distutils.core import setup, Extension

setup(name='gmf',
      version='1.0',
      description='Patch extraction tools for coronary angiograms',
      author='Fernando Cervantes',
      author_email='iie.fercer@gmail.com',
      
      ext_modules=[Extension('gmf', ['include/gmf.c'],
                             define_macros=[('BUILDING_PYTHON_MODULE',''), ('NDEBUG',)],
                             include_dirs=[np.get_include(), 'usr/local/include'],
                             libraries=['fftw3', 'm'],
                             )],
      )
