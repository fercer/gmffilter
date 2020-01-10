#!/usr/bin/env python

import numpy as np
from distutils.core import setup, Extension
import os

setup(name='gmf',
      version='1.0',
      description='Patch extraction tools for coronary angiograms',
      author='Fernando Cervantes',
      author_email='iie.fercer@gmail.com',
      
      ext_modules=[Extension('gmf', ['include/gmf.c'],
                             define_macros=[('BUILDING_PYTHON_MODULE',''), ('NDEBUG',)],
                             library_dirs=[os.environ['FFTW_LIBS_PATH']],
                             include_dirs=[np.get_include(), os.environ['FFTW_INCLUDE_PATH']],
                             libraries=['fftw3'],
                             )],
      )
