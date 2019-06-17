import os

from distutils.core import setup, Extension

gmf = Extension('gmf', 
                define_macros=[('BUILDING_PYTHON_MODULE',)],
                include_dirs=['C:\\fftw\\include', 'C:\\Users\\cimat\\Anaconda3\\envs\\ptenv\\Lib\\site-packages\\numpy\\core\\include', 'C:\\Users\\cimat\\Anaconda3\\envs\\ptenv\\include'],
                library_dirs=['C:\\fftw\\lib', 'C:\\Users\\cimat\\Anaconda3\\envs\\ptenv\\Lib\\site-packages\\numpy\\core\\lib', 'C:\\Users\\cimat\\Anaconda3\\envs\\ptenv\\libs'],
                libraries=['fftw3', 'npymath', 'C:\\Users\\cimat\\Anaconda3\\envs\\ptenv\\libs\\_tkinter', 'C:\\Users\\cimat\\Anaconda3\\envs\\ptenv\\libs\\python3', 'C:\\Users\\cimat\\Anaconda3\\envs\\ptenv\\libs\\python37'],
                sources = ['.\\include\\gmf.c'])

setup(name='gmf',
        version='1.0',
        description='Gaussian matchde filter for numpy arrays',
        author='Fernando Cervantes-Sanchez',
        author_email='iie.fercer@gmail.com',
        ext_modules = [gmf])
