from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy



setup(
    ext_modules=cythonize("summer_persistence_analysis.pyx"),
    include_dirs=[numpy.get_include()]
)


# python setup_cython.py build_ext --inplace
