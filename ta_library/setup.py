from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

examples_extension = Extension(
    name="etalib",
    sources=["etalib.pyx"],
    libraries=["eta"],
    library_dirs=["lib"],
    include_dirs=["lib", numpy.get_include()]
)
setup(
    name="etalib",
    ext_modules=cythonize([examples_extension]),
    include_dirs=[numpy.get_include()]
)