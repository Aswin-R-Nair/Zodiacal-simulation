from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy

# setup(
#     ext_modules=[
#         Extension("cSim2", ["cSim2.c"],
#                   include_dirs=[numpy.get_include()]),
#     ],
# )


setup(
    ext_modules=cythonize("Sim2_cmod.pyx"),
    include_dirs=[numpy.get_include()]
)