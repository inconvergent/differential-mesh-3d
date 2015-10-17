#!/usr/bin/python

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext
import numpy

_extra = ['-fopenmp', '-O3', '-ffast-math']


extensions = [
  Extension('mesh3d',
            sources = ['./src/mesh3d.pyx'],
            extra_compile_args = _extra,
            extra_link_args = ['-fopenmp']
    ),
  Extension('differentialMesh3d',
            sources = ['./src/differentialMesh3d.pyx'],
            extra_compile_args = _extra,
            extra_link_args = ['-fopenmp']
  )
]

setup(
  name = "differential-mesh-3d",
  version = '0.0.1',
  author = '@inconvergent',
  install_requires = ['numpy>=1.8.2', 'cython>=0.20.0'],
  license = 'MIT',
  cmdclass={'build_ext' : build_ext},
  include_dirs = [numpy.get_include()],
  ext_modules = cythonize(extensions,include_path = [numpy.get_include()])
)
