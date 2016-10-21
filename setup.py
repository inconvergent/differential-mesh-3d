#!/usr/bin/python3

try:
  from setuptools import setup
  from setuptools.extension import Extension
except Exception:
  from distutils.core import setup
  from distutils.extension import Extension

from Cython.Build import cythonize
from Cython.Distutils import build_ext

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
    version = '0.1.0',
    author = '@inconvergent',
    install_requires = ['numpy', 'cython'],
    license = 'MIT',
    cmdclass={'build_ext' : build_ext},
    ext_modules = cythonize(extensions)
    )
