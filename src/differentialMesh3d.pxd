# -*- coding: utf-8 -*-
# cython: profile=True

from __future__ import division

cimport mesh3d

from libc.stdlib cimport malloc, free
from zonemap3d cimport Zonemap3d


cdef class DifferentialMesh3d(mesh3d.Mesh3d):

  cdef double nearl
  cdef double farl

  cdef double *DX
  cdef double *DY
  cdef double *DZ
  cdef double *DI

  ## FUNCTIONS

  cdef long __reject(
    self,
    long v1,
    double *dx,
    double *dy,
    double *dz,
    double stp,
    long *vertices,
    double *dst,
    long num,
    long *connected,
    long cnum
  ) nogil

  cdef long __attract(
    self,
    long v1,
    double *dx,
    double *dy,
    double *dz,
    double stp,
    long *vertices,
    long num,
  ) nogil

  cdef long __smooth_intensity(
    self,
    long v1,
    double alpha,
    double *old,
    double *new,
    long *vertices,
    long num
  ) nogil

  ## EXTERNAL

  cpdef long position_noise(
    self,
    double[:,:] a,
    long scale_intensity
  )

  cpdef long optimize_position(
    self,
    double reject_stp,
    double attract_stp,
    double diminish,
    double alpha,
    long scale_intensity
  )

