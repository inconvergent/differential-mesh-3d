# -*- coding: utf-8 -*-
# cython: profile=True

from __future__ import division

cimport mesh3d
cimport numpy as np

from libc.stdlib cimport malloc, free
from zonemap3d cimport Zonemap3d


cdef class DifferentialMesh3d(mesh3d.Mesh3d):

  cdef double nearl

  cdef double farl

  cdef double *DX

  cdef double *DY

  cdef double *DZ

  ## FUNCTIONS

  cdef long __edge_vertex_force(self, long he1, long v1, double scale) nogil

  cdef long __triangle_force(self, double scale) nogil

  cdef long __reject(self, long v, double stp, long *vertices, double *dst) nogil

  cdef long __attract(self, long v1, double stp)

  cdef long __unfold(self, double scale) nogil

  cdef long __smooth_intensity(self, double alpha) nogil

  ## EXTERNAL

  cpdef long smooth_intensity(self, double alpha)

  cpdef long position_noise(
    self,
    np.ndarray[double, mode="c",ndim=2] a,
    long scale_intensity
  )

  cpdef long optimize_position(
    self,
    double reject_stp,
    double triangle_stp,
    double attract_stp,
    double unfold_stp,
    double cohesion_stp,
    long itt,
    long scale_intensity
  )

