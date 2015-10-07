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

  cdef long num_sources

  cdef double source_rad

  cdef double *DX

  cdef double *DY

  cdef double *DZ

  cdef double *SX

  cdef double *SY

  cdef double *SZ

  cdef Zonemap3d source_zonemap

  ## FUNCTIONS

  cdef long __edge_vertex_force(self, long he1, long v1, double scale) nogil

  cdef long __triangle_force(self, double scale) nogil

  cdef long __reject(self, double scale) nogil

  cdef long __attract(self, double scale) nogil

  cdef long __unfold(self, double scale) #nogil

  cdef long __find_nearby_sources(self) nogil

  cdef long __smooth_intensity(self, double alpha) nogil

  ## EXTERNAL

  cpdef long initialize_sources(self, list sources, double source_rad)

  cpdef long find_nearby_sources(self)

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

