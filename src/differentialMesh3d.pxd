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

  cdef int num_sources

  cdef double source_rad

  cdef double *DX

  cdef double *DY

  cdef double *DZ

  cdef double *SX

  cdef double *SY

  cdef double *SZ

  cdef Zonemap3d source_zonemap

  ## FUNCTIONS

  cdef int __edge_vertex_force(self, int he1, int v1, double scale) nogil

  cdef int __triangle_force(self, double scale) nogil

  cdef int __reject(self, double scale) nogil

  cdef int __attract(self, double scale) nogil

  cdef int __unfold(self, double scale) nogil

  cdef int __find_nearby_sources(self) nogil

  cdef int __smooth_intensity(self) nogil

  ## EXTERNAL

  cpdef int initialize_sources(self, list sources, double source_rad)

  cpdef int find_nearby_sources(self)

  cpdef int smooth_intensity(self)

  cpdef int position_noise(self, np.ndarray[double, mode="c",ndim=2] a, int scale_intensity)

  cpdef int optimize_position(self, double step, int itt, int scale_intensity)

