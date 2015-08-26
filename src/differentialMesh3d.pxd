# -*- coding: utf-8 -*-
# cython: profile=True

from __future__ import division

cimport mesh3d
from libc.stdlib cimport malloc, free
from zonemap3d cimport Zonemap3d


cdef class DifferentialMesh3d(mesh3d.Mesh3d):

  cdef float nearl

  cdef float farl

  cdef int num_sources

  cdef float source_rad

  cdef float *DX

  cdef float *DY

  cdef float *DZ

  cdef float *SX

  cdef float *SY

  cdef float *SZ

  #cdef Zonemap3d source_zonemap

  ## FUNCTIONS

  cdef int __reject(self, float scale) nogil

  cdef int __attract(self, float scale) nogil

  ## EXTERNAL

  cpdef int optimize_position(self, float step, int itt=*)

