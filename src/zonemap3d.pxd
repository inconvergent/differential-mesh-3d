# -*- coding: utf-8 -*-
# cython: profile=True

from __future__ import division

cdef struct s_Z:
  long i
  long size
  long count
  long *ZV

ctypedef s_Z sZ

cdef class Zonemap3d:

  cdef long vnum

  cdef long vsize

  cdef long nz

  cdef long total_zones

  cdef long greatest_zone_size

  ## ARRAYS

  cdef double *X

  cdef double *Y

  cdef double *Z

  cdef long *VZ

  cdef sZ **ZONES

  ## FUNCTIONS

  cdef void __init_zones(self) nogil

  cdef long __add_vertex(self, long v1) nogil

  cdef long __del_vertex(self, long v1) nogil

  cdef long __add_v_to_zone(self, long z1, long v1) nogil

  cdef long __extend_zv_of_zone(self, sZ *zone) nogil

  cdef long __remove_v_from_zone(self, long zone, long v1) nogil

  cdef long __get_z(self, double x, double y, double z) nogil

  cdef long __update_v(self, long v1) nogil

  cdef long __sphere_vertices(self, double x, double y, double z, double rad, long *vertices) nogil

  cdef long __sphere_vertices_dst(self, double x, double y, double z, double rad, long *vertices, double *dst) nogil

  cdef long __sphere_is_free(self, double x, double y, double z, double rad) nogil

  cdef long __sphere_is_free_ignore(self, double x, double y, double z, long v, double rad) nogil

  cdef long __get_max_sphere_count(self) nogil

  cdef long __get_greatest_zone_size(self) nogil

  cdef void __assign_xyz_arrays(self, double *x, double *y, double *z) nogil

  ## INFO

  cpdef list _perftest(self, long nmax, long num_points, long num_lookup)

  cpdef long add_vertex(self, long v1)

  cpdef long del_vertex(self, long v1)

  cpdef long update_v(self, long v1)

  cpdef long sphere_is_free(self, double x, double y, double z, double rad)

  cpdef long get_max_sphere_count(self)

  cpdef long get_greatest_zone_size(self)

  cpdef long get_vnum(self)

  cpdef list get_zone_info_dicts(self)

