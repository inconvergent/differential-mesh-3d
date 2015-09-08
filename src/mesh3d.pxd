# -*- coding: utf-8 -*-
# cython: profile=True

cimport numpy as np
from zonemap3d cimport Zonemap3d

cdef struct s_HE:
  int id # half edge id
  int gen # generation of edge
  int first # vertex 1
  int last # vertex 2
  int twin # twin edge
  int next # next edge id
  int face # adjacent face id

ctypedef s_HE sHE


cdef class Mesh3d:

  cdef int nmax

  cdef int vnum

  cdef int henum

  cdef int fnum

  cdef int nz

  cdef double zonewidth

  ## ARRAYS

  cdef double *X

  cdef double *Y

  cdef double *Z

  cdef double *I

  cdef sHE *HE

  cdef int *VHE

  cdef int *FHE

  cdef Zonemap3d zonemap

  ## FUNCTIONS

  ## INTERNAL

  cdef int __valid_new_vertex(self, double x, double y, double z) nogil

  cdef int __new_vertex(self, double x, double y, double z) nogil

  cdef int __new_edge(self, int first, int last) nogil

  cdef int __new_edge_from_edge(self, int he1, int last) nogil

  cdef int __new_face(self, int he1) nogil

  cdef void __set_face_of_three_edges(self, int face, int he1, int he2, int he3) nogil

  cdef void __set_gen_of_three_edges(self, int gen, int he1, int he2, int he3) nogil

  cdef void __set_edge_of_face(self, int face, int he1) nogil

  cdef void __set_mutual_twins(self, int he1, int he2) nogil

  cdef int __is_surface_edge(self, int t1) nogil

  cdef int __next_surface(self, int he1, int direction) nogil

  cdef int __edge_duplicate_test(self, int he1, int a, int b) nogil

  cdef int __flip_edge(self, int he1, double limit) nogil

  cdef int __split_edge(self, int he1) nogil

  cdef int __set_next_of_triangle(self, int he1, int he2, int he3) nogil

  cdef int __split_internal_edge(self, int he1) nogil

  cdef int __split_surface_edge(self, int he1) nogil

  cdef int __split_all_longest_triangle_edges(self, double limit) nogil

  cdef double __get_edge_length(self, int he1) nogil

  cdef int __edge_integrity(self, int he1) nogil

  cdef int __triangle_integrity(self, int face) nogil

  cdef int __safe_vertex_positions(self, double limit) nogil

  cdef double __get_edge_intensity(self, int he1) nogil

  cdef void __set_vertex_intensity(self, int v1, double i) nogil

  cdef void __set_edge_intensity(self, int he1, double i) nogil

  ## EXTERNAL

  cpdef int edge_integrity(self, int he1)

  cpdef int triangle_integrity(self, int face)

  cpdef int safe_vertex_positions(self, double limit)

  cpdef list new_faces_in_ngon(self, double x1, double y1, double z1, double rad, int num)

  cpdef int next_surface(self, int he1, int direction)

  cpdef int flip_edge(self, int he1, double limit)

  cpdef int split_edge(self, int he1)

  cpdef int optimize_edges(self, double split_limit, double flip_limit)

  ## GET DATA

  cpdef int np_get_vertices(self, np.ndarray[double, mode="c",ndim=2] a)

  cpdef int np_get_triangles_vertices(self, np.ndarray[long, mode="c",ndim=2] a)

  cpdef int np_get_triangles_gen(self, np.ndarray[long, mode="c",ndim=1] a)

  cpdef double get_edge_intensity(self, int he1)

  cpdef int set_vertex_intensity(self, int v1, double i)

  cpdef int set_edge_intensity(self, int he1, double i)

  cpdef double get_triangle_intensity(self, int f1)

  cpdef list get_triangle_edges(self, int f1)

  ## INFO

  cpdef int initiate_faces(self, list vertices, list faces)

  cpdef int diminish_all_vertex_intensity(self, double d)

  cpdef int is_surface_edge(self, int t1)

  cpdef double get_edge_length(self, int he1)

  cpdef dict get_edge_dict(self, int he1)

  cpdef list get_triangle_dicts(self, int f)

  cpdef int get_vnum(self)

  cpdef int get_henum(self)

  cpdef int get_fnum(self)

