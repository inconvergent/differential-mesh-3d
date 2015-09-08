# -*- coding: utf-8 -*-
# cython: profile=True

from __future__ import division

cimport cython
from libc.stdlib cimport malloc, free

from libc.math cimport sqrt
from libc.math cimport pow as cpow
from libc.math cimport cos
from libc.math cimport sin
from libc.math cimport fabs
from libc.math cimport M_PI

from helpers cimport int_array_init
from helpers cimport double_array_init
from helpers cimport vcross

import numpy as np
cimport numpy as np
cimport cython

#cimport zonemap3d


cdef double TWOPI = M_PI*2


def dict_list_add(dict d, k, v):

  if d.has_key(k):
    d[k].append(v)
  else:
    d[k] = [v]


cdef class Mesh3d:
  """
  a triangular mesh optimized for operations like splitting edges in two and
  appending triangles to mesh surface. the structure also supports rebalancing
  by splitting long edges and flipping triangle edges. this code assumes that
  all faces are triangles, and everything will fail horribly if that is not the
  case.

  we utilize the half-edge data structure to store edges. this means that all
  edges have a two associated vertices (first, last), an associated face, and a
  twin edge.

  all faces are rotated counterclockwise.

  the entire mesh must exist within the unit square.

  """

  def __init__(self, int nmax, double zonewidth):
    """
    initialize triangular mesh.

    - nmax is the maximal number of vertices/edges/triangles. this space is
    reserved upon instantiation

    ### TODO: dynamically allocate this space


    """

    self.nmax = nmax

    self.vnum = 0

    self.henum = 0

    self.fnum = 0

    self.zonewidth = zonewidth

    self.nz = int(1.0 /zonewidth)

    if self.nz<3:
      self.nz = 1
      self.zonewidth = 1.0

    self.zonemap = Zonemap3d(self.nz)
    self.zonemap.__assign_xyz_arrays(self.X, self.Y, self.Z)

    print('nmax: {:d}'.format(nmax))
    print('number of zones: {:d}'.format(self.nz))
    print('zonewidth: {:f}'.format(zonewidth))

    return

  def __cinit__(self, int nmax, *arg, **args):

    self.X = <double *>malloc(nmax*sizeof(double))
    double_array_init(self.X,nmax,0.)

    self.Y = <double *>malloc(nmax*sizeof(double))
    double_array_init(self.Y,nmax,0.)

    self.Z = <double *>malloc(nmax*sizeof(double))
    double_array_init(self.Z,nmax,0.)

    self.I = <double *>malloc(nmax*sizeof(double))
    double_array_init(self.I,nmax,0.)

    self.HE = <sHE *>malloc(nmax*sizeof(sHE))

    self.VHE = <int *>malloc(nmax*sizeof(int))
    int_array_init(self.VHE,nmax,-1)

    self.FHE = <int *>malloc(nmax*sizeof(int))
    int_array_init(self.FHE,nmax,-1)

    return

  def __dealloc__(self):

    free(self.X)

    free(self.Y)

    free(self.Z)

    free(self.HE)

    free(self.VHE)

    free(self.FHE)

    return

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cdef int __valid_new_vertex(self, double x, double y, double z) :#nogil:
    """
    check that x,y is within unit square
    """

    if x<0. or x>1.:

      return -1

    if y<0. or y>1.:

      return -1

    if z<0. or z>1.:

      return -1

    return 1

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cdef int __safe_vertex_positions(self, double limit) :#nogil:
    """
    check that all vertices are within limit of unit square boundary
    """

    cdef int vnum = self.vnum
    cdef int i

    for i in xrange(vnum):

      if self.X[i]<limit or self.X[i]>1.-limit:

        return -1

      if self.Y[i]<limit or self.Y[i]>1.-limit:

        return -1

      if self.Z[i]<limit or self.Z[i]>1.-limit:

        return -1

    return 1

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cdef double __get_edge_intensity(self, int he1) :#nogil:

    return (self.I[self.HE[he1].first]+self.I[self.HE[he1].last])*0.5

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cdef void __set_vertex_intensity(self, int v1, double i) :#nogil:

    self.I[v1] = i

    return

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cdef void __set_edge_intensity(self, int he1, double i) :#nogil:

    self.I[self.HE[he1].first] = i
    self.I[self.HE[he1].last] = i

    return

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cdef int __new_vertex(self, double x, double y, double z) :#nogil:
    """
    adds a vertex x,y. returns id of new vertex
    """

    if self.__valid_new_vertex(x,y,z)<0:
      return -1

    cdef int vnum = self.vnum

    self.X[vnum] = x
    self.Y[vnum] = y
    self.Z[vnum] = z

    self.zonemap.__add_vertex(vnum)

    self.vnum += 1
    return vnum

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cdef int __new_edge(self, int first, int last) :#nogil:
    """
    adds edge from vertex, first, to vertex, last. returns id of new edge
    """

    cdef int henum = self.henum

    cdef sHE he

    he.id = henum
    he.gen = -1
    he.first = first
    he.last = last
    he.face = -1
    he.twin = -1
    he.next = -1

    self.HE[henum] = he

    self.henum += 1
    return henum

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cdef int __new_edge_from_edge(self, int he1, int last) :#nogil:
    """
    appends edge to end of edge, he, to vertex, last. returns id of new edge
    """

    cdef int henum = self.henum

    cdef int first = self.HE[he1].last
    self.HE[he1].next = henum

    cdef sHE he

    he.id = henum
    he.gen = -1
    he.first = first
    he.last = last
    he.face = -1
    he.twin = -1
    he.next = -1

    self.HE[henum] = he

    self.henum += 1
    return henum

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cdef int __new_face(self, int he1) :#nogil:
    """
    creates a face and associates it with edge, he. returns id of new face
    """

    cdef int fnum = self.fnum

    self.FHE[fnum] = he1

    self.fnum += 1
    return fnum

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cdef void __set_face_of_three_edges(self, int f, int he1, int he2, int he3) :#nogil:

    self.HE[he1].face = f
    self.HE[he2].face = f
    self.HE[he3].face = f

    return

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cdef void __set_gen_of_three_edges(self, int gen, int he1, int he2, int he3) :#nogil:

    self.HE[he1].gen = gen
    self.HE[he2].gen = gen
    self.HE[he3].gen = gen

    return

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cdef void __set_edge_of_face(self, int f, int he1) :#nogil:

    self.FHE[f] = he1

    return

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cdef void __set_mutual_twins(self, int he1, int he2) :#nogil:

    self.HE[he1].twin = he2
    self.HE[he2].twin = he1

    return

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cdef int __is_surface_edge(self, int he) :#nogil:

    if self.HE[he].twin<0:
      return 1

    return -1

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cdef int __next_surface(self, int he1, int direction) :#nogil:

    cdef int he2 = -1000
    cdef int twin = -2000

    if direction>0:

      ## counter clockwise / backward surface

      he2 = self.HE[self.HE[he1].next].next
      twin = self.HE[he2].twin

      if self.__edge_integrity(he2)<0:

        return -1

      while twin>-1:
        he2 = self.HE[self.HE[twin].next].next
        twin = self.HE[he2].twin

        if self.__edge_integrity(he2)<0:

          return -1

        if twin == he1:

          return -1

      return he2

    else:

      ## clockwise / forward surface

      he2 = self.HE[he1].next
      twin = self.HE[he2].twin

      while twin>-1:
        he2 = self.HE[twin].next
        twin = self.HE[he2].twin

        if self.__edge_integrity(he2)<0:

          return -1

        if twin == he1:
          return -2

      return he2

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cdef int __edge_duplicate_test(self, int he1, int a, int b) :#nogil:
    """
    assuming he1 is connected to either vertex a or b.
    test if there is an edge ab or ba connected to he1.

    returns -1 if such an edge exists.
    otherwise it returns 1
    it returns -2 if there was another error.
    """

    if self.__edge_integrity(he1)<0:

      return -2

    cdef int he2
    cdef int twin

    ## counter-clockwise / backward

    he2 = self.HE[self.HE[he1].next].next
    twin = self.HE[he2].twin

    while twin>-1 and twin!=he1:

      he2 = self.HE[self.HE[twin].next].next
      twin = self.HE[he2].twin

      if (self.HE[he2].first == a and self.HE[he2].last == b) or\
          (self.HE[he2].last == a and self.HE[he2].first == b):

        return -1

    ## clockwise / forward

    he2 = self.HE[he1].next
    twin = self.HE[he2].twin

    while twin>-1 and twin!=he1:

      he2 = self.HE[twin].next
      twin = self.HE[he2].twin

      if (self.HE[he2].first == a and self.HE[he2].last == b) or\
          (self.HE[he2].last == a and self.HE[he2].first == b):

        return -1

    return 1

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cdef int __flip_edge(self, int he1, double limit) :#nogil:

    if self.__is_surface_edge(he1)>0:

      return -1

    cdef int the1 = self.HE[he1].twin
    cdef int f1 = self.HE[he1].face
    cdef int f2 = self.HE[the1].face

    cdef int bc1 = self.HE[he1].next
    cdef int ad1 = self.HE[the1].next
    cdef int db1 = self.HE[ad1].next
    cdef int ca1 = self.HE[bc1].next

    cdef int a = self.HE[he1].first
    cdef int b = self.HE[the1].first
    cdef int c = self.HE[ca1].first
    cdef int d = self.HE[db1].first

    cdef double ablen = sqrt(cpow(self.X[a]-self.X[b],2)+
                            cpow(self.Y[a]-self.Y[b],2)+
                            cpow(self.Z[a]-self.Z[b],2))

    cdef double dclen = sqrt(cpow(self.X[d]-self.X[c],2)+
                            cpow(self.Y[d]-self.Y[c],2)+
                            cpow(self.Z[d]-self.Z[c],2))

    if self.__edge_duplicate_test(db1, c,d)!=1:

      return -1

    if self.__edge_duplicate_test(ca1, c,d)!=1:

      return -1

    if ablen<dclen:

      return -1

    if dclen<limit or ablen<limit:

      return -1

    self.__set_mutual_twins(he1, the1)

    self.__set_next_of_triangle(he1,db1,bc1)
    self.__set_next_of_triangle(the1,ca1,ad1)

    self.__set_face_of_three_edges(f1, he1, db1, bc1)
    self.__set_face_of_three_edges(f2, the1, ca1, ad1)

    cdef int gen = self.HE[he1].gen
    if gen>self.HE[the1].gen:
      gen = self.HE[the1].gen

    self.__set_gen_of_three_edges(gen, he1, db1, bc1)
    self.__set_gen_of_three_edges(gen, the1, ca1, ad1)

    self.__set_edge_of_face(f1,he1) #
    self.__set_edge_of_face(f2,the1)

    self.HE[he1].first = c
    self.HE[he1].last = d
    self.HE[the1].first = d
    self.HE[the1].last = c

    return 1

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cdef int __set_next_of_triangle(self, int he1, int he2, int he3) :#nogil:

    self.HE[he1].next = he2
    self.HE[he2].next = he3
    self.HE[he3].next = he1

    return 1

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cdef int __split_edge(self, int he1) :#nogil:

    if self.__is_surface_edge(he1)>0:

      return self.__split_surface_edge(he1)

    else:

      return self.__split_internal_edge(he1)

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cdef int __split_internal_edge(self, int he1) :#nogil:

    cdef sHE he = self.HE[he1]
    cdef sHE the = self.HE[he.twin]

    cdef int the1 = he.twin
    cdef int bc1 = he.next
    cdef sHE bc = self.HE[bc1]

    cdef int ad1 = the.next
    cdef sHE ad = self.HE[ad1]

    cdef int db1 = ad.next
    cdef int ca1 = bc.next

    cdef int f1 = he.face
    cdef int f2 = the.face

    cdef int a = he.first
    cdef int b = he.last
    cdef int c = bc.last
    cdef int d = ad.last

    if a == b or a == c or a == d or b == c or b == d or c == d:

      return -1

    cdef double x = (self.X[a] + self.X[b])*0.5
    cdef double y = (self.Y[a] + self.Y[b])*0.5
    cdef double z = (self.Z[a] + self.Z[b])*0.5

    cdef int e = self.__new_vertex(x,y,z)

    cdef int de1 = self.__new_edge(d,e)
    cdef int ec1 = self.__new_edge(e,c)
    cdef int ed1 = self.__new_edge(e,d)
    cdef int ce1 = self.__new_edge(c,e)

    cdef int be1 = self.__new_edge(b,e)
    cdef int eb1 = self.__new_edge(e,b)

    self.HE[the1].first = e
    self.HE[he1].last = e

    self.__set_mutual_twins(de1,ed1)
    self.__set_mutual_twins(be1,eb1)
    self.__set_mutual_twins(ce1,ec1)
    self.__set_mutual_twins(he1,the1) #

    self.__set_next_of_triangle(ad1, de1, the1)
    self.__set_next_of_triangle(be1, ed1, db1)
    self.__set_next_of_triangle(eb1, bc1, ce1)
    self.__set_next_of_triangle(ec1, ca1, he1)

    cdef int f3 = self.__new_face(bc1)
    cdef int f4 = self.__new_face(db1)

    self.__set_edge_of_face(f1, he1)
    self.__set_edge_of_face(f2, the1)

    self.__set_face_of_three_edges(f1, ca1, he1, ec1)
    self.__set_face_of_three_edges(f2, ad1, de1, the1)
    self.__set_face_of_three_edges(f3, eb1, ce1, bc1)
    self.__set_face_of_three_edges(f4, ed1, db1, be1)

    self.__set_gen_of_three_edges(he.gen, ca1, he1, ec1)
    self.__set_gen_of_three_edges(he.gen, ad1, de1, the1)
    self.__set_gen_of_three_edges(he.gen, eb1, ce1, bc1)
    self.__set_gen_of_three_edges(he.gen, ed1, db1, be1)

    return 1

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cdef int __split_surface_edge(self, int he1) :#nogil:

    cdef sHE he = self.HE[he1]

    cdef int bc1 = he.next
    cdef sHE bc = self.HE[bc1]

    cdef int ca1 = bc.next
    cdef sHE ca = self.HE[ca1]

    cdef int f1 = he.face

    cdef int a = he.first
    cdef int b = he.last
    cdef int c = bc.last

    if a == b or a == c or b == c:

      return -1

    cdef double x = (self.X[a] + self.X[b])*0.5
    cdef double y = (self.Y[a] + self.Y[b])*0.5
    cdef double z = (self.Z[a] + self.Z[b])*0.5

    cdef int e = self.__new_vertex(x,y,z)

    cdef int ec1 = self.__new_edge(e,c)
    cdef int ce1 = self.__new_edge(c,e)

    cdef int eb1 = self.__new_edge(e,b)

    self.HE[he1].last = e

    self.__set_vertex_intensity(e, (self.I[he.first] + self.I[he.last])*0.5)

    self.__set_mutual_twins(ce1,ec1)

    self.__set_next_of_triangle(eb1, bc1, ce1)
    self.__set_next_of_triangle(ec1, ca1, he1)

    cdef int f3 = self.__new_face(bc1)

    self.__set_edge_of_face(f1, he1)

    self.__set_face_of_three_edges(f1, ca1, he1, ec1)
    self.__set_face_of_three_edges(f3, eb1, ce1, bc1)

    self.__set_gen_of_three_edges(he.gen, ca1, he1, ec1)
    self.__set_gen_of_three_edges(he.gen, eb1, ce1, bc1)

    return 1

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cdef int __split_all_longest_triangle_edges(self, double limit) :#nogil:

    cdef int fnum = self.fnum
    cdef int f

    cdef int he1
    cdef int he2
    cdef int he3

    cdef double l1
    cdef double l2
    cdef double l3

    for f in reversed(xrange(fnum)):

      he1 = self.FHE[f]
      he2 = self.HE[he1].next
      he3 = self.HE[he2].next

      l1 = self.__get_edge_length(he1)

      if l1<=0.0:

        continue

      l2 = self.__get_edge_length(he2)

      if l2<=0.0:

        continue

      l3 = self.__get_edge_length(he3)

      if l3<=0.0:

        continue

      if l1>limit and l1>l2 and l1>l3:

        self.__split_edge(he1)

      elif l2>limit and l2>l1 and l2>l3:

        self.__split_edge(he2)

      elif l3>limit and l3>l1 and l3>l2:

        self.__split_edge(he3)

    return 1

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cdef double __get_edge_length(self, int he) :#nogil:

    cdef int first = self.HE[he].first
    cdef int last = self.HE[he].last
    cdef double dx = self.X[first] - self.X[last]
    cdef double dy = self.Y[first] - self.Y[last]
    cdef double dz = self.Z[first] - self.Z[last]

    return sqrt(dx*dx+dy*dy+dz*dz)

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cdef int __edge_integrity(self, int he1) :#nogil:

    if he1<0:

      return -100

    cdef sHE he = self.HE[he1]
    cdef sHE the
    cdef int thenext

    if he.next == he1:

      return -1

    if he.last == he.first:

      return -2

    if he.twin>-1:

      if he.twin == he1:

        return -4

      the = self.HE[he.twin]
      thenext = the.next

      if self.HE[thenext].last == self.HE[he.next].last:

        return -3

      if the.first == he.first or the.last == he.last:

        return -5

    return 1

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cdef int __triangle_integrity(self, int face) :#nogil:
    """
    1: triangle is ok
    -2: edges are not linked
    -3: duplicate vertex in edges
    """

    cdef int he1 = self.FHE[face]
    cdef int he2 = self.HE[he1].next
    cdef int he3 = self.HE[he2].next

    if self.HE[he1].last == self.HE[he1].first or\
      self.HE[he2].last == self.HE[he2].first or\
      self.HE[he3].last == self.HE[he3].first:

      return -3

    if self.HE[he1].last != self.HE[he2].first or\
      self.HE[he2].last != self.HE[he3].first or\
      self.HE[he3].last != self.HE[he1].first:

      return -2

    return 1

  ##@cython.wraparound(False)
  ##@cython.boundscheck(False)
  ##@cython.nonecheck(False)
  cpdef int initiate_faces(self, list vertices, list faces):

    cdef double x
    cdef double y
    cdef double z

    cdef int he1
    cdef int he2
    cdef int he3

    cdef int v1
    cdef int v2
    cdef int v3

    cdef tuple k1
    cdef tuple k2
    cdef tuple k3

    cdef int f1

    cdef dict facemap = {}
    facemap = {}

    cdef list vv

    i = 0

    for x,y,z in vertices:

      self.__new_vertex(x,y,z)

    for v1,v2,v3 in faces:

      he1 = self.__new_edge(v1,v2)
      he2 = self.__new_edge(v2,v3)
      he3 = self.__new_edge(v3,v1)

      self.__set_next_of_triangle(he1,he2,he3)
      f1 = self.__new_face(he1)
      self.__set_face_of_three_edges(f1,he1,he2,he3)
      self.__set_gen_of_three_edges(0,he1,he2,he3)

      k1 = tuple(sorted([v1, v2]))
      k2 = tuple(sorted([v2, v3]))
      k3 = tuple(sorted([v3, v1]))

      dict_list_add(facemap, k1, he1)
      dict_list_add(facemap, k2, he2)
      dict_list_add(facemap, k3, he3)

    for k1,vv in facemap.iteritems():

      if len(vv)>1:

        he1 = <int>vv[0]
        he2 = <int>vv[1]

        self.__set_mutual_twins(he1,he2)

    return 1

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cpdef int edge_integrity(self, int face):

    return self.__edge_integrity(face)

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cpdef int triangle_integrity(self, int face):

    return self.__triangle_integrity(face)

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cpdef int safe_vertex_positions(self, double limit):

    return self.__safe_vertex_positions(limit)

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cpdef list new_faces_in_ngon(self, double x1, double y1, double z1, double rad, int num):

    cdef list vertices = []
    cdef list edges = []
    cdef list outside_edges = []

    cdef double the = 0.0
    cdef double thediff = TWOPI/num
    cdef double x
    cdef double y
    cdef double z
    cdef int f
    cdef int i
    cdef int first
    cdef int last
    cdef int first_edge
    cdef int outside_edge
    cdef int edge

    cdef int vmid = self.__new_vertex(x1,y1,z1)

    for i in xrange(num):

      x = x1 + cos(the)*rad
      y = y1 + sin(the)*rad
      z = z1
      vertices.append(self.__new_vertex(x,y,z))
      the += thediff

    for i in xrange(num):

      first = i
      last = i+1

      if i>=num-1:
        first = i
        last = 0

      first_edge = self.__new_edge(vmid, vertices[first])
      outside_edge = self.__new_edge(vertices[first], vertices[last])
      edge = self.__new_edge(vertices[last], vmid)

      edges.append(first_edge)
      edges.append(edge)
      outside_edges.append(outside_edge)

      self.__set_next_of_triangle(first_edge, outside_edge, edge)

      f = self.__new_face(first_edge)
      self.__set_face_of_three_edges(f, first_edge, edge, outside_edge)
      self.__set_edge_of_face(f, first_edge)
      self.__set_gen_of_three_edges(0, first_edge, edge, outside_edge)

    for i in xrange(num-1):

      self.__set_mutual_twins(edges[2*i+1],edges[2*i+2])

    self.__set_mutual_twins(edges[0],edges[2*num-1])

    return outside_edges

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cpdef int next_surface(self, int he, int direction):

    return self.__next_surface(he, direction)

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cpdef int flip_edge(self, int he1, double limit):

    return self.__flip_edge(he1, limit)

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cpdef int split_edge(self, int he):

    return self.__split_edge(he)

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cpdef int optimize_edges(self, double split_limit, double flip_limit):

    cdef int he1

    ## TODO: maybe try to do this in one run?

    for he1 in reversed(xrange(self.henum)):

      self.__flip_edge(he1, flip_limit)

    self.__split_all_longest_triangle_edges(split_limit)

    return self.henum

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cpdef int np_get_vertices(self, np.ndarray[double, mode="c",ndim=2] a):
    """
    """

    cdef int v
    cdef int vnum = self.vnum

    for v in xrange(vnum):

      a[v, 0] = self.X[v]
      a[v, 1] = self.Y[v]
      a[v, 2] = self.Z[v]

    return vnum

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cpdef int np_get_triangles_vertices(self, np.ndarray[long, mode="c",ndim=2] a):
    """
    """

    cdef int fnum = self.fnum
    cdef int f
    cdef int next
    cdef int* vertices = [-1,-1,-1]

    cdef sHE he

    for f in xrange(fnum):

      he = self.HE[self.FHE[f]]
      next = he.next
      vertices[0] = he.first

      he = self.HE[next]
      next = he.next
      vertices[1] = he.first

      he = self.HE[next]
      vertices[2] = he.first

      a[f, 0] = vertices[0]
      a[f, 1] = vertices[1]
      a[f, 2] = vertices[2]

    return fnum

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cpdef int np_get_triangles_gen(self, np.ndarray[long, mode="c",ndim=1] a):
    """
    assigns triangle generations to a
    returns number of triangles
    """

    cdef int f
    for f in xrange(self.fnum):

      a[f] = self.HE[self.FHE[f]].gen

    return self.fnum

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cpdef double get_edge_intensity(self, int he1):

    return self.__get_edge_intensity(he1)

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cpdef int set_vertex_intensity(self, int v1, double i):

    self.__set_vertex_intensity(v1, i)

    return 1

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cpdef int set_edge_intensity(self, int he1, double i):

    self.__set_edge_intensity(he1, i)

    return 1

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cpdef double get_triangle_intensity(self, int f1):

    cdef int he1 = self.FHE[f1]
    cdef int he2 = self.HE[he1].next

    return (self.I[self.HE[he1].first] +
      self.I[self.HE[he1].last] +
      self.I[self.HE[he2].last])/3.0

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cpdef list get_triangle_edges(self, int f1):
    """
    """

    cdef int e1 = self.FHE[f1]
    cdef int e2 = self.HE[e1].next
    cdef int e3 = self.HE[e2].next

    return [e1,e2,e3]

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cpdef int diminish_all_vertex_intensity(self, double d):

    cdef int v

    for v in xrange(self.vnum):

      self.I[v] *= d

    return 1

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cpdef int is_surface_edge(self, int t1):
    """
    returns 1 if edge is on the surface of the mesh. that is, the edge has no
    twin edge
    """

    return self.__is_surface_edge(t1)

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cpdef double get_edge_length(self, int he1):

    return self.__get_edge_length(he1)

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cpdef dict get_edge_dict(self, int he1):
    """
    get some debug info about the edge
    """

    return {
      'id': self.HE[he1].id,
      'gen': self.HE[he1].gen,
      'first': self.HE[he1].first,
      'last': self.HE[he1].last,
      'twin': self.HE[he1].twin,
      'face': self.HE[he1].face,
      'next': self.HE[he1].next
    }

  #@cython.wraparound(False)
  #@cython.boundscheck(False)
  #@cython.nonecheck(False)
  cpdef list get_triangle_dicts(self, int f):
    """
    get some debug info about the triangle
    """

    cdef int he1 = self.FHE[f]
    cdef int he2 = self.HE[he1].next
    cdef int he3 = self.HE[he2].next

    return [
      (he1, self.get_edge_dict(he1)),
      (he2, self.get_edge_dict(he2)),
      (he3, self.get_edge_dict(he3)),
    ]

  #@cython.nonecheck(False)
  cpdef int get_vnum(self):
    """
    number of vertices
    """
    return self.vnum

  #@cython.nonecheck(False)
  cpdef int get_henum(self):
    """
    number of (half) edges
    """
    return self.henum

  #@cython.nonecheck(False)
  cpdef int get_fnum(self):
    """
    number of faces
    """
    return self.fnum

