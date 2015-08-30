# -*- coding: utf-8 -*-
# cython: profile=True

from __future__ import division

cimport cython
cimport mesh3d

from zonemap3d cimport Zonemap3d

from cython.parallel import parallel, prange

from libc.math cimport sqrt
from libc.math cimport cos
from libc.math cimport sin

from helpers cimport float_array_init
from helpers cimport int_array_init
from helpers cimport vcross

import numpy as np
cimport numpy as np


cdef class DifferentialMesh3d(mesh3d.Mesh3d):

  def __init__(self, int nmax, float zonewidth, float nearl, float farl):

    mesh3d.Mesh3d.__init__(self, nmax, zonewidth)

    """
    - nearl is the closest comfortable distance between two vertices.

    - farl is the distance beyond which disconnected vertices will ignore
    each other
    """

    self.nearl = nearl

    self.farl = farl

    self.num_sources = 0

    #self.source_zonemap = Zonemap3d(self.nz)
    #self.source_zonemap.__assign_xyz_arrays(self.SX, self.SY, self.SZ)

    print('nearl: {:f}'.format(nearl))
    print('farl: {:f}'.format(farl))

    return

  def __cinit__(self, int nmax, *arg, **args):

    self.DX = <float *>malloc(nmax*sizeof(float))

    self.DY = <float *>malloc(nmax*sizeof(float))

    self.DZ = <float *>malloc(nmax*sizeof(float))

    self.SX = <float *>malloc(nmax*sizeof(float))

    self.SY = <float *>malloc(nmax*sizeof(float))

    self.SZ = <float *>malloc(nmax*sizeof(float))

    return

  def __dealloc__(self):

    free(self.DX)

    free(self.DY)

    free(self.DZ)

    free(self.SX)

    free(self.SY)

    free(self.SZ)

    return

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  @cython.cdivision(True)
  cdef int __reject(self, float scale) nogil:
    """
    all vertices will move away from all neighboring (closer than farl)
    vertices
    """

    cdef float farl = self.farl
    cdef int vnum = self.vnum

    cdef int v
    cdef int k
    cdef int neigh

    cdef float x
    cdef float y
    cdef float z
    cdef float dx
    cdef float dy
    cdef float dz
    cdef float nrm
    cdef float force

    cdef float resx = 0.
    cdef float resy = 0.
    cdef float resz = 0.

    cdef int asize = self.zonemap.__get_greatest_zone_size()*9
    cdef int *vertices
    cdef int neighbor_num

    vertices = <int *>malloc(asize*sizeof(int))

    for v in xrange(vnum):

      x = self.X[v]
      y = self.Y[v]
      z = self.Z[v]

      neighbor_num = self.zonemap.__sphere_vertices(x, y, z, farl, vertices)

      resx = 0.
      resy = 0.
      resz = 0.

      for k in range(neighbor_num):

        neigh = vertices[k]

        if neigh == v:

          continue

        dx = x-self.X[neigh]
        dy = y-self.Y[neigh]
        dz = z-self.Z[neigh]
        nrm = sqrt(dx*dx+dy*dy+dz*dz)

        if nrm>farl or nrm<=0.0:

          continue

        force = (farl-nrm)/nrm

        resx += dx*force
        resy += dy*force
        resz += dz*force

      self.DX[v] += resx
      self.DY[v] += resy
      self.DZ[v] += resz

    free(vertices)

    return 1


  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  @cython.cdivision(True)
  cdef int __attract(self, float scale) nogil:
    """
    vertices will move towards all connected vertices further away than
    nearl
    """

    cdef int v1
    cdef int v2
    cdef int k

    cdef float nearl = self.nearl

    cdef float dx
    cdef float dy
    cdef float dz
    cdef float nrm

    cdef float s

    for k in xrange(self.henum):

      v1 = self.HE[k].first
      v2 = self.HE[k].last

      dx = self.X[v2]-self.X[v1]
      dy = self.Y[v2]-self.Y[v1]
      dz = self.Z[v2]-self.Z[v1]
      nrm = sqrt(dx*dx+dy*dy+dz*dz)

      if self.HE[k].twin>-1:

        # internal edge. has two opposing half edges
        # half the force because it is applied twice
        s = scale*0.5

        if nrm>nearl and nrm>0.0:

          ## attract
          self.DX[v1] += dx/nrm*s
          self.DY[v1] += dy/nrm*s
          self.DZ[v1] += dz/nrm*s

        elif nrm<nearl*0.5 and nrm>0.0:

          ## reject
          self.DX[v1] -= dx/nrm*s
          self.DY[v1] -= dy/nrm*s
          self.DZ[v1] -= dz/nrm*s

      else:

        # surface edge has one half edge, and they are all rotated the same way
        # half the force because we apply it twice.
        s = scale*0.5

        if nrm>nearl and nrm>0.0:

          ## attract
          self.DX[v1] += dx/nrm*s
          self.DY[v1] += dy/nrm*s
          self.DZ[v1] += dz/nrm*s
          # and the other vertex
          self.DX[v2] -= dx/nrm*s
          self.DY[v2] -= dy/nrm*s
          self.DZ[v2] -= dz/nrm*s

        elif nrm<nearl*0.5 and nrm>0.0:

          ## reject

          self.DX[v1] -= dx/nrm*s
          self.DY[v1] -= dy/nrm*s
          self.DZ[v1] -= dz/nrm*s
          # and the other vertex
          self.DX[v2] += dx/nrm*s
          self.DY[v2] += dy/nrm*s
          self.DZ[v2] += dz/nrm*s

    return 1

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cpdef int optimize_position(self, float step, int itt=1):

    cdef int v
    cdef int i

    cdef float reject_scale = 1.0
    cdef float scale = 0.1

    for i in xrange(itt):

      float_array_init(self.DX, self.vnum, 0.)
      float_array_init(self.DY, self.vnum, 0.)
      float_array_init(self.DZ, self.vnum, 0.)

      self.__reject(reject_scale)
      self.__attract(scale)

      for v in range(self.vnum):

        self.X[v] += self.DX[v]*step
        self.Y[v] += self.DY[v]*step
        self.Z[v] += self.DZ[v]*step

    return 1

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cpdef int position_noise(self, np.ndarray[double, mode="c",ndim=2] a):

    cdef int v

    for v in xrange(self.vnum):

      self.X[v] += a[v,0]
      self.Y[v] += a[v,1]
      self.Z[v] += a[v,2]

    return 1

