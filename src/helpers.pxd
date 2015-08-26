# -*- coding: utf-8 -*-

cimport cython

from libc.math cimport sqrt
from libc.math cimport pow

cdef extern from "stdlib.h":
  void qsort(void *base, int nmemb, int size,
       int(*compar)(const void *, const void *)) nogil

cdef inline void int_array_init(int *a, int n, int v):
  """
  initialize integer array a of length n with integer value v
  """
  cdef int i
  for i in xrange(n):
    a[i] = v
  return

cdef inline void float_array_init(float *a, int n, float v):
  """
  initialize float array a of length n with float value v
  """
  cdef int i
  for i in xrange(n):
    a[i] = v
  return

#cdef inline int compare(const void *aa, const void *bb):
  #"""
  #compare function used with qsort
  #"""

  #cdef int a = (<int *>(aa))[0]
  #cdef int b = (<int *>(bb))[0]

  #if a < b:
    #return -1
  #elif a > b:
    #return 1
  #else:
    #return 0

cdef inline float xcross(float x1, float y1,
                         float x2, float y2) nogil:
  """
  cross [x1,y1] and [x2,y2]
  """
  
  return x1*y2-y1*x2

cdef inline float vcross(float x1, float y1,
                         float x2, float y2, float x3, float y3) nogil:
  """
  cross of v12 and v13
  """
  
  cdef float xa = x2-x1
  cdef float ya = y2-y1
  cdef float xb = x3-x1
  cdef float yb = y3-y1
  return xa*yb-ya*xb

cdef inline float cross(float x1, float y1, float x2, float y2, 
                        float x3, float y3, float x4, float y4) nogil:
  """
  cross of v12 and v34
  """

  cdef float xa = x2-x1
  cdef float ya = y2-y1
  cdef float xb = x4-x3
  cdef float yb = y4-y3
  return xa*yb-ya*xb
