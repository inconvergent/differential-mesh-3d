# -*- coding: utf-8 -*-

cimport cython

#from libc.math cimport sqrt
#from libc.math cimport pow
#from libc.math cimport fabs

#cdef extern from "stdlib.h":
  #void qsort(void *base, long nmemb, long size,
       #<long>(*compar)(const void *, const void *)) nogil

cdef inline void long_array_init(long *a, long n, long v) nogil:
  cdef long i
  for i in xrange(n):
    a[i] = v
  return

cdef inline void float_array_init(float *a, long n, float v) nogil:
  cdef long i
  for i in xrange(n):
    a[i] = v
  return

cdef inline void double_array_init(double *a, long n, double v) nogil:
  cdef long i
  for i in xrange(n):
    a[i] = v
  return

#cdef inline long positive_nearly_zero(double a):

  #if a<=1.0e-9:
    #return 1
  #else:
    #return -1

#cdef inline long compare(const void *aa, const void *bb):
  #"""
  #compare function used with qsort
  #"""

  #cdef long a = (<long *>(aa))[0]
  #cdef long b = (<long *>(bb))[0]

  #if a < b:
    #return -1
  #elif a > b:
    #return 1
  #else:
    #return 0

#cdef inline double xcross(double x1, double y1,
                         #double x2, double y2) nogil:
  #"""
  #cross [x1,y1] and [x2,y2]
  #"""

  #return x1*y2-y1*x2

#cdef inline double vcross(double x1, double y1,
                         #double x2, double y2, double x3, double y3) nogil:
  #"""
  #cross of v12 and v13
  #"""

  #cdef double xa = x2-x1
  #cdef double ya = y2-y1
  #cdef double xb = x3-x1
  #cdef double yb = y3-y1
  #return xa*yb-ya*xb

#cdef inline double cross(double x1, double y1, double x2, double y2,
                        #double x3, double y3, double x4, double y4) nogil:
  #"""
  #cross of v12 and v34
  #"""

  #cdef double xa = x2-x1
  #cdef double ya = y2-y1
  #cdef double xb = x4-x3
  #cdef double yb = y4-y3
  #return xa*yb-ya*xb

