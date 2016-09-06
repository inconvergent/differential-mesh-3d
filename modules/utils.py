# -*- coding: utf-8 -*-


def get_surface_vertices(dm):

  res = []

  for he in range(dm.get_henum()):
    e = dm.is_surface_edge(he)
    if e>0:
      d = dm.get_edge_dict(he)
      res.append(d['first'])
      res.append(d['last'])

  return list(set(res))

def get_seed_selector(dm, t, sr):
  from numpy import array
  from numpy import arange
  from numpy.random import random

  if t == 'surface':

    def f():
      vertices = array(get_surface_vertices(dm))
      rm = (random(size=len(vertices))<sr).nonzero()[0]
      if len(rm)<1:
        return array([])
      return vertices[rm]

  elif t == 'random':

    def f():
      vn = dm.get_vnum()
      vertices = arange(vn)
      rm = (random(size=vn)<sr).nonzero()[0]
      if len(rm)<1:
        return array([])
      return vertices[rm]

  return f

