# -*- coding: utf-8 -*-

def make_info_str(args):
  s = ''
  for k in vars(args):
    s += '# ' + str(k) + ': ' + str(getattr(args,k)) + '\n'
  return s


def print_stats(steps,dm, meta=False):
  from time import strftime
  from time import time

  if isinstance(meta, str):
    meta = ' | {:s}'.format(meta)
  else:
    meta = ''

  print(
      '{:s} | stp: {:d} sec: {:.2f} v: {:d} e: {:d} f: {:d}{:s}'
      .format(
          strftime('%d/%m/%y %H:%M:%S'),
          steps,
          time()-dm.get_start_time(),
          dm.get_vnum(),
          dm.get_henum(),
          dm.get_fnum(),
          meta
          )
      )

  return

def get_exporter(dm, fn, nmax):
  from numpy import zeros
  from .geometry import move_scale
  from iutils.ioOBJ import export

  np_verts = zeros((nmax, 3), 'float')
  np_tris = zeros((nmax, 3), 'int')
  np_int = zeros(nmax, 'float')

  def e():
    vnum = dm.np_get_vertices(np_verts)
    tnum = dm.np_get_triangles_vertices(np_tris)
    dm.np_get_triangles_intensity(np_int)
    move_scale(np_verts[:vnum, :], s=1000)
    export(
        'thing_mesh',
        fn.name(),
        verts=np_verts[:vnum, :],
        tris=np_tris[:tnum, :]
        )
  return e


def get_surface_vertices(dm):
  res = []

  for he in range(dm.get_henum()):
    e = dm.is_surface_edge(he)
    if e > 0:
      d = dm.get_edge_dict(he)
      res.append(d['first'])
      res.append(d['last'])

  return list(set(res))

def get_seed_selector(dm, t, sr):
  from numpy import array
  from numpy import arange
  from numpy import ones
  from numpy.random import random

  if sr >= 1:
    get_mask = lambda n, sr: ones(n, 'bool')
  else:
    get_mask = lambda n, sr: (random(size=n) < sr).nonzero()[0]

  if t == 'surface':
    def f():
      vertices = array(get_surface_vertices(dm))
      rm = get_mask(len(vertices), sr)
      if len(rm) < 1:
        return array([])
      return vertices[rm]

  elif t == 'random':
    def f():
      vn = dm.get_vnum()
      vertices = arange(vn)
      rm = get_mask(len(vertices), sr)
      if len(rm) < 1:
        return array([])
      return vertices[rm]
  else:
    raise ValueError('use "surface" or "random".')

  return f

