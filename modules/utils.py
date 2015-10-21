#!/usr/bin/python
# -*- coding: utf-8 -*-

def random_unit_vec(num, scale):

  from numpy.random import normal
  from numpy.linalg import norm
  from numpy import reshape

  rnd = normal(size=(num,3))
  d = norm(rnd,axis=1)
  rnd[:] /= reshape(d, (num,1))

  return rnd*scale

def load_obj(
  fn,
  sx=[1.0,1.0,1.0],
  mx=[0,5,0.5,0.5]
):

  from codecs import open
  from numpy import row_stack

  vertices = []
  faces = []

  with open(fn, 'r', encoding='utf8') as f:

    for l in f:
      if l.startswith('#'):
        continue

      values = l.split()
      if not values:
        continue
      if values[0] == 'v':
        vertices.append([float(v) for v in values[1:]])

      if values[0] == 'f':
        face = [int(v.split('//')[0])-1 for v in values[1:]]
        faces.append(face)

  np_vertices = row_stack(vertices)

  xmax = np_vertices[:,0].max()
  xmin = np_vertices[:,0].min()
  ymax = np_vertices[:,1].max()
  ymin = np_vertices[:,1].min()
  zmax = np_vertices[:,2].max()
  zmin = np_vertices[:,2].min()
  dx = xmax - xmin
  dy = ymax - ymin
  dz = zmax - zmin

  print('original')
  print('x min max, {:0.8f} {:0.8f}, dst: {:0.8f}'.format(xmin,xmax,dx))
  print('y min max, {:0.8f} {:0.8f}, dst: {:0.8f}'.format(ymin,ymax,dy))
  print('z min max, {:0.8f} {:0.8f}, dst: {:0.8f}'.format(zmin,zmax,dz))

  np_vertices /= max([dx,dy,dz])

  np_vertices[:,0] *= sx[0]
  np_vertices[:,1] *= sx[1]
  np_vertices[:,2] *= sx[2]

  np_vertices[:,0] += mx[0]
  np_vertices[:,1] += mx[1]
  np_vertices[:,2] += mx[2]

  xmax = np_vertices[:,0].max()
  xmin = np_vertices[:,0].min()
  ymax = np_vertices[:,1].max()
  ymin = np_vertices[:,1].min()
  zmax = np_vertices[:,2].max()
  zmin = np_vertices[:,2].min()
  dx = xmax - xmin
  dy = ymax - ymin
  dz = zmax - zmin

  print('rescaled')
  print('x min max, {:0.8f} {:0.8f}, dst: {:0.8f}'.format(xmin,xmax,dx))
  print('y min max, {:0.8f} {:0.8f}, dst: {:0.8f}'.format(ymin,ymax,dy))
  print('z min max, {:0.8f} {:0.8f}, dst: {:0.8f}'.format(zmin,zmax,dz))

  return {
    'faces': faces,
    'vertices': [list(row) for row in np_vertices]
  }

def export_obj(dm, obj_name, fn, write_intensity=False, meta=False):

  from numpy import zeros
  from codecs import open
  from time import time

  vnum = dm.get_vnum()
  fnum = dm.get_fnum()
  henum = dm.get_henum()
  np_verts = zeros((vnum,3),'float')
  np_tris = zeros((fnum,3),'int')

  runtime = time()-dm.get_start_time()

  dm.np_get_vertices(np_verts)
  dm.np_get_triangles_vertices(np_tris)

  intensity = None

  if write_intensity:
    intensity = zeros(vnum,'double')
    dm.get_vertices_intensity(intensity)

  print('storing mesh ...')
  print('num vertices: {:d}, num triangles: {:d}'.format(vnum, fnum))

  with open(fn, 'wb', encoding='utf8') as f:

    if meta:
      f.write('# meta:\n')
      f.write(meta+'\n')

    f.write('# info:\n')

    f.write('# vnum: {:d}\n# henum: {:d}\n# fnum: {:d}\n# runtime: {:f}\n\n'
      .format(vnum, fnum, henum, runtime))

    f.write('o {:s}\n'.format(obj_name))

    for v in np_verts[:vnum,:]:
      f.write('v {:f} {:f} {:f}\n'.format(*v))

    f.write('s off\n')

    for t in np_tris[:fnum,:]:
      t += 1
      f.write('f {:d} {:d} {:d}\n'.format(*t))

  if write_intensity:

    with open(fn+'.x', 'wb', encoding='utf8') as f:

      f.write('o {:s}\n'.format(obj_name))

      for i in intensity[:vnum]:
        f.write('c {:f} {:f} {:f}\n'.format(*[i]*3))

    print('done.')

def get_surface_vertices(dm):

  res = []

  for he in xrange(dm.get_henum()):
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

