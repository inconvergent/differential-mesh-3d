#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function


NMAX = 10e7
ITT = 150000
OPT_ITT = 1

STP = 1.0e-5
NEARL = 0.0068
H = NEARL*1.5

FARL = 0.03

EXPORT_ITT = 1000


def random_unit_vec(num, scale):

  from numpy.random import normal
  from numpy.linalg import norm
  from numpy import reshape

  rnd = normal(size=(num,3))
  d = norm(rnd,axis=1)
  rnd[:] /= reshape(d, (num,1))

  return rnd*scale

def load(dm,fn):

  from json import load
  from codecs import open

  with open(fn, 'rb', encoding='utf8') as f:

    data = load(f)

  vertices = data['vertices']
  faces = data['triangles']

  dm.initiate_faces(vertices, faces)

def export(dm,fn):

  from numpy import zeros
  from json import dump
  from codecs import open

  np_verts = zeros((NMAX,3),'float')
  np_tris = zeros((NMAX,3),'int')

  vnum = dm.np_get_vertices(np_verts)
  tnum = dm.np_get_triangles_vertices(np_tris)

  data = {
    'vertices': list([list(row) for row in np_verts[:vnum,:]]),
    'triangles': list([list(row) for row in np_tris[:tnum,:]])
  }

  print('storing mesh ...')
  print('num vertices: {:d}, num triangles: {:d}'.format(vnum, tnum))

  with open(fn, 'wb', encoding='utf8') as f:

    dump(data, f)

    print('done.')


def main():

  from differentialMesh3d import DifferentialMesh3d
  from time import time
  from modules.helpers import print_stats

  fn_in = './data/base.json'
  fn_out = './res/res'

  DM = DifferentialMesh3d(NMAX, 2*FARL, NEARL, FARL)

  load(DM, fn_in)

  for i in xrange(ITT):

    try:

      t1 = time()

      DM.optimize_position(STP, OPT_ITT)

      vnum = DM.get_vnum()

      noise = random_unit_vec(vnum, 4.3*STP)
      DM.position_noise(noise)

      DM.optimize_edges(H, STP)

      print_stats(i, time()-t1, DM)

      if i%EXPORT_ITT==0:

        ## TODO: write obj directly. it's not that hard.
        fn = '{:s}_{:06d}.json'.format(fn_out, i)
        export(DM, fn)
        print('wrote: ' + fn)

    except KeyboardInterrupt:

      break

  export(DM, fn_out+'_final.json')


if __name__ == '__main__' :

    main()

