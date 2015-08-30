#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

#from numpy import pi
from numpy import sqrt
#from numpy import zeros
#from numpy import cos
#from numpy import sin
#from numpy import power
#from numpy import floor
#from numpy.random import random

NMAX = 10e7

STP = 0.000001

ITT = 100000
OPT_ITT = 1

NEARL = 0.025
H = 0.0273
FARL = 0.1

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

  from numpy.random import random

  fn_in = './res/base.json'
  fn_out = './res/res.json'

  DM = DifferentialMesh3d(NMAX, 2*FARL, NEARL, FARL)

  load(DM, fn_in)

  for i in xrange(ITT):

    try:

      t1 = time()

      DM.optimize_position(STP, OPT_ITT)

      vnum = DM.get_vnum()

      for v in xrange(vnum):

        if DM.is_surface_edge(v)>0:
          print('surface ',v)
          raise KeyboardInterrupt

      noise = (1.0-2*random(size=(vnum,3)))*STP*0.1
      DM.position_noise(noise)

      DM.optimize_edges(H, 10000)

      print_stats(i, time()-t1, DM)

    except KeyboardInterrupt:

      break

  export(DM, fn_out)


if __name__ == '__main__' :

    main()

