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
SIZE = 10000
ONE = 1./SIZE

RAD = 3*ONE
H = sqrt(3.)*RAD
NEARL = 2*RAD
FARL = RAD*20

OPT_STP = ONE
OPT_ITT = 1

MID = 0.5

LINEWIDTH = 1*ONE

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

  fn_in = './res/sphere.json'
  fn_out = './res/dat.json'

  DM = DifferentialMesh3d(NMAX, 2*FARL, NEARL, FARL)

  load(DM, fn_in)

  for i in xrange(72):

    t1 = time()

    DM.optimize_position(OPT_STP, OPT_ITT)

    DM.optimize_edges(H*1.8, NEARL*0.5)

    print_stats(i, time()-t1, DM)

  export(DM, fn_out)


if __name__ == '__main__' :

    main()

