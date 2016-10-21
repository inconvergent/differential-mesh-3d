#!/usr/bin/python3
# -*- coding: utf-8 -*-

from time import time
from numpy import zeros

NMAX = 10000000

PROCS = 4
STP = 1.0e-6

NEARL = 0.001
FARL = 0.004
ZONEWIDTH = 0.004

REJECT = 1.0*STP
ATTRACT = 0.9*STP

DIMINISH = 0.99
SMOOTH = 0.05

SPLIT_LIMIT = NEARL*0.9
FLIP_LIMIT = NEARL*0.5

OBJ = './data/cyl.obj'

SEEDTYPE = 'surface'
SEEDRATIO = 1.0
SEED_FREQ = 2

SCALE = 0.006

EXPORT = 1000
STAT = 100
DRAW = 100

from OpenGL.GL import glNewList
from OpenGL.GL import glBegin
from OpenGL.GL import glEnd
from OpenGL.GL import glEndList
from OpenGL.GL import glColor3f
from OpenGL.GL import glNormal3f
from OpenGL.GL import glColor3d
from OpenGL.GL import glVertex3f
from OpenGL.GL import GL_COMPILE
from OpenGL.GL import GL_TRIANGLES

from numpy import cross
from numpy import array
from numpy import mean
from numpy import concatenate
from numpy.linalg import norm

from numpy.random import random

def show_geometry(np_verts, np_tris, np_int, tnum, vnum):
  glNewList(1, GL_COMPILE)
  glBegin(GL_TRIANGLES)
  glColor3f(1, 1, 1)

  coords = np_verts[np_tris[:tnum,:],:]
  mi = coords[:, :, 0].min(axis=0)
  ma = coords[:, :, 0].max(axis=0)

  for f, vv in enumerate(coords):
    v1 = vv[2, :] - vv[0, :]
    v2 = vv[1, :] - vv[0, :]
    n = cross(v2, v1).squeeze()
    n /= norm(n)
    glNormal3f(*n)
    c = np_int[f]
    glColor3f(c, c, c)
    glVertex3f(*vv[0, :].squeeze())
    glVertex3f(*vv[1, :].squeeze())
    glVertex3f(*vv[2, :].squeeze())

  glEnd()
  glEndList()

  return concatenate((mi,ma))

def move_scale(verts, s=1.0):
  mm = mean(verts, axis=0)
  verts[:, :] -= mm
  verts[:, :] *= s


def main():

  from differentialMesh3d import DifferentialMesh3d
  from modules.utils import print_stats
  from modules.utils import get_seed_selector
  from iutils.ioOBJ import export

  from iutils.ioOBJ import load_move_scale as load_obj
  from iutils.random import random_unit_vec

  from fn import Fn
  fn = Fn(prefix='./res/', postfix='.obj')


  DM = DifferentialMesh3d(
      nmax=NMAX,
      zonewidth=ZONEWIDTH,
      nearl=NEARL,
      farl=FARL,
      procs=PROCS
      )

  data = load_obj(
      OBJ,
      s = SCALE,
      mx = [0.5]*3
      )
  info = DM.initiate_faces(list(data['vertices']), list(data['faces']))
  if info['minedge']<NEARL:
    return

  seed_selector = get_seed_selector(DM, SEEDTYPE, SEEDRATIO)

  noise = random_unit_vec(DM.get_vnum(), STP*100.)
  DM.position_noise(noise, scale_intensity=-1)

  DM.optimize_edges(SPLIT_LIMIT, FLIP_LIMIT)

  for he in range(DM.get_henum()):
    DM.set_edge_intensity(he, 1.0)

  np_verts = zeros((NMAX,3),'float')
  np_tris = zeros((NMAX,3),'int')
  np_int = zeros(NMAX,'float')

  def do_export():
    vnum = DM.np_get_vertices(np_verts)
    tnum = DM.np_get_triangles_vertices(np_tris)
    DM.np_get_triangles_intensity(np_int)
    move_scale(np_verts[:vnum,:], s=1000)
    export(
        'thing_mesh',
        fn.name(),
        verts=np_verts[:vnum,:],
        tris=np_tris[:tnum,:]
        )

  def geometry_generator():
    seeds = seed_selector()
    i = 0
    k = 0
    box = array((-1,-1,-1,1,1), 'float')
    while True:
      i += 1
      for _ in range(50):
        k += 1
        DM.optimize_position(
            REJECT,
            ATTRACT,
            DIMINISH,
            SMOOTH,
            scale_intensity=1
            )

        DM.optimize_edges(SPLIT_LIMIT, FLIP_LIMIT)

        if len(seeds)>0:
          DM.set_vertices_intensity(seeds, 1.0)

        if k%SEED_FREQ == 0:
          seeds = seed_selector()

      print_stats(k, DM, meta='alive v: {:d}'.format(len(seeds)))
      vnum = DM.np_get_vertices(np_verts)
      tnum = DM.np_get_triangles_vertices(np_tris)
      DM.np_get_triangles_intensity(np_int)
      move_scale(np_verts[:vnum,:])
      box = show_geometry(np_verts, np_tris, np_int, tnum, vnum)

      yield box

  from view3d import View3d
  v3d = View3d(1000)
  v3d.start(geometry_generator)

  do_export()


if __name__ == '__main__' :
  main()

