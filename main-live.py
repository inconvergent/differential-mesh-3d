#!/usr/bin/python3
# -*- coding: utf-8 -*-

from numpy import zeros

PROCS = 4
STP = 1.0e-6
NEARL = 0.003
FARL = 0.05
REJECT = 0.5*STP
ATTRACT = 0.9*STP
UNFOLD = 0.1*STP
TRIANGLE = 0.01*STP
DIMINISH = 0.99
SMOOTH = 0.08
STAT = 100
SPLIT_LIMIT = NEARL*1.2
FLIP_LIMIT = NEARL*0.5
SEED_FREQ = 500
NMAX = 10000000
OBJ = './data/square.obj'
SEEDTYPE = 'random'
SEEDRATIO = 1.0
SCALE = 0.02

from OpenGL.GL import glNewList
from OpenGL.GL import glBegin
from OpenGL.GL import glEnd
from OpenGL.GL import glEndList
from OpenGL.GL import glColor3f
from OpenGL.GL import glNormal3f
from OpenGL.GL import glVertex3f
from OpenGL.GL import GL_COMPILE
from OpenGL.GL import GL_TRIANGLES

from numpy.random import random
from numpy import cross
from numpy.linalg import norm

def show_geometry(np_verts, np_tris, tnum, vnum):
  glNewList(1, GL_COMPILE)
  glBegin(GL_TRIANGLES)
  glColor3f(1, 1, 1)

  for vv in np_verts[np_tris[:tnum,:],:]:
    v1 = vv[2, :] - vv[0, :]
    v2 = vv[1, :] - vv[0, :]
    n = cross(v2, v1).squeeze()
    n /= norm(n)
    glNormal3f(*n)
    glVertex3f(*vv[0, :].squeeze())
    glVertex3f(*vv[1, :].squeeze())
    glVertex3f(*vv[2, :].squeeze())

  glEnd()
  glEndList()


def main():

  from differentialMesh3d import DifferentialMesh3d
  from modules.utils import print_stats
  from modules.utils import get_seed_selector

  from iutils.ioOBJ import load_move_scale as load_obj
  from iutils.random import random_unit_vec


  DM = DifferentialMesh3d(
      nmax = NMAX,
      zonewidth = FARL,
      nearl = NEARL,
      farl = FARL,
      procs = PROCS
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

  noise = random_unit_vec(DM.get_vnum(), STP*1000.)
  DM.position_noise(noise, scale_intensity=-1)

  DM.optimize_edges(SPLIT_LIMIT, FLIP_LIMIT)

  for he in range(DM.get_henum()):
    DM.set_edge_intensity(he, 1.0)

  def geometry_generator():
    seeds = seed_selector()
    i = 0
    np_verts = zeros((NMAX,3),'float')
    np_tris = zeros((NMAX,3),'int')
    while True:
      for k in range(100):
        i += 1
        DM.optimize_position(
            REJECT,
            ATTRACT,
            UNFOLD,
            TRIANGLE,
            DIMINISH,
            SMOOTH,
            scale_intensity=1
            )

        DM.optimize_edges(SPLIT_LIMIT, FLIP_LIMIT)

        if len(seeds)>0:
          DM.set_vertices_intensity(seeds, 1.0)

        if i%SEED_FREQ == 0:
          seeds = seed_selector()

        if i%STAT==0:
          print_stats(i, DM, meta='alive v: {:d}'.format(len(seeds)))

        num_verts = DM.get_vnum()
        vnum = DM.np_get_vertices(np_verts)
        tnum = DM.np_get_triangles_vertices(np_tris)

        show_geometry(np_verts, np_tris, tnum, vnum)
      yield i

  from view3d import View3d
  v3d = View3d(1000)
  v3d.start(geometry_generator)


if __name__ == '__main__' :
  main()

