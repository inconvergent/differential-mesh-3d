# -*- coding: utf-8 -*-

from OpenGL.GL import glNewList
from OpenGL.GL import glBegin
from OpenGL.GL import glEnd
from OpenGL.GL import glEndList
from OpenGL.GL import glColor3f
from OpenGL.GL import glNormal3f
from OpenGL.GL import glVertex3f
from OpenGL.GL import GL_COMPILE
from OpenGL.GL import GL_TRIANGLES

from numpy import cross
from numpy import mean
from numpy import concatenate
from numpy.linalg import norm

def show_geometry(np_verts, np_tris, np_int, tnum, vnum):
  glNewList(1, GL_COMPILE)
  glBegin(GL_TRIANGLES)
  glColor3f(1, 1, 1)

  coords = np_verts[np_tris[:tnum, :], :]
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
