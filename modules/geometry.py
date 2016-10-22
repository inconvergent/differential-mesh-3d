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
from numpy import zeros
from numpy import mean
from numpy import concatenate
from numpy.linalg import norm

def get_show_geometry(dm, nmax):

  np_verts = zeros((nmax, 3), 'float')
  np_tris = zeros((nmax, 3), 'int')
  np_int = zeros(nmax, 'float')

  def show_geometry():
    glNewList(1, GL_COMPILE)
    glBegin(GL_TRIANGLES)
    glColor3f(1, 1, 1)

    vnum = dm.np_get_vertices(np_verts)
    tnum = dm.np_get_triangles_vertices(np_tris)
    dm.np_get_triangles_intensity(np_int)
    move_scale(np_verts[:vnum, :])

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
    return concatenate((mi, ma))
  return show_geometry

def move_scale(verts, s=1.0):
  mm = mean(verts, axis=0)
  verts[:, :] -= mm
  verts[:, :] *= s
