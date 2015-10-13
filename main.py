#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

from modules.utils import export_obj
from modules.utils import load_obj
from modules.utils import random_unit_vec
from modules.utils import get_surface_edges

PROCS = 4

NMAX = int(10e6)
ITT = 100000
OPT_ITT = 1

NEARL = 0.003
H = NEARL*1.2

FARL = 0.03

FLIP_LIMIT = NEARL*0.5

EXPORT_ITT = 1000
STAT_ITT = 10

SCALE = [0.009]*3
MOVE = [0.5]*3


STP = 1.0e-7
REJECT_STP = STP
ATTRACT_STP = STP*0.1
UNFOLD_STP = STP*0.01



def main(argv):

  from differentialMesh3d import DifferentialMesh3d
  from time import time
  from modules.helpers import print_stats
  from numpy import unique
  from numpy import array
  from numpy.random import random
  from numpy.random import randint

  name = argv[0]
  fn_obj = './data/base.obj'
  fn_out = './res/{:s}'.format(name)

  DM = DifferentialMesh3d(NMAX, FARL, NEARL, FARL, PROCS)

  data = load_obj(
    fn_obj,
    sx = SCALE,
    mx = MOVE
  )
  info = DM.initiate_faces(data['vertices'], data['faces'])
  if info['minedge']<NEARL:
    return

  noise = random_unit_vec(DM.get_vnum(), STP*1000.)
  DM.position_noise(noise, scale_intensity=-1)

  alive_vertices = list(l for l in set(get_surface_edges(DM)))

  DM.optimize_edges(H, FLIP_LIMIT)

  for he in xrange(DM.get_henum()):
    DM.set_edge_intensity(he, 1.0)

  for i in xrange(ITT):

    try:

      t1 = time()

      DM.optimize_position(
        REJECT_STP,
        ATTRACT_STP,
        UNFOLD_STP,
        OPT_ITT,
        scale_intensity=1
      )

      DM.optimize_edges(H, FLIP_LIMIT)

      DM.diminish_all_vertex_intensity(0.99)

      # TODO: this is for testing
      alive_vertices = list(l for l in set(get_surface_edges(DM)))

      #if i%100 == 0:
        #alive_vertices = list(l for l in set(get_surface_edges(DM)) if random()<1)
        #alive_vertices = set(randint(DM.get_vnum(), size=DM.get_vnum()))
      #print('number of alive vertices: {:d}'.format(len(alive_vertices)))

      if len(alive_vertices)>0:
        DM.set_vertices_intensity(array([v for v in alive_vertices]), 1.0)

      DM.smooth_intensity(0.08)

      if i%STAT_ITT==0:
        print_stats(i, time()-t1, DM)

      if i%EXPORT_ITT==0:
        fn = '{:s}_{:08d}.obj'.format(fn_out, i)
        export_obj(DM, 'thing_mesh', fn, write_intensity=False)

    except KeyboardInterrupt:

      break


if __name__ == '__main__' :

  import sys

  argv = sys.argv

  if False:

    import pstats, cProfile
    fn = './profile/profile'
    cProfile.run('main(argv[1:])',fn)
    p = pstats.Stats(fn)
    p.strip_dirs().sort_stats('cumulative').print_stats()

  else:

    main(argv[1:])

