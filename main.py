#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

from modules.utils import export_obj
from modules.utils import load_obj
from modules.utils import random_unit_vec
from modules.utils import get_surface_edges


NMAX = int(10e7)
ITT = int(10e7)
OPT_ITT = 1

NEARL = 0.0028
H = NEARL*1.2

FARL = 0.05

FLIP_LIMIT = NEARL*0.5

EXPORT_ITT = 1000
STAT_ITT = 1


#STP = 1.0e-6
STP = 1.0e-7
REJECT_STP = STP*1.0
TRIANGLE_STP = STP*0.1
ATTRACT_STP = STP*0.1
UNFOLD_STP = STP*0.2
COHESION_STP = STP*0.



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

  DM = DifferentialMesh3d(NMAX, FARL, NEARL, FARL)

  data = load_obj(
    fn_obj,
    sx = [0.01]*3,
    mx = [0.5]*3
  )
  DM.initiate_faces(data['vertices'], data['faces'])

  noise = random_unit_vec(DM.get_vnum(), 1.0e-5)
  DM.position_noise(noise, scale_intensity=-1)

  #alive_vertices = set(randint(DM.get_vnum(), size=DM.get_vnum()))
  #print(alive_vertices)

  alive_vertices = list(l for l in set(get_surface_edges(DM)) if random()<0.5)
  print(alive_vertices)

  DM.optimize_edges(H, FLIP_LIMIT)

  for he in xrange(DM.get_henum()):
    DM.set_edge_intensity(he, 1.0)

  for i in xrange(ITT):

    try:

      t1 = time()

      DM.optimize_position(
        REJECT_STP,
        TRIANGLE_STP,
        ATTRACT_STP,
        UNFOLD_STP,
        COHESION_STP,
        OPT_ITT,
        scale_intensity=1
      )

      DM.optimize_edges(H, FLIP_LIMIT)

      DM.diminish_all_vertex_intensity(0.99)

      if i%500 == 0:
        alive_vertices = list(l for l in set(get_surface_edges(DM)) if random()<0.3)
        print(alive_vertices)

      if len(alive_vertices)>0:
        DM.set_vertices_intensity(array([v for v in alive_vertices]), 1.0)
      #for he in unique((random(DM.get_henum())<0.009).nonzero()[0]):
        #DM.add_edge_intensity(he, 0.05)

      DM.smooth_intensity(0.1)

      if i%STAT_ITT==0:
        print_stats(i, time()-t1, DM)

      if i%EXPORT_ITT==0:
        fn = '{:s}_{:06d}.obj'.format(fn_out, i)
        export_obj(DM, 'thing_mesh', fn,NMAX)

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

