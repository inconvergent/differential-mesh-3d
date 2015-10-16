#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

from modules.utils import export_obj
from modules.utils import load_obj
from modules.utils import random_unit_vec
from modules.utils import get_surface_edges

OPT_ITT = 1
MOVE = [0.5]*3


def main(args):

  from differentialMesh3d import DifferentialMesh3d
  from time import time
  from modules.helpers import print_stats
  from numpy import unique
  from numpy import array
  from numpy.random import random
  from numpy.random import randint


  h = args.nearl*1.2
  flip_limit = args.nearl*1.2

  DM = DifferentialMesh3d(args.nmax, args.farl, args.nearl, args.farl, args.procs)

  data = load_obj(
    args.obj,
    sx = [args.scale]*3,
    mx = MOVE
  )
  info = DM.initiate_faces(data['vertices'], data['faces'])
  if info['minedge']<args.nearl:
    return

  noise = random_unit_vec(DM.get_vnum(), args.stp*1000.)
  DM.position_noise(noise, scale_intensity=-1)

  alive_vertices = list(l for l in set(get_surface_edges(DM)))

  DM.optimize_edges(h, flip_limit)

  for he in xrange(DM.get_henum()):
    DM.set_edge_intensity(he, 1.0)

  for i in xrange(args.itt):

    try:

      t1 = time()

      DM.optimize_position(
        args.reject*args.stp,
        args.attract*args.stp,
        args.unfold*args.stp,
        OPT_ITT,
        scale_intensity=1
      )

      DM.optimize_edges(h, flip_limit)

      DM.diminish_all_vertex_intensity(0.99)

      if i%100 == 0:
        alive_vertices = list(l for l in set(get_surface_edges(DM)) if random()<0.7)
        print('number of alive vertices: {:d}'.format(len(alive_vertices)))

      if len(alive_vertices)>0:
        DM.set_vertices_intensity(array([v for v in alive_vertices]), 1.0)

      DM.smooth_intensity(0.08)

      if i%args.stat==0:
        print_stats(i, time()-t1, DM)

      if i%args.export==0:
        fn = '{:s}_{:08d}.obj'.format(args.out, i)
        export_obj(DM, 'thing_mesh', fn, write_intensity=False)

    except KeyboardInterrupt:

      break


if __name__ == '__main__' :

  from modules.helpers import get_args

  args = get_args()
  print(args)

  if args.profile:

    import pstats, cProfile
    cProfile.run('main(args)','./profile/profile')
    p = pstats.Stats('./profile/profile')
    p.strip_dirs().sort_stats('cumulative').print_stats()

  else:

    main(args)

