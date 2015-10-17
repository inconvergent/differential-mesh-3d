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
  from numpy import array
  from numpy.random import random


  reject = args.reject*args.stp
  attract = args.attract*args.stp
  unfold = args.unfold*args.stp
  triangle = args.triangle*args.stp
  diminish = args.diminish
  smooth = args.smooth
  stat = args.stat
  export = args.export
  out = args.out
  split_limit = args.nearl*1.2
  flip_limit = args.nearl*0.5
  vnum_max = args.vnum

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

  alive_vertices = get_surface_edges(DM)

  print_stats(-2, 0.0, DM)
  DM.optimize_edges(split_limit, flip_limit)
  print_stats(-1, 0.0, DM)

  for he in xrange(DM.get_henum()):
    DM.set_edge_intensity(he, 1.0)

  for i in xrange(args.itt):

    try:

      t1 = time()

      DM.optimize_position(reject, attract, unfold, triangle, OPT_ITT, scale_intensity=1)

      DM.optimize_edges(split_limit, flip_limit)

      DM.diminish_all_vertex_intensity(diminish)

      if i%100 == 0:
        alive_vertices = [l for l in get_surface_edges(DM) if random()<0.95]
        print('number of alive vertices: {:d}'.format(len(alive_vertices)))

      if len(alive_vertices)>0:
        DM.set_vertices_intensity(array(alive_vertices), 1.0)

      DM.smooth_intensity(smooth)

      if i%stat==0:
        print_stats(i, time()-t1, DM)

      if i%export==0:
        fn = '{:s}_{:08d}.obj'.format(out, i)
        export_obj(DM, 'thing_mesh', fn, write_intensity=False)

      if DM.get_vnum()>vnum_max:
        return

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

