#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

from modules.utils import export_obj
from modules.utils import load_obj
from modules.utils import random_unit_vec
from modules.utils import get_surface_edges


MOVE = [0.5]*3


def main(args):

  from differentialMesh3d import DifferentialMesh3d
  from modules.helpers import print_stats
  from numpy import array
  from numpy.random import random
  from modules.helpers import make_info_str


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
  seed_ratio = args.seedRatio
  seed_freq = args.seedFreq
  vnum_max = args.vnum

  DM = DifferentialMesh3d(
    nmax = args.nmax,
    zonewidth = args.farl,
    nearl = args.nearl,
    farl = args.farl,
    procs = args.procs
  )

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

  alive_vertices = array(get_surface_edges(DM))

  DM.optimize_edges(split_limit, flip_limit)

  for he in xrange(DM.get_henum()):
    DM.set_edge_intensity(he, 1.0)

  for i in xrange(args.itt):

    try:

      DM.optimize_position(
        reject,
        attract,
        unfold,
        triangle,
        diminish,
        smooth,
        scale_intensity=1
      )

      DM.optimize_edges(split_limit, flip_limit)

      ## TODO: make this more configurable.
      ## this only works for open surfaces.

      if i%seed_freq == 0:
        alive_vertices = array([l for l in get_surface_edges(DM) if random()<seed_ratio])
        print('number of alive vertices: {:d}'.format(len(alive_vertices)))

      if len(alive_vertices)>0:
        DM.set_vertices_intensity(alive_vertices, 1.0)

      if i%stat==0:
        print_stats(i, DM)

      if i%export==0:
        export_obj(
          DM,
          'thing_mesh',
          '{:s}_{:012d}.obj'.format(out, i),
          write_intensity=False,
          meta=make_info_str(args)
        )

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

