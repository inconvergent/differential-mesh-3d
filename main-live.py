#!/usr/bin/python3
# -*- coding: utf-8 -*-

from numpy import array

NMAX = 1000000

PROCS = 4
STP = 5.0e-5

NEARL = 0.01
FARL = 0.08
ZONEWIDTH = 0.08

REJECT = 0.1*STP
ATTRACT = 0.9*STP

DIMINISH = 0.99
SMOOTH = 0.1

SPLIT_LIMIT = NEARL*1.2
FLIP_LIMIT = NEARL*0.5

OBJ = './data/cyl.obj'
SCALE = 0.06

SEEDTYPE = 'surface'
SEEDRATIO = 1.0
SEED_FREQ = 10

SPEEDUP = 5



def main():

  from differentialMesh3d import DifferentialMesh3d
  from modules.ioOBJ import load_move_scale as load_obj
  from modules.random import random_unit_vec
  from modules.geometry import get_show_geometry
  from modules.utils import get_seed_selector
  from modules.utils import print_stats
  from modules.utils import get_exporter

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
      s=SCALE,
      mx=[0.5] * 3
      )
  info = DM.initiate_faces(list(data['vertices']), list(data['faces']))

  if info['min_edge'] < NEARL:
    print('all edges are too short. try scaling up initial size.')
    return

  seed_selector = get_seed_selector(DM, SEEDTYPE, SEEDRATIO)

  noise = random_unit_vec(DM.get_vnum(), STP*100.)
  DM.position_noise(noise, scale_intensity=-1)

  DM.optimize_edges(SPLIT_LIMIT, FLIP_LIMIT)

  for he in range(DM.get_henum()):
    DM.set_edge_intensity(he, 1.0)

  show_geometry = get_show_geometry(DM, NMAX)

  def geometry_generator():
    seeds = seed_selector()
    i = 0
    k = 0
    box = array((-1, -1, -1, 1, 1), 'float')
    while True:
      i += 1
      for _ in range(SPEEDUP):
        k += 1
        DM.optimize_position(
            REJECT,
            ATTRACT,
            DIMINISH,
            SMOOTH,
            scale_intensity=1
            )

        DM.optimize_edges(SPLIT_LIMIT, FLIP_LIMIT)

        if len(seeds) > 0:
          DM.set_vertices_intensity(seeds, 1.0)

        if k%SEED_FREQ == 0:
          seeds = seed_selector()

      print_stats(k, DM, meta='alive v: {:d}'.format(len(seeds)))

      box = show_geometry()

      yield box

  from view3d import View3d
  v3d = View3d(
      size=512,
      initial_scale=1.3,
      autorotate=True,
      save_img=True
      )
  v3d.start(geometry_generator)

  export = get_exporter(DM, fn, NMAX)
  export()


if __name__ == '__main__' :
  main()

