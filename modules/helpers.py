#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

def get_args():

  import argparse

  parser = argparse.ArgumentParser()
  parser.add_argument(
    '--obj',
    type=str,
    default='./data/base.obj'
  )
  parser.add_argument(
    '--procs',
    type=int,
    default=4,
    help='number of processors.'
  )
  parser.add_argument(
    '--profile',
    type=bool,
    default=False,
    help='run in profiler mode?'
  )
  parser.add_argument(
    '--nearl',
    type=float,
    default=0.003
  )
  parser.add_argument(
    '--farl',
    type=float,
    default=0.03
  )
  parser.add_argument(
    '--stp',
    type=float,
    default=1.0e-7
  )
  parser.add_argument(
    '--reject',
    type=float,
    default=1.0
  )
  parser.add_argument(
    '--attract',
    type=float,
    default=0.3
  )
  parser.add_argument(
    '--unfold',
    type=float,
    default=0.1
  )
  parser.add_argument(
    '--nmax',
    type=int,
    default=1000000
  )
  parser.add_argument(
    '--itt',
    type=int,
    default=10000000000
  )
  parser.add_argument(
    '--stat',
    type=int,
    default=100
  )
  parser.add_argument(
    '--export',
    type=int,
    default=1000
  )
  parser.add_argument(
    '--out',
    type=str,
    default='./res/res'
  )
  parser.add_argument(
    '--scale',
    type=float,
    default=0.009
  )

  return parser.parse_args()


def print_stats(steps,t_diff,dm):

  import time

  s = '{:s} | steps: {:d} time: {:.5f} vnum: {:d} henum: {:d} fnum: {:d}'.format(
    time.strftime('%d/%m/%y %H:%M:%S'),
    steps,
    t_diff,
    dm.get_vnum(),
    dm.get_henum(),
    dm.get_fnum()
  )

  print(s)

  return

