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
    default=0.05
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
    '--triangle',
    type=float,
    default=0.1
  )
  parser.add_argument(
    '--diminish',
    type=float,
    default=0.99
  )
  parser.add_argument(
    '--smooth',
    type=float,
    default=0.08
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
    '--vnum',
    type=int,
    default=10000000000
  )
  parser.add_argument(
    '--stat',
    type=int,
    default=100
  )
  parser.add_argument(
    '--seedFreq',
    type=int,
    default=100
  )
  parser.add_argument(
    '--seedRatio',
    type=float,
    default=1.0
  )
  parser.add_argument(
    '--seedType',
    type=str,
    default='random'
  )
  parser.add_argument(
    '--exportLeap',
    type=int,
    default=100
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

def make_info_str(args):
  s = ''
  for k in vars(args):
    s += '# ' + str(k) + ': ' + str(getattr(args,k)) + '\n'
  return s


def print_stats(steps,dm, meta=False):

  from time import strftime
  from time import time

  if isinstance(meta, str):
    meta = ' | {:s}'.format(meta)
  else:
    meta = ''

  print(
    '{:s} | stp: {:d} sec: {:.2f} v: {:d} e: {:d} f: {:d}{:s}'
      .format(
      strftime('%d/%m/%y %H:%M:%S'),
      steps,
      time()-dm.get_start_time(),
      dm.get_vnum(),
      dm.get_henum(),
      dm.get_fnum(),
      meta
    )
  )

  return

