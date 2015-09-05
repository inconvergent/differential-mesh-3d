#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function


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

