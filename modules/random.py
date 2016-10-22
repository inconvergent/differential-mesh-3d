# -*- coding: utf-8 -*-

from numpy.linalg import norm
from numpy import reshape


def random_unit_vec(num, scale):
  from numpy.random import normal

  rnd = normal(size=(num, 3))
  d = norm(rnd, axis=1)
  rnd[:] /= reshape(d, (num, 1))
  return rnd*scale

