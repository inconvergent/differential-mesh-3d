# -*- coding: utf-8 -*-

from sys import path
path.append('.')

import bpy


def main(argv):

  # from modules.blender_utils import Obj
  from dddUtils.blender import Obj

  name = argv[0]

  print('importing: ' + name)

  O = Obj(name, 'a')
  # O.get_vertex_color()
  O.smooth()
  O.set_smooth_shade()
  # O.move_rescale([-0.5]*3, 100)
  O.move_rescale(set_pivot=[0.5,-0.5,0.5], pos=[0,0,0], scale=100)
  O.apply_mat()

  prefix = name.split('_')[0]

  bpy.ops.wm.save_as_mainfile(filepath='{:s}.blend'.format(prefix))


if __name__ == '__main__':

  import sys
  argv = sys.argv
  argv = argv[argv.index("--") + 1:]
  main(argv)

