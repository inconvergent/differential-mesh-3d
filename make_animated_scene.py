# -*- coding: utf-8 -*-

from sys import path
path.append('.')

import bpy


def main(argv):

  from glob import glob
  from os import chdir
  from modules.blender_utils import Obj

  name = argv[0]
  dirname = './res/'

  objs = []

  count = 0

  chdir(dirname)
  for fn in sorted(glob('{:s}_*.obj'.format(name))):

    print('importing: ' + fn)

    O = Obj(fn, 'a')
    O.get_vertex_color()
    O.smooth(1)
    O.move_rescale([-0.5]*3, 100)
    O.animate_vis(count, count+1)
    O.apply_mat()
    objs.append(O)

    count += 1

  bpy.data.scenes['Scene'].frame_current = 1
  bpy.data.scenes['Scene'].frame_end = count-1

  #chdir('..')
  bpy.ops.wm.save_as_mainfile(filepath='./scene_ani_{:s}.blend'.format(name))


if __name__ == '__main__':

  import sys
  argv = sys.argv
  argv = argv[argv.index("--") + 1:]
  main(argv)

