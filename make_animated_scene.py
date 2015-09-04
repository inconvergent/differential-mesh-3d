import bpy


class Obj(object):

  def __init__(self, fn, obj_name):

    self.obj_name = obj_name
    self.obj = self.__import(fn)

    return

  def __import(self, fn):

    bpy.ops.object.select_all(action='DESELECT')

    bpy.ops.import_scene.obj(
      filepath=fn,
      use_smooth_groups=False,
      use_edges=True,
    )

    obj = bpy.context.selected_objects[0]

    return obj

  def smooth(self, levels):

    bpy.context.scene.objects.active = self.obj
    #bpy.context.scene.objects.selected = self.obj

    bpy.ops.object.modifier_add(type='SUBSURF')
    self.obj.modifiers['Subsurf'].levels = levels
    self.obj.modifiers['Subsurf'].render_levels = levels

    bpy.ops.object.shade_smooth()

  def __set_vis(self, frame, vis=True):

    bpy.context.scene.objects.active = self.obj

    bpy.data.scenes['Scene'].frame_current = frame
    bpy.context.active_object.hide = not vis
    bpy.context.active_object.hide_render = not vis

    bpy.context.active_object.keyframe_insert(
      data_path="hide",
      index=-1,
      frame=frame
    )
    bpy.context.active_object.keyframe_insert(
      data_path="hide_render",
      index=-1,
      frame=frame
    )

  def animate_vis(self, ain, aout):

    self.__set_vis(0, False)
    self.__set_vis(ain, True)
    self.__set_vis(aout, False)

  def apply_mat(self):

    mat = bpy.data.materials["Material"]
    self.obj.data.materials.append(mat)


def main():

  from time import time
  import glob, os

  dirname = './res/'

  objs = []

  count = 0

  os.chdir(dirname)
  for fn in sorted(glob.glob('res_*.obj')):

    print('importing: ' + fn)

    t1 = time()

    O = Obj(fn,'a')
    O.smooth(1)
    O.animate_vis(count, count+1)
    O.apply_mat()
    objs.append(O)

    count += 1

    print('\ntime:',time()-t1,'\n\n')

  print('asdf')

  os.chdir('..')
  bpy.ops.wm.save_as_mainfile(filepath='./final.blend')


if __name__ == '__main__':

  main()

