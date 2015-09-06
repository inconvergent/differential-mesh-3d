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

  def move_rescale(self, pos, scale):

    bpy.ops.object.origin_set(type='GEOMETRY_ORIGIN')
    obj = self.obj

    sx,sy,sz = obj.scale
    sx *= scale
    sy *= scale
    sz *= scale

    obj.scale = ((sx,sy,sz))

  def spheres(self):

    scn = bpy.context.scene
    obj = self.obj
    world = obj.matrix_world
    mat = bpy.data.materials["Material"]

    levels = 1

    bpy.ops.surface.primitive_nurbs_surface_sphere_add(
      radius = 0.2,
      location = (0,0,0)
    )
    sphere = bpy.context.active_object
    sphere.data.materials.append(mat)

    bpy.ops.object.modifier_add(type='SUBSURF')
    sphere.modifiers['Subsurf'].levels = levels
    sphere.modifiers['Subsurf'].render_levels = levels

    mesh = sphere.data

    for i,v in enumerate(obj.data.polygons):
      loc = world*v.center
      ob = bpy.data.objects.new('one', mesh)
      ob.location = loc
      scn.objects.link(ob)
      print(i)


  #def __set_vis(self, frame, vis=True):

    #bpy.context.scene.objects.active = self.obj

    #bpy.data.scenes['Scene'].frame_current = frame
    #bpy.context.active_object.hide = not vis
    #bpy.context.active_object.hide_render = not vis

    #bpy.context.active_object.keyframe_insert(
      #data_path="hide",
      #index=-1,
      #frame=frame
    #)
    #bpy.context.active_object.keyframe_insert(
      #data_path="hide_render",
      #index=-1,
      #frame=frame
    #)

  #def animate_vis(self, ain, aout):

    #self.__set_vis(0, False)
    #self.__set_vis(ain, True)
    #self.__set_vis(aout, False)

  #def apply_mat(self):

    #mat = bpy.data.materials["Material"]
    #self.obj.data.materials.append(mat)


def main(argv):

  from time import time

  fn = argv[0]

  print('importing: ' + fn)

  t1 = time()

  O = Obj(fn,'a')
  O.move_rescale([-0.5]*3, 100)
  O.spheres()

  print('\ntime:',time()-t1,'\n\n')

  bpy.ops.wm.save_as_mainfile(filepath='./scene_sphere.blend')


if __name__ == '__main__':

  import sys
  argv = sys.argv
  argv = argv[argv.index("--") + 1:]
  main(argv)

