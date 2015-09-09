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

  def del_mesh(self):

    bpy.ops.object.select_all(action='DESELECT')
    self.obj.select=True
    bpy.ops.object.delete()

  def spheres(self):

    from numpy import array
    from numpy import row_stack
    from numpy import diff
    from numpy import mean
    from numpy.linalg import norm
    from collections import defaultdict

    base_scale = 0.45

    scn = bpy.context.scene
    obj = self.obj
    world = obj.matrix_world
    mat = bpy.data.materials["Material"]


    bpy.ops.surface.primitive_nurbs_surface_sphere_add(
      radius = 1,
      location = (0,0,0)
    )
    sphere = bpy.context.active_object
    sphere.data.materials.append(mat)
    bpy.context.active_object.hide = True
    bpy.context.active_object.hide_render = True

    mesh = sphere.data

    vertex_map = defaultdict(list)
    vertices = obj.data.vertices
    edges = obj.data.edges
    vnum = len(vertices)
    enum = len(edges)
    vert_co = row_stack([world*v.co for v in vertices])

    print('\ncalculating sizes:\n')

    for i,e in enumerate(edges):

      vv = array([v for v in e.vertices],'int')
      dx = diff(vert_co[vv,:],axis=0).squeeze()
      #dx = diff(row_stack([ \for v in vv]),axis=0).squeeze()
      dd = norm(dx)

      vertex_map[vv[0]].append(dd)
      vertex_map[vv[1]].append(dd)

      if i%100==0:
        print(i, enum)

    vertex_weight = {k:mean(d) for k,d in vertex_map.items()}

    print('\nplacing objects:\n')

    for i,(v,rad) in enumerate(vertex_weight.items()):

      o = bpy.data.objects.new('one', mesh)
      #o.location = world*vertices[v].co
      o.location = vert_co[v,:]

      scale = rad*base_scale

      sx,sy,sz = o.scale
      sx *= scale
      sy *= scale
      sz *= scale

      o.scale = ((sx,sy,sz))
      scn.objects.link(o)

      if i%100==0:
        print(i, vnum)



def main(argv):

  from time import time

  fn = argv[0]

  print('importing: ' + fn)

  t1 = time()

  O = Obj(fn,'a')
  O.move_rescale([-0.5]*3, 100)
  O.spheres()
  O.del_mesh()

  print('\ntime:',time()-t1,'\n\n')

  bpy.ops.wm.save_as_mainfile(filepath='./scene_sphere.blend')


if __name__ == '__main__':

  import sys
  argv = sys.argv
  argv = argv[argv.index("--") + 1:]
  main(argv)

