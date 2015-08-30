import bpy
import bmesh


class Mesh(object):

  def __init__(self,in_fn):

    obj_name = 'thing'

    self.obj_name = obj_name

    bpy.ops.object.select_pattern(pattern=obj_name)

    return

  def __get_bmesh(self):

    bpy.data.objects[self.obj_name].select = True
    self.obj = bpy.context.active_object
    bpy.ops.object.mode_set(mode='EDIT')
    bm = bmesh.from_edit_mesh(self.obj.data)

    return bm

  def __to_mesh(self):

    bpy.ops.object.mode_set(mode='OBJECT')

    return

  def build(self, fn):

    from json import dump
    from codecs import open


    world = bpy.context.active_object.matrix_world


    bm = self.__get_bmesh()

    vertices = []
    triangles = []

    for v in bm.verts:

      g_co = world*v.co

      vertices.append([g_co.x, g_co.y, g_co.z])

    for p in bm.faces:

      face = [loop.vert.index for loop in p.loops]
      triangles.append(face)

    data = {
      'vertices': vertices,
      'triangles': triangles
    }

    print('exporting mesh ...')
    print('num vertices: {:d}, num triangles: {:d}'.format(
      len(vertices), len(triangles)))

    with open(fn, 'wb', encoding='utf8') as f:

      dump(data, f)

      print('done.')

def main():

  from time import time

  fn_in = './res/sphere.blend'
  fn_out = './res/base.json'

  t1 = time()

  LM = Mesh(fn_in)

  LM.build(fn_out)

  print('\ntime:',time()-t1,'\n\n')

  return


if __name__ == '__main__':

  main()

