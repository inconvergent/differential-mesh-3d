import bpy
import bmesh


class Mesh(object):

  def __init__(self,in_fn):

    names = ['Cube', 'Icosphere', 'Sphere']
    for name in names:
      try:
        bpy.data.objects[name].select = True
        bpy.ops.object.delete()
      except Exception:
        pass

    obj_name = 'thing'
    self.obj_name = obj_name

    mesh = bpy.data.meshes.new('mesh')
    obj = bpy.data.objects.new(obj_name, mesh)

    scn = bpy.context.scene
    scn.objects.link(obj)
    scn.objects.active = obj
    obj.select = True

    self.obj_name = obj_name
    self.obj = obj

    self.data = self.__load_from_file(in_fn)

    return

  def __load_from_file(self,fn):

    from json import load
    from codecs import open

    with open(fn,'rb', encoding='utf8') as f:

      data = load(f)

      for k,v in data.items():
        print(k, len(v))

    return data

  def __get_bmesh(self):

    bpy.data.objects[self.obj_name].select = True
    self.obj = bpy.context.active_object
    bpy.ops.object.mode_set(mode='EDIT')
    bm = bmesh.from_edit_mesh(self.obj.data)

    return bm

  def __to_mesh(self):

    bpy.ops.object.mode_set(mode='OBJECT')

    return

  def rescale(self, scale):

    obj_name = self.obj_name

    bpy.ops.object.origin_set(type='GEOMETRY_ORIGIN')

    o = bpy.data.objects[obj_name]

    x,y,z = o.scale

    x *= scale
    y *= scale
    z *= scale

    o.scale = ((x,y,z))
    #bpy.data.objects[obj_name].dimensions.z *= scale


  def delete_from_scene(self):

      bpy.data.objects[self.obj_name].select = True
      bpy.ops.object.delete()


  def build(self):

    bm = self.__get_bmesh()

    vertices = []
    faces = []

    for v in self.data['vertices']:

      vertices.append(bm.verts.new(v))

    for f in self.data['triangles']:

      faces.append(bm.faces.new([vertices[v] for v in f]))

    self.__to_mesh()

  def save_blend(self,fn):

    bpy.ops.wm.save_as_mainfile(filepath=fn)

    return

  def save_obj(self,fn):

    bpy.ops.export_scene.obj(
      filepath=fn,
      use_smooth_groups=False,
      use_triangles=True,
      use_edges=True,
      use_normals=False,
      use_materials=False
    )

    return

def main():

  from time import time
  import glob, os

  dirname = './res/'

  os.chdir(dirname)
  for fn in sorted(glob.glob('res_*.json')):

    print('exporting: ' + fn)

    t1 = time()

    LM = Mesh(fn)
    LM.build()
    LM.rescale(100)

    fn_out = fn[:-4] + 'obj'

    LM.save_obj(fn_out)

    LM.delete_from_scene()

    print('\ntime:',time()-t1,'\n\n')

  return


if __name__ == '__main__':

  main()

