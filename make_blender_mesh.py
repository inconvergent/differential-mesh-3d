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

  def build(self):

    bm = self.__get_bmesh()

    vertices = []
    faces = []

    for v in self.data['vertices']:

      vertices.append(bm.verts.new(v))

    for f in self.data['triangles']:

      faces.append(bm.faces.new([vertices[v] for v in f]))

    self.__to_mesh()

  def save(self,fn):

    bpy.ops.wm.save_as_mainfile(filepath=fn)

    return

def main():

  from time import time

  fn_in = './res/dat.json'
  fn_out = './res/res.blend'

  t1 = time()

  LM = Mesh(fn_in)

  LM.build()

  LM.save(fn_out)

  print('\ntime:',time()-t1,'\n\n')

  return


if __name__ == '__main__':

  #import argparse

  #parser = argparse.ArgumentParser()
  #parser.add_argument(
    #'-i',
    #'--input',
    #help='input pickle file',
    #default='./res/in.pkl'
  #)

  #parser.add_argument(
    #'-o',
    #'--output',
    #help='output blender file',
    #default='./res/res.blend'
  #)
  #args = parser.parse_args()

  main()

