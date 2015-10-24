# -*- coding: utf-8 -*-

import bpy

class Obj(object):

  def __init__(self, fn, obj_name):

    self.fn = fn
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

  def get_vertex_color(self):

    from mathutils import Color

    colors = []

    try:

      with open(self.fn+'.x', 'r', encoding='utf8') as f:

        for l in f:
          if l.startswith('#'):
            continue

          values = l.split()
          if not values:
            continue

          if values[0] == 'c':
            c = [float(v) for v in values[1:]]
            colors.append(c)

    except FileNotFoundError:
      return

    mesh = self.obj.data

    if not mesh.vertex_colors:
      mesh.vertex_colors.new()

    #print(mesh.vertex_colors)
    col = mesh.vertex_colors.active

    print(col)



    num = len(colors)

    numv = len(self.obj.data.polygons)

    #for i,c in enumerate(colors):

      #col.data[i].color = c

    i = 0

    print(Color(([0.1,0.1,0.1])))

    for poly in self.obj.data.polygons:
      loop = poly.loop_indices
      verts = poly.vertices
      #print([v for v in poly.loop_indices])
      #print([v for v in poly.vertices])
      for idx,v in zip(loop,verts):
        col.data[idx].color = Color(colors[v])
        i += 1

    print(num, numv, len(col.data), i)

    #for i,(rgb) in enumerate(row_stack(colors)):
      #col.data[i].color = [float(i)/num]*3

    #print(colors)


  def move_rescale(self, pos, scale):

    #obj_name = self.obj_name
    bpy.ops.object.origin_set(type='GEOMETRY_ORIGIN')

    #o = bpy.data.objects[obj_name]
    obj = self.obj

    sx,sy,sz = obj.scale

    sx *= scale
    sy *= scale
    sz *= scale

    obj.scale = ((sx,sy,sz))
    #bpy.data.objects[obj_name].dimensions.z *= scale

  def smooth(self, levels):

    bpy.context.scene.objects.active = self.obj
    #bpy.context.scene.objects.selected = self.obj

    bpy.ops.object.modifier_add(type='SUBSURF')
    self.obj.modifiers['Subsurf'].levels = 1
    self.obj.modifiers['Subsurf'].render_levels = 2

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

