#!/bin/bash

set -e

here=$(pwd)
blender="/opt/blender-2.75ax64/blender"

name="$1"

mesh="$here/export_base_geometry.py"

"$blender"  "$here/data/sphere.blend" -b -P "$mesh"  

