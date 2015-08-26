#!/bin/bash

set -e

here=$(pwd)
blender="/opt/blender-2.75ax64/blender"

name="$1"

mesh="$here/make_blender_mesh.py"

"$blender"  "$here/res/empty.blend" -b -P "$mesh"  
#echo "please specify a name for the result.";
#exit 1;

