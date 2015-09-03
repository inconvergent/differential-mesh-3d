#!/bin/bash

set -e

here=$(pwd)
blender="/opt/blender-2.75ax64/blender"

name="$1"

run_script="$here/export_res_geometry.py"

rm -f ./res/res*.obj

"$blender"  "$here/data/empty.blend" -b -P "$run_script"

