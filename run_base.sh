#!/bin/bash

set -e

here=$(pwd)
blender="/opt/blender-2.75ax64/blender"

name="$1"

run_script="$here/export_base_geometry.py"

rm -f ./data/base.json

"$blender"  "$here/data/sphere.blend" -b -P "$run_script"

