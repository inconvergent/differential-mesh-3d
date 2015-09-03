#!/bin/bash

set -e

here=$(pwd)
blender="/opt/blender-2.75ax64/blender"

name="$1"

run_script="$here/make_animated_scene.py"

"$blender"  "$here/data/scene.blend" -b -P "$run_script"

