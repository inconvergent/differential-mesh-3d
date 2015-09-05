#!/bin/bash

set -e

here=$(pwd)

rm -f ./res/res*.obj

./main_sources.py

