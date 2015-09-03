#!/bin/bash

set -e

here=$(pwd)

rm -f ./res/res*.json

./main.py

