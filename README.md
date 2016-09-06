# Differential Mesh 3d

![img](/img/img.png?raw=true "img")


Differential Mesh is an extension of the Hyphae
(https://github.com/inconvergent/hyphae) algorithm. I started working on it
with the intention of mimicking the growth of (certain types of) lichen.

This three dimensional extension is something else entirely ...

## Prerequisites

In order for this code to run you must first download and install:

*    `zonemap3d`: https://github.com/inconvergent/zonemap-3d
*    `iutils`: https://github.com/inconvergent/iutils
*    `fn`: https://github.com/inconvergent/fn-python3 (only used to generate file
       names, you can skip this by changing file names in `main.py`)

## Other Dependencies

The code also depends on:

*    `numpy`
*    `scipy`
*    `cython`
*    `python-cairo` (do not install with pip, this generally does not work)

## Running it on Linux (Ubuntu)

To install the libraries locally, run `install`. I have only tested this
code in Ubuntu 14.04 LTS, but my guess is that it should work on most other
platforms as well.

Run using the `main-run` script (asdf is a prefix for the generated files and
can be chosen freely):

    ./main-run asdf

