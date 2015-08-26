# Differential Mesh 3d

Differential Mesh is an extension of the Hyphae
(https://github.com/inconvergent/hyphae) algorithm. I started working on it
with the intention of mimicking the growth of (certain types of) lichen. The
results, to me, are equal parts crystal-like and biological.

## Prerequisites

In order for this code to run you must first download and install these two
repositories:

*    `zonemap3d`: https://github.com/inconvergent/zonemap-3d

## Other Dependencies

The code also depends on:

*    `numpy`
*    `scipy`
*    `cython`
*    `python-cairo` (do not install with pip, this generally does not work)

## Running it on Linux (Ubuntu)

To install the libraries locally, run `make`. I have only tested this code in
Ubuntu 14.04 LTS, but my guess is that it should work on most other platforms
platforms as well.  However i know that the scripted install in `make` will not
work in Windows

## Running it on Windows?

The code will probably work just fine under Windows, but I'm not sure how to
install it. (Let me know if you get it working!)

## Similar code

If you find this alorithm insteresting you might also want to check out:

*    https://github.com/inconvergent/differential-line
*    https://github.com/inconvergent/differential-mesh

-----------
http://inconvergent.net

