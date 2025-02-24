#### Time-stamp: <run.sh  clemensh  2005-11-09 13:12 EST>

set -x

## Run conversion like this:
##./convert-to-dx GaAs_QD_1band_cube.nd_rAtom_0_ascii GaAs_QD_1band_cube.nd_wf_1_ascii wf.dx
## Nemo writes DX native format now.  You need wf_* and dx_shape.

## Write tiff files.
dx -script dx-multiple-isosurfaces-boundary-axes.net   # writes image.tiff
dx -script dx-multiple-isosurfaces-boundary-axes-x.net # writes image-x.tiff
dx -script dx-multiple-isosurfaces-boundary-axes-y.net # writes image-y.tiff
dx -script dx-multiple-isosurfaces-boundary-axes-z.net # writes image-z.tiff

#### local variables:
#### end:
