#### Time-stamp: <example.sh  clemensh  2005-11-08 16:01 EST>

set -x

## Run conversion like this:
##./convert-to-dx GaAs_QD_1band_cube.nd_rAtom_0_ascii GaAs_QD_1band_cube.nd_wf_1_ascii wf.dx

## Write tiff files.
dx -script dx-isosurface-axes.net   # writes image.tiff
dx -script dx-isosurface-axes-x.net # writes image-x.tiff
dx -script dx-isosurface-axes-y.net # writes image-y.tiff
dx -script dx-isosurface-axes-z.net # writes image-z.tiff

#### local variables:
#### end:
