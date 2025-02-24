#!/bin/tcsh
ln -s ../nemo3d-x86_64-linux_intel9_64_mpi_intelfast.ex fmtdat.ex
./fmtdat.ex -a2 silicon4.nd_Ek
awk '$4 > 0 {print $4}' silicon4.nd_Ek_ascii > Ek

