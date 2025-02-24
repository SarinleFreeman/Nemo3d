#! /bin/sh

# *****************************************************************************
# The Jet Propulsion Laboratory (JPL) NanoElectronicMOdeling-3D package.
# Copyright (C) 2002 California Institute of Technology (Caltech)
#
# This application is free software, which you can redistribute and/or modify
# under the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation; either version 2.1 of the License, or (at your
# option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library; see the file COPYING. If not, write to the
# Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA  02111-1307  USA
#
# For additional information, please contact
#   Gerhard Klimeck (gekco@jpl.nasa.gov)
#   Fabiano Oyafuso (fabiano@jpl.nasa.gov)
#
# Written by:  Chris Bowen
#              Gerhard Klimeck
#              Fabiano Oyafuso
#              Seungwon Lee
#              Olga Lazarenkova
#              Hook Hua
#
# This product includes software developed by the Apache Software Foundation
# (http://www.apache.org/).
#
# *****************************************************************************
# $Header: /repo/NEMO_3D/nemo3d_setup.sh,v 1.13 2007/09/19 20:53:24 stevenclark Exp $
# *****************************************************************************


# *****************************************************************************
# After untarring NEMO3D.tar.gz, run this script to setup the third party
# software, symbolic links, and patches. Then build NEMO3D before running.
# *****************************************************************************


# MACOSX's default tar has a limitation, so we use gnutar
if [ `uname` = "Darwin" ]
then
	TAR="gnutar"
else
	TAR="tar"
fi


# Untarball all the *.tar.gz files in build/dist/
echo "#### Untaring ARPACK into the NEMO3D directory"
gzip -dc build/dist/arpack96.tar.gz | $TAR -xvf -

echo "#### Applying ARPACK patch"
gzip -dc build/dist/patch.tar.gz | $TAR -xvf -

echo "#### Untaring PARPACK into the NEMO3D directory"
gzip -dc build/dist/parpack96.tar.gz | $TAR -xvf -

echo "#### Applying PARPACK patch"
gzip -dc build/dist/ppatch.tar.gz | $TAR -xvf -

echo "#### Untaring LAPACK into the NEMO3D directory"
gzip -dc build/dist/lapack3.tar.gz | $TAR -xvf -

echo "#### Untaring XERCES-C++ into the NEMO3D directory"
gzip -dc build/dist/xerces-c-src2_1_0.tar.gz | $TAR -xvf -


# Apply patch fixes
echo "#### Patching incomplete XERCES-C++ 2.1.0 implementation on IRIX platforms"
cp -f build/patches/IRIXPlatformUtils.cpp xerces-c-src2_1_0/src/xercesc/util/Platforms/IRIX/

echo "#### Patching a bug in XERCES-C++ 2.1.0 that occurs with Intel compilers"
cp -f build/patches/Makefile.in xerces-c-src2_1_0/src/xercesc/util/Compilers/
echo "#### Patching a bug in XERCES-C++ 2.1.0 that occurs with CrayXt platforms"
cp -f build/patches/LinuxPlatformUtils.cpp xerces-c-src2_1_0/src/xercesc/util/Platforms/Linux

echo "#### Patching a bug with GNU GCC in MACOSX with PARPACK code with maximum line spacing."
cp -f build/patches/pznaup2.f           ARPACK/PARPACK/SRC/MPI/



# Create symbolic links needed by "eigen"
echo "#### Creating symbolic links in eigen to LAPACK and nml..."
cd eigen
ln -s ../LAPACK LAPACK
ln -s ../nml nml
cd ..


echo "======================================================================"
echo "Edit 'NEMO3D/build/make.inc' before running 'NEMO3D/build/Makefile'"
echo "Make sure env variable 'XERCESCROOT' is properly set before building."
echo "  OR use this command:"
echo "  env XERCESCROOT=\"`pwd`/../xerces-c-src2_1_0\" make"
echo "======================================================================"
echo ""
