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
# $Header: /repo/NEMO_3D/CVS_nemo3d_setup.sh,v 1.8 2006/12/18 18:01:41 gekco Exp $
# *****************************************************************************


# Check the existance of $USER login variable so we can use it to log into the CVS repository to checkout the files.
if [ -z $USER ]
then
    echo ""
    echo "NOTICE: The \$USER variable is not defined. It is used to log into the CVS repository.\n"
    echo ""
    echo "Usage: env USER=<your_login_on_max.ecn.purdue.edu> ./nemo3d_setup.sh"
    echo ""
    exit 1
else
    echo "#### Getting files from CVS..."
    cvs -d :ext:$USER@max.ecn.purdue.edu:/repo co build nml eigen XMLtoC nemo3d viz PostNEMO3D TRACEMIN
#    cvs -d :ext:$USER@max.ecn.purdue.edu:/repo co build
#    cvs -d :ext:$USER@max.ecn.purdue.edu:/repo co nml
#    cvs -d :ext:$USER@max.ecn.purdue.edu:/repo co eigen
#    cvs -d :ext:$USER@max.ecn.purdue.edu:/repo co XMLtoC
#    cvs -d :ext:$USER@max.ecn.purdue.edu:/repo co nemo3d
#    cvs -d :ext:$USER@max.ecn.purdue.edu:/repo co viz
#    cvs -d :ext:$USER@max.ecn.purdue.edu:/repo co PostNEMO3D
fi


# After checking out NEMO3D components from CVS, run the "nemo3d_setup.sh" script to
# setup the third party packages (Xerces C++, LAPACK, BLAS, ARPACK, PARPACK),
# symbolic links, and patches.
./nemo3d_setup.sh

