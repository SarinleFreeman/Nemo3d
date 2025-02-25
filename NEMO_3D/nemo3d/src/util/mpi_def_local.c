/*****************************************************************************
The Jet Propulsion Laboratory (JPL) NanoElectronicMOdeling-3D package.
Copyright (C) 2002 California Institute of Technology (Caltech)

This application is free software, which you can redistribute and/or modify
under the terms of the GNU Lesser General Public License as published by the
Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this library; see the file COPYING. If not, write to the
Free Software Foundation, Inc.,
59 Temple Place, Suite 330,
Boston, MA  02111-1307  USA

For additional information, please contact
  Gerhard Klimeck (gekco@jpl.nasa.gov)
  Fabiano Oyafuso (fabiano@jpl.nasa.gov)

Written by:  Chris Bowen
             Gerhard Klimeck
             Fabiano Oyafuso
             Seungwon Lee
             Olga Lazarenkova
             Hook Hua

This product includes software developed by the Apache Software Foundation
(http://www.apache.org/).

*****************************************************************************
$Header: /repo/nemo3d/src/util/mpi_def_local.c,v 1.6 2003/10/08 16:18:46 hook
Exp $
*****************************************************************************/

#include "mpi_def_local.h"

#ifdef MPI3d
#include "run3d_mpi.h"
#else
#include "mpi_fake.h"
#endif

int mpi_n3d_id = 0;
int mpi_n3d_numprocs = 1;
int mpi_n3d_masterid = 0;

// -------------------------------------------------------
// stuff that may be used later
// -------------------------------------------------------

mpiCommClass mpiComm;

mpiCommClass::mpiCommClass() {
  this->mpi_ID = 0;
  this->mpi_Nproc = 1;
  this->mpi_Root = 0;
}

mpiCommClass::mpiCommClass(bool periodicX) {
  mpi_Root = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_Nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_ID);

  mpi_left = (mpi_ID ? mpi_ID - 1 : (periodicX ? mpi_Nproc - 1 : -1));
  mpi_right = (mpi_ID < mpi_Nproc - 1 ? mpi_ID + 1 : (periodicX ? 0 : -1));
}
