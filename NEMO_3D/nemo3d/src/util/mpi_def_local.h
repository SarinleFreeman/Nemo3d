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
$Header: /repo/nemo3d/src/util/mpi_def_local.h,v 1.5 2003/10/08 16:18:47 hook
Exp $
*****************************************************************************/

#ifndef MPI_DEF_LOCAL_H
#define MPI_DEF_LOCAL_H

/* Global mpi variable definition */

extern int mpi_n3d_id;
extern int mpi_n3d_masterid;
extern int mpi_n3d_numprocs;

// -------------------------------------------------------
// stuff that may be used later
// -------------------------------------------------------

class mpiCommClass {
  int mpi_Root;
  int mpi_ID;
  int mpi_Nproc;

  int mpi_left, mpi_right;

public:
  mpiCommClass();

  mpiCommClass(bool);
};

extern mpiCommClass mpiComm;

#endif /* MPI_DEF_LOCAL_H */
