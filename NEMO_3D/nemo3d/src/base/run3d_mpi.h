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
$Header: /repo/nemo3d/src/base/run3d_mpi.h,v 1.7 2007/09/19 20:53:26 stevenclark
Exp $
*****************************************************************************/

#ifndef RUN3D_MPI_H
#define RUN3D_MPI_H

#include "is_def.h"

#if (defined MPI3d && !defined FAKE_MPI)

#ifdef UNSET_SEEK
#undef SEEK_SET
#undef SEEK_CUR
#undef SEEK_END
#endif
#include <mpi.h>
#define Allreduce_MPI_sp(pa, pb, c, d, e, f) MPI_Allreduce(pa, pb, c, d, e, f)
#ifdef UNSET_SEEK
/* #include <stdio.h> */
/* The possibilities for the third argument to `fseek'.
   These values should not be changed.  */
#ifndef SEEK_SET
#define SEEK_SET 0 /* Seek from beginning of file.  */
#endif
#ifndef SEEK_CUR
#define SEEK_CUR 1 /* Seek from current position.  */
#endif
#ifndef SEEK_END
#define SEEK_END 2 /* Seek from end of file.  */
#endif
#endif

#else

#include "mpi_fake.h"
#ifdef NoGnuStatementExpre
void Allreduce_MPI_sp(double *pa, double *pb, int count, MPI_Datatype datatype,
                      MPI_Op Operation, MPI_Comm comm);
#else
#define Allreduce_MPI_sp(pa, pb, c, d, e, f)                                   \
  ({                                                                           \
    int iLOC;                                                                  \
    for (iLOC = 0; iLOC < (c); iLOC++)                                         \
      (pb)[iLOC] = (pa)[iLOC];                                                 \
  })
#endif

#endif /* (defined MPI3d && !defined FAKE_MPI) */

/* #ifdef MPI3d */
_IS_DEF MPI_Status mpi_n3d_status;
/* #endif */

_IS_DEF int mpi_n3d_numprocs;
_IS_DEF int mpi_n3d_id;
_IS_DEF int mpi_n3d_masterid;

#ifdef MPI_TIMING
_IS_DEF double start_time;
_IS_DEF double start_time2;
_IS_DEF double start_time3;
_IS_DEF double time_lanczos_iter;
_IS_DEF double time_lanczos_iter_all;
_IS_DEF double time_lanczos_setup_once;
_IS_DEF double time_barrier;
_IS_DEF double time_bcast;
#include "MPI_Timing.h"
#endif /* MPI_TIMING */

#define MPI_rc_distr 23001

#endif /* RUN3D_MPI_H */
