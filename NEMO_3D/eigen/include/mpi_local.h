#ifndef	MPI_LOCAL_H
#define	MPI_LOCAL_H

/*****************************************************************************
The JPL Parallel Eigensolvers package.
Copyright (C) 2002 California Institute of Technology (Caltech)

This file is part of
The Jet Propulsion Laboratory (JPL) Parallel Eigensolvers package.

This library is free software, which you can redistribute and/or modify
under the terms of the GNU Lesser General Public License as published by the
Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this library; see the file COPYING. If not, write to the
Free Software Foundation, Inc.,
59 Temple Place, Suite 330,
Boston, MA  02111-1307  USA

For additional information, please contact
  Gerhard Klimeck (gekco@jpl.nasa.gov)
  Fabiano Oyafuso (fabiano@jpl.nasa.gov)
  E. Robert Tisdale (E.Robert.Tisdale@jpl.nasa.gov)

Written by: Chris Bowen
            Gerhard Klimeck
            Fabiano Oyafuso
            E. Robert Tisdale

*****************************************************************************
$Header: /repo/eigen/include/mpi_local.h,v 1.8 2007/09/19 20:53:25 stevenclark Exp $
*****************************************************************************/


#if (defined MPI3d && !defined FAKE_MPI)

#ifdef UNSET_SEEK
#undef SEEK_SET
#undef SEEK_CUR
#undef SEEK_END
#endif
#include <mpi.h>
#define Allreduce_MPI_sp(pa,pb,c,d,e,f) MPI_Allreduce(pa,pb,c,d,e,f)
#ifdef UNSET_SEEK
/* #include <stdio.h> */
/* The possibilities for the third argument to `fseek'.
   These values should not be changed.  */
#ifndef SEEK_SET
#define SEEK_SET	0	/* Seek from beginning of file.  */
#endif
#ifndef SEEK_CUR
#define SEEK_CUR	1	/* Seek from current position.  */
#endif
#ifndef SEEK_END
#define SEEK_END	2	/* Seek from end of file.  */
#endif
#endif

#else /* not (defined MPI3d && !defined FAKE_MPI) */

#include "mpi_fake.h"

#ifdef NoGnuStatementExpre
void Allreduce_MPI_sp(double *pa,double *pb, int count, MPI_Datatype datatype, MPI_Op Operation, MPI_Comm comm);
#else
#define Allreduce_MPI_sp(pa,pb,c,d,e,f) \
        ({int iLOC; for(iLOC=0; iLOC<(c); iLOC++) (pb)[iLOC]=(pa)[iLOC];})
#endif



#endif /* (defined MPI3d && !defined FAKE_MPI) */

extern int mpi_n3d_id, mpi_n3d_masterid, mpi_n3d_numprocs;

#ifdef  MPI_TIMING

extern double start_time;
extern double start_time2;
extern double start_time3;
extern double time_lanczos_iter;
extern double time_lanczos_iter_all;
extern double time_lanczos_setup_once;
extern double time_barrier;
extern double time_bcast;


/* NOTE: "NEMO3D/nemo3d/src/base/MPI_Timing.h" checks to see if this header's
   MPI_LOCAL_H has been defined. If not, it defines these MPI tic-tocs. */

#ifndef MPI_TIME_INIT
#define MPI_TIME_INIT(T)        double T = 0.0
#endif

#ifndef MPI_TIME_TIC
#define MPI_TIME_TIC(T)         double T = MPI_Wtime()
#endif

#ifndef MPI_TIC
#define MPI_TIC(T)              T = MPI_Wtime()
#endif

#ifndef MPI_TOC
#define MPI_TOC(T,T0)           T += (MPI_Wtime() - T0)
#endif

#else /* MPI_TIMING not defined */

#ifndef MPI_TIME_INIT
#define MPI_TIME_INIT(T)
#endif

#ifndef MPI_TIME_TIC
#define MPI_TIME_TIC(T)         double T = 0.0
#endif

#ifndef MPI_TIC
#define MPI_TIC(T)
#endif

#ifndef MPI_TOC
#define MPI_TOC(T,T0)
#endif

#endif /* MPI_TIMING */


#endif /* MPI_LOCAL_H */

