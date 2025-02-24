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
$Header: /repo/eigen/include/mpi_fake.h,v 1.4 2003/10/08 16:15:46 hook Exp $
*****************************************************************************/


/*  These are the definitions to use our fake version of MPI in mpi_stub.c 
 *  Guaranteed only to make PGAPack compile and link without MPI, but should
 *  work from user C and fortran programs.
 *
 *  Some of these are from mpich's mpi.h, others are custom.
 *
 *  Author: Brian P. Walenz
 */

#ifndef MPI_FAKE_H
#define MPI_FAKE_H 1

#if (defined MPI3d && !defined FAKE_MPI)
#include <mpi.h>
#else

#include <stdio.h>

typedef void *  MPI_Comm;
typedef void *  MPI_Datatype;
typedef long    MPI_Aint;
typedef int     MPI_Op;

typedef struct {
    int     MPI_SOURCE;
    int     MPI_TAG;
    int     MPI_ERROR;
} MPI_Status;



#define MPI_BYTE            (void *)NULL
#define MPI_CHAR            (void *)NULL
#define MPI_DOUBLE          (void *)NULL
#define MPI_DOUBLE_COMPLEX  (void *)NULL
#define MPI_FLOAT           (void *)NULL
#define MPI_INT             (void *)NULL
#define MPI_INTEGER        (void *)NULL
#define MPI_LONG            (void *)NULL
#define MPI_LONG_DOUBLE     (void *)NULL
#define MPI_PACKED          (void *)NULL
#define MPI_SHORT           (void *)NULL
#define MPI_UNSIGNED_CHAR   (void *)NULL
#define MPI_UNSIGNED        (void *)NULL
#define MPI_UNSIGNED_LONG   (void *)NULL
#define MPI_UNSIGNED_SHORT  (void *)NULL

#define MPI_COMM_WORLD      (void *)NULL
#define MPI_COMM_SELF       (void *)NULL

#define MPI_BOTTOM          (void *)0

#define MPI_PROC_NULL       (-1)
#define MPI_ANY_SOURCE      (-2)
#define MPI_ANY_TAG         (-1)
#define MPI_SUM             (MPI_Op)(102)

#define MPI_MIN            (0)  /* type MPI_Op */
#define MPI_MAX            (0)  /* type MPI_Op */

/*  Declare prototypes for the MPI functions.  */
int MPI_Address(void *, MPI_Aint *);
int MPI_Bcast(void *, int, MPI_Datatype, int, MPI_Comm);
int MPI_Barrier(MPI_Comm);
int MPI_Comm_dup(MPI_Comm, MPI_Comm *);
int MPI_Comm_free(MPI_Comm *);
int MPI_Comm_rank(MPI_Comm, int *);
int MPI_Comm_size(MPI_Comm, int *);
int MPI_Finalize(void);
int MPI_Init(int *, char ***);
int MPI_Initialized(int *);
int MPI_Probe(int, int, MPI_Comm, MPI_Status *);
int MPI_Send(void *, int, MPI_Datatype, int, int, MPI_Comm);
int MPI_Recv(void *, int, MPI_Datatype, int, int, MPI_Comm, MPI_Status *);
int MPI_Sendrecv(void *, int, MPI_Datatype, int, int, void *, int,
	         MPI_Datatype, int, int, MPI_Comm, MPI_Status *);
int MPI_Type_commit(MPI_Datatype *);
int MPI_Type_free(MPI_Datatype *);
int MPI_Type_struct(int, int *, MPI_Aint *, MPI_Datatype *, MPI_Datatype *);
int MPI_Allreduce(void* buf, void *variable, int count, MPI_Datatype datatype, MPI_Op Operation, MPI_Comm comm);
int MPI_Get_count(MPI_Status *, MPI_Datatype, int *);
double MPI_Wtime(void);

/* #if defined(FORTRANCAP) */
#if defined(FORTRAN_ALLCAPS)
#define mpi_address_      MPI_ADDRESS
#define mpi_bcast_        MPI_BCAST
#define mpi_comm_dup_     MPI_COMM_DUP
#define mpi_comm_free_    MPI_COMM_FREE
#define mpi_comm_rank_    MPI_COMM_RANK
#define mpi_comm_size_    MPI_COMM_SIZE
#define mpi_finalize_     MPI_FINALIZE
#define mpi_init_         MPI_INIT
#define mpi_initialized_  MPI_INITIALIZED
#define mpi_probe_        MPI_PROBE
#define mpi_send_         MPI_SEND
#define mpi_recv_         MPI_RECV
#define mpi_sendrecv_     MPI_SENDRECV
#define mpi_type_commit_  MPI_TYPE_COMMIT
#define mpi_type_free_    MPI_TYPE_FREE
#define mpi_type_struct_  MPI_TYPE_STRUCT
/* #elif !defined(FORTRANUNDERSCORE) */
#elif !defined(FORTRAN_UNDERSCORE)
#define mpi_address_      mpi_address
#define mpi_bcast_        mpi_bcast
#define mpi_comm_dup_     mpi_comm_dup
#define mpi_comm_free_    mpi_comm_free
#define mpi_comm_rank_    mpi_comm_rank
#define mpi_comm_size_    mpi_comm_size
#define mpi_finalize_     mpi_finalize
#define mpi_init_         mpi_init
#define mpi_initialized_  mpi_initialized
#define mpi_probe_        mpi_probe
#define mpi_send_         mpi_send
#define mpi_recv_         mpi_recv
#define mpi_sendrecv_     mpi_sendrecv
#define mpi_type_commit_  mpi_type_commit
#define mpi_type_free_    mpi_type_free
#define mpi_type_struct_  mpi_type_struct
#endif

void mpi_address_(void **location, MPI_Aint *address);
void mpi_bcast_(void **n, int *com, MPI_Datatype *dt, int *r, MPI_Comm *c, int *ie);
void mpi_comm_dup_(MPI_Comm *comm, MPI_Comm **newcomm, int *ie);
void mpi_comm_free_(MPI_Comm **comm, int *ie);
void mpi_comm_rank_(MPI_Comm *comm, int *rank, int *ie);
void mpi_comm_size_(MPI_Comm *comm, int *size, int *ie);
void mpi_finalize_(int *ie);
void mpi_init_(int *ie);
void mpi_initialized_(int *flag, int *ie);
void mpi_probe_(int *source, int *tag, MPI_Comm *comm, MPI_Status *status,
		int *ie);
void mpi_send_(void *buf, int *count, MPI_Datatype *datatype, int *dest,
	       int *tag, MPI_Comm *comm, int *ie);
void mpi_recv_(void *buf, int *count, MPI_Datatype *datatype, int *source,
	       int *tag, MPI_Comm *comm, MPI_Status *status, int *ie);
void mpi_sendrecv_(void *sendbuf, int *sendcount, MPI_Datatype *sendtype,
		   int *dest, int *sendtag, void *recvbuf, int *recvcount,
		   MPI_Datatype *recvtype, int *source, int *recvtag,
		   MPI_Comm *comm, MPI_Status *status, int *ie);
void mpi_type_commit_(MPI_Datatype **datatype, int *ie);
void mpi_type_free_(MPI_Datatype **datatype, int *ie);
void mpi_type_struct_(int *count, int **array_of_blocklengths,
		      MPI_Aint **array_of_displacements,
		      MPI_Datatype **array_of_types, MPI_Datatype **newtype, int *ie);

#endif /* MPI3d */

#endif /* MPI_FAKE_H */

