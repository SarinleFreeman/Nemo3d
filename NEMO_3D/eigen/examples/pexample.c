/*****************************************************************************
The <namespace>.
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
$Header: /repo/eigen/examples/pexample.c,v 1.3 2003/10/08 16:15:39 hook Exp $
*****************************************************************************/

#include<lanczos.h>
#include<mpi.h>

/* Link to PARPACK Fortran 77 [simple] matrix-vector multiplier.	*/
/* See pzsmvmul.f for a description.					*/
#ifdef	FORTRAN_UNDERSCORE
#define av            av_
#endif/*FORTRAN_UNDERSCORE	*/

#ifdef	FORTRAN_UNDERSCORE2
#define av            av__
#endif/*FORTRAN_UNDERSCORE2	*/

#ifdef	FORTRAN_ALLCAPS
#define av            AV
#endif/*FORTRAN_ALLCAPS		*/

#ifdef  __cplusplus
extern "C" {
#endif/*__cplusplus		*/
void (av)(const int* pC, const int* pL, const int* pN,
    nml_dcscalar* pB, const nml_dcscalar* pV, nml_dcscalar* pW);
#ifdef	__cplusplus
  }
#endif/*__cplusplus		*/

void matmul(const int* argument[], nml_dcscalar* pW, const nml_dcscalar* pV) {
  const int	communicator	= *argument[0];
  const int	processors	= *argument[2];
  const int	blocks		= *argument[3];
  const int	nx		= processors*blocks;
  const int	length		= nx*nx/processors;
  nml_dcvector	*pB		= nml_dcv_new(length);
  (av)(&communicator, &length, &nx, pB, pV, pW);
  nml_dcv_delete(pB);
  }

/* Diagonal matrix-vector multiplier.					*/
/*
void matmul(const int* argument[], nml_dcscalar* pW, const nml_dcscalar* pV) {
  const int	rank		= *argument[1];
  const	int	processors	= *argument[2];
  const int	blocks		= *argument[3];
  const	int	nx		= processors*blocks;
  const	int	n		= nx*nx;
  const
  nml_extent	length		= n/processors;
  nml_offset	j		= 0;
  for (j = 0; j < length; ++j) {
    double	d_j = (n - rank*length) - j;
    pW[j].r = d_j*pV[j].r;
    pW[j].i = d_j*pV[j].i;
    }
  }
*/
double start_time;
double start_time2;
double start_time3;
double time_lanczos_iter;
double time_lanczos_iter_all;
double time_lanczos_setup_once;
double time_barrier;
double time_bcast;


int mpi_n3d_id	 = 0;
int mpi_n3d_masterid = 0;
int mpi_n3d_numprocs = 0;

int
main(int argc, char* argv[]) {
  const int	communicator = MPI_COMM_WORLD;
  const int	initialized = (MPI_SUCCESS == MPI_Init(&argc, &argv));
  const int	processors
    = (MPI_SUCCESS == MPI_Comm_size(communicator, &mpi_n3d_numprocs))?
    mpi_n3d_numprocs: mpi_n3d_numprocs;
  const int	rank
    = (MPI_SUCCESS == MPI_Comm_rank(communicator, &mpi_n3d_id))? mpi_n3d_id: mpi_n3d_id;
  const int	blocks	   = (1 < argc)? atoi(argv[1]): 2;
  const int	nx = processors*blocks;	/* block matrix size		*/
  const	int	n = nx*nx;		/* 	problem size		*/
  const int	values = (2 < argc)? atoi(argv[2]): 18;
  const int	maxitr = (3 < argc)? atoi(argv[3]): 1024;

  int	verbose = 0;		/* shutup verbose diagnostics		*/
  FILE	*fp_log = stdout;	/* direct diagnostics to standard output*/
  FILE	*fp_instant = fopen("eig_instant", "a");
  FILE	*fp_trace   = fopen("eig_trace",   "a");
  const	int*	argument[] = {&communicator, &rank, &processors, &blocks};
				/* argument list for matmul		*/
  const double	tolerance = 1.0e-8;	/* convergence tolerance	*/
  const double	resolution = 100.0*tolerance;
				/* eigenvalue separation minimum	*/
  const double	emax = 8.0;	/* maximum [eigenvalue] energy		*/
  const double	emin = 4.0;	/* minimum [eigenvalue] energy		*/
  const
  nml_extent	imax = maxitr;	/* maximum number of iterations + 1	*/
  const
  int ConvCheckSkipRate  = 10;	/* Convergence Check Skip Rate		*/
/*int ConvCheckSkipRate  = 100;*//* Convergence Check Skip Rate		*/
  const
/*int ConvCheckStartIter = 100;*//* Convergence Check Start Iteration	*/
  int ConvCheckStartIter = 10;	/* Convergence Check Start Iteration	*/
  const
  nml_extent requested = values;/* number of eigenvalues requested	*/
  const
  nml_extent	length = n/processors;	/* segment length		*/
  nml_dcvector	*q_n1	= nml_dcv_fill(nml_dcv_new(length),
				       nml_dcmplx(0.0, 0.0));
  nml_dcvector	*q_n	= nml_dcv_new(length);
  nml_dcvector	*r_n	= nml_dcv_fill(nml_dcv_new(length),
				       nml_dcmplx(1.0, 0.0));
  nml_extent	iterations = 0;	/* actual number of iterations taken	*/
  nml_d3bands	*T = nml_d3_new(imax);		/* tridiagonal matrix	*/
  nml_dvector*	value = nml_dv_new(imax);	/* eigenvalue vector	*/

  nml_extent	eigenvalues  = eigenvaluesLanczos(value, T, &iterations,
	r_n, q_n, q_n1, length, requested,
	ConvCheckStartIter, ConvCheckSkipRate, imax,
	emin, emax, tolerance, resolution, matmul, argument,
	fp_trace, fp_instant, fp_log, verbose);

  if (initialized)
    MPI_Finalize();

  if (mpi_n3d_masterid == mpi_n3d_id) {
    fprintf(stdout, "value =");
    { nml_offset	eigenvalue = 0;
      for (eigenvalue = 0; eigenvalue < eigenvalues; ++eigenvalue) {
	if (0 == eigenvalue%6)
	  fprintf(stdout, "\n");
        fprintf(stdout, "%12.9f", value[eigenvalue]);
	}
      }
    fprintf(stdout, "\n\n");
    fprintf(stdout, "pexample\n\n");
    fprintf(stdout, "Size of the matrix is   %d\n", n);
    fprintf(stdout, "The number of processors is   %d\n", processors);
    fprintf(stdout,
	"The number of Ritz values requested is   %d\n", requested);
    fprintf(stdout,
	"The number of converged Ritz values is   %d\n", eigenvalues);
    fprintf(stdout,
	"The number of Lanczos update iterations taken is   %d\n", iterations);
    fprintf(stdout, "The maximum number of iterations is   %d\n", imax);
    fprintf(stdout, "The convergence criterion is   %g\n", tolerance);
    }

  /*
  if (mpi_n3d_masterid == mpi_n3d_id) {
    nml_offset	iteration = 0;
    for (iteration = 0; iteration < iterations; ++iteration)
      fprintf(stdout, "%i\t%f\t%f\n",
	iteration, T[0][iteration], T[1][iteration]);
    }
  */

  nml_dv_delete(value);	/* eigenvalue vector	*/
  nml_d3_delete(T);	/* tridiagonal matrix	*/
  nml_dcv_delete(r_n);	/* r_{n} = q_{n+1}*beta_{n}			*/
  nml_dcv_delete(q_n);	/* current  complex Krylov vector q_{n}		*/
  nml_dcv_delete(q_n1);	/* previous complex Krylov vector q_{n-1}	*/
  
  return 0;
  }

