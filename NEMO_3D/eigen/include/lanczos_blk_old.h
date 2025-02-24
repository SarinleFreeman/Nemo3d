#ifndef	LANCZOS_BLK_H
#define	LANCZOS_BLK_H 1

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
$Header: /repo/eigen/include/lanczos_blk_old.h,v 1.1 2006/10/11 18:45:49 rrahman Exp $
*****************************************************************************/


#include "io.h"
#include "lanczos.h"
#include "mpi_local.h"
#include "nml_dcvector.h"
#include "nml_ivector.h"
#include "nml_dcmatrix.h"
#include "nml_d3bands.h"
#include "matrix.h"
#include<float.h>

int solve_pentamatrix_vector(
    int iteration,
    nml_dvector *aa,
    nml_dcvector *b,
    nml_dvector *g,
    nml_dcvector *x,
    nml_dcvector *rhs,
    nml_dscalar shift
    );
/*  Find a few eigenvalues in a complex hermitian matrix
    using Lanczos' algorithm.						*/
nml_extent eigenvaluesBlockLanczos(	/* number of eigenvalues actually found	*/
    nml_dvector	 *value,	/*   out real eigenvalue vector		*/
    nml_dvector   *ALPHA,
    nml_dvector   *GAMMA,
    nml_dcvector   *BETA,
    nml_extent	 *pIterations,	/*   out actual number of iterations	*/
    nml_dcvector *r_n,		/* inout r_{n} = q_{n+1}*beta_{n}	*/
    nml_dcvector *r2_n,		/* inout r_{n} = q_{n+1}*beta_{n}	*/
    nml_dcvector *q_n,		/* inout current  complex Lanczos vector*/
    nml_dcvector *q_n1,		/* inout previous complex Lanczos vector*/
    nml_extent	 length,	/* in extent of vectors r_n, q_n & q_n1	*/
    nml_extent	 requested1,	/* in requested number of eigenvalues	*/
    nml_extent	 requested2,	/* in requested number of eigenvalues	*/
    int	ConvCheckStartIter,	/* in convergence check start iteration	*/
    int	ConvCheckSkipRate,	/* in convergence check skip rate	*/
    nml_extent	 imax,		/* in maximum number of iterations + 1	*/
    nml_dscalar	 emin1,		/* in eigenvalue search lower bound	*/
    nml_dscalar	 emax1,		/* in eigenvalue search upper bound	*/
    nml_dscalar	 emin2,		/* in eigenvalue search lower bound	*/
    nml_dscalar	 emax2,		/* in eigenvalue search upper bound	*/
    nml_dscalar	 tolerance,	/* in eigenvalue convergence tolerance	*/
    nml_dscalar	 resolution,	/* in eigenvalue separation minimum	*/
    void (*matmul)(const int**, nml_dcscalar*, const nml_dcscalar*),
    const int*	 argument[],	/* matrix-vector multiply argument list	*/
    FILE	 *fp_trace,	/* in   trace log text file pointer	*/
    FILE	 *fp_log,	/* in message log text file pointer	*/
    int		 verbose	/* in verbose diagnostic messages	*/
    );

/*  Compute the eigenvectors corresponding to the eigenvalues found
    in a complex hermitian matrix using Lanczos' algorithm.		*/
nml_extent eigenvectorsBlockLanczos(	/* number of eigenvalues actually found	*/
    nml_dcvector *cvalue,	/*   out complex eigenvalues		*/
    nml_dcmatrix *vector,	/*   out complex eigenvectors		*/
    const
    nml_dvector  *value,	/* in    real    eigenvalues		*/
    nml_extent	 eigenvalues,	/* in	 number of eigenvalues		*/
    nml_extent	 iterations,	/* in    number of iterations		*/
    nml_dvector   *ALPHA,
    nml_dvector   *GAMMA,
    nml_dcvector   *BETA,
    nml_dcscalar *r_n,		/* inout r_{n} = q_{n+1}*beta_{n}	*/
    nml_dcvector *r2_n,		/* inout r_{n} = q_{n+1}*beta_{n}	*/
    nml_dcscalar *q_n,		/* inout current  complex Lanczos vector*/
    nml_dcscalar *q_n1,		/* inout previous complex Lanczos vector*/
    nml_extent	 length,	/* in extent of vectors r_n, q_n & q_n1	*/
    void (*matmul)(const int**, nml_dcscalar*, const nml_dcscalar*),
    const int*	 argument[],	/* matrix-vector multiply argument list	*/
    nml_dscalar  tolerance,     /* in eigenvalue convergence tolerance  */
    FILE	 *fp_trace,	/* in   trace log text file pointer	*/
    FILE	 *fp_log,	/* in message log text file pointer	*/
    int		 verbose	/* in verbose diagnostic messages	*/
    );

/* Deprecated -- Calculates requested eigenvalues and eigenvectors
   of the Hamiltonian in the energy range [emin:emax].
   imax is the maximum number of iterations + 1.			*/
int eigensolve_blk_lanczos(		/* number of eigenvalues actually found	*/
    nml_dvector	 **pValue,	/* out real eigenvalue vector		*/
    int	ConvCheckStartIter,	/* in convergence check start iteration	*/
    int	ConvCheckSkipRate,	/* in convergence check skip rate	*/
    int		 imax,		/* in maximum number of iterations + 1	*/
    int		 requested1,	/* in requested number of eigenvalues	*/
    int		 requested2,	/* in requested number of eigenvalues	*/
    nml_extent	 length,	/* in extent of vectors r, q and q_n1	*/
    int*	 argument[],	/* matrix-vector multiply argument list	*/
    nml_dscalar	 emin1,		/* in eigenvalue search lower bound	*/
    nml_dscalar	 emax1,		/* in eigenvalue search upper bound	*/
    nml_dscalar	 emin2,		/* in eigenvalue search lower bound	*/
    nml_dscalar	 emax2,		/* in eigenvalue search upper bound	*/
    nml_dscalar	 tolerance,	/* in eigenvalue convergence tolerance	*/
    int		 verbose,	/* in verbose diagnostics messages	*/
    const char	 *fn_trace,	/* in eigenvalue log text file name	*/
    const char	 *fn_pendiag,	/* in pendiagonal matrix binary file name */
    const char	 *fn_eigval,	/* in eigenvalue  vector binary file name */
    nml_dcvector *r_n,          /* inout r_{n} = q_{n+1}*beta_{n}       */
    const nml_dcvector *r_0,    /* in    initial guess vector           */
    nml_dcvector *q_n,          /* inout current  complex Lanczos vector*/
    nml_dcvector *q_n1,         /* inout previous complex Lanczos vector*/
    nml_dcvector **pCValue,     /*   out complex eigenvalue  vector     */
    nml_dcmatrix **pCVector,    /*   out complex eigenvector matrix     */
    void (*)(const int**, nml_dcscalar*, const nml_dcscalar*),
    int		 restoreState	/* in read pendiagonal matrix from file	*/
    );

#endif /*LANCZOS_BLK_H		*/

