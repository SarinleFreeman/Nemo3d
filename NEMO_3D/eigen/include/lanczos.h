#ifndef	LANCZOS_H
#define	LANCZOS_H 1

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
$Header: /repo/eigen/include/lanczos.h,v 1.10 2007/09/19 20:53:24 stevenclark Exp $
*****************************************************************************/


#include "io.h"
#include "mpi_local.h"
#include "nml_dcvector.h"
#include "nml_ivector.h"
#include "nml_dcmatrix.h"
#include "nml_d3bands.h"
#include "matrix.h"
#include<float.h>

/*! \file lanczos.h
 *  \brief Calculates requested eigenvalues and eigenvectors
 *  of the Hamiltonian in the energy range [emin:emax] using Lanczos
 *  eigensolver algorithm
 */

#define MTRX_SCALE 1.0  /* this was HAMILTONIAN_SCALE_VALUE */
/*! \define Matrix Scale value, set to 1.0
 */

#ifdef	FORTRAN_UNDERSCORE
#define dstebz            dstebz_
#endif/*FORTRAN_UNDERSCORE	*/

#ifdef	FORTRAN_UNDERSCORE2
#define dstebz            dstebz__
#endif/*FORTRAN_UNDERSCORE2	*/

#ifdef	FORTRAN_ALLCAPS
#define dstebz            DSTEBZ
#endif/*FORTRAN_ALLCAPS		*/

#ifdef	__cplusplus
extern "C" {
#endif/*__cplusplus	*/
void (dstebz)(const char*, const char*, const int*,
    const double*, const double*, const int*, const int*,
    const double*, const nml_dvector*, const nml_dvector*, int*, int*,
    nml_dvector*, nml_ivector*, nml_ivector*, nml_dvector*, nml_ivector*, int*);
#ifdef	__cplusplus
  }
#endif/*__cplusplus	*/

#if (defined MPI3d && !defined FAKE_MPI)
#include <mpi.h>
#define Allreduce_MPI_sp(pa,pb,c,d,e,f) MPI_Allreduce(pa,pb,c,d,e,f)
#else /* not (defined MPI3d && !defined FAKE_MPI) */
#ifdef NoGnuStatementExpre
void Allreduce_MPI_sp(double *pa,double *pb, int count, MPI_Datatype datatype, MPI_Op Operation, MPI_Comm comm);
#else
#define Allreduce_MPI_sp(pa,pb,c,d,e,f) \
        ({int iLOC; for(iLOC=0; iLOC<(c); iLOC++) (pb)[iLOC]=(pa)[iLOC];})
#endif
#endif

void rsort(
    nml_dvector *v, 
    nml_extent n   
    ) ;
/*!< This function is deprecated. See rsort_twovec */

void rsort_twovec( 
		  double *v,
		  int *second_vec , 
		  int n
		  ) ;
/*!< Sort the vector of eigenvalues in ascending order
 * @param v INOUT: Eigenvalue array
 * @param second_vec INOUT: Iteration number for each eigenvaluess
 * @param n IN: Number of eigenvalues
 * @return void
 */

nml_dvector* copyNmlDvector(nml_dvector* pv, nml_extent vOffsetStart,
    const nml_dvector* pw, nml_extent wOffsetStart, nml_extent copyLength) ;
/*!< Copies nml_dvector
 * @param pv OUT: Destination vector, where copied vector is stored
 * @param vOffsetStart IN: Offset of destination vector
 * @param pw IN: Source vector, the vector to be copied
 * @param wOffsetStart IN: Offset of source vector
 * @param copyLength IN: Copy length
 * @return nml_dvector(not considered)
 */

nml_d3bands* copyNmlD3bands(nml_d3bands* pT, nml_extent tOffsetStart,
    const nml_d3bands* pU, nml_extent uOffsetStart, nml_extent copyLength) ;
/*!< Copies nml_d3bands(nml_dscalar)
 * @param pT OUT: Destination vector, where copied vector is stored
 * @param tOffsetStart IN: Offset of destination vector
 * @param pU IN: Source vector, the vector to be copied
 * @param uOffsetStart IN: Offset of source vector
 * @param copyLength IN: Copy length
 * @param nml_d3bands(not considered)
 */

nml_dscalar InverseIteration( 
    nml_d3bands *h,
    nml_dvector *x,
    nml_dvector *y,
    nml_dscalar lambda,
    FILE        *fp_log 
    ) ;
/*!< This function is deprecated. See originalInverseIteration */

int originalInverseIteration( 
    nml_d3bands *T,     
    nml_dvector *x,     
    nml_dvector *y,     
    nml_dscalar lambda  
    ) ;
/*!< Return refined eigenvalue estimate using inverse iteration
 * @param T INOUT: Real symmetric tridiagonal matrix
 * @param x INOUT: Workspace vector
 * @param y INOUT: Eigenvector
 * @param lambda IN: Real eigenvalue estimate
 * @return 0
 */

void find_converged_eigval(
    int         iterations, 
    const
    nml_d3bands *T,         
    nml_dvector *eigv_conv,  
    nml_ivector * convIterationNumbers ,
    int         *pEigenvalues,
    int         *n_nosuccess, 
    int         *n_old,    
    nml_dvector *eigv_old, 
    FILE        *fp_trace, 
    int         verbose, 
    nml_dscalar emin,     
    nml_dscalar emax,     
    nml_dscalar tolerance, 
    nml_dscalar resolution, 
    FILE        *fp_log      
    ) ;
/*!< Find the list of converged eigenvalues between check intervals
 * @param iterations IN: Number of current Lanczos iterations
 * @param T IN: Real symmetric tridiagonal matrix
 * @param eigv_conv INOUT: List of actual converged eigenvalues
 * @param convIterationNumbers INOUT: Iteration number when the eigenvalues converged
 * @param pEigenvalues OUT: Number of converged eigenvalues
 * @param n_nosuccess OUT: Number of success
 * @param n_old INOUT: Number of converged eigenvalues on previous check
 * @param eigv_old INOUT: List of eigenvalues on previous check
 * @param fp_trace IN: Pointer to the trace file
 * @param verbose IN: Diagnostic messages (1: enabled, 0: disabled)
 * @param emin IN: Lower bound of eigenvalue range
 * @param emax IN: Upper bound of eigenvalue range
 * @param tolerance IN: Convergence tolerance
 * @param resolution IN: Minimum spacing between eigenvalues 
 * @param fp_log IN: Pointer to the log file
 * @return void
 */

void saveTridiagonalMatrix(
    const
    nml_d3bands *T,     
    int iterations,    
    int processors,   
    FILE *fp_tridiag 
    ) ;
/*!< Save Tridiagonal Matrix
 * @param T IN: Real symmetric tridiagonal matrix to be saved
 * @param iterations IN: Number of final Lanczos iterations
 * @param processors IN: Number of processors
 * @param fp_tridiag INOUT: Pointer to the binary file where the tridiagonal matrix is stored
 * @return void
 */

void restoreTridiagonalMatrix(
    nml_d3bands **pT,  
    nml_extent
    *pIterations,   
    nml_extent
    *pProcessors,  
    FILE *fp_tridiag
    ) ;
/*!< Restore Tridiagonal Matrix
 * @param pT OUT: Real symmetric tridiagonal matrix to be saved
 * @param pIterations OUT: Number of final Lanczos iterations
 * @param pProcessors INOUT: Number of processors
 * @param fp_tridiag INOUT: Pointer to the binary file where the tridiagonal matrix is stored
 * @return void
 */

void saveEigenvalues(
    const
    nml_dvector *value,
    int eigenvalues,
    FILE *fp_eigval   
    ) ;
/*!< Save eigenvalues to file
 * @param value IN: Real eigenvalue vector
 * @param eigenvalues IN: Number of eigenvalues
 * @param fp_eigval INOUT: Pointer to the binary file where the real eigenvalues are stored
 * @return void
 */

void restoreEigenvalues(
    nml_dvector **pValue,
    nml_extent
    *pEigenvalues,
    FILE *fp_eigval
    ) ;
/*!< Restore eigenvalues to file
 * @param pValue OUT: Real eigenvalue vector
 * @param pEigenvalues OUT: Number of eigenvalues
 * @param fp_eigval INOUT: Pointer to the binary file where the real eigenvalues are stored
 * @return void
 */

int solve_rtridiag_ii(
    nml_dvector	*x,
    const
    nml_d3bands	*T,
    const
    nml_dvector	*y	
    );
/*!< Solve x*T' = y for x
 * @param x OUT: Real vector
 * @param T IN: Real symmetric tridiagonal matrix
 * @param y IN: real vector
 * @return 0
 */

nml_extent eigenvaluesLanczos(
    nml_dvector	 *value,
    nml_ivector * convIterationNumbers ,
    nml_d3bands	 *T,
    nml_extent	 *pIterations,
    nml_dcvector *r_n,	
    nml_dcvector *q_n,
    nml_dcvector *q_n1,	
    nml_extent	 length,	
    nml_extent	 requested1,
    nml_extent	 requested2,
    int	ConvCheckStartIter,
    int	ConvCheckSkipRate,
    nml_extent	 imax,	
    nml_dscalar	 emin1,
    nml_dscalar	 emax1,
    nml_dscalar	 emin2,
    nml_dscalar	 emax2,	
    nml_dscalar	 tolerance,
    nml_dscalar	 resolution,
    void (*matmul)(const int**, nml_dcscalar*, const nml_dcscalar*),
    const int*	 argument[],
    FILE	 *fp_trace,
    FILE	 *fp_log,
    int		 verbose	
    );
/*!< Find a few eigenvalues in a complex Hermitian matrix using Lanczos algorithm
 * @param value OUT: Real eigenvalue vector
 * @param convIterationNumbers OUT: Iteration numbers for each converged eigenvalues
 * @param T OUT: Symmetric tridiagonal matrix
 * @param pIterations OUT: Final number of iterations
 * @param r_n INOUT: r_{n} = q_{n+1} * beta_{n}
 * @param q_n INOUT: Current complex Lanczos vector
 * @param q_n1 INOUT: Previous complex Lanczos vector
 * @param length IN: Length of r_n, q_n & q_n1
 * @param requested1 IN: Requested number of eigenvalues in valence band
 * @param requested2 IN: Requested number of eigenvalues in conduction band
 * @param ConvCheckStartIter IN: Iteration number for first convergence check
 * @param ConvCheckSkipRate IN: Convergence check interval
 * @param imax IN: Maximum number of iterations + 1
 * @param emin1 IN: Minimum energy range in valence band
 * @param emax1 IN: Maximum energy range in valence band
 * @param emin2 IN: Minimum energy range in conduction band
 * @param emax2 IN: Maximum energy range in conduction band
 * @param tolerance IN: Convergence tolerance
 * @param resolution IN: Minimum seperation between eigenvalues
 * @param matmul IN: Function pointer to matrix-vector multiply routine
 * @param argument IN: Matrix-vector multiply routine argument list
 * @param fp_trace IN: File pointer to trace log file
 * @param fp_log IN: File pointer tot log file
 * @param verbose IN: Turn on diagnostic message (1: enabled, 0: disabled)
 * @return Number of eigenvalues found
 */

nml_extent eigenvectorsLanczos(
    nml_dcvector *cvalue,
    nml_dcmatrix *vector,
    const
    nml_dvector  *value,
    const
    nml_ivector * convIterationNumbers ,
    nml_extent	 eigenvalues,
    nml_extent	 iterations,
    const
    nml_d3bands	 *T,	
    nml_dcscalar *r_n,
    nml_dcscalar *q_n,
    nml_dcscalar *q_n1,
    nml_extent	 length,
    void (*matmul)(const int**, nml_dcscalar*, const nml_dcscalar*),
    const int*	 argument[],
    nml_dscalar	 tolerance,
    FILE	 *fp_trace,
    FILE	 *fp_log,
    int		 verbose
    );
/*!< Compute the eigenvectors corrersponding to the eigenvalues found in a complex
 * Hermitian matrix using Lanczos algorithm
 * @param cvalue OUT: Complex eigenvalue vector
 * @param vector OUT: Complex set of eigenvectors
 * @param value IN: Real eigenvalues
 * @param convIterationNumbers OUT: Iteration numbers for each converged eigenvalues
 * @param eigenvalues IN: Number of eigenvalues
 * @param iterations IN: Number of iterations
 * @param T IN: Symmetric tridiagonal matrix
 * @param r_n INOUT: r_{n} = q_{n+1} * beta_{n}
 * @param q_n INOUT: Current complex Lanczos vector
 * @param q_n1 INOUT: Previous complex Lanczos vector
 * @param length IN: Length of r_n, q_n & q_n1
 * @param matmul IN: Function pointer to matrix-vector multiply routine
 * @param argument IN: Matrix-vector multiply routine argument list
 * @param tolerance IN: Convergence tolerance
 * @param fp_trace IN: File pointer to trace log file
 * @param fp_log IN: File pointer tot log file
 * @param verbose IN: Turn on diagnostic message (1: enabled, 0: disabled)
 * @return Number of eigenvalues found
 */

int eigensolve_lanczos(		
    nml_dvector	 **pValue,
    int	ConvCheckStartIter,
    int	ConvCheckSkipRate,
    int		 imax,	
    int		 requested1,
    int		 requested2,
    nml_extent	 length,	
    int*	 argument[],	
    nml_dscalar	 emin1,
    nml_dscalar	 emax1,
    nml_dscalar	 emin2,
    nml_dscalar	 emax2,	
    nml_dscalar	 tolerance,
    int		 verbose,	
    const char	 *fn_trace,	
    const char	 *fn_tridiag,
    const char	 *fn_eigval,
    nml_dcvector *r_n,     
    const nml_dcvector *r_0,
    nml_dcvector *q_n,   
    nml_dcvector *q_n1,   
    nml_dcvector **pCValue,
    nml_dcmatrix **pCVector,
    void (*)(const int**, nml_dcscalar*, const nml_dcscalar*),
    int		 restoreState
    );
/*!< Calculates requested eigenvalues and eigenvectors of the Hamiltonian in the energy range [emin:emax]
 * @param pValue OUT: Real eigenvalue vector
 * @param ConvCheckStartIter IN: Iteration number for first convergence check
 * @param ConvCheckSkipRate IN: Convergence check interval
 * @param imax IN: Maximum number of iterations + 1
 * @param requested1 IN: Requested number of eigenvalues in valence band
 * @param requested2 IN: Requested number of eigenvalues in conduction band
 * @param length IN: Length of r_n, q_n & q_n1
 * @param argument IN: Matrix-vector multiply routine argument list
 * @param emin1 IN: Minimum energy range in valence band
 * @param emax1 IN: Maximum energy range in valence band
 * @param emin2 IN: Minimum energy range in conduction band
 * @param emax2 IN: Maximum energy range in conduction band
 * @param tolerance IN: Convergence tolerance
 * @param verbose IN: Turn on diagnostic message (1: enabled, 0: disabled)
 * @param fn_trace IN: File pointer to trace log file
 * @param fn_tridiag IN: The name of the file where tridiagonal matrix is stored
 * @param fn_eigval IN: The name of the file where eigenvalue vector is stored
 * @param r_n INOUT: r_{n} = q_{n+1} * beta_{n}
 * @param r_0 IN: Initial guess vector
 * @param q_n INOUT: Current complex Lanczos vector
 * @param q_n1 INOUT: Previous complex Lanczos vector
 * @param pCValue OUT: Complex eigenvalue vector
 * @param pCVector OUT: Complex eigenvector matrix
 * @param matmul IN: Function pointer to matrix-vector multiply routine
 * @param restoreState IN: Read tridiagonal matrix from file (1: enabled, 0: disabled)
 * @return Number of eigenvalues found
 */

#endif /*LANCZOS_H		*/

