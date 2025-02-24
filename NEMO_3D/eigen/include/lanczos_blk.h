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
$Header: /repo/eigen/include/lanczos_blk.h,v 1.6 2007/07/25 20:52:52 lee509 Exp $
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
#include<vector>

using namespace std;



/*! \file lanczos_blk.h
 *  \brief Calculates requested eigenvalues and eigenvectors 
 *  of the Hamiltonian in the energy range [emin:emax] using Block Lanczos
 *  eigensolver algorithm
 */  

/* imax is the maximum number of iterations + 1.*/

//#define BLANCZ_MATMULT 
/*! \define This definition utilizes matrix-matrix multiplication (Y=HX).
 *  This definition is turned on by default. If one intends to use
 *  matrix-vector multiplication (y=Hx), please comment out.
 */

int eigensolve_blk_lanczos(	
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
    const char	 *fn_pendiag,
    const char	 *fn_eigval,	
    nml_dcvector *r_n,         
    const nml_dcvector *r_0,  
    nml_dcvector *q_n,       
    nml_dcvector *q_n1,     
    nml_dcvector **pCValue,
    nml_dcmatrix **pCVector,

	#ifndef BLANCZ_MATMULT 
    void (*)(const int**, nml_dcscalar*, const nml_dcscalar*),
	#else
    void (*)(const int**, nml_dcmatrix*, const nml_dcmatrix*),
	#endif

    int		 restoreState,
    int p,
    double deflation_tolerance,
    int eig_dump,
    int num_vb,
    int num_cb,
    int max_no_success,
    int save_BT_band,
    int save_BT_SMF,
    int read_BT_band,
    int deflation_on,
    int conv_method,
    int conv_vec_method		
    );
/*!< Calls Block Lanczos eigensolver routine, used by lanczos_driver in eigsys3d_par.c
 * @param pValue OUT: Real eigenvalue vector
 * @param ConvCheckStartIter IN: First iteration count to check convergence 
 * @param ConvCheckSkipRate IN: Convergence check skip rate
 * @param imax IN: Maximum number of iterations + 1
 * @param requested1 IN: Requested number of eigenvalues in valence band
 * @param requested2 IN: Requested number of eigenvalues in conduction band
 * @param length IN: The length of vectors used in this algorithm (nloc * basis function)
 * @param argument IN: Matrix-matrix multiply or matrix-vector multiply argument list
 * @param emin1 IN: Minimum energy range in valence band
 * @param emax1 IN: Maximum energy range in valence band
 * @param emin2 IN: Minimum energy range in conduction band
 * @param emax2 IN: Maximum energy range in conduction band
 * @param tolerance IN: Convergence tolerance
 * @param verbose IN: Print out debug messages (1: enable 0: disable)
 * @param fn_trace IN: Binary file name for eigenvalue logs
 * @param fn_pendiag IN: Binary file name for tridiagonal matrix
 * @param fn_eigval IN: Binary file name for eigenvalue vector
 * @param r_n INOUT: r_{n} = q_{n+1} * beta_{n}
 * @param r_0 IN: Randomly generated vector for initial guess
 * @param q_n INOUT: Current complex Lanczos vector
 * @param q_n1 INOUT: Previous complex Lanczos vector
 * @param pCValue OUT: Current complex eigenvalue vector
 * @param pCVector OUT: Current complex eigenvector matrix
 * @param matmul IN: Matrix-vector or matrix matrix multiply routine (Function pointer)
 * @param restoreState IN: Read Tridiagonal matrix from file (1: enable, 0: disable)
 * @param p IN: Block size
 * @param deflation_tolerance IN: Tolerance when deflation mode is enabled
 * @param eig_dump IN: Dump eigenvalues (1: enable, 0: disable)
 * @param num_vb IN: Number of requested eigenvalues in valence band
 * @param num_cb IN: Number of requested eigenvalues in conduction band
 * @param max_no_success IN: Maximum number of iterations to abort without any eigenvalues
 * @param save_BT_band IN: Save block tridiagonal matrix in band matrix format (1: enable, 0: disable)
 * @param save_BT_SMF IN: Save block tridiagonal matrix in sparse matrix format (1: enable, 0: disable)
 * @param read_BT_band IN: Read block tridiagonal matrix in band matrix format (1: enable, 0: disable)
 * @param deflation_on IN: Set deflation mode (1: enable, 0: disable)
 * @param conv_method IN: Select convergence method (1: eigenvalues,  2: eigenvectors)
 * @param conv_vec_method IN: If conv_method = 2, 1 for ALL components and 2 for selected components
 * @return Number of eigenvalues found
 */


#endif 	

