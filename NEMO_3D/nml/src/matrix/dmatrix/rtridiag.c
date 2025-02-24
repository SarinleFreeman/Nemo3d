/*****************************************************************************
Copyright (C) 2002 California Institute of Technology (Caltech)

This file is part of
The NanoElectronic MOdeling (NEMO) Math Library.

This library is free software which you can redistribute and/or modify
under the terms of the GNU Library General Public License
as published by the Free Software Foundation;
either version 2, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this library; see the file COPYING. If not, write to the
Free Software Foundation, Inc.,
59 Temple Place, Suite 330,
Boston, MA  02111-1307  USA

Written by E. Robert Tisdale

*****************************************************************************
$Header: /repo/nml/src/matrix/dmatrix/rtridiag.c,v 1.2 2003/10/08 16:19:43 hook Exp $
*****************************************************************************/


#include "rtridiag.h"

/* Allocate memory for a real tridiagonal matrix.			*/
rtridiag Rtridiag(int dim) {
  return nml_d3_new(dim);
  }

/* Unless T already exists and has dimension dim,
   free T and allocate a new T with dimension dim.			*/
rtridiag check_rtridiag(rtridiag pT, int dim) {
  if ((NULL == pT) || (nml_d3_extent(pT) != (nml_extent)dim)) {
    nml_d3_resize(&pT, dim);
    }
  return pT;
  }

/* Allocates a new real tridiagonal matrix Rtridiag "T",
 * copies the source rtridiag "S" onto target rtridiag "T"
 * then returns "T".  
 */
rtridiag cp_rtridiag(rtridiag pT, rtridiag pS) {
  pT = check_rtridiag(pT, nml_d3_extent(pS));
  nml_dv_vcpy(pT[-1], pS[-1]);
  nml_dv_vcpy(pT[ 0], pS[ 0]);
  nml_dv_vcpy(pT[+1], pS[+1]);
  return pT;
  }

/* Deallocates memory for a real tridiagonal matrix.			*/
void rm_rtridiag(rtridiag* ppT) {
  if (NULL != ppT) {
    nml_d3_delete(*ppT);
    *ppT = NULL;
    }
  }

/* Multiplies a tridiagonal real matrix with a real vector.		*/
nml_dvector* mul_rtridiag(nml_dvector* py,
    const rtridiag pT, const nml_dvector* px) {		/* y <-- Tx 	*/
  return nml_dv_v3dot(py, px, pT);
  }

/* Solves the system "Tx = y"  for x
   where T is a real tridiagonal system.				*/
nml_dvector* solve_rtridiag(nml_dvector* px,
    const rtridiag pT, const nml_dvector* py) {
  return nml_dv_solveTridiagonal(px, pT, py);
  }

