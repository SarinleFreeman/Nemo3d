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
$Header: /repo/nml/src/matrix/dcmatrix/ctridiag.c,v 1.2 2003/10/08 16:19:38 hook Exp $
*****************************************************************************/

#include "ctridiag.h"

/*
 * Allocate memory for a complex tridiagonal matrix.
 */

ctridiag Ctridiag(int dim) {
  return nml_dc3_new(dim);
  }
/* Unless T already exists and has dimension dim,
   free T and allocate a new T with dimension dim,			*/
ctridiag check_ctridiag(ctridiag pT, int dim) {
  if ((NULL == pT) || (nml_dc3_extent(pT) != (nml_extent)dim)) {
    nml_dc3_resize(&pT, dim);
    }
  return pT;
  }

/* Allocates a new complex tridiagonal matrix Ctridiag "T",
   copies the source ctridiag "S" onto target ctridiag "T"
   then returns "T".							*/
ctridiag cp_ctridiag(ctridiag pT, ctridiag pS) {
  pT = check_ctridiag(pT, nml_dc3_extent(pS));
  nml_dcv_vcpy(pT[-1], pS[-1]);
  nml_dcv_vcpy(pT[ 0], pS[ 0]);
  nml_dcv_vcpy(pT[+1], pS[+1]);
  return pT;
  }

/* Deallocates memory for the tridiagonal complex matrix.		*/
void rm_ctridiag(ctridiag* ppT) {
  if (NULL != ppT) {
    nml_dc3_delete(*ppT);
    *ppT = NULL;
    }
  }

/* Multiplies a complex tridiagonal with a complex vector.		*/
nml_dcvector* mul_ctridiag(nml_dcvector* py,
    const ctridiag pT, const nml_dcvector* px) {	/* y <-- Tx	*/
  return nml_dcv_v3dot(py, px, pT);
  }

/* Solves the system "Tx = y"  for x
   where T is a complex tridiagonal system.				*/
nml_dcvector* solve_ctridiag(nml_dcvector* px,
    const ctridiag pT, const nml_dcvector* py) {
  return nml_dcv_solveTridiagonal(px, pT, py);
  }

