/*****************************************************************************
Copyright (C) 2002-2004 California Institute of Technology (Caltech)

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

Written by Marek J. Korkusinski
*****************************************************************************/

#include "ktridiag.h"

/*
 * Allocate memory for a complex tridiagonal matrix.
 */

ktridiag Ktridiag(int dim) {
  return nml_fc3_new(dim);
  }
/* Unless T already exists and has dimension dim,
   free T and allocate a new T with dimension dim,			*/
ktridiag check_ktridiag(ktridiag pT, int dim) {
  if ((NULL == pT) || (nml_fc3_extent(pT) != (nml_extent)dim)) {
    nml_fc3_resize(&pT, dim);
    }
  return pT;
  }

/* Allocates a new complex tridiagonal matrix Ctridiag "T",
   copies the source ctridiag "S" onto target ctridiag "T"
   then returns "T".							*/
ktridiag cp_ktridiag(ktridiag pT, ktridiag pS) {
  pT = check_ktridiag(pT, nml_fc3_extent(pS));
  nml_fcv_vcpy(pT[-1], pS[-1]);
  nml_fcv_vcpy(pT[ 0], pS[ 0]);
  nml_fcv_vcpy(pT[+1], pS[+1]);
  return pT;
  }

/* Deallocates memory for the tridiagonal complex matrix.		*/
void rm_ktridiag(ktridiag* ppT) {
  if (NULL != ppT) {
    nml_fc3_delete(*ppT);
    *ppT = NULL;
    }
  }

/* Multiplies a complex tridiagonal with a complex vector.		*/
nml_fcvector* mul_ktridiag(nml_fcvector* py,
    const ktridiag pT, const nml_fcvector* px) {	/* y <-- Tx	*/
  return nml_fcv_v3dot(py, px, pT);
  }

/* Solves the system "Tx = y"  for x
   where T is a complex tridiagonal system.				*/
nml_fcvector* solve_ktridiag(nml_fcvector* px,
    const ktridiag pT, const nml_fcvector* py) {
  return nml_fcv_solveTridiagonal(px, pT, py);
  }

