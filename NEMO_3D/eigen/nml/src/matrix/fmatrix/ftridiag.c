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


#include "ftridiag.h"

/* Allocate memory for a real tridiagonal matrix.			*/
ftridiag Ftridiag(int dim) {
  return nml_f3_new(dim);
  }

/* Unless T already exists and has dimension dim,
   free T and allocate a new T with dimension dim.			*/
ftridiag check_ftridiag(ftridiag pT, int dim) {
  if ((NULL == pT) || (nml_f3_extent(pT) != (nml_extent)dim)) {
    nml_f3_resize(&pT, dim);
    }
  return pT;
  }

/* Allocates a new real tridiagonal matrix Rtridiag "T",
 * copies the source rtridiag "S" onto target rtridiag "T"
 * then returns "T".  
 */
ftridiag cp_ftridiag(ftridiag pT, ftridiag pS) {
  pT = check_ftridiag(pT, nml_f3_extent(pS));
  nml_fv_vcpy(pT[-1], pS[-1]);
  nml_fv_vcpy(pT[ 0], pS[ 0]);
  nml_fv_vcpy(pT[+1], pS[+1]);
  return pT;
  }

/* Deallocates memory for a real tridiagonal matrix.			*/
void rm_ftridiag(ftridiag* ppT) {
  if (NULL != ppT) {
    nml_f3_delete(*ppT);
    *ppT = NULL;
    }
  }

/* Multiplies a tridiagonal real matrix with a real vector.		*/
nml_fvector* mul_ftridiag(nml_fvector* py,
    const ftridiag pT, const nml_fvector* px) {		/* y <-- Tx 	*/
  return nml_fv_v3dot(py, px, pT);
  }

/* Solves the system "Tx = y"  for x
   where T is a real tridiagonal system.				*/
nml_fvector* solve_ftridiag(nml_fvector* px,
    const ftridiag pT, const nml_fvector* py) {
  return nml_fv_solveTridiagonal(px, pT, py);
  }

