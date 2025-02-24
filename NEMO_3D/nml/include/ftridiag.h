#ifndef FTRIDIAG_H
#define FTRIDIAG_H 1

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

#include "nml_f3bands.h"
typedef nml_f3bands* ftridiag;

#include "matrix.h"
#define tddim(t) *int_ptr((t-2))


/*  element numbering convention:
|t[ 0][0] t[ 1][0]                           | 
|t[-1][1] t[ 0][1] t[ 1][1]                  |
|         t[-1][2] t[ 0][2] t[ 1][2]         |
|                     ...      ...      ...  |
|                           t[-1][m] t[ 0][m]|,
where m = dim - 1.							*/

ftridiag Ftridiag(int n);
ftridiag Ftridiag_clone(const ftridiag);
ftridiag tridiag_copy(ftridiag p, const ftridiag t);

ftridiag check_ftridiag(ftridiag t, int n);
/* If t exists and has dimension n, returns t,
   otherwise allocates new t with dim n,
   frees old t if t != NULL.						*/

ftridiag cp_ftridiag(ftridiag t, ftridiag s);
/* copies s into t and returns t.					*/

void rm_ftridiag(ftridiag* t_ptr);

nml_fvector* solve_ftridiag(nml_fvector* x, const ftridiag t, const nml_fvector* y);
/* solves tx = y, t tridiagonal. t and y are preserved.			*/

nml_fvector* mul_ftridiag(nml_fvector* y, ftridiag t, const nml_fvector* x);
/* evaluates y = tx.							*/

#endif /*FTRIDIAG_H */

