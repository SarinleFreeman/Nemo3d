#ifndef RTRIDIAG_H
#define RTRIDIAG_H 1

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
$Header: /repo/nml/include/rtridiag.h,v 1.3 2003/10/08 16:19:23 hook Exp $
*****************************************************************************/

#include "nml_d3bands.h"
typedef nml_d3bands* rtridiag;

#include "matrix.h"
#define tddim(t) *int_ptr((t-2))


/*  element numbering convention:
|t[ 0][0] t[ 1][0]                           | 
|t[-1][1] t[ 0][1] t[ 1][1]                  |
|         t[-1][2] t[ 0][2] t[ 1][2]         |
|                     ...      ...      ...  |
|                           t[-1][m] t[ 0][m]|,
where m = dim - 1.							*/

rtridiag Rtridiag(int n);
rtridiag Rtridiag_clone(const rtridiag);
rtridiag tridiag_copy(rtridiag p, const rtridiag t);

rtridiag check_rtridiag(rtridiag t, int n);
/* If t exists and has dimension n, returns t,
   otherwise allocates new t with dim n,
   frees old t if t != NULL.						*/

rtridiag cp_rtridiag(rtridiag t, rtridiag s);
/* copies s into t and returns t.					*/

void rm_rtridiag(rtridiag* t_ptr);

nml_dvector* solve_rtridiag(nml_dvector* x, const rtridiag t, const nml_dvector* y);
/* solves tx = y, t tridiagonal. t and y are preserved.			*/

nml_dvector* mul_rtridiag(nml_dvector* y, rtridiag t, const nml_dvector* x);
/* evaluates y = tx.							*/

#endif /*RTRIDIAG_H */

