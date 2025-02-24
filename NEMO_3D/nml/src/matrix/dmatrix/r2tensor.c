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
$Header: /repo/nml/src/matrix/dmatrix/r2tensor.c,v 1.3 2003/10/08 16:19:41 hook Exp $
*****************************************************************************/


#include "nml_d2DBase.h"
#undef NML_INLINE
#undef NML_MACROS
/* Include r2tensor type definition.					*/
#include "r2tensor.h"

/* Constructor								*/
r2tensor	R2tensor(int lb2, int ub2, int lb1, int ub1) {
  nml_d2DBase*	pp = nml_d2D_new(lb2, ub2, lb1, ub1);
  if (NULL != pp)
    return pp - lb2;
  else
    return NULL;
  }

/* Destructor								*/
void		rm_r2tensor(r2tensor* ppM, int lb2) {
  if (NULL != ppM) {
    nml_d2D_delete(*ppM + lb2);
    *ppM = NULL;
    }
  }

/* Reconstructor							*/
r2tensor	check_r2tensor(r2tensor pM,
    int lb2, int ub2, int lb1, int ub1) {
  if (NULL != pM) {
    nml_d2DBase*	p2 = pM + lb2;
    if (nml_d2D_lower2(p2) == lb2) {
      if ((nml_d2D_upper2(p2) != ub2)
       || (nml_d2D_lower1(p2) != lb1)
       || (nml_d2D_upper1(p2) != ub1)) {
	p2 = nml_d2D_resize(&p2, lb2, ub2, lb1, ub1);
	pM = (NULL == p2)? NULL: p2 - lb2;
	}
      }
    else {
      nml_message(
"In function check_r2tensor(r2tensor, int, int, int, int):\n"
"second argument does not match r2tensor lower bound.");
      }
    }
  else {		/* (NULL == pM)		*/
    pM = R2tensor(lb2, ub2, lb1, ub1);
    }
  return pM;
  }
