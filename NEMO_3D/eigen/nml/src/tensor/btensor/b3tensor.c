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


#include "nml_b3DBase.h"
#undef NML_INLINE
#undef NML_MACROS
/* Include b3tensor type definition.					*/
#include "b3tensor.h"

/* Constructor								*/
b3tensor	B3tensor(int lb3, int ub3, int lb2, int ub2, int lb1, int ub1) {
  nml_b3DBase*	ppp = nml_b3D_new(lb3, ub3, lb2, ub2, lb1, ub1);
  if (NULL != ppp)
    ppp -= lb3;
  return ppp;
  }

/* Destructor								*/
void		rm_b3tensor(b3tensor* ppT, int lb3) {
  if (NULL != ppT) {
    nml_b3D_delete(*ppT + lb3);
    *ppT = NULL;
    }
  }

/* Reconstructor							*/
b3tensor	check_b3tensor(b3tensor pT,
    int lb3, int ub3, int lb2, int ub2, int lb1, int ub1) {
  if (NULL != pT) {
    nml_b3DBase*	p3 = pT + lb3;
    if (nml_b3D_lower2(p3) == lb3) {
      if ((nml_b3D_upper3(p3) != ub3)
       || (nml_b3D_lower2(p3) != lb2)
       || (nml_b3D_upper2(p3) != ub2)
       || (nml_b3D_lower1(p3) != lb1)
       || (nml_b3D_upper1(p3) != ub1)) {
	p3 = nml_b3D_resize(&p3, ub3, lb3, lb2, ub2, lb1, ub1);
	pT = (NULL == p3)? NULL: p3 - lb3;
	}
      }
    else {
      nml_message(
"In function check_b3tensor(b3tensor, int, int, int, int, int, int):\n"
"second argument does not match b3tensor lower bound.");
      }
    }
  else {		/* (NULL == pT)		*/
    pT = nml_b3D_new(lb3, ub3, lb2, ub2, lb1, ub1);
    }
  return pT;
  }

