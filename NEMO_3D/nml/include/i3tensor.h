#ifndef I3TENSOR_H
#define I3TENSOR_H 1

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
$Header: /repo/nml/include/i3tensor.h,v 1.3 2003/10/08 16:19:15 hook Exp $
*****************************************************************************/


#include "tensor.h"
/* Include nml_i3DBase type definition.					*/
#include "nml_i3DBase.h"

typedef nml_i3DBase* i3tensor;

#ifdef	NML_INLINE
/* Constructor								*/
inline static
i3tensor	I3tensor(int lb3, int ub3, int lb2, int ub2, int lb1, int ub1) {
  nml_i3DBase*	ppp = nml_i3D_new(lb3, ub3, lb2, ub2, lb1, ub1);
  if (NULL != ppp)
    ppp -= lb3;
  return ppp;
  }
/* Destructor								*/
inline static
void		rm_i3tensor(i3tensor* ppT, int lb3) {
  if (NULL != ppT) {
    nml_i3D_delete(*ppT + lb3);
    *ppT = NULL;
    }
  }
/* Reconstructor							*/
inline static
i3tensor	check_i3tensor(i3tensor pT,
    int lb3, int ub3, int lb2, int ub2, int lb1, int ub1) {
  if (NULL != pT) {
    nml_i3DBase*	p3 = pT + lb3;
    if (nml_i3D_lower2(p3) == lb3) {
      if ((nml_i3D_upper3(p3) != ub3)
       || (nml_i3D_lower2(p3) != lb2)
       || (nml_i3D_upper2(p3) != ub2)
       || (nml_i3D_lower1(p3) != lb1)
       || (nml_i3D_upper1(p3) != ub1)) {
	p3 = nml_i3D_resize(&p3, ub3, lb3, lb2, ub2, lb1, ub1);
	pT = (NULL == p3)? NULL: p3 - lb3;
	}
      }
    else {
      nml_message(
"In function check_i3tensor(i3tensor, int, int, int, int, int, int):\n"
"second argument does not match i3tensor lower bound.");
      }
    }
  else {		/* (NULL == pT)		*/
    pT = nml_i3D_new(lb3, ub3, lb2, ub2, lb1, ub1);
    }
  return pT;
  }
#else /*NML_INLINE	*/
/* Constructor								*/
i3tensor	I3tensor(int, int, int, int, int, int);
/* Destructor								*/
void		rm_i3tensor(i3tensor*, int);
/* Reconstructor							*/
i3tensor	check_i3tensor(i3tensor, int, int, int, int, int, int);
#endif/*NML_INLINE	*/

#ifdef	NML_MACROS
#endif/*NML_MACROS	*/

#endif /* I3TENSOR_H */
