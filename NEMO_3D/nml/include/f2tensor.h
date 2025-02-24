#ifndef F2TENSOR_H
#define F2TENSOR_H 1

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


#include "tensor.h"
/* Include nml_f2DBase type definition.					*/
#include "nml_f2DBase.h"

typedef nml_f2DBase* f2tensor;

#ifdef	NML_INLINE
/* Constructor								*/
inline static
f2tensor	F2tensor(int lb2, int ub2, int lb1, int ub1) {
  nml_f2DBase*	pp = nml_f2D_new(lb2, ub2, lb1, ub1);
  if (NULL != pp)
    return pp - lb2;
  else
    return NULL;
  }
/* Destructor								*/
inline static
void		rm_f2tensor(f2tensor* ppM, int lb2) {
  if (NULL != ppM) {
    nml_f2D_delete(*ppM + lb2);
    *ppM = NULL;
    }
  }
/* Reconstructor							*/
inline static
f2tensor	check_f2tensor(f2tensor pM,
    int lb2, int ub2, int lb1, int ub1) {
  if (NULL != pM) {
    nml_f2DBase*	p2 = pM + lb2;
    if (nml_f2D_lower2(p2) == lb2) {
      if ((nml_f2D_upper2(p2) != ub2)
       || (nml_f2D_lower1(p2) != lb1)
       || (nml_f2D_upper1(p2) != ub1)) {
	p2 = nml_f2D_resize(&p2, lb2, ub2, lb1, ub1);
	pM = (NULL == p2)? NULL: p2 - lb2;
	}
      }
    else {
      nml_message(
"In function check_f2tensor(f2tensor, int, int, int, int):\n"
"second argument does not match f2tensor lower bound.");
      }
    }
  else {		/* (NULL == pM)		*/
    pM = F2tensor(lb2, ub2, lb1, ub1);
    }
  return pM;
  }
#else /*NML_INLINE	*/
/* Constructor								*/
f2tensor	F2tensor(int, int, int, int);
/* Destructor								*/
void		rm_f2tensor(f2tensor*, int);
/* Reconstructor							*/
f2tensor	check_f2tensor(f2tensor, int, int, int, int);
#endif/*NML_INLINE	*/

#ifdef	NML_MACROS
#endif/*NML_MACROS	*/

#endif /* F2TENSOR_H */
