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
$Header: /repo/nml/src/matrix/imatrix/imatrix.c,v 1.3 2003/10/08 16:19:46 hook Exp $
*****************************************************************************/

#include "nml_imatrix.h"
#undef NML_INLINE
#undef NML_MACROS
#include "imatrix.h"

imatrix	Imatrix(int n, int m) {
  return nml_im_new(n, m);
  }
#define Imatrix(n, m) nml_im_new((n), (m))

void	rm_imatrix(imatrix* ppM) {
  if (NULL != ppM) {
    nml_im_delete(*ppM);
    *ppM = NULL;
    }
  }
#define rm_imatrix(ppM) \
do { \
  if (NULL != (ppM)) { \
    nml_im_delete(*(ppM)); \
    *(ppM) = NULL; \
    } \
  } while(0)

imatrix	check_imatrix(imatrix pM, int n, int m) {
  if (NULL == pM)
    pM = nml_im_new(n, m);
  else {		/* (NULL != pM)	*/
    if ((nml_im_extent2(pM) != (nml_extent)n)
     || (nml_im_extent1(pM) != (nml_extent)m))
      pM = nml_im_resize(&pM, n, m);
    }
  return pM;
  }

