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
$Header: /repo/nml/src/vector/dcvector/cvector.c,v 1.3 2003/10/08 16:20:08 hook Exp $
*****************************************************************************/


#include "nml_dcvector.h"
#undef NML_INLINE
#undef NML_MACROS
#include "cvector.h"

/* Constructor								*/
cvectr	Cvectr(int n) {
  return nml_dcv_new(n);
  }
#define	Cvectr(n) nml_dcv_new(n)

/* Destructor								*/
void	rm_cvectr(cvectr* ppv) {
  if (NULL != ppv) {
    nml_dcv_delete(*ppv);
    *ppv = NULL;
    }
  }
#define rm_cvectr(ppv) \
do { \
  if (NULL != (ppv)) { \
    nml_dcv_delete(*(ppv)); \
    *(ppv) = NULL; \
    } \
  } while(0)

/* Reconstructor							*/
cvectr	check_cvectr(cvectr pv, int n) {
  if (NULL == pv)
    pv = nml_dcv_new(n);
  else {		/* (NULL != pv)	*/
    if (nml_dcv_extent(pv) != (nml_extent)n)
      pv = nml_dcv_resize(&pv, n);
    }
  return pv;
  }

