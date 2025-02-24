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


#include "nml_fcvector.h"
#undef NML_INLINE
#undef NML_MACROS
#include "kvector.h"

/* Constructor								*/
kvectr	Kvectr(int n) {
  return nml_fcv_new(n);
  }
#define	Kvectr(n) nml_fcv_new(n)

/* Destructor								*/
void	rm_kvectr(kvectr* ppv) {
  if (NULL != ppv) {
    nml_fcv_delete(*ppv);
    *ppv = NULL;
    }
  }
#define rm_kvectr(ppv) \
do { \
  if (NULL != (ppv)) { \
    nml_fcv_delete(*(ppv)); \
    *(ppv) = NULL; \
    } \
  } while(0)

/* Reconstructor							*/
kvectr	check_kvectr(kvectr pv, int n) {
  if (NULL == pv)
    pv = nml_fcv_new(n);
  else {		/* (NULL != pv)	*/
    if (nml_fcv_extent(pv) != (nml_extent)n)
      pv = nml_fcv_resize(&pv, n);
    }
  return pv;
  }

