#ifndef _nml_bmatrix_h
#define _nml_bmatrix_h 1

/*****************************************************************************
The NEMO Math Library.
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
$Header: /repo/nml/src/matrix/matrix.hP,v 1.2 2003/10/08 16:19:35 hook Exp $
*****************************************************************************/



/* Include vector type definitions.					*/
#include<nml_bvector.h>

typedef nml_bscalar*	nml_bmatrix;

#ifdef	NML_INLINE
/* Functions								*/
inline static/*
The NEMO Math Library.
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

You should have received a copy of the GNU Library General Public License
along with this library.  If not, write to the Free Software Foundation,
Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

Written by E. Robert Tisdale
*/
nml_extent		(nml_bm_extent1)(const nml_bmatrix* pM) {
  return *((const nml_extent*)(pM - 1));
  }
inline static
nml_extent		(nml_bm_extent2)(const nml_bmatrix* pM) {
  return *((const nml_extent*)(pM - 2));
  }

/* Constructor								*/
inline static
nml_bmatrix*	(nml_bm_new)(nml_extent m, nml_extent n) {
  nml_bscalar**	pp
    = (nml_bscalar**)nml_calloc(3+m, sizeof(nml_bscalar*));
  if (NULL == pp) {
    nml_message("In function nml_bm_new(nml_extent, nml_extent):\n"
      "unable to allocate memory.");
    }
  else {		/* (NULL != pp)		*/
    *pp = NULL;		     ++pp;
    *((nml_extent*)pp)	= m; ++pp;
    *((nml_extent*)pp)	= n; ++pp;
    if (0 < m) {
      if (0 < n) {
	pp[0] = (nml_bscalar*)nml_bv_new(n*m);
	if (NULL == pp[0]) {
	  nml_offset	i = 1;
          for (i = 1; i < m; ++i)
	    pp[i] = NULL;
	  nml_message("In function nml_bm_new(nml_extent, nml_extent):\n"
	    "unable to allocate memory.");
	  }
	else {		/* (NULL != pp[0])	*/
	  nml_offset	i = 1;
	  for (i = 1; i < m; ++i)
	    pp[i] = pp[i-1] + n;
	  }
	}
      else {		/* (n <= 0)		*/
	nml_offset	i = 0;
        for (i = 0; i < m; ++i)
	  pp[i] = NULL;
	}
      }
    }
  return pp;
  }

/* Destructor								*/
inline static
void			(nml_bm_delete)(nml_bmatrix* pM) {
  if (NULL != pM) {
    nml_bv_delete(pM[0]);
    nml_free(pM - 3);
    }
  }

/* Reconstructor							*/
inline static
nml_bmatrix*	(nml_bm_resize)(nml_bmatrix* *ppM,
    nml_extent m, nml_extent n) {
  if (NULL != ppM) {
    (nml_bm_delete)(*ppM);
    *ppM = (nml_bm_new)(m, n);
    return *ppM;
    }
  else {
    nml_message(
"In function nml_bm_resize(nml_bmatrix**, nml_extent, nml_extent):\n"
"invalid pointer to nml_bmatrix*.");
    return NULL;
    }
  }
#else /*NML_INLINE	*/
/* Functions								*/
nml_extent		(nml_bm_extent1)(const nml_bmatrix*);
nml_extent		(nml_bm_extent2)(const nml_bmatrix*);
/* Constructor								*/
nml_bmatrix*	(nml_bm_new)(nml_extent, nml_extent);
/* Destructor								*/
void			(nml_bm_delete)(nml_bmatrix*);
/* Reconstructor							*/
nml_bmatrix*	(nml_bm_resize)(nml_bmatrix**,
  nml_extent, nml_extent);
#endif/*NML_INLINE	*/

int			(nml_bm_fprintf)(FILE*,
    const nml_bmatrix*, int, unsigned int, nml_fmtflags, unsigned int);

#ifdef	NML_MACROS
#define nml_bm_extent1(pM) (*((const nml_extent*)((pM) - 1)))
#define nml_bm_extent2(pM) (*((const nml_extent*)((pM) - 2)))
#endif/*NML_MACROS	*/

#endif /* _nml_bmatrix_h */
