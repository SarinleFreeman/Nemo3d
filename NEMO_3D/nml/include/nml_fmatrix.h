#ifndef _nml_fmatrix_h
#define _nml_fmatrix_h 1

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
#include<nml_fvector.h>

typedef nml_fscalar*	nml_fmatrix;

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
nml_extent		(nml_fm_extent1)(const nml_fmatrix* pM) {
  return *((const nml_extent*)(pM - 1));
  }
inline static
nml_extent		(nml_fm_extent2)(const nml_fmatrix* pM) {
  return *((const nml_extent*)(pM - 2));
  }

/* Constructor								*/
inline static
nml_fmatrix*	(nml_fm_new)(nml_extent m, nml_extent n) {
  nml_fscalar**	pp
    = (nml_fscalar**)nml_calloc(3+m, sizeof(nml_fscalar*));
  if (NULL == pp) {
    nml_message("In function nml_fm_new(nml_extent, nml_extent):\n"
      "unable to allocate memory.");
    }
  else {		/* (NULL != pp)		*/
    *pp = NULL;		     ++pp;
    *((nml_extent*)pp)	= m; ++pp;
    *((nml_extent*)pp)	= n; ++pp;
    if (0 < m) {
      if (0 < n) {
	pp[0] = (nml_fscalar*)nml_fv_new(n*m);
	if (NULL == pp[0]) {
	  nml_offset	i = 1;
          for (i = 1; i < m; ++i)
	    pp[i] = NULL;
	  nml_message("In function nml_fm_new(nml_extent, nml_extent):\n"
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
void			(nml_fm_delete)(nml_fmatrix* pM) {
  if (NULL != pM) {
    nml_fv_delete(pM[0]);
    nml_free(pM - 3);
    }
  }

/* Reconstructor							*/
inline static
nml_fmatrix*	(nml_fm_resize)(nml_fmatrix* *ppM,
    nml_extent m, nml_extent n) {
  if (NULL != ppM) {
    (nml_fm_delete)(*ppM);
    *ppM = (nml_fm_new)(m, n);
    return *ppM;
    }
  else {
    nml_message(
"In function nml_fm_resize(nml_fmatrix**, nml_extent, nml_extent):\n"
"invalid pointer to nml_fmatrix*.");
    return NULL;
    }
  }
#else /*NML_INLINE	*/
/* Functions								*/
nml_extent		(nml_fm_extent1)(const nml_fmatrix*);
nml_extent		(nml_fm_extent2)(const nml_fmatrix*);
/* Constructor								*/
nml_fmatrix*	(nml_fm_new)(nml_extent, nml_extent);
/* Destructor								*/
void			(nml_fm_delete)(nml_fmatrix*);
/* Reconstructor							*/
nml_fmatrix*	(nml_fm_resize)(nml_fmatrix**,
  nml_extent, nml_extent);
#endif/*NML_INLINE	*/

int			(nml_fm_fprintf)(FILE*,
    const nml_fmatrix*, int, unsigned int, nml_fmtflags, unsigned int);

#ifdef	NML_MACROS
#define nml_fm_extent1(pM) (*((const nml_extent*)((pM) - 1)))
#define nml_fm_extent2(pM) (*((const nml_extent*)((pM) - 2)))
#endif/*NML_MACROS	*/

#endif /* _nml_fmatrix_h */
