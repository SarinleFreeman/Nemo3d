#ifndef _nml_i2DBase_h
#define _nml_i2DBase_h 1

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
$Header: /repo/nml/src/matrix/2DBase.hP,v 1.2 2003/10/08 16:19:31 hook Exp $
*****************************************************************************/


/* Include scalar type definitions.					*/
#include<nml_iscalar.h>

typedef nml_iscalar*	nml_i2DBase;

#ifdef	NML_INLINE
/* Functions								*/
inline static
nml_index		(nml_i2D_lower2)(const nml_i2DBase* pM) {
  return *((const nml_index*)(pM - 1));
  }
inline static
nml_index		(nml_i2D_upper2)(const nml_i2DBase* pM) {
  return *((const nml_index*)(pM - 2));
  }
inline static
nml_extent		(nml_i2D_extent2)(const nml_i2DBase* pM) {
  return 1+nml_i2D_upper2(pM)-nml_i2D_lower2(pM);
  }
inline static
nml_index		(nml_i2D_lower1)(const nml_i2DBase* pM) {
  return *((const nml_index*)(pM - 3));
  }
inline static
nml_index		(nml_i2D_upper1)(const nml_i2DBase* pM) {
  return *((const nml_index*)(pM - 4));
  }
inline static
nml_extent		(nml_i2D_extent1)(const nml_i2DBase* pM) {
  return 1+nml_i2D_upper1(pM)-nml_i2D_lower1(pM);
  }

/* Constructor								*/
inline static
nml_i2DBase*	(nml_i2D_new)(
    nml_index lb2, nml_index ub2, nml_index lb1, nml_index ub1) {
  nml_iscalar**	pp = NULL;
  nml_index		m = 1 + ub2 - lb2;
  nml_index		n = 1 + ub1 - lb1;
  if (m < 0) {
    nml_message(
"In function nml_i2D_new(nml_index, nml_index, nml_index, nml_index):\n"
"extent2 is negative.");
    }
  else
  if (n < 0) {
    nml_message(
"In function nml_i2D_new(nml_index, nml_index, nml_index, nml_index):\n"
"extent1 is negative.");
    }
  else {		/* (0 <= m) && (0 <= n)	*/
    pp = (nml_iscalar**)nml_calloc(4 + m, sizeof(nml_iscalar*));
    if (NULL != pp) {
      *((nml_index*)pp) = ub1; ++pp;
      *((nml_index*)pp) = lb1; ++pp;
      *((nml_index*)pp) = ub2; ++pp;
      *((nml_index*)pp) = lb2; ++pp;
      if (0 < m) {
	if (0 < n) {
	  nml_iscalar*	p
	    = (nml_iscalar*)nml_calloc(m*n, sizeof(nml_iscalar));
	  if (NULL != p) {
	    pp[0] = p - lb1;
	    { nml_index	i = 1;
	      for (i = 1; i < m; ++i)
		pp[i] = pp[i-1] + n;
	      }
	    }
	  else {	/* (NULL == p)		*/
	    int	i = 1;
            for (i = 1; i < m; ++i)
	      pp[i] = NULL;
	    nml_message(
"In function nml_i2D_new(nml_index, nml_index, nml_index, nml_index):\n"
"unable to allocate memory.");
	    }
	  }
	else {		/* (n <= 0)		*/
	  int	i = 0;
          for (i = 0; i < m; ++i)
	    pp[i] = NULL;
	  }
	}
      }
    else {		/* (NULL == pp)		*/
      nml_message(
"In function nml_i2D_new(nml_index, nml_index, nml_index, nml_index):\n"
"unable to allocate memory.");
      }
    }
  return pp;
  }

/* Destructor								*/
inline static
void			(nml_i2D_delete)(nml_i2DBase* pM) {
  if (NULL != pM) {
    if (NULL != pM[0]) {
      nml_index	lb1 = nml_i2D_lower1(pM);
      nml_free(pM[0] + lb1);
      }
    nml_free(pM - 4);
    }
  }

/* Reconstructor							*/
inline static
nml_i2DBase*	(nml_i2D_resize)(nml_i2DBase* *ppM,
    nml_index lb2, nml_index ub2, nml_index lb1, nml_index ub1) {
  if (NULL != ppM) {
    (nml_i2D_delete)(*ppM);
    *ppM = (nml_i2D_new)(lb2, ub2, lb1, ub1);
    return *ppM;
    }
  else {
    nml_message(
"In function nml_i2D_resize(nml_i2DBase**, nml_index, nml_index,\n"
"  nml_index, nml_index): invalid pointer to nml_i2DBase*.");
    return NULL;
    }
  }
#else /*NML_INLINE	*/
/* Functions								*/
nml_index		(nml_i2D_lower2)(const nml_i2DBase*);
nml_index		(nml_i2D_upper2)(const nml_i2DBase*);
nml_extent		(nml_i2D_extent2)(const nml_i2DBase*);

nml_index		(nml_i2D_lower1)(const nml_i2DBase*);
nml_index		(nml_i2D_upper1)(const nml_i2DBase*);
nml_extent		(nml_i2D_extent1)(const nml_i2DBase*);
/* Constructor								*/
nml_i2DBase*	(nml_i2D_new)(
    nml_index, nml_index, nml_index, nml_index);
/* Destructor								*/
void			(nml_i2D_delete)(nml_i2DBase*);
/* Reconstructor							*/
nml_i2DBase*	(nml_i2D_resize)(nml_i2DBase**,
    nml_index, nml_index, nml_index, nml_index);
#endif/*NML_INLINE	*/

int			(nml_i2D_fprintf)(FILE*,
    const nml_i2DBase*, int, unsigned int, nml_fmtflags, unsigned int);

#ifdef	NML_MACROS
#define nml_i2D_lower2(pM) (*((const nml_index*)((pM) - 1)))
#define nml_i2D_upper2(pM) (*((const nml_index*)((pM) - 2)))
#define nml_i2D_extent2(pM) (1+nml_i2D_upper2(pM)-nml_i2D_lower2(pM))
#define nml_i2D_lower1(pM) (*((const nml_index*)((pM) - 3)))
#define nml_i2D_upper1(pM) (*((const nml_index*)((pM) - 4)))
#define nml_i2D_extent1(pM) (1+nml_i2D_upper1(pM)-nml_i2D_lower1(pM))
#endif/*NML_MACROS	*/

#endif /* _nml_i2DBase_h */
