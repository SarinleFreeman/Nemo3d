#ifndef _nml_fc3DBase_h
#define _nml_fc3DBase_h 1

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
$Header: /repo/nml/src/tensor/3DBase.hP,v 1.2 2003/10/08 16:19:54 hook Exp $
*****************************************************************************/



/* Include float precision real 3DBase type definitions.		*/
#include<nml_f3DBase.h>

/* Include scalar type definitions.					*/
#include<nml_fcscalar.h>

typedef nml_fcscalar**	nml_fc3DBase;

#ifdef	NML_INLINE
/* Functions								*/
inline static
nml_index		(nml_fc3D_lower3)(const nml_fc3DBase* pT) {
  return *((const nml_index*)(pT - 1));
  }
inline static
nml_index		(nml_fc3D_upper3)(const nml_fc3DBase* pT) {
  return *((const nml_index*)(pT - 2));
  }
inline static
nml_extent		(nml_fc3D_extent3)(const nml_fc3DBase* pT) {
  return 1+nml_fc3D_upper3(pT)-nml_fc3D_lower3(pT);
  }
inline static
nml_index		(nml_fc3D_lower2)(const nml_fc3DBase* pT) {
  return *((const nml_index*)(pT - 3));
  }
inline static
nml_index		(nml_fc3D_upper2)(const nml_fc3DBase* pT) {
  return *((const nml_index*)(pT - 4));
  }
inline static
nml_extent		(nml_fc3D_extent2)(const nml_fc3DBase* pT) {
  return 1+nml_fc3D_upper2(pT)-nml_fc3D_lower2(pT);
  }
inline static
nml_index		(nml_fc3D_lower1)(const nml_fc3DBase* pT) {
  return *((const nml_index*)(pT - 5));
  }
inline static
nml_index		(nml_fc3D_upper1)(const nml_fc3DBase* pT) {
  return *((const nml_index*)(pT - 6));
  }
inline static
nml_extent		(nml_fc3D_extent1)(const nml_fc3DBase* pT) {
  return 1+nml_fc3D_upper1(pT)-nml_fc3D_lower1(pT);
  }

/* Constructor								*/
inline static
nml_fc3DBase*	(nml_fc3D_new)(nml_index lb3, nml_index ub3,
    nml_index lb2, nml_index ub2, nml_index lb1, nml_index ub1) {
  nml_fcscalar***	ppp = NULL;
  nml_index		l = 1 + ub3 - lb3;
  nml_index		m = 1 + ub2 - lb2;
  nml_index		n = 1 + ub1 - lb1;
  if (l < 0) {
    nml_message(
"In function nml_fc3D_new(nml_index, nml_index, nml_index, nml_index,\n"
"  nml_index, nml_index): extent3 is negative.");
    }
  else
  if (m < 0) {
    nml_message(
"In function nml_fc3D_new(nml_index, nml_index, nml_index, nml_index,\n"
"  nml_index, nml_index): extent2 is negative.");
    }
  else
  if (n < 0) {
    nml_message(
"In function nml_fc3D_new(nml_index, nml_index, nml_index, nml_index,\n"
"  nml_index, nml_index): extent1 is negative.");
    }
  else {		/* (0 <= l) && (0 <= m) && (0 <= n)	*/
    ppp = (nml_fcscalar***)nml_calloc(6 + l, sizeof(nml_fcscalar**));
    if (NULL != ppp) {
      *((nml_index*)ppp) = ub1; ++ppp;
      *((nml_index*)ppp) = lb1; ++ppp;
      *((nml_index*)ppp) = ub2; ++ppp;
      *((nml_index*)ppp) = lb2; ++ppp;
      *((nml_index*)ppp) = ub3; ++ppp;
      *((nml_index*)ppp) = lb3; ++ppp;
      if (0 < l) {
	if (0 < m) {	/* (0 < l) && (0 < m)	*/
	  nml_fcscalar** pp
	    = (nml_fcscalar**)nml_calloc(l*m, sizeof(nml_fcscalar*));
	  if (NULL != pp) {
	    ppp[0] = pp - lb2;
	    { nml_index	h = 1;
	      for (h = 1; h < l; ++h)
		ppp[h] = ppp[h-1] + m;
	      }
	    if (0 < n) {
	      nml_fcscalar*	p
		= (nml_fcscalar*)nml_calloc(l*m*n, sizeof(nml_fcscalar));
	      if (NULL != p) {
		ppp[0][lb2] = p - lb1;
		{ nml_index	h = 1;
		  for (h = 1; h < l; ++h)
		    ppp[h][lb2] = ppp[h-1][lb2] + m*n;
		  for (h = 0; h < l; ++h) {
		    nml_index	i = lb2 + 1;
		    for (i = lb2 + 1; i <= ub2; ++i)
		      ppp[h][i] = ppp[h][i-1] + n;
		    }
		  }
		}

	      else {	/* (NULL == p)		*/
		nml_index	h = 0;
		for (h = 0; h < l; ++h) {
		  nml_index	i = lb2;
		  for (i = lb2; i <= ub2; ++i)
		    ppp[h][i] = NULL;
		  }
		nml_message(
"In function nml_fc3D_new(nml_index, nml_index, nml_index, nml_index,\n"
"  nml_index, nml_index): unable to allocate memory.");
		}
	      }
	    else {	/* (n <= 0)		*/
	      nml_index	h = 0;
	      for (h = 0; h < l; ++h) {
		nml_index	i = lb2;
		for (i = lb2; i <= ub2; ++i)
		  ppp[h][i] = NULL;
		}
	      }
	    }
	  else {	/* (NULL == pp)		*/
	    nml_index	h = 0;
	    for (h = 0; h < l; ++h)
	      ppp[h] = NULL;
	    nml_message(
"In function nml_fc3D_new(nml_index, nml_index, nml_index, nml_index,\n"
"  nml_index, nml_index): unable to allocate memory.");
	    }
	  }
	else {		/* (m <= 0)		*/
	  nml_index	h = 0;
	  for (h = 0; h < l; ++h)
	    ppp[h] = NULL;
	  }
	}
      }
    else {		/* (NULL == ppp)	*/
      nml_message(
"In function nml_fc3D_new(nml_index, nml_index, nml_index, nml_index,\n"
"  nml_index, nml_index): unable to allocate memory.");
      }
    }
  return ppp;
  }

/* Destructor								*/
inline static
void			(nml_fc3D_delete)(nml_fc3DBase* pT) {
  if (NULL != pT) {
    if (NULL != pT[0]) {
      nml_index	lb2 = nml_fc3D_lower2(pT);
      if (NULL != pT[0][lb2]) {
	nml_index	lb1 = nml_fc3D_lower1(pT);
	nml_free(pT[0][lb2] + lb1);
	}
      nml_free(pT[0] + lb2);
      }
    nml_free(pT - 6);
    }
  }

/* Reconstructor							*/
inline static
nml_fc3DBase*	(nml_fc3D_resize)(nml_fc3DBase* *ppT,
    nml_index lb3, nml_index ub3, nml_index lb2, nml_index ub2,
    nml_index lb1, nml_index ub1) {
  if (NULL != ppT) {
    (nml_fc3D_delete)(*ppT);
    *ppT = (nml_fc3D_new)(lb3, ub3, lb2, ub2, lb1, ub1);
    return *ppT;
    }
  else {
    nml_message(
"In function nml_fc3D_resize(nml_fc3DBase**,\n"
"  nml_index, nml_index, nml_index, nml_index, nml_index, nml_index):\n"
"invalid pointer to nml_fc3DBase*.");
    return NULL;
    }
  }
#else /*NML_INLINE	*/
/* Functions								*/
nml_index		(nml_fc3D_lower3)(const nml_fc3DBase*);
nml_index		(nml_fc3D_upper3)(const nml_fc3DBase*);
nml_extent		(nml_fc3D_extent3)(const nml_fc3DBase*);

nml_index		(nml_fc3D_lower2)(const nml_fc3DBase*);
nml_index		(nml_fc3D_upper2)(const nml_fc3DBase*);
nml_extent		(nml_fc3D_extent2)(const nml_fc3DBase*);

nml_index		(nml_fc3D_lower1)(const nml_fc3DBase*);
nml_index		(nml_fc3D_upper1)(const nml_fc3DBase*);
nml_extent		(nml_fc3D_extent1)(const nml_fc3DBase*);
/* Constructor								*/
nml_fc3DBase*	(nml_fc3D_new)(
    nml_index, nml_index, nml_index, nml_index, nml_index, nml_index);
/* Destructor								*/
void			(nml_fc3D_delete)(nml_fc3DBase*);
/* Reconstructor							*/
nml_fc3DBase*	(nml_fc3D_resize)(nml_fc3DBase**,
    nml_index, nml_index, nml_index, nml_index, nml_index, nml_index);
#endif/*NML_INLINE	*/

int			(nml_fc3D_fprintf)(FILE*,
    const nml_fc3DBase*, int, unsigned int, nml_fmtflags, unsigned int);

#ifdef	NML_MACROS
#define nml_fc3D_lower3(pT) (*((const nml_index*)((pT) - 1)))
#define nml_fc3D_upper3(pT) (*((const nml_index*)((pT) - 2)))
#define nml_fc3D_extent3(pT) (1+nml_fc3D_upper3(pT)-nml_fc3D_lower3(pT))
#define nml_fc3D_lower2(pT) (*((const nml_index*)((pT) - 3)))
#define nml_fc3D_upper2(pT) (*((const nml_index*)((pT) - 4)))
#define nml_fc3D_extent2(pT) (1+nml_fc3D_upper2(pT)-nml_fc3D_lower2(pT))
#define nml_fc3D_lower1(pT) (*((const nml_index*)((pT) - 5)))
#define nml_fc3D_upper1(pT) (*((const nml_index*)((pT) - 6)))
#define nml_fc3D_extent1(pT) (1+nml_fc3D_upper1(pT)-nml_fc3D_lower1(pT))
#endif/*NML_MACROS	*/

#endif /* _nml_fc3DBase_h */
