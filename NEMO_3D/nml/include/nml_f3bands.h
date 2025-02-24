#ifndef _nml_f3bands_h
#define _nml_f3bands_h 1

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
$Header: /repo/nml/src/matrix/3bands.hP,v 1.2 2003/10/08 16:19:33 hook Exp $
*****************************************************************************/

/*
subscripting convention:
|T[ 0][0]   T[+1][0]                                                         |
|T[-1][1]   T[ 0][1]   T[+1][1]                                              |
|           T[-1][2]   T[ 0][2]   T[+1][2]                                   |
|                          .          .          .                           |
|                                     .          .          .                |
|                                            T[-1][m-1] T[ 0][m-1] T[+1][m-2]|
|                                                       T[-1][m-1] T[ 0][m-1]|
where nml_extent m = nml_f3_extent(pT);				*/

/* Include vector type definitions.					*/
#include<nml_fvector.h>

typedef nml_fscalar*	nml_f3bands;

#ifdef	NML_INLINE
/* Functions								*/
inline static
nml_extent		(nml_f3_extent)(const nml_f3bands* pT) {
  return *((const nml_extent*)(pT - 2));
  }

inline static
nml_f3bands*	(nml_f3_fill)(
    nml_f3bands* pT, const nml_fscalar x) {
  nml_fv_fill(pT[+1], x);
  nml_fv_fill(pT[ 0], x);
  nml_fv_fill(pT[-1], x);
  return pT;
  }

inline static
nml_f3bands*	(nml_f3_smul)(
    nml_f3bands* pT, const nml_fscalar x) {
  nml_fv_smul(pT[+1], x);
  nml_fv_smul(pT[ 0], x);
  nml_fv_smul(pT[-1], x);
  return pT;
  }

inline static
nml_f3bands*	(nml_f3_sdiv)(
    nml_f3bands* pT, const nml_fscalar x) {
  nml_fv_sdiv(pT[+1], x);
  nml_fv_sdiv(pT[ 0], x);
  nml_fv_sdiv(pT[-1], x);
  return pT;
  }

inline static
nml_f3bands*	(nml_f3_sadd)(
    nml_f3bands* pT, const nml_fscalar x) {
  nml_fv_sadd(pT[+1], x);
  nml_fv_sadd(pT[ 0], x);
  nml_fv_sadd(pT[-1], x);
  return pT;
  }

inline static
nml_f3bands*	(nml_f3_ssub)(
    nml_f3bands* pT, const nml_fscalar x) {
  nml_fv_ssub(pT[+1], x);
  nml_fv_ssub(pT[ 0], x);
  nml_fv_ssub(pT[-1], x);
  return pT;
  }

inline static
nml_f3bands*	(nml_f3_3cpy)(
    nml_f3bands* pT, const nml_f3bands* pU) {
  nml_fv_vcpy(pT[+1], pU[+1]);
  nml_fv_vcpy(pT[ 0], pU[ 0]);
  nml_fv_vcpy(pT[-1], pU[-1]);
  return pT;
  }

/* Constructor								*/
inline static
nml_f3bands*	(nml_f3_new)(nml_extent n) {
  nml_fscalar**	pT
    = (nml_fscalar**)nml_calloc(4*sizeof(nml_fscalar*));
  if (NULL == pT) {
    nml_message("In function nml_f3_new(nml_extent):\n"
      "unable to allocate memory.");
    }
  else {		/* (NULL != pT)		*/
    *((nml_extent*)pT) = n;	++pT;
    *pT = nml_fv_new(n);	++pT;	/* T[-1]			*/
    *pT = nml_fv_new(n);	++pT;	/* T[ 0]			*/
    *pT = nml_fv_new(n);	--pT;	/* T[+1]			*/
    }
  return pT;
  }
inline static
nml_f3bands*	(nml_f3_clone)(const nml_f3bands* pU) {
  return nml_f3_3cpy(nml_f3_new(nml_f3_extent(pU)), pU);
  }

/* Destructor								*/
inline static
void			(nml_f3_delete)(nml_f3bands* pT) {
  if (NULL != pT) {
    ++pT; nml_fv_delete(*pT);	/* T[+1]			*/
    --pT; nml_fv_delete(*pT);	/* T[ 0]			*/
    --pT; nml_fv_delete(*pT);	/* T[-1]			*/
    --pT; nml_free(pT);
    }
  }

/* Reconstructor							*/
inline static
nml_f3bands*	(nml_f3_resize)(nml_f3bands* *ppT,
    nml_extent n) {
  if (NULL != ppT) {
    (nml_f3_delete)(*ppT);
    *ppT = (nml_f3_new)(n);
    return *ppT;
    }
  else {
    nml_message(
	"In function nml_f3_resize(nml_f3bands**, nml_extent):\n"
	"invalid pointer to nml_f3bands*.");
    return NULL;
    }
  }
#else /*NML_INLINE	*/
/* Functions								*/
nml_extent		(nml_f3_extent)(const nml_f3bands*);
nml_f3bands*	(nml_f3_fill)(
    nml_f3bands*, const nml_fscalar);
nml_f3bands*	(nml_f3_smul)(
    nml_f3bands*, const nml_fscalar);
nml_f3bands*	(nml_f3_sdiv)(
    nml_f3bands*, const nml_fscalar);
nml_f3bands*	(nml_f3_sadd)(
    nml_f3bands*, const nml_fscalar);
nml_f3bands*	(nml_f3_ssub)(
    nml_f3bands*, const nml_fscalar);
nml_f3bands*	(nml_f3_3cpy)(
    nml_f3bands*, const nml_f3bands*);
nml_fvector*	(nml_fv_v3dot)(
    nml_fvector*, const nml_fvector*, const nml_f3bands*);
nml_fvector*	(nml_fv_solveTridiagonal)(
    nml_fvector*, const nml_f3bands*, const nml_fvector*);
/* Constructor								*/
nml_f3bands*	(nml_f3_new)(nml_extent);
nml_f3bands*	(nml_f3_clone)(const nml_f3bands* pU);
/* Destructor								*/
void			(nml_f3_delete)(nml_f3bands*);
/* Reconstructor							*/
nml_f3bands*	(nml_f3_resize)(nml_f3bands**,
			    nml_extent);
#endif/*NML_INLINE	*/

int			(nml_f3_fprintf)(FILE*,
    const nml_f3bands*, int, unsigned int, nml_fmtflags, unsigned int);

#ifdef	NML_MACROS
#define nml_f3_extent(pT) (*((const nml_extent*)((pT) - 2)))
#define nml_f3_clone(pT) \
	nml_f3_3cpy(nml_f3_new(nml_f3_extent(pT)), (pT))
#endif/*NML_MACROS	*/

#endif /* _nml_f3bands_h */
