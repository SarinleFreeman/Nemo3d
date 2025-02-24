#ifndef _nml_fc3bands_h
#define _nml_fc3bands_h 1

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
where nml_extent m = nml_fc3_extent(pT);				*/

/* Include float precision real 3bands type definitions.		*/
#include<nml_f3bands.h>

/* Include vector type definitions.					*/
#include<nml_fcvector.h>

typedef nml_fcscalar*	nml_fc3bands;

#ifdef	NML_INLINE
/* Functions								*/
inline static
nml_extent		(nml_fc3_extent)(const nml_fc3bands* pT) {
  return *((const nml_extent*)(pT - 2));
  }

inline static
nml_fc3bands*	(nml_fc3_fill)(
    nml_fc3bands* pT, const nml_fcscalar x) {
  nml_fcv_fill(pT[+1], x);
  nml_fcv_fill(pT[ 0], x);
  nml_fcv_fill(pT[-1], x);
  return pT;
  }

inline static
nml_fc3bands*	(nml_fc3_smul)(
    nml_fc3bands* pT, const nml_fcscalar x) {
  nml_fcv_smul(pT[+1], x);
  nml_fcv_smul(pT[ 0], x);
  nml_fcv_smul(pT[-1], x);
  return pT;
  }

inline static
nml_fc3bands*	(nml_fc3_sdiv)(
    nml_fc3bands* pT, const nml_fcscalar x) {
  nml_fcv_sdiv(pT[+1], x);
  nml_fcv_sdiv(pT[ 0], x);
  nml_fcv_sdiv(pT[-1], x);
  return pT;
  }

inline static
nml_fc3bands*	(nml_fc3_sadd)(
    nml_fc3bands* pT, const nml_fcscalar x) {
  nml_fcv_sadd(pT[+1], x);
  nml_fcv_sadd(pT[ 0], x);
  nml_fcv_sadd(pT[-1], x);
  return pT;
  }

inline static
nml_fc3bands*	(nml_fc3_ssub)(
    nml_fc3bands* pT, const nml_fcscalar x) {
  nml_fcv_ssub(pT[+1], x);
  nml_fcv_ssub(pT[ 0], x);
  nml_fcv_ssub(pT[-1], x);
  return pT;
  }

inline static
nml_fc3bands*	(nml_fc3_3cpy)(
    nml_fc3bands* pT, const nml_fc3bands* pU) {
  nml_fcv_vcpy(pT[+1], pU[+1]);
  nml_fcv_vcpy(pT[ 0], pU[ 0]);
  nml_fcv_vcpy(pT[-1], pU[-1]);
  return pT;
  }

/* Constructor								*/
inline static
nml_fc3bands*	(nml_fc3_new)(nml_extent n) {
  nml_fcscalar**	pT
    = (nml_fcscalar**)nml_calloc(4*sizeof(nml_fcscalar*));
  if (NULL == pT) {
    nml_message("In function nml_fc3_new(nml_extent):\n"
      "unable to allocate memory.");
    }
  else {		/* (NULL != pT)		*/
    *((nml_extent*)pT) = n;	++pT;
    *pT = nml_fcv_new(n);	++pT;	/* T[-1]			*/
    *pT = nml_fcv_new(n);	++pT;	/* T[ 0]			*/
    *pT = nml_fcv_new(n);	--pT;	/* T[+1]			*/
    }
  return pT;
  }
inline static
nml_fc3bands*	(nml_fc3_clone)(const nml_fc3bands* pU) {
  return nml_fc3_3cpy(nml_fc3_new(nml_fc3_extent(pU)), pU);
  }

/* Destructor								*/
inline static
void			(nml_fc3_delete)(nml_fc3bands* pT) {
  if (NULL != pT) {
    ++pT; nml_fcv_delete(*pT);	/* T[+1]			*/
    --pT; nml_fcv_delete(*pT);	/* T[ 0]			*/
    --pT; nml_fcv_delete(*pT);	/* T[-1]			*/
    --pT; nml_free(pT);
    }
  }

/* Reconstructor							*/
inline static
nml_fc3bands*	(nml_fc3_resize)(nml_fc3bands* *ppT,
    nml_extent n) {
  if (NULL != ppT) {
    (nml_fc3_delete)(*ppT);
    *ppT = (nml_fc3_new)(n);
    return *ppT;
    }
  else {
    nml_message(
	"In function nml_fc3_resize(nml_fc3bands**, nml_extent):\n"
	"invalid pointer to nml_fc3bands*.");
    return NULL;
    }
  }
#else /*NML_INLINE	*/
/* Functions								*/
nml_extent		(nml_fc3_extent)(const nml_fc3bands*);
nml_fc3bands*	(nml_fc3_fill)(
    nml_fc3bands*, const nml_fcscalar);
nml_fc3bands*	(nml_fc3_smul)(
    nml_fc3bands*, const nml_fcscalar);
nml_fc3bands*	(nml_fc3_sdiv)(
    nml_fc3bands*, const nml_fcscalar);
nml_fc3bands*	(nml_fc3_sadd)(
    nml_fc3bands*, const nml_fcscalar);
nml_fc3bands*	(nml_fc3_ssub)(
    nml_fc3bands*, const nml_fcscalar);
nml_fc3bands*	(nml_fc3_3cpy)(
    nml_fc3bands*, const nml_fc3bands*);
nml_fcvector*	(nml_fcv_v3dot)(
    nml_fcvector*, const nml_fcvector*, const nml_fc3bands*);
nml_fcvector*	(nml_fcv_solveTridiagonal)(
    nml_fcvector*, const nml_fc3bands*, const nml_fcvector*);
/* Constructor								*/
nml_fc3bands*	(nml_fc3_new)(nml_extent);
nml_fc3bands*	(nml_fc3_clone)(const nml_fc3bands* pU);
/* Destructor								*/
void			(nml_fc3_delete)(nml_fc3bands*);
/* Reconstructor							*/
nml_fc3bands*	(nml_fc3_resize)(nml_fc3bands**,
			    nml_extent);
#endif/*NML_INLINE	*/

int			(nml_fc3_fprintf)(FILE*,
    const nml_fc3bands*, int, unsigned int, nml_fmtflags, unsigned int);

#ifdef	NML_MACROS
#define nml_fc3_extent(pT) (*((const nml_extent*)((pT) - 2)))
#define nml_fc3_clone(pT) \
	nml_fc3_3cpy(nml_fc3_new(nml_fc3_extent(pT)), (pT))
#endif/*NML_MACROS	*/

#endif /* _nml_fc3bands_h */
