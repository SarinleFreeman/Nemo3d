#ifndef _nml_dc3bands_h
#define _nml_dc3bands_h 1

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
where nml_extent m = nml_dc3_extent(pT);				*/

/* Include double precision real 3bands type definitions.		*/
#include<nml_d3bands.h>

/* Include vector type definitions.					*/
#include<nml_dcvector.h>

typedef nml_dcscalar*	nml_dc3bands;

#ifdef	NML_INLINE
/* Functions								*/
inline static
nml_extent		(nml_dc3_extent)(const nml_dc3bands* pT) {
  return *((const nml_extent*)(pT - 2));
  }

inline static
nml_dc3bands*	(nml_dc3_fill)(
    nml_dc3bands* pT, const nml_dcscalar x) {
  nml_dcv_fill(pT[+1], x);
  nml_dcv_fill(pT[ 0], x);
  nml_dcv_fill(pT[-1], x);
  return pT;
  }

inline static
nml_dc3bands*	(nml_dc3_smul)(
    nml_dc3bands* pT, const nml_dcscalar x) {
  nml_dcv_smul(pT[+1], x);
  nml_dcv_smul(pT[ 0], x);
  nml_dcv_smul(pT[-1], x);
  return pT;
  }

inline static
nml_dc3bands*	(nml_dc3_sdiv)(
    nml_dc3bands* pT, const nml_dcscalar x) {
  nml_dcv_sdiv(pT[+1], x);
  nml_dcv_sdiv(pT[ 0], x);
  nml_dcv_sdiv(pT[-1], x);
  return pT;
  }

inline static
nml_dc3bands*	(nml_dc3_sadd)(
    nml_dc3bands* pT, const nml_dcscalar x) {
  nml_dcv_sadd(pT[+1], x);
  nml_dcv_sadd(pT[ 0], x);
  nml_dcv_sadd(pT[-1], x);
  return pT;
  }

inline static
nml_dc3bands*	(nml_dc3_ssub)(
    nml_dc3bands* pT, const nml_dcscalar x) {
  nml_dcv_ssub(pT[+1], x);
  nml_dcv_ssub(pT[ 0], x);
  nml_dcv_ssub(pT[-1], x);
  return pT;
  }

inline static
nml_dc3bands*	(nml_dc3_3cpy)(
    nml_dc3bands* pT, const nml_dc3bands* pU) {
  nml_dcv_vcpy(pT[+1], pU[+1]);
  nml_dcv_vcpy(pT[ 0], pU[ 0]);
  nml_dcv_vcpy(pT[-1], pU[-1]);
  return pT;
  }

/* Constructor								*/
inline static
nml_dc3bands*	(nml_dc3_new)(nml_extent n) {
  nml_dcscalar**	pT
    = (nml_dcscalar**)nml_calloc(4*sizeof(nml_dcscalar*));
  if (NULL == pT) {
    nml_message("In function nml_dc3_new(nml_extent):\n"
      "unable to allocate memory.");
    }
  else {		/* (NULL != pT)		*/
    *((nml_extent*)pT) = n;	++pT;
    *pT = nml_dcv_new(n);	++pT;	/* T[-1]			*/
    *pT = nml_dcv_new(n);	++pT;	/* T[ 0]			*/
    *pT = nml_dcv_new(n);	--pT;	/* T[+1]			*/
    }
  return pT;
  }
inline static
nml_dc3bands*	(nml_dc3_clone)(const nml_dc3bands* pU) {
  return nml_dc3_3cpy(nml_dc3_new(nml_dc3_extent(pU)), pU);
  }

/* Destructor								*/
inline static
void			(nml_dc3_delete)(nml_dc3bands* pT) {
  if (NULL != pT) {
    ++pT; nml_dcv_delete(*pT);	/* T[+1]			*/
    --pT; nml_dcv_delete(*pT);	/* T[ 0]			*/
    --pT; nml_dcv_delete(*pT);	/* T[-1]			*/
    --pT; nml_free(pT);
    }
  }

/* Reconstructor							*/
inline static
nml_dc3bands*	(nml_dc3_resize)(nml_dc3bands* *ppT,
    nml_extent n) {
  if (NULL != ppT) {
    (nml_dc3_delete)(*ppT);
    *ppT = (nml_dc3_new)(n);
    return *ppT;
    }
  else {
    nml_message(
	"In function nml_dc3_resize(nml_dc3bands**, nml_extent):\n"
	"invalid pointer to nml_dc3bands*.");
    return NULL;
    }
  }
#else /*NML_INLINE	*/
/* Functions								*/
nml_extent		(nml_dc3_extent)(const nml_dc3bands*);
nml_dc3bands*	(nml_dc3_fill)(
    nml_dc3bands*, const nml_dcscalar);
nml_dc3bands*	(nml_dc3_smul)(
    nml_dc3bands*, const nml_dcscalar);
nml_dc3bands*	(nml_dc3_sdiv)(
    nml_dc3bands*, const nml_dcscalar);
nml_dc3bands*	(nml_dc3_sadd)(
    nml_dc3bands*, const nml_dcscalar);
nml_dc3bands*	(nml_dc3_ssub)(
    nml_dc3bands*, const nml_dcscalar);
nml_dc3bands*	(nml_dc3_3cpy)(
    nml_dc3bands*, const nml_dc3bands*);
nml_dcvector*	(nml_dcv_v3dot)(
    nml_dcvector*, const nml_dcvector*, const nml_dc3bands*);
nml_dcvector*	(nml_dcv_solveTridiagonal)(
    nml_dcvector*, const nml_dc3bands*, const nml_dcvector*);
/* Constructor								*/
nml_dc3bands*	(nml_dc3_new)(nml_extent);
nml_dc3bands*	(nml_dc3_clone)(const nml_dc3bands* pU);
/* Destructor								*/
void			(nml_dc3_delete)(nml_dc3bands*);
/* Reconstructor							*/
nml_dc3bands*	(nml_dc3_resize)(nml_dc3bands**,
			    nml_extent);
#endif/*NML_INLINE	*/

int			(nml_dc3_fprintf)(FILE*,
    const nml_dc3bands*, int, unsigned int, nml_fmtflags, unsigned int);

#ifdef	NML_MACROS
#define nml_dc3_extent(pT) (*((const nml_extent*)((pT) - 2)))
#define nml_dc3_clone(pT) \
	nml_dc3_3cpy(nml_dc3_new(nml_dc3_extent(pT)), (pT))
#endif/*NML_MACROS	*/

#endif /* _nml_dc3bands_h */
