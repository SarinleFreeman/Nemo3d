#ifndef _nml_ivector_h
#define _nml_ivector_h 1

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
Modified by Marek J. Korkusinski, November 2004
*****************************************************************************
$Header: /repo/nml/src/vector/vector.hP,v 1.6 2004/12/02 12:11:55 marek Exp $
*****************************************************************************/


#include<stdlib.h>

/* Include scalar type definitions.					*/
#include<nml_iscalar.h>

typedef nml_iscalar	nml_ivector;

#define EXT_MEM_HANDLING

#ifdef EXT_MEM_HANDLING
typedef struct 
{
   nml_extent noElements ;
   double start ;
} nml_idarray ;

static nml_extent nml_ioffsetLength = ( (nml_extent) 
                        ((char *) &((nml_idarray *) 0)->start)) ;
#endif /*   EXT_MEM_HANDLING   */



#ifdef	NML_INLINE
/* Functions								*/
inline static
nml_extent		(nml_iv_extent)(const nml_ivector* pv) {
#ifdef EXT_MEM_HANDLING
   return(  
            ( (nml_idarray *)
                     (    ((char *)pv) - nml_ioffsetLength     )
            ) -> noElements
         ) ;
#else   /* EXT_MEM_HANDLING */
  return *((const nml_extent*)(pv - 1));
#endif   /*   EXT_MEM_HANDLING   */
  }

inline static
nml_ivector*	(nml_iv_fill)(nml_ivector* pv,
    nml_iscalar x) {
  nml_extent		n = nml_iv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    pv[j] = x;
  return pv;
  }

inline static
nml_ivector*	(nml_iv_smul)(nml_ivector* pv,
    nml_iscalar x) {
  nml_extent		n = nml_iv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_iMUL(pv[j], pv[j], x);
  return pv;
  }

inline static
nml_ivector*	(nml_iv_sdiv)(nml_ivector* pv,
    nml_iscalar x) {
  nml_extent		n = nml_iv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_iDIV(pv[j], pv[j], x);
  return pv;
  }

inline static
nml_ivector*	(nml_iv_sadd)(nml_ivector* pv,
    nml_iscalar x) {
  nml_extent		n = nml_iv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_iSUB(pv[j], pv[j], x);
  return pv;
  }

inline static
nml_ivector*	(nml_iv_ssub)(nml_ivector* pv,
    nml_iscalar x) {
  nml_extent		n = nml_iv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_iSUB(pv[j], pv[j], x);
  return pv;
  }

inline static
nml_iscalar		(nml_iv_sminmag)(const nml_ivector* pv) {
  nml_extent		n = nml_iv_extent(pv);
  if (0 < n) {
    nml_iscalar	y = NML_iABS(pv[0]);
    nml_offset		j = 0;
    for (j = 0; j < n; ++j) {
      const
      nml_iscalar	t = NML_iABS(pv[j]);
      if (t < y) y = t;
      }
    return y;
    }
#ifdef	NML_ERRANT
  nml_message("In function nml_iv_sminmag(const nml_ivector*),\n"
      "zero vector extent.");
#endif/*NML_ERRANT		*/
  return (nml_iscalar)0;
  }

inline static
nml_iscalar		(nml_iv_smaxmag)(const nml_ivector* pv) {
  nml_extent		n = nml_iv_extent(pv);
  if (0 < n) {
    nml_iscalar	y = NML_iABS(pv[0]);
    nml_offset		j = 0;
    for (j = 1; j < n; ++j) {
      const
      nml_iscalar	t = NML_iABS(pv[j]);
      if (y < t) y = t;
      }
    return y;
    }
#ifdef	NML_ERRANT
  nml_message("In function nml_iv_smaxmag(const nml_ivector*),\n"
      "zero vector extent.");
#endif/*NML_ERRANT		*/
  return (nml_iscalar)0;
  }

inline static
nml_iscalar		(nml_iv_ssqr)(const nml_ivector* pv) {
  nml_iscalar		t = 0;
  nml_extent		n = nml_iv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    t += NML_iSQR(pv[j]);
  return t;
  }

inline static
nml_ivector*	(nml_iv_vcpy)(nml_ivector* pv,
    const nml_ivector* pw) {
  nml_extent		n = nml_iv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    pv[j] = pw[j];
  return pv;
  }

/* Constructor								*/

#ifdef EXT_MEM_HANDLING
inline static
nml_ivector*	(nml_iv_new)(nml_extent n) {
   nml_idarray * p ;
   nml_extent * pointerCopy ;

   p = ( nml_idarray * ) nml_calloc(nml_ioffsetLength +
                                       n * sizeof(nml_iscalar) , 1 ) ;
  if (NULL == p) {
    nml_message("In function nml_iv_new(nml_extent):\n"
      "unable to allocate memory.");
       }
  else {
      pointerCopy = (nml_extent * ) p ;
      *pointerCopy = n ;
       } 
   return( (nml_ivector*) &p->start  ) ;
  }
#else  /*   EXT_MEM_HANDLING   */
inline static
nml_ivector*	(nml_iv_new)(nml_extent n) {
     nml_iscalar*	p
       = (nml_iscalar*)nml_calloc(1+n, sizeof(nml_iscalar));
  if (NULL == p) {
    nml_message("In function nml_iv_new(nml_extent):\n"
      "unable to allocate memory.");
    }
  else {
       *((nml_extent*)p)	= n; ++p;
    }
  return p;
  }
#endif  /*   EXT_MEM_HANDLING   */

inline static
nml_ivector*	(nml_iv_clone)(const nml_ivector* pv) {
  nml_extent		n = nml_iv_extent(pv);
  return nml_iv_vcpy(nml_iv_new(n), pv);
  }

/* Destructor								*/
#ifdef EXT_MEM_HANDLING
inline static
void			(nml_iv_delete)(nml_ivector* pv) {
   nml_idarray * pointerToDarray ;
	if ( NULL != pv ) {
   pointerToDarray = 
   ( nml_idarray * ) ( (char*) pv - nml_ioffsetLength ) ;
   nml_free( pointerToDarray ) ;
	}
  }
#else  /*   EXT_MEM_HANDLING   */
inline static
void			(nml_iv_delete)(nml_ivector* pv) {
  nml_free(pv - 1);
  }
#endif  /*   EXT_MEM_HANDLING   */

/* Reconstructor							*/
inline static
nml_ivector*	(nml_iv_resize)(nml_ivector* *ppv, nml_extent n) {
  if (NULL != ppv) {
    (nml_iv_delete)(*ppv);
    *ppv = (nml_iv_new)(n);
    return *ppv;
    }
  else {
    nml_message(
"In function nml_iv_resize(nml_ivector**, nml_extent):\n"
"invalid pointer to nml_ivector*.");
    return NULL;
    }
  }

#else /*NML_INLINE	*/
/* Functions								*/
nml_extent		(nml_iv_extent)(const nml_ivector*);
nml_ivector*	(nml_iv_fill)(nml_ivector*, nml_iscalar);
nml_ivector*	(nml_iv_smul)(nml_ivector*, nml_iscalar);
nml_ivector*	(nml_iv_sdiv)(nml_ivector*, nml_iscalar);
nml_ivector*	(nml_iv_sadd)(nml_ivector*, nml_iscalar);
nml_ivector*	(nml_iv_ssub)(nml_ivector*, nml_iscalar);
nml_iscalar           (nml_iv_sminmag)(const nml_ivector*);
nml_iscalar           (nml_iv_smaxmag)(const nml_ivector*);
nml_iscalar		(nml_iv_ssqr)(const nml_ivector*);
nml_ivector*	(nml_iv_vcpy)(
			nml_ivector*, const nml_ivector* pw);
/* Constructor								*/
nml_ivector*	(nml_iv_new)(nml_extent);
nml_ivector*	(nml_iv_clone)(const nml_ivector*);
/* Destructor								*/
void			(nml_iv_delete)(nml_ivector*);
/* Reconstructor							*/
nml_ivector*	(nml_iv_resize)(nml_ivector**, nml_extent);
#endif/*NML_INLINE	*/

int			(nml_iv_fprintf)(FILE*,
    const nml_ivector*, int, unsigned int, nml_fmtflags, unsigned int);

#ifdef	NML_MACROS

#ifdef EXT_MEM_HANDLING
#define nml_iv_extent(pv) (((nml_idarray *) \
   (((char *)pv) - nml_ioffsetLength)) -> noElements )
#else   /* EXT_MEM_HANDLING */
#define nml_iv_extent(pv) (*((const nml_extent*)((pv) - 1)))
#endif   /*   EXT_MEM_HANDLING   */

#define nml_iv_clone(pv) \
nml_iv_vcpy(nml_iv_new(nml_iv_extent(pv)), (pv))
#endif/*NML_MACROS	*/

#endif/* _nml_ivector_h */
