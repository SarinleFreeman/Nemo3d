#ifndef _nml_fvector_h
#define _nml_fvector_h 1

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
#include<nml_fscalar.h>

typedef nml_fscalar	nml_fvector;

#define EXT_MEM_HANDLING

#ifdef EXT_MEM_HANDLING
typedef struct 
{
   nml_extent noElements ;
   double start ;
} nml_fdarray ;

static nml_extent nml_foffsetLength = ( (nml_extent) 
                        ((char *) &((nml_fdarray *) 0)->start)) ;
#endif /*   EXT_MEM_HANDLING   */



#ifdef	NML_INLINE
/* Functions								*/
inline static
nml_extent		(nml_fv_extent)(const nml_fvector* pv) {
#ifdef EXT_MEM_HANDLING
   return(  
            ( (nml_fdarray *)
                     (    ((char *)pv) - nml_foffsetLength     )
            ) -> noElements
         ) ;
#else   /* EXT_MEM_HANDLING */
  return *((const nml_extent*)(pv - 1));
#endif   /*   EXT_MEM_HANDLING   */
  }

inline static
nml_fvector*	(nml_fv_fill)(nml_fvector* pv,
    nml_fscalar x) {
  nml_extent		n = nml_fv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    pv[j] = x;
  return pv;
  }

inline static
nml_fvector*	(nml_fv_smul)(nml_fvector* pv,
    nml_fscalar x) {
  nml_extent		n = nml_fv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_FMUL(pv[j], pv[j], x);
  return pv;
  }

inline static
nml_fvector*	(nml_fv_sdiv)(nml_fvector* pv,
    nml_fscalar x) {
  nml_extent		n = nml_fv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_FDIV(pv[j], pv[j], x);
  return pv;
  }

inline static
nml_fvector*	(nml_fv_sadd)(nml_fvector* pv,
    nml_fscalar x) {
  nml_extent		n = nml_fv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_FSUB(pv[j], pv[j], x);
  return pv;
  }

inline static
nml_fvector*	(nml_fv_ssub)(nml_fvector* pv,
    nml_fscalar x) {
  nml_extent		n = nml_fv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_FSUB(pv[j], pv[j], x);
  return pv;
  }

inline static
nml_fscalar		(nml_fv_sminmag)(const nml_fvector* pv) {
  nml_extent		n = nml_fv_extent(pv);
  if (0 < n) {
    nml_fscalar	y = NML_FABS(pv[0]);
    nml_offset		j = 0;
    for (j = 0; j < n; ++j) {
      const
      nml_fscalar	t = NML_FABS(pv[j]);
      if (t < y) y = t;
      }
    return y;
    }
#ifdef	NML_ERRANT
  nml_message("In function nml_fv_sminmag(const nml_fvector*),\n"
      "zero vector extent.");
#endif/*NML_ERRANT		*/
  return (nml_fscalar)0;
  }

inline static
nml_fscalar		(nml_fv_smaxmag)(const nml_fvector* pv) {
  nml_extent		n = nml_fv_extent(pv);
  if (0 < n) {
    nml_fscalar	y = NML_FABS(pv[0]);
    nml_offset		j = 0;
    for (j = 1; j < n; ++j) {
      const
      nml_fscalar	t = NML_FABS(pv[j]);
      if (y < t) y = t;
      }
    return y;
    }
#ifdef	NML_ERRANT
  nml_message("In function nml_fv_smaxmag(const nml_fvector*),\n"
      "zero vector extent.");
#endif/*NML_ERRANT		*/
  return (nml_fscalar)0;
  }

inline static
nml_fscalar		(nml_fv_ssqr)(const nml_fvector* pv) {
  nml_fscalar		t = 0;
  nml_extent		n = nml_fv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    t += NML_FSQR(pv[j]);
  return t;
  }

inline static
nml_fscalar		(nml_fv_sabs)(const nml_fvector* pv) {
  return sqrt(nml_fv_ssqr(pv));
  }

inline static
nml_fvector*	(nml_fv_vcpy)(nml_fvector* pv,
    const nml_fvector* pw) {
  nml_extent		n = nml_fv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    pv[j] = pw[j];
  return pv;
  }

/* Constructor								*/

#ifdef EXT_MEM_HANDLING
inline static
nml_fvector*	(nml_fv_new)(nml_extent n) {
   nml_fdarray * p ;
   nml_extent * pointerCopy ;

   p = ( nml_fdarray * ) nml_calloc(nml_foffsetLength +
                                       n * sizeof(nml_fscalar) , 1 ) ;
  if (NULL == p) {
    nml_message("In function nml_fv_new(nml_extent):\n"
      "unable to allocate memory.");
       }
  else {
      pointerCopy = (nml_extent * ) p ;
      *pointerCopy = n ;
       } 
   return( (nml_fvector*) &p->start  ) ;
  }
#else  /*   EXT_MEM_HANDLING   */
inline static
nml_fvector*	(nml_fv_new)(nml_extent n) {
     nml_fscalar*	p
       = (nml_fscalar*)nml_calloc(1+n, sizeof(nml_fscalar));
  if (NULL == p) {
    nml_message("In function nml_fv_new(nml_extent):\n"
      "unable to allocate memory.");
    }
  else {
       *((nml_extent*)p)	= n; ++p;
    }
  return p;
  }
#endif  /*   EXT_MEM_HANDLING   */

inline static
nml_fvector*	(nml_fv_clone)(const nml_fvector* pv) {
  nml_extent		n = nml_fv_extent(pv);
  return nml_fv_vcpy(nml_fv_new(n), pv);
  }

/* Destructor								*/
#ifdef EXT_MEM_HANDLING
inline static
void			(nml_fv_delete)(nml_fvector* pv) {
   nml_fdarray * pointerToDarray ;
	if ( NULL != pv ) {
   pointerToDarray = 
   ( nml_fdarray * ) ( (char*) pv - nml_foffsetLength ) ;
   nml_free( pointerToDarray ) ;
	}
  }
#else  /*   EXT_MEM_HANDLING   */
inline static
void			(nml_fv_delete)(nml_fvector* pv) {
  nml_free(pv - 1);
  }
#endif  /*   EXT_MEM_HANDLING   */

/* Reconstructor							*/
inline static
nml_fvector*	(nml_fv_resize)(nml_fvector* *ppv, nml_extent n) {
  if (NULL != ppv) {
    (nml_fv_delete)(*ppv);
    *ppv = (nml_fv_new)(n);
    return *ppv;
    }
  else {
    nml_message(
"In function nml_fv_resize(nml_fvector**, nml_extent):\n"
"invalid pointer to nml_fvector*.");
    return NULL;
    }
  }

#else /*NML_INLINE	*/
/* Functions								*/
nml_extent		(nml_fv_extent)(const nml_fvector*);
nml_fvector*	(nml_fv_fill)(nml_fvector*, nml_fscalar);
nml_fvector*	(nml_fv_smul)(nml_fvector*, nml_fscalar);
nml_fvector*	(nml_fv_sdiv)(nml_fvector*, nml_fscalar);
nml_fvector*	(nml_fv_sadd)(nml_fvector*, nml_fscalar);
nml_fvector*	(nml_fv_ssub)(nml_fvector*, nml_fscalar);
nml_fscalar           (nml_fv_sminmag)(const nml_fvector*);
nml_fscalar           (nml_fv_smaxmag)(const nml_fvector*);
nml_fscalar		(nml_fv_ssqr)(const nml_fvector*);
nml_fscalar		(nml_fv_sabs)(const nml_fvector*);
nml_fvector*	(nml_fv_vcpy)(
			nml_fvector*, const nml_fvector* pw);
/* Constructor								*/
nml_fvector*	(nml_fv_new)(nml_extent);
nml_fvector*	(nml_fv_clone)(const nml_fvector*);
/* Destructor								*/
void			(nml_fv_delete)(nml_fvector*);
/* Reconstructor							*/
nml_fvector*	(nml_fv_resize)(nml_fvector**, nml_extent);
#endif/*NML_INLINE	*/

int			(nml_fv_fprintf)(FILE*,
    const nml_fvector*, int, unsigned int, nml_fmtflags, unsigned int);

#ifdef	NML_MACROS

#ifdef EXT_MEM_HANDLING
#define nml_fv_extent(pv) (((nml_fdarray *) \
   (((char *)pv) - nml_foffsetLength)) -> noElements )
#else   /* EXT_MEM_HANDLING */
#define nml_fv_extent(pv) (*((const nml_extent*)((pv) - 1)))
#endif   /*   EXT_MEM_HANDLING   */

#define nml_fv_sabs(pv) sqrt(nml_fv_ssqr(pv));
#define nml_fv_clone(pv) \
nml_fv_vcpy(nml_fv_new(nml_fv_extent(pv)), (pv))
#endif/*NML_MACROS	*/

#endif/* _nml_fvector_h */
