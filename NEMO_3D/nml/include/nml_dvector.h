#ifndef _nml_dvector_h
#define _nml_dvector_h 1

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
#include<nml_dscalar.h>

typedef nml_dscalar	nml_dvector;

#define EXT_MEM_HANDLING

#ifdef EXT_MEM_HANDLING
typedef struct 
{
   nml_extent noElements ;
   double start ;
} nml_ddarray ;

static nml_extent nml_doffsetLength = ( (nml_extent) 
                        ((char *) &((nml_ddarray *) 0)->start)) ;
#endif /*   EXT_MEM_HANDLING   */



#ifdef	NML_INLINE
/* Functions								*/
inline static
nml_extent		(nml_dv_extent)(const nml_dvector* pv) {
#ifdef EXT_MEM_HANDLING
   return(  
            ( (nml_ddarray *)
                     (    ((char *)pv) - nml_doffsetLength     )
            ) -> noElements
         ) ;
#else   /* EXT_MEM_HANDLING */
  return *((const nml_extent*)(pv - 1));
#endif   /*   EXT_MEM_HANDLING   */
  }

inline static
nml_dvector*	(nml_dv_fill)(nml_dvector* pv,
    nml_dscalar x) {
  nml_extent		n = nml_dv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    pv[j] = x;
  return pv;
  }

inline static
nml_dvector*	(nml_dv_smul)(nml_dvector* pv,
    nml_dscalar x) {
  nml_extent		n = nml_dv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_DMUL(pv[j], pv[j], x);
  return pv;
  }

inline static
nml_dvector*	(nml_dv_sdiv)(nml_dvector* pv,
    nml_dscalar x) {
  nml_extent		n = nml_dv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_DDIV(pv[j], pv[j], x);
  return pv;
  }

inline static
nml_dvector*	(nml_dv_sadd)(nml_dvector* pv,
    nml_dscalar x) {
  nml_extent		n = nml_dv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_DSUB(pv[j], pv[j], x);
  return pv;
  }

inline static
nml_dvector*	(nml_dv_ssub)(nml_dvector* pv,
    nml_dscalar x) {
  nml_extent		n = nml_dv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_DSUB(pv[j], pv[j], x);
  return pv;
  }

inline static
nml_dscalar		(nml_dv_sminmag)(const nml_dvector* pv) {
  nml_extent		n = nml_dv_extent(pv);
  if (0 < n) {
    nml_dscalar	y = NML_DABS(pv[0]);
    nml_offset		j = 0;
    for (j = 0; j < n; ++j) {
      const
      nml_dscalar	t = NML_DABS(pv[j]);
      if (t < y) y = t;
      }
    return y;
    }
#ifdef	NML_ERRANT
  nml_message("In function nml_dv_sminmag(const nml_dvector*),\n"
      "zero vector extent.");
#endif/*NML_ERRANT		*/
  return (nml_dscalar)0;
  }

inline static
nml_dscalar		(nml_dv_smaxmag)(const nml_dvector* pv) {
  nml_extent		n = nml_dv_extent(pv);
  if (0 < n) {
    nml_dscalar	y = NML_DABS(pv[0]);
    nml_offset		j = 0;
    for (j = 1; j < n; ++j) {
      const
      nml_dscalar	t = NML_DABS(pv[j]);
      if (y < t) y = t;
      }
    return y;
    }
#ifdef	NML_ERRANT
  nml_message("In function nml_dv_smaxmag(const nml_dvector*),\n"
      "zero vector extent.");
#endif/*NML_ERRANT		*/
  return (nml_dscalar)0;
  }

inline static
nml_dscalar		(nml_dv_ssqr)(const nml_dvector* pv) {
  nml_dscalar		t = 0;
  nml_extent		n = nml_dv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    t += NML_DSQR(pv[j]);
  return t;
  }

inline static
nml_dscalar		(nml_dv_sabs)(const nml_dvector* pv) {
  return sqrt(nml_dv_ssqr(pv));
  }

inline static
nml_dvector*	(nml_dv_vcpy)(nml_dvector* pv,
    const nml_dvector* pw) {
  nml_extent		n = nml_dv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    pv[j] = pw[j];
  return pv;
  }

/* Constructor								*/

#ifdef EXT_MEM_HANDLING
inline static
nml_dvector*	(nml_dv_new)(nml_extent n) {
   nml_ddarray * p ;
   nml_extent * pointerCopy ;

   p = ( nml_ddarray * ) nml_calloc(nml_doffsetLength +
                                       n * sizeof(nml_dscalar) , 1 ) ;
  if (NULL == p) {
    nml_message("In function nml_dv_new(nml_extent):\n"
      "unable to allocate memory.");
       }
  else {
      pointerCopy = (nml_extent * ) p ;
      *pointerCopy = n ;
       } 
   return( (nml_dvector*) &p->start  ) ;
  }
#else  /*   EXT_MEM_HANDLING   */
inline static
nml_dvector*	(nml_dv_new)(nml_extent n) {
     nml_dscalar*	p
       = (nml_dscalar*)nml_calloc(1+n, sizeof(nml_dscalar));
  if (NULL == p) {
    nml_message("In function nml_dv_new(nml_extent):\n"
      "unable to allocate memory.");
    }
  else {
       *((nml_extent*)p)	= n; ++p;
    }
  return p;
  }
#endif  /*   EXT_MEM_HANDLING   */

inline static
nml_dvector*	(nml_dv_clone)(const nml_dvector* pv) {
  nml_extent		n = nml_dv_extent(pv);
  return nml_dv_vcpy(nml_dv_new(n), pv);
  }

/* Destructor								*/
#ifdef EXT_MEM_HANDLING
inline static
void			(nml_dv_delete)(nml_dvector* pv) {
   nml_ddarray * pointerToDarray ;
	if ( NULL != pv ) {
   pointerToDarray = 
   ( nml_ddarray * ) ( (char*) pv - nml_doffsetLength ) ;
   nml_free( pointerToDarray ) ;
	}
  }
#else  /*   EXT_MEM_HANDLING   */
inline static
void			(nml_dv_delete)(nml_dvector* pv) {
  nml_free(pv - 1);
  }
#endif  /*   EXT_MEM_HANDLING   */

/* Reconstructor							*/
inline static
nml_dvector*	(nml_dv_resize)(nml_dvector* *ppv, nml_extent n) {
  if (NULL != ppv) {
    (nml_dv_delete)(*ppv);
    *ppv = (nml_dv_new)(n);
    return *ppv;
    }
  else {
    nml_message(
"In function nml_dv_resize(nml_dvector**, nml_extent):\n"
"invalid pointer to nml_dvector*.");
    return NULL;
    }
  }

#else /*NML_INLINE	*/
/* Functions								*/
nml_extent		(nml_dv_extent)(const nml_dvector*);
nml_dvector*	(nml_dv_fill)(nml_dvector*, nml_dscalar);
nml_dvector*	(nml_dv_smul)(nml_dvector*, nml_dscalar);
nml_dvector*	(nml_dv_sdiv)(nml_dvector*, nml_dscalar);
nml_dvector*	(nml_dv_sadd)(nml_dvector*, nml_dscalar);
nml_dvector*	(nml_dv_ssub)(nml_dvector*, nml_dscalar);
nml_dscalar           (nml_dv_sminmag)(const nml_dvector*);
nml_dscalar           (nml_dv_smaxmag)(const nml_dvector*);
nml_dscalar		(nml_dv_ssqr)(const nml_dvector*);
nml_dscalar		(nml_dv_sabs)(const nml_dvector*);
nml_dvector*	(nml_dv_vcpy)(
			nml_dvector*, const nml_dvector* pw);
/* Constructor								*/
nml_dvector*	(nml_dv_new)(nml_extent);
nml_dvector*	(nml_dv_clone)(const nml_dvector*);
/* Destructor								*/
void			(nml_dv_delete)(nml_dvector*);
/* Reconstructor							*/
nml_dvector*	(nml_dv_resize)(nml_dvector**, nml_extent);
#endif/*NML_INLINE	*/

int			(nml_dv_fprintf)(FILE*,
    const nml_dvector*, int, unsigned int, nml_fmtflags, unsigned int);

#ifdef	NML_MACROS

#ifdef EXT_MEM_HANDLING
#define nml_dv_extent(pv) (((nml_ddarray *) \
   (((char *)pv) - nml_doffsetLength)) -> noElements )
#else   /* EXT_MEM_HANDLING */
#define nml_dv_extent(pv) (*((const nml_extent*)((pv) - 1)))
#endif   /*   EXT_MEM_HANDLING   */

#define nml_dv_sabs(pv) sqrt(nml_dv_ssqr(pv));
#define nml_dv_clone(pv) \
nml_dv_vcpy(nml_dv_new(nml_dv_extent(pv)), (pv))
#endif/*NML_MACROS	*/

#endif/* _nml_dvector_h */
