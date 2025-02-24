#ifndef _nml_s1DBase_h
#define _nml_s1DBase_h 1

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
$Header: /repo/nml/src/vector/1DBase.hP,v 1.6 2004/12/02 12:11:55 marek Exp $
*****************************************************************************/


#include<stdlib.h>

/* Include scalar type definitions.					*/
#include<nml_sscalar.h>

typedef nml_sscalar	nml_s1DBase;

#define EXT_MEM_HANDLING

#ifdef EXT_MEM_HANDLING
typedef struct 
{
   nml_index upper ;
   nml_index lower ;
   double start ;
} nml_s1D_darray ;

static nml_extent nml_s1D_offsetLength = ( (nml_extent) 
                        ((char *) &((nml_s1D_darray *) 0)->start)) ;
#endif /*   EXT_MEM_HANDLING   */


#ifdef	NML_INLINE
/* Functions								*/
#ifdef EXT_MEM_HANDLING
inline static
nml_index		(nml_s1D_lower)(const nml_s1DBase* pv) {
   return(  
            ( (nml_s1D_darray *)
                     (    ((char *)pv) - nml_s1D_offsetLength     )
            ) -> lower
         ) ;
  }
inline static
nml_index		(nml_s1D_upper)(const nml_s1DBase* pv) {
   return(  
            ( (nml_s1D_darray *)
                     (    ((char *)pv) - nml_s1D_offsetLength     )
            ) -> upper
         ) ;
  }

#else /*   EXT_MEM_HANDLING   */

inline static
nml_index		(nml_s1D_lower)(const nml_s1DBase* pv) {
  return *((const nml_index*)(pv - 1));
  }
inline static
nml_index		(nml_s1D_upper)(const nml_s1DBase* pv) {
  return *((const nml_index*)(pv - 2));
  }
#endif  /*   EXT_MEM_HANDLING   */

inline static
nml_extent		(nml_s1D_extent)(const nml_s1DBase* pv) {
  return 1+nml_s1D_upper(pv)-nml_s1D_lower(pv);
  }

/* Constructor								*/
#ifdef EXT_MEM_HANDLING
inline static
nml_s1DBase*	(nml_s1D_new)(nml_index lb, nml_index ub) {
   nml_s1D_darray * p = NULL ;
   nml_index * pointerCopy ;
   nml_index		n = 1 + ub - lb;
  if (n < 0) {
    nml_message("In function nml_s1D_new(nml_index. nml_index):\n"
      "extent is negative.");
    }
  else {		/* (0 <= n)		*/
    p = ( nml_s1D_darray * ) nml_calloc(nml_s1D_offsetLength +
                                       n * sizeof(nml_sscalar) , 1 ) ;
    if (NULL == p) {
      nml_message("In function nml_s1D_new(nml_index. nml_index):\n"
	"unable to allocate memory.");
      }
    else {		/* (NULL != p)		*/
      pointerCopy = ( nml_index * ) p ;
      *pointerCopy = ub ; 
      ++pointerCopy ;
      *pointerCopy = lb ;
      }
    }
  return( (nml_sscalar*)  &p->start ) ;
  }
#else /*   EXT_MEM_HANDLING   */
inline static
nml_s1DBase*	(nml_s1D_new)(nml_index lb, nml_index ub) {
  nml_sscalar*	p = NULL;
  nml_index		n = 1 + ub - lb;
  if (n < 0) {
    nml_message("In function nml_s1D_new(nml_index. nml_index):\n"
      "extent is negative.");
    }
  else {		/* (0 <= n)		*/
    p = (nml_sscalar*)nml_calloc(2 + n, sizeof(nml_sscalar));
    if (NULL == p) {
      nml_message("In function nml_s1D_new(nml_index. nml_index):\n"
	"unable to allocate memory.");
      }
    else {		/* (NULL != p)		*/
      *((nml_index*)p)	= ub; ++p;
      *((nml_index*)p)	= lb; ++p;
      }
    }
  return p;
  }
#endif  /*   EXT_MEM_HANDLING   */

/* Destructor								*/
#ifdef EXT_MEM_HANDLING
inline static
void			(nml_s1D_delete)(nml_s1DBase* pv) {
   nml_s1D_darray * pointerToDarray ;
	if (NULL != pv ) {
   pointerToDarray = 
   ( nml_s1D_darray * ) ( (char*) pv - nml_s1D_offsetLength ) ;
   nml_free( pointerToDarray ) ;
	}
  }
#else /*   EXT_MEM_HANDLING   */
inline static
void			(nml_s1D_delete)(nml_s1DBase* pv) {
  nml_free(pv - 2);
  }
#endif  /*   EXT_MEM_HANDLING   */

/* Reconstructor							*/
inline static
nml_s1DBase*	(nml_s1D_resize)(nml_s1DBase* *ppv,
    nml_index lb, nml_index ub) {
  if (NULL != ppv) {
    (nml_s1D_delete)(*ppv);
    *ppv = (nml_s1D_new)(lb, ub);
    return *ppv;
    }
  else {		/* (NULL == ppv)	*/
    nml_message(
"In function nml_s1D_resize(nml_s1DBase**, nml_index, nml_index):\n"
"invalid pointer to nml_s1DBase*.");
    return NULL;
    }
  }
#else /*NML_INLINE	*/
/* Functions								*/
nml_index		(nml_s1D_lower)(const nml_s1DBase*);
nml_index		(nml_s1D_upper)(const nml_s1DBase*);
nml_extent		(nml_s1D_extent)(const nml_s1DBase*);
/* Constructor								*/
nml_s1DBase*	(nml_s1D_new)(nml_index, nml_index);
/* Destructor								*/
void			(nml_s1D_delete)(nml_s1DBase*);
/* Reconstructor							*/
nml_s1DBase*	(nml_s1D_resize)(nml_s1DBase**,
    nml_index, nml_index);
#endif/*NML_INLINE	*/

int			(nml_s1D_fprintf)(FILE*,
    const nml_s1DBase*, int, unsigned int, nml_fmtflags, unsigned int);

#ifdef	NML_MACROS

#ifdef EXT_MEM_HANDLING
#define nml_s1D_lower(pv) (((nml_s1D_darray *) \
     (((char *)pv) - nml_s1D_offsetLength )  ) -> lower  )
#define nml_s1D_upper(pv) (((nml_s1D_darray *) \
     (((char *)pv) - nml_s1D_offsetLength )  ) -> upper  )
#else /*   EXT_MEM_HANDLING   */
#define nml_s1D_lower(pv) (*((const nml_index*)((pv) - 1)))
#define nml_s1D_upper(pv) (*((const nml_index*)((pv) - 2)))
#endif /*   EXT_MEM_HANDLING   */

#define nml_s1D_extent(pv) (1+nml_s1D_upper(pv)-nml_s1D_lower(pv))
#endif/*NML_MACROS	*/

#endif/* _nml_s1DBase_h */
