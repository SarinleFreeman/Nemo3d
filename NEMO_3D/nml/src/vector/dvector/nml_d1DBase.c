
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
$Header: /repo/nml/src/vector/1DBase.cP,v 1.4 2004/11/26 22:56:42 marek Exp $
*****************************************************************************/


#undef NML_INLINE
#undef NML_MACROS
#include<nml_d1DBase.h>

/* Functions								*/

#ifdef EXT_MEM_HANDLING
nml_index		(nml_d1D_lower)(const nml_d1DBase* pv) {
   return(  
            ( (nml_d1D_darray *)
                     (    ((char *)pv) - nml_d1D_offsetLength     )
            ) -> lower
         ) ;
  }
#define nml_d1D_lower(pv) (((nml_d1D_darray *) \
     (((char *)pv) - nml_d1D_offsetLength )  ) -> lower  )

nml_index		(nml_d1D_upper)(const nml_d1DBase* pv) {
   return(  
            ( (nml_d1D_darray *)
                     (    ((char *)pv) - nml_d1D_offsetLength     )
            ) -> upper
         ) ;
  }
#define nml_d1D_upper(pv) (((nml_d1D_darray *) \
     (((char *)pv) - nml_d1D_offsetLength )  ) -> upper  )

#else /*   EXT_MEM_HANDLING   */

nml_index		(nml_d1D_lower)(const nml_d1DBase* pv) {
  return *((const nml_index*)(pv - 1));
  }
#define nml_d1D_lower(pv) (*((const nml_index*)((pv) - 1)))
nml_index		(nml_d1D_upper)(const nml_d1DBase* pv) {
  return *((const nml_index*)(pv - 2));
  }
#define nml_d1D_upper(pv) (*((const nml_index*)((pv) - 2)))
#endif  /*   EXT_MEM_HANDLING   */

nml_extent		(nml_d1D_extent)(const nml_d1DBase* pv) {
  return 1+nml_d1D_upper(pv)-nml_d1D_lower(pv);
  }
#define nml_d1D_extent(pv) (1+nml_d1D_upper(pv)-nml_d1D_lower(pv))

int
(nml_d1D_fprintf)(FILE* stream, const nml_d1DBase* pv, 
  int width, unsigned int precision, nml_fmtflags flags, unsigned int columns) {
  int		characters;
  int		total = 0;
  const
  nml_extent	n = nml_d1D_extent(pv);
  nml_offset	j = 0;
  columns = (0 < columns)? columns: 4;
  for (j = 0; j < n; ++j) {
    if (0 < j) {
      characters
	= (0 == j%columns)? fprintf(stream, "\n"): fprintf(stream, " ");
      if (0 <= characters) { total += characters; } else { return -1; }
      }
    characters = nml_dfprintf(
      stream, pv[j], width, precision, flags);
    if (0 <= characters) { total += characters; } else { return -1; }
    }
  characters = fprintf(stream, "\n");
  if (0 <= characters) { total += characters; } else { return -1; }
  return total;
  }



/* Constructor								*/
#ifdef EXT_MEM_HANDLING
nml_d1DBase*	(nml_d1D_new)(nml_index lb, nml_index ub) {
   nml_d1D_darray * p = NULL ;
   nml_index * pointerCopy ;
   nml_index		n = 1 + ub - lb;
  if (n < 0) {
    nml_message("In function nml_d1D_new(nml_index. nml_index):\n"
      "extent is negative.");
    }
  else {		/* (0 <= n)		*/
    p = ( nml_d1D_darray * ) nml_calloc(nml_d1D_offsetLength +
                                       n * sizeof(nml_dscalar) , 1 ) ;
    if (NULL == p) {
      nml_message("In function nml_d1D_new(nml_index. nml_index):\n"
	"unable to allocate memory.");
      }
    else {		/* (NULL != p)		*/
      pointerCopy = ( nml_index * ) p ;
      *pointerCopy = ub ; 
      ++pointerCopy ;
      *pointerCopy = lb ;
      }
    }
  return( (nml_dscalar*)  &p->start ) ;
  }
#else /*   EXT_MEM_HANDLING   */

nml_d1DBase*	(nml_d1D_new)(nml_index lb, nml_index ub) {
  nml_dscalar*	p = NULL;
  nml_index		n = 1 + ub - lb;
  if (n < 0) {
    nml_message("In function nml_d1D_new(nml_index, nml_index):\n"
      "extent is negative.");
    }
  else {		/* (0 <= n)		*/
    p = (nml_dscalar*)nml_calloc(2+n, sizeof(nml_dscalar));
    if (NULL == p) {
      nml_message("In function nml_d1D_new(nml_index, nml_index):\n"
	"unable to allocate memory.");
      }
    else {		/* (NULL != p)		*/
      *((nml_index*)p)	= ub; ++p;
      *((nml_index*)p)	= lb; ++p;
      }
    }
  return p;
  }
#endif /*   EXT_MEM_HANDLING   */

/* Destructor								*/

#ifdef EXT_MEM_HANDLING
void			(nml_d1D_delete)(nml_d1DBase* pv) {
   nml_d1D_darray * pointerToDarray ;
	if (NULL != pv) {
   pointerToDarray = 
   ( nml_d1D_darray * ) ( (char*) pv - nml_d1D_offsetLength ) ;
   nml_free( pointerToDarray ) ;
	}
  }
#else /*   EXT_MEM_HANDLING   */
void			(nml_d1D_delete)(nml_d1DBase* pv) {
  nml_free(pv - 2);
  }
#endif  /*   EXT_MEM_HANDLING   */


/* Reconstructor							*/
nml_d1DBase*	(nml_d1D_resize)(nml_d1DBase* *ppv,
    nml_index lb, nml_index ub) {
  if (NULL != ppv) {
    (nml_d1D_delete)(*ppv);
    *ppv = (nml_d1D_new)(lb, ub);
    return *ppv;
    }
  else {		/* (NULL == ppv)	*/
    nml_message(
"In function nml_d1D_resize(nml_d1DBase**, nml_index, nml_index):\n"
"invalid pointer to nml_d1DBase*.");
    return NULL;
    }
  }

