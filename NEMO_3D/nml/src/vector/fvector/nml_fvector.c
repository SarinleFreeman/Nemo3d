
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
$Header: /repo/nml/src/vector/vector.cP,v 1.4 2004/11/26 22:56:42 marek Exp $
*****************************************************************************/



#undef NML_INLINE
#undef NML_MACROS
#include<nml_fvector.h>

/* Functions								*/
#ifdef EXT_MEM_HANDLING
nml_extent		(nml_fv_extent)(const nml_fvector* pv) {
   return(  
            ( (nml_fdarray *)
                     (    ((char *)pv) - nml_foffsetLength     )
            ) -> noElements
         ) ;
  }
#define nml_fv_extent(pv) (((nml_fdarray *) \
   (((char *)pv) - nml_foffsetLength)) -> noElements )
#else   /* EXT_MEM_HANDLING */
nml_extent		(nml_fv_extent)(const nml_fvector* pv) {
  return *((const nml_extent*)(pv - 1));
  }
#define nml_fv_extent(pv) (*((const nml_extent*)((pv) - 1)))
#endif   /* EXT_MEM_HANDLING */

nml_fvector*	(nml_fv_fill)(nml_fvector* pv,
    nml_fscalar x) {
  nml_extent		n = nml_fv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    pv[j] = x;
  return pv;
  }

nml_fvector*	(nml_fv_smul)(nml_fvector* pv,
    nml_fscalar x) {
  nml_extent		n = nml_fv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_FMUL(pv[j], pv[j], x);
  return pv;
  }

nml_fvector*	(nml_fv_sdiv)(nml_fvector* pv,
    nml_fscalar x) {
  nml_extent		n = nml_fv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_FDIV(pv[j], pv[j], x);
  return pv;
  }

nml_fvector*	(nml_fv_sadd)(nml_fvector* pv,
    nml_fscalar x) {
  nml_extent		n = nml_fv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_FSUB(pv[j], pv[j], x);
  return pv;
  }

nml_fvector*	(nml_fv_ssub)(nml_fvector* pv,
    nml_fscalar x) {
  nml_extent		n = nml_fv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_FSUB(pv[j], pv[j], x);
  return pv;
  }

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
#ifdef  NML_ERRANT
  nml_message("In function nml_fv_sminmag(const nml_fvector*),\n"
      "zero vector extent.");
#endif/*NML_ERRANT              */
  return (nml_fscalar)0;
  }

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
#ifdef  NML_ERRANT
  nml_message("In function nml_fv_smaxmag(const nml_fvector*),\n"
      "zero vector extent.");
#endif/*NML_ERRANT              */
  return (nml_fscalar)0;
  }

nml_fscalar		(nml_fv_ssqr)(const nml_fvector* pv) {
  nml_fscalar		t = 0;
  nml_extent		n = nml_fv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    t += NML_FSQR(pv[j]);
  return t;
  }

nml_fscalar		(nml_fv_sabs)(const nml_fvector* pv) {
  return sqrt(nml_fv_ssqr(pv));
  }

nml_fvector*	(nml_fv_vcpy)(nml_fvector* pv,
    const nml_fvector* pw) {
  nml_extent		n = nml_fv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    pv[j] = pw[j];
  return pv;
  }


int
(nml_fv_fprintf)(FILE* stream, const nml_fvector* pv, 
  int width, unsigned int precision, nml_fmtflags flags, unsigned int columns) {
  int		characters;
  int		total = 0;
  const
  nml_extent	n = nml_fv_extent(pv);
  nml_offset	j = 0;
  columns = (0 < columns)? columns: 4;
  for (j = 0; j < n; ++j) {
    if (0 < j) {
      characters
	= (0 == j%columns)? fprintf(stream, "\n"): fprintf(stream, " ");
      if (0 <= characters) { total += characters; } else { return -1; }
      }
    characters = nml_ffprintf(
      stream, pv[j], width, precision, flags);
    if (0 <= characters) { total += characters; } else { return -1; }
    }
  characters = fprintf(stream, "\n");
  if (0 <= characters) { total += characters; } else { return -1; }
  return total;
  }

/* Constructor								*/

#ifdef EXT_MEM_HANDLING
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

nml_fvector*	(nml_fv_clone)(const nml_fvector* pv) {
  nml_extent		n = nml_fv_extent(pv);
  return nml_fv_vcpy(nml_fv_new(n), pv);
  }

/* Destructor								*/
#ifdef EXT_MEM_HANDLING
void			(nml_fv_delete)(nml_fvector* pv) {
   nml_fdarray * pointerToDarray ;
	if (NULL != pv ) {
   pointerToDarray = 
   ( nml_fdarray * ) ( (char*) pv - nml_foffsetLength ) ;
   nml_free( pointerToDarray ) ;
	}
  }
#else  /*   EXT_MEM_HANDLING   */
void			(nml_fv_delete)(nml_fvector* pv) {
  if (NULL != pv)
  nml_free(pv - 1);
  }
#endif  /*   EXT_MEM_HANDLING   */



/* Reconstructor							*/
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

