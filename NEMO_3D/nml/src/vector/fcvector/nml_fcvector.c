
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
#include<nml_fcvector.h>

/* Functions								*/
#ifdef EXT_MEM_HANDLING
nml_extent		(nml_fcv_extent)(const nml_fcvector* pv) {
   return(  
            ( (nml_fcdarray *)
                     (    ((char *)pv) - nml_fcoffsetLength     )
            ) -> noElements
         ) ;
  }
#define nml_fcv_extent(pv) (((nml_fcdarray *) \
   (((char *)pv) - nml_fcoffsetLength)) -> noElements )
#else   /* EXT_MEM_HANDLING */
nml_extent		(nml_fcv_extent)(const nml_fcvector* pv) {
  return *((const nml_extent*)(pv - 1));
  }
#define nml_fcv_extent(pv) (*((const nml_extent*)((pv) - 1)))
#endif   /* EXT_MEM_HANDLING */

nml_fcvector*	(nml_fcv_fill)(nml_fcvector* pv,
    nml_fcscalar x) {
  nml_extent		n = nml_fcv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    pv[j] = x;
  return pv;
  }

nml_fcvector*	(nml_fcv_smul)(nml_fcvector* pv,
    nml_fcscalar x) {
  nml_extent		n = nml_fcv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_FCMUL(pv[j], pv[j], x);
  return pv;
  }

nml_fcvector*	(nml_fcv_sdiv)(nml_fcvector* pv,
    nml_fcscalar x) {
  nml_extent		n = nml_fcv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_FCDIV(pv[j], pv[j], x);
  return pv;
  }

nml_fcvector*	(nml_fcv_sadd)(nml_fcvector* pv,
    nml_fcscalar x) {
  nml_extent		n = nml_fcv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_FCSUB(pv[j], pv[j], x);
  return pv;
  }

nml_fcvector*	(nml_fcv_ssub)(nml_fcvector* pv,
    nml_fcscalar x) {
  nml_extent		n = nml_fcv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_FCSUB(pv[j], pv[j], x);
  return pv;
  }

nml_fscalar		(nml_fcv_sminmag)(const nml_fcvector* pv) {
  nml_extent		n = nml_fcv_extent(pv);
  if (0 < n) {
    nml_fscalar	y = NML_FCABS(pv[0]);
    nml_offset		j = 0;
    for (j = 0; j < n; ++j) {
      const
      nml_fscalar	t = NML_FCABS(pv[j]);
      if (t < y) y = t;
      }
    return y;
    }
#ifdef  NML_ERRANT
  nml_message("In function nml_fcv_sminmag(const nml_fcvector*),\n"
      "zero vector extent.");
#endif/*NML_ERRANT              */
  return (nml_fscalar)0;
  }

nml_fscalar		(nml_fcv_smaxmag)(const nml_fcvector* pv) {
  nml_extent		n = nml_fcv_extent(pv);
  if (0 < n) {
    nml_fscalar	y = NML_FCABS(pv[0]);
    nml_offset		j = 0;
    for (j = 1; j < n; ++j) {
      const
      nml_fscalar	t = NML_FCABS(pv[j]);
      if (y < t) y = t;
      }
    return y;
    }
#ifdef  NML_ERRANT
  nml_message("In function nml_fcv_smaxmag(const nml_fcvector*),\n"
      "zero vector extent.");
#endif/*NML_ERRANT              */
  return (nml_fscalar)0;
  }

nml_fscalar		(nml_fcv_ssqr)(const nml_fcvector* pv) {
  nml_fscalar		t = 0;
  nml_extent		n = nml_fcv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    t += NML_FCSQR(pv[j]);
  return t;
  }

nml_fscalar		(nml_fcv_sabs)(const nml_fcvector* pv) {
  return sqrt(nml_fcv_ssqr(pv));
  }

nml_fcvector*	(nml_fcv_vcpy)(nml_fcvector* pv,
    const nml_fcvector* pw) {
  nml_extent		n = nml_fcv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    pv[j] = pw[j];
  return pv;
  }


int
(nml_fcv_fprintf)(FILE* stream, const nml_fcvector* pv, 
  int width, unsigned int precision, nml_fmtflags flags, unsigned int columns) {
  int		characters;
  int		total = 0;
  const
  nml_extent	n = nml_fcv_extent(pv);
  nml_offset	j = 0;
  columns = (0 < columns)? columns: 4;
  for (j = 0; j < n; ++j) {
    if (0 < j) {
      characters
	= (0 == j%columns)? fprintf(stream, "\n"): fprintf(stream, " ");
      if (0 <= characters) { total += characters; } else { return -1; }
      }
    characters = nml_fcfprintf(
      stream, pv[j], width, precision, flags);
    if (0 <= characters) { total += characters; } else { return -1; }
    }
  characters = fprintf(stream, "\n");
  if (0 <= characters) { total += characters; } else { return -1; }
  return total;
  }

/* Constructor								*/

#ifdef EXT_MEM_HANDLING
nml_fcvector*	(nml_fcv_new)(nml_extent n) {
   nml_fcdarray * p ;
   nml_extent * pointerCopy ;

   p = ( nml_fcdarray * ) nml_calloc(nml_fcoffsetLength +
                                       n * sizeof(nml_fcscalar) , 1 ) ;
  if (NULL == p) {
    nml_message("In function nml_fcv_new(nml_extent):\n"
      "unable to allocate memory.");
       }
  else {
      pointerCopy = (nml_extent * ) p ;
      *pointerCopy = n ;
       } 
   return( (nml_fcvector*) &p->start  ) ;
  }
#else  /*   EXT_MEM_HANDLING   */
nml_fcvector*	(nml_fcv_new)(nml_extent n) {
     nml_fcscalar*	p
       = (nml_fcscalar*)nml_calloc(1+n, sizeof(nml_fcscalar));
  if (NULL == p) {
    nml_message("In function nml_fcv_new(nml_extent):\n"
      "unable to allocate memory.");
    }
  else {
       *((nml_extent*)p)	= n; ++p;
    }
  return p;
  }
#endif  /*   EXT_MEM_HANDLING   */

nml_fcvector*	(nml_fcv_clone)(const nml_fcvector* pv) {
  nml_extent		n = nml_fcv_extent(pv);
  return nml_fcv_vcpy(nml_fcv_new(n), pv);
  }

/* Destructor								*/
#ifdef EXT_MEM_HANDLING
void			(nml_fcv_delete)(nml_fcvector* pv) {
   nml_fcdarray * pointerToDarray ;
	if (NULL != pv ) {
   pointerToDarray = 
   ( nml_fcdarray * ) ( (char*) pv - nml_fcoffsetLength ) ;
   nml_free( pointerToDarray ) ;
	}
  }
#else  /*   EXT_MEM_HANDLING   */
void			(nml_fcv_delete)(nml_fcvector* pv) {
  if (NULL != pv)
  nml_free(pv - 1);
  }
#endif  /*   EXT_MEM_HANDLING   */



/* Reconstructor							*/
nml_fcvector*	(nml_fcv_resize)(nml_fcvector* *ppv, nml_extent n) {
  if (NULL != ppv) {
    (nml_fcv_delete)(*ppv);
    *ppv = (nml_fcv_new)(n);
    return *ppv;
    }
  else {
    nml_message(
"In function nml_fcv_resize(nml_fcvector**, nml_extent):\n"
"invalid pointer to nml_fcvector*.");
    return NULL;
    }
  }

