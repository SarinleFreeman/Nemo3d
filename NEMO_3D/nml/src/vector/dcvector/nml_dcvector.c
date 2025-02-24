
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
#include<nml_dcvector.h>

/* Functions								*/
#ifdef EXT_MEM_HANDLING
nml_extent		(nml_dcv_extent)(const nml_dcvector* pv) {
   return(  
            ( (nml_dcdarray *)
                     (    ((char *)pv) - nml_dcoffsetLength     )
            ) -> noElements
         ) ;
  }
#define nml_dcv_extent(pv) (((nml_dcdarray *) \
   (((char *)pv) - nml_dcoffsetLength)) -> noElements )
#else   /* EXT_MEM_HANDLING */
nml_extent		(nml_dcv_extent)(const nml_dcvector* pv) {
  return *((const nml_extent*)(pv - 1));
  }
#define nml_dcv_extent(pv) (*((const nml_extent*)((pv) - 1)))
#endif   /* EXT_MEM_HANDLING */

nml_dcvector*	(nml_dcv_fill)(nml_dcvector* pv,
    nml_dcscalar x) {
  nml_extent		n = nml_dcv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    pv[j] = x;
  return pv;
  }

nml_dcvector*	(nml_dcv_smul)(nml_dcvector* pv,
    nml_dcscalar x) {
  nml_extent		n = nml_dcv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_DCMUL(pv[j], pv[j], x);
  return pv;
  }

nml_dcvector*	(nml_dcv_sdiv)(nml_dcvector* pv,
    nml_dcscalar x) {
  nml_extent		n = nml_dcv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_DCDIV(pv[j], pv[j], x);
  return pv;
  }

nml_dcvector*	(nml_dcv_sadd)(nml_dcvector* pv,
    nml_dcscalar x) {
  nml_extent		n = nml_dcv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_DCSUB(pv[j], pv[j], x);
  return pv;
  }

nml_dcvector*	(nml_dcv_ssub)(nml_dcvector* pv,
    nml_dcscalar x) {
  nml_extent		n = nml_dcv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_DCSUB(pv[j], pv[j], x);
  return pv;
  }

nml_dscalar		(nml_dcv_sminmag)(const nml_dcvector* pv) {
  nml_extent		n = nml_dcv_extent(pv);
  if (0 < n) {
    nml_dscalar	y = NML_DCABS(pv[0]);
    nml_offset		j = 0;
    for (j = 0; j < n; ++j) {
      const
      nml_dscalar	t = NML_DCABS(pv[j]);
      if (t < y) y = t;
      }
    return y;
    }
#ifdef  NML_ERRANT
  nml_message("In function nml_dcv_sminmag(const nml_dcvector*),\n"
      "zero vector extent.");
#endif/*NML_ERRANT              */
  return (nml_dscalar)0;
  }

nml_dscalar		(nml_dcv_smaxmag)(const nml_dcvector* pv) {
  nml_extent		n = nml_dcv_extent(pv);
  if (0 < n) {
    nml_dscalar	y = NML_DCABS(pv[0]);
    nml_offset		j = 0;
    for (j = 1; j < n; ++j) {
      const
      nml_dscalar	t = NML_DCABS(pv[j]);
      if (y < t) y = t;
      }
    return y;
    }
#ifdef  NML_ERRANT
  nml_message("In function nml_dcv_smaxmag(const nml_dcvector*),\n"
      "zero vector extent.");
#endif/*NML_ERRANT              */
  return (nml_dscalar)0;
  }

nml_dscalar		(nml_dcv_ssqr)(const nml_dcvector* pv) {
  nml_dscalar		t = 0;
  nml_extent		n = nml_dcv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    t += NML_DCSQR(pv[j]);
  return t;
  }

nml_dscalar		(nml_dcv_sabs)(const nml_dcvector* pv) {
  return sqrt(nml_dcv_ssqr(pv));
  }

nml_dcvector*	(nml_dcv_vcpy)(nml_dcvector* pv,
    const nml_dcvector* pw) {
  nml_extent		n = nml_dcv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    pv[j] = pw[j];
  return pv;
  }


int
(nml_dcv_fprintf)(FILE* stream, const nml_dcvector* pv, 
  int width, unsigned int precision, nml_fmtflags flags, unsigned int columns) {
  int		characters;
  int		total = 0;
  const
  nml_extent	n = nml_dcv_extent(pv);
  nml_offset	j = 0;
  columns = (0 < columns)? columns: 4;
  for (j = 0; j < n; ++j) {
    if (0 < j) {
      characters
	= (0 == j%columns)? fprintf(stream, "\n"): fprintf(stream, " ");
      if (0 <= characters) { total += characters; } else { return -1; }
      }
    characters = nml_dcfprintf(
      stream, pv[j], width, precision, flags);
    if (0 <= characters) { total += characters; } else { return -1; }
    }
  characters = fprintf(stream, "\n");
  if (0 <= characters) { total += characters; } else { return -1; }
  return total;
  }

/* Constructor								*/

#ifdef EXT_MEM_HANDLING
nml_dcvector*	(nml_dcv_new)(nml_extent n) {
   nml_dcdarray * p ;
   nml_extent * pointerCopy ;

   p = ( nml_dcdarray * ) nml_calloc(nml_dcoffsetLength +
                                       n * sizeof(nml_dcscalar) , 1 ) ;
  if (NULL == p) {
    nml_message("In function nml_dcv_new(nml_extent):\n"
      "unable to allocate memory.");
       }
  else {
      pointerCopy = (nml_extent * ) p ;
      *pointerCopy = n ;
       } 
   return( (nml_dcvector*) &p->start  ) ;
  }
#else  /*   EXT_MEM_HANDLING   */
nml_dcvector*	(nml_dcv_new)(nml_extent n) {
     nml_dcscalar*	p
       = (nml_dcscalar*)nml_calloc(1+n, sizeof(nml_dcscalar));
  if (NULL == p) {
    nml_message("In function nml_dcv_new(nml_extent):\n"
      "unable to allocate memory.");
    }
  else {
       *((nml_extent*)p)	= n; ++p;
    }
  return p;
  }
#endif  /*   EXT_MEM_HANDLING   */

nml_dcvector*	(nml_dcv_clone)(const nml_dcvector* pv) {
  nml_extent		n = nml_dcv_extent(pv);
  return nml_dcv_vcpy(nml_dcv_new(n), pv);
  }

/* Destructor								*/
#ifdef EXT_MEM_HANDLING
void			(nml_dcv_delete)(nml_dcvector* pv) {
   nml_dcdarray * pointerToDarray ;
	if (NULL != pv ) {
   pointerToDarray = 
   ( nml_dcdarray * ) ( (char*) pv - nml_dcoffsetLength ) ;
   nml_free( pointerToDarray ) ;
	}
  }
#else  /*   EXT_MEM_HANDLING   */
void			(nml_dcv_delete)(nml_dcvector* pv) {
  if (NULL != pv)
  nml_free(pv - 1);
  }
#endif  /*   EXT_MEM_HANDLING   */



/* Reconstructor							*/
nml_dcvector*	(nml_dcv_resize)(nml_dcvector* *ppv, nml_extent n) {
  if (NULL != ppv) {
    (nml_dcv_delete)(*ppv);
    *ppv = (nml_dcv_new)(n);
    return *ppv;
    }
  else {
    nml_message(
"In function nml_dcv_resize(nml_dcvector**, nml_extent):\n"
"invalid pointer to nml_dcvector*.");
    return NULL;
    }
  }

