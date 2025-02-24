
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
#include<nml_bvector.h>

/* Functions								*/
#ifdef EXT_MEM_HANDLING
nml_extent		(nml_bv_extent)(const nml_bvector* pv) {
   return(  
            ( (nml_bdarray *)
                     (    ((char *)pv) - nml_boffsetLength     )
            ) -> noElements
         ) ;
  }
#define nml_bv_extent(pv) (((nml_bdarray *) \
   (((char *)pv) - nml_boffsetLength)) -> noElements )
#else   /* EXT_MEM_HANDLING */
nml_extent		(nml_bv_extent)(const nml_bvector* pv) {
  return *((const nml_extent*)(pv - 1));
  }
#define nml_bv_extent(pv) (*((const nml_extent*)((pv) - 1)))
#endif   /* EXT_MEM_HANDLING */

nml_bvector*	(nml_bv_fill)(nml_bvector* pv,
    nml_bscalar x) {
  nml_extent		n = nml_bv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    pv[j] = x;
  return pv;
  }

nml_bvector*	(nml_bv_smul)(nml_bvector* pv,
    nml_bscalar x) {
  nml_extent		n = nml_bv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_bMUL(pv[j], pv[j], x);
  return pv;
  }

nml_bvector*	(nml_bv_sdiv)(nml_bvector* pv,
    nml_bscalar x) {
  nml_extent		n = nml_bv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_bDIV(pv[j], pv[j], x);
  return pv;
  }

nml_bvector*	(nml_bv_sadd)(nml_bvector* pv,
    nml_bscalar x) {
  nml_extent		n = nml_bv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_bSUB(pv[j], pv[j], x);
  return pv;
  }

nml_bvector*	(nml_bv_ssub)(nml_bvector* pv,
    nml_bscalar x) {
  nml_extent		n = nml_bv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_bSUB(pv[j], pv[j], x);
  return pv;
  }

nml_bscalar		(nml_bv_sminmag)(const nml_bvector* pv) {
  nml_extent		n = nml_bv_extent(pv);
  if (0 < n) {
    nml_bscalar	y = NML_bABS(pv[0]);
    nml_offset		j = 0;
    for (j = 0; j < n; ++j) {
      const
      nml_bscalar	t = NML_bABS(pv[j]);
      if (t < y) y = t;
      }
    return y;
    }
#ifdef  NML_ERRANT
  nml_message("In function nml_bv_sminmag(const nml_bvector*),\n"
      "zero vector extent.");
#endif/*NML_ERRANT              */
  return (nml_bscalar)0;
  }

nml_bscalar		(nml_bv_smaxmag)(const nml_bvector* pv) {
  nml_extent		n = nml_bv_extent(pv);
  if (0 < n) {
    nml_bscalar	y = NML_bABS(pv[0]);
    nml_offset		j = 0;
    for (j = 1; j < n; ++j) {
      const
      nml_bscalar	t = NML_bABS(pv[j]);
      if (y < t) y = t;
      }
    return y;
    }
#ifdef  NML_ERRANT
  nml_message("In function nml_bv_smaxmag(const nml_bvector*),\n"
      "zero vector extent.");
#endif/*NML_ERRANT              */
  return (nml_bscalar)0;
  }

nml_bscalar		(nml_bv_ssqr)(const nml_bvector* pv) {
  nml_bscalar		t = 0;
  nml_extent		n = nml_bv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    t += NML_bSQR(pv[j]);
  return t;
  }

nml_bvector*	(nml_bv_vcpy)(nml_bvector* pv,
    const nml_bvector* pw) {
  nml_extent		n = nml_bv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    pv[j] = pw[j];
  return pv;
  }


int
(nml_bv_fprintf)(FILE* stream, const nml_bvector* pv, 
  int width, unsigned int precision, nml_fmtflags flags, unsigned int columns) {
  int		characters;
  int		total = 0;
  const
  nml_extent	n = nml_bv_extent(pv);
  nml_offset	j = 0;
  columns = (0 < columns)? columns: 4;
  for (j = 0; j < n; ++j) {
    if (0 < j) {
      characters
	= (0 == j%columns)? fprintf(stream, "\n"): fprintf(stream, " ");
      if (0 <= characters) { total += characters; } else { return -1; }
      }
    characters = nml_bfprintf(
      stream, pv[j], width, precision, flags);
    if (0 <= characters) { total += characters; } else { return -1; }
    }
  characters = fprintf(stream, "\n");
  if (0 <= characters) { total += characters; } else { return -1; }
  return total;
  }

/* Constructor								*/

#ifdef EXT_MEM_HANDLING
nml_bvector*	(nml_bv_new)(nml_extent n) {
   nml_bdarray * p ;
   nml_extent * pointerCopy ;

   p = ( nml_bdarray * ) nml_calloc(nml_boffsetLength +
                                       n * sizeof(nml_bscalar) , 1 ) ;
  if (NULL == p) {
    nml_message("In function nml_bv_new(nml_extent):\n"
      "unable to allocate memory.");
       }
  else {
      pointerCopy = (nml_extent * ) p ;
      *pointerCopy = n ;
       } 
   return( (nml_bvector*) &p->start  ) ;
  }
#else  /*   EXT_MEM_HANDLING   */
nml_bvector*	(nml_bv_new)(nml_extent n) {
     nml_bscalar*	p
       = (nml_bscalar*)nml_calloc(1+n, sizeof(nml_bscalar));
  if (NULL == p) {
    nml_message("In function nml_bv_new(nml_extent):\n"
      "unable to allocate memory.");
    }
  else {
       *((nml_extent*)p)	= n; ++p;
    }
  return p;
  }
#endif  /*   EXT_MEM_HANDLING   */

nml_bvector*	(nml_bv_clone)(const nml_bvector* pv) {
  nml_extent		n = nml_bv_extent(pv);
  return nml_bv_vcpy(nml_bv_new(n), pv);
  }

/* Destructor								*/
#ifdef EXT_MEM_HANDLING
void			(nml_bv_delete)(nml_bvector* pv) {
   nml_bdarray * pointerToDarray ;
	if (NULL != pv ) {
   pointerToDarray = 
   ( nml_bdarray * ) ( (char*) pv - nml_boffsetLength ) ;
   nml_free( pointerToDarray ) ;
	}
  }
#else  /*   EXT_MEM_HANDLING   */
void			(nml_bv_delete)(nml_bvector* pv) {
  if (NULL != pv)
  nml_free(pv - 1);
  }
#endif  /*   EXT_MEM_HANDLING   */



/* Reconstructor							*/
nml_bvector*	(nml_bv_resize)(nml_bvector* *ppv, nml_extent n) {
  if (NULL != ppv) {
    (nml_bv_delete)(*ppv);
    *ppv = (nml_bv_new)(n);
    return *ppv;
    }
  else {
    nml_message(
"In function nml_bv_resize(nml_bvector**, nml_extent):\n"
"invalid pointer to nml_bvector*.");
    return NULL;
    }
  }

