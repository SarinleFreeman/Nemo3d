
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
#include<nml_svector.h>

/* Functions								*/
#ifdef EXT_MEM_HANDLING
nml_extent		(nml_sv_extent)(const nml_svector* pv) {
   return(  
            ( (nml_sdarray *)
                     (    ((char *)pv) - nml_soffsetLength     )
            ) -> noElements
         ) ;
  }
#define nml_sv_extent(pv) (((nml_sdarray *) \
   (((char *)pv) - nml_soffsetLength)) -> noElements )
#else   /* EXT_MEM_HANDLING */
nml_extent		(nml_sv_extent)(const nml_svector* pv) {
  return *((const nml_extent*)(pv - 1));
  }
#define nml_sv_extent(pv) (*((const nml_extent*)((pv) - 1)))
#endif   /* EXT_MEM_HANDLING */

nml_svector*	(nml_sv_fill)(nml_svector* pv,
    nml_sscalar x) {
  nml_extent		n = nml_sv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    pv[j] = x;
  return pv;
  }

nml_svector*	(nml_sv_smul)(nml_svector* pv,
    nml_sscalar x) {
  nml_extent		n = nml_sv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_SMUL(pv[j], pv[j], x);
  return pv;
  }

nml_svector*	(nml_sv_sdiv)(nml_svector* pv,
    nml_sscalar x) {
  nml_extent		n = nml_sv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_SDIV(pv[j], pv[j], x);
  return pv;
  }

nml_svector*	(nml_sv_sadd)(nml_svector* pv,
    nml_sscalar x) {
  nml_extent		n = nml_sv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_SSUB(pv[j], pv[j], x);
  return pv;
  }

nml_svector*	(nml_sv_ssub)(nml_svector* pv,
    nml_sscalar x) {
  nml_extent		n = nml_sv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_SSUB(pv[j], pv[j], x);
  return pv;
  }

nml_sscalar		(nml_sv_sminmag)(const nml_svector* pv) {
  nml_extent		n = nml_sv_extent(pv);
  if (0 < n) {
    nml_sscalar	y = NML_SABS(pv[0]);
    nml_offset		j = 0;
    for (j = 0; j < n; ++j) {
      const
      nml_sscalar	t = NML_SABS(pv[j]);
      if (t < y) y = t;
      }
    return y;
    }
#ifdef  NML_ERRANT
  nml_message("In function nml_sv_sminmag(const nml_svector*),\n"
      "zero vector extent.");
#endif/*NML_ERRANT              */
  return (nml_sscalar)0;
  }

nml_sscalar		(nml_sv_smaxmag)(const nml_svector* pv) {
  nml_extent		n = nml_sv_extent(pv);
  if (0 < n) {
    nml_sscalar	y = NML_SABS(pv[0]);
    nml_offset		j = 0;
    for (j = 1; j < n; ++j) {
      const
      nml_sscalar	t = NML_SABS(pv[j]);
      if (y < t) y = t;
      }
    return y;
    }
#ifdef  NML_ERRANT
  nml_message("In function nml_sv_smaxmag(const nml_svector*),\n"
      "zero vector extent.");
#endif/*NML_ERRANT              */
  return (nml_sscalar)0;
  }

nml_sscalar		(nml_sv_ssqr)(const nml_svector* pv) {
  nml_sscalar		t = 0;
  nml_extent		n = nml_sv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    t += NML_SSQR(pv[j]);
  return t;
  }

nml_svector*	(nml_sv_vcpy)(nml_svector* pv,
    const nml_svector* pw) {
  nml_extent		n = nml_sv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    pv[j] = pw[j];
  return pv;
  }


int
(nml_sv_fprintf)(FILE* stream, const nml_svector* pv, 
  int width, unsigned int precision, nml_fmtflags flags, unsigned int columns) {
  int		characters;
  int		total = 0;
  const
  nml_extent	n = nml_sv_extent(pv);
  nml_offset	j = 0;
  columns = (0 < columns)? columns: 4;
  for (j = 0; j < n; ++j) {
    if (0 < j) {
      characters
	= (0 == j%columns)? fprintf(stream, "\n"): fprintf(stream, " ");
      if (0 <= characters) { total += characters; } else { return -1; }
      }
    characters = nml_sfprintf(
      stream, pv[j], width, precision, flags);
    if (0 <= characters) { total += characters; } else { return -1; }
    }
  characters = fprintf(stream, "\n");
  if (0 <= characters) { total += characters; } else { return -1; }
  return total;
  }

/* Constructor								*/

#ifdef EXT_MEM_HANDLING
nml_svector*	(nml_sv_new)(nml_extent n) {
   nml_sdarray * p ;
   nml_extent * pointerCopy ;

   p = ( nml_sdarray * ) nml_calloc(nml_soffsetLength +
                                       n * sizeof(nml_sscalar) , 1 ) ;
  if (NULL == p) {
    nml_message("In function nml_sv_new(nml_extent):\n"
      "unable to allocate memory.");
       }
  else {
      pointerCopy = (nml_extent * ) p ;
      *pointerCopy = n ;
       } 
   return( (nml_svector*) &p->start  ) ;
  }
#else  /*   EXT_MEM_HANDLING   */
nml_svector*	(nml_sv_new)(nml_extent n) {
     nml_sscalar*	p
       = (nml_sscalar*)nml_calloc(1+n, sizeof(nml_sscalar));
  if (NULL == p) {
    nml_message("In function nml_sv_new(nml_extent):\n"
      "unable to allocate memory.");
    }
  else {
       *((nml_extent*)p)	= n; ++p;
    }
  return p;
  }
#endif  /*   EXT_MEM_HANDLING   */

nml_svector*	(nml_sv_clone)(const nml_svector* pv) {
  nml_extent		n = nml_sv_extent(pv);
  return nml_sv_vcpy(nml_sv_new(n), pv);
  }

/* Destructor								*/
#ifdef EXT_MEM_HANDLING
void			(nml_sv_delete)(nml_svector* pv) {
   nml_sdarray * pointerToDarray ;
	if (NULL != pv ) {
   pointerToDarray = 
   ( nml_sdarray * ) ( (char*) pv - nml_soffsetLength ) ;
   nml_free( pointerToDarray ) ;
	}
  }
#else  /*   EXT_MEM_HANDLING   */
void			(nml_sv_delete)(nml_svector* pv) {
  if (NULL != pv)
  nml_free(pv - 1);
  }
#endif  /*   EXT_MEM_HANDLING   */



/* Reconstructor							*/
nml_svector*	(nml_sv_resize)(nml_svector* *ppv, nml_extent n) {
  if (NULL != ppv) {
    (nml_sv_delete)(*ppv);
    *ppv = (nml_sv_new)(n);
    return *ppv;
    }
  else {
    nml_message(
"In function nml_sv_resize(nml_svector**, nml_extent):\n"
"invalid pointer to nml_svector*.");
    return NULL;
    }
  }

