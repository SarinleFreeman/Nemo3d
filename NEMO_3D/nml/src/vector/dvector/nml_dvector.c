
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
#include<nml_dvector.h>

/* Functions								*/
#ifdef EXT_MEM_HANDLING
nml_extent		(nml_dv_extent)(const nml_dvector* pv) {
   return(  
            ( (nml_ddarray *)
                     (    ((char *)pv) - nml_doffsetLength     )
            ) -> noElements
         ) ;
  }
#define nml_dv_extent(pv) (((nml_ddarray *) \
   (((char *)pv) - nml_doffsetLength)) -> noElements )
#else   /* EXT_MEM_HANDLING */
nml_extent		(nml_dv_extent)(const nml_dvector* pv) {
  return *((const nml_extent*)(pv - 1));
  }
#define nml_dv_extent(pv) (*((const nml_extent*)((pv) - 1)))
#endif   /* EXT_MEM_HANDLING */

nml_dvector*	(nml_dv_fill)(nml_dvector* pv,
    nml_dscalar x) {
  nml_extent		n = nml_dv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    pv[j] = x;
  return pv;
  }

nml_dvector*	(nml_dv_smul)(nml_dvector* pv,
    nml_dscalar x) {
  nml_extent		n = nml_dv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_DMUL(pv[j], pv[j], x);
  return pv;
  }

nml_dvector*	(nml_dv_sdiv)(nml_dvector* pv,
    nml_dscalar x) {
  nml_extent		n = nml_dv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_DDIV(pv[j], pv[j], x);
  return pv;
  }

nml_dvector*	(nml_dv_sadd)(nml_dvector* pv,
    nml_dscalar x) {
  nml_extent		n = nml_dv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_DSUB(pv[j], pv[j], x);
  return pv;
  }

nml_dvector*	(nml_dv_ssub)(nml_dvector* pv,
    nml_dscalar x) {
  nml_extent		n = nml_dv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_DSUB(pv[j], pv[j], x);
  return pv;
  }

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
#ifdef  NML_ERRANT
  nml_message("In function nml_dv_sminmag(const nml_dvector*),\n"
      "zero vector extent.");
#endif/*NML_ERRANT              */
  return (nml_dscalar)0;
  }

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
#ifdef  NML_ERRANT
  nml_message("In function nml_dv_smaxmag(const nml_dvector*),\n"
      "zero vector extent.");
#endif/*NML_ERRANT              */
  return (nml_dscalar)0;
  }

nml_dscalar		(nml_dv_ssqr)(const nml_dvector* pv) {
  nml_dscalar		t = 0;
  nml_extent		n = nml_dv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    t += NML_DSQR(pv[j]);
  return t;
  }

nml_dscalar		(nml_dv_sabs)(const nml_dvector* pv) {
  return sqrt(nml_dv_ssqr(pv));
  }

nml_dvector*	(nml_dv_vcpy)(nml_dvector* pv,
    const nml_dvector* pw) {
  nml_extent		n = nml_dv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    pv[j] = pw[j];
  return pv;
  }


int
(nml_dv_fprintf)(FILE* stream, const nml_dvector* pv, 
  int width, unsigned int precision, nml_fmtflags flags, unsigned int columns) {
  int		characters;
  int		total = 0;
  const
  nml_extent	n = nml_dv_extent(pv);
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

nml_dvector*	(nml_dv_clone)(const nml_dvector* pv) {
  nml_extent		n = nml_dv_extent(pv);
  return nml_dv_vcpy(nml_dv_new(n), pv);
  }

/* Destructor								*/
#ifdef EXT_MEM_HANDLING
void			(nml_dv_delete)(nml_dvector* pv) {
   nml_ddarray * pointerToDarray ;
	if (NULL != pv ) {
   pointerToDarray = 
   ( nml_ddarray * ) ( (char*) pv - nml_doffsetLength ) ;
   nml_free( pointerToDarray ) ;
	}
  }
#else  /*   EXT_MEM_HANDLING   */
void			(nml_dv_delete)(nml_dvector* pv) {
  if (NULL != pv)
  nml_free(pv - 1);
  }
#endif  /*   EXT_MEM_HANDLING   */



/* Reconstructor							*/
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

