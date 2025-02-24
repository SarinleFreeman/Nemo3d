
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
#include<nml_ivector.h>

/* Functions								*/
#ifdef EXT_MEM_HANDLING
nml_extent		(nml_iv_extent)(const nml_ivector* pv) {
   return(  
            ( (nml_idarray *)
                     (    ((char *)pv) - nml_ioffsetLength     )
            ) -> noElements
         ) ;
  }
#define nml_iv_extent(pv) (((nml_idarray *) \
   (((char *)pv) - nml_ioffsetLength)) -> noElements )
#else   /* EXT_MEM_HANDLING */
nml_extent		(nml_iv_extent)(const nml_ivector* pv) {
  return *((const nml_extent*)(pv - 1));
  }
#define nml_iv_extent(pv) (*((const nml_extent*)((pv) - 1)))
#endif   /* EXT_MEM_HANDLING */

nml_ivector*	(nml_iv_fill)(nml_ivector* pv,
    nml_iscalar x) {
  nml_extent		n = nml_iv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    pv[j] = x;
  return pv;
  }

nml_ivector*	(nml_iv_smul)(nml_ivector* pv,
    nml_iscalar x) {
  nml_extent		n = nml_iv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_iMUL(pv[j], pv[j], x);
  return pv;
  }

nml_ivector*	(nml_iv_sdiv)(nml_ivector* pv,
    nml_iscalar x) {
  nml_extent		n = nml_iv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_iDIV(pv[j], pv[j], x);
  return pv;
  }

nml_ivector*	(nml_iv_sadd)(nml_ivector* pv,
    nml_iscalar x) {
  nml_extent		n = nml_iv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_iSUB(pv[j], pv[j], x);
  return pv;
  }

nml_ivector*	(nml_iv_ssub)(nml_ivector* pv,
    nml_iscalar x) {
  nml_extent		n = nml_iv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    NML_iSUB(pv[j], pv[j], x);
  return pv;
  }

nml_iscalar		(nml_iv_sminmag)(const nml_ivector* pv) {
  nml_extent		n = nml_iv_extent(pv);
  if (0 < n) {
    nml_iscalar	y = NML_iABS(pv[0]);
    nml_offset		j = 0;
    for (j = 0; j < n; ++j) {
      const
      nml_iscalar	t = NML_iABS(pv[j]);
      if (t < y) y = t;
      }
    return y;
    }
#ifdef  NML_ERRANT
  nml_message("In function nml_iv_sminmag(const nml_ivector*),\n"
      "zero vector extent.");
#endif/*NML_ERRANT              */
  return (nml_iscalar)0;
  }

nml_iscalar		(nml_iv_smaxmag)(const nml_ivector* pv) {
  nml_extent		n = nml_iv_extent(pv);
  if (0 < n) {
    nml_iscalar	y = NML_iABS(pv[0]);
    nml_offset		j = 0;
    for (j = 1; j < n; ++j) {
      const
      nml_iscalar	t = NML_iABS(pv[j]);
      if (y < t) y = t;
      }
    return y;
    }
#ifdef  NML_ERRANT
  nml_message("In function nml_iv_smaxmag(const nml_ivector*),\n"
      "zero vector extent.");
#endif/*NML_ERRANT              */
  return (nml_iscalar)0;
  }

nml_iscalar		(nml_iv_ssqr)(const nml_ivector* pv) {
  nml_iscalar		t = 0;
  nml_extent		n = nml_iv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    t += NML_iSQR(pv[j]);
  return t;
  }

nml_ivector*	(nml_iv_vcpy)(nml_ivector* pv,
    const nml_ivector* pw) {
  nml_extent		n = nml_iv_extent(pv);
  nml_offset		j = 0;
  for (j = 0; j < n; ++j)
    pv[j] = pw[j];
  return pv;
  }


int
(nml_iv_fprintf)(FILE* stream, const nml_ivector* pv, 
  int width, unsigned int precision, nml_fmtflags flags, unsigned int columns) {
  int		characters;
  int		total = 0;
  const
  nml_extent	n = nml_iv_extent(pv);
  nml_offset	j = 0;
  columns = (0 < columns)? columns: 4;
  for (j = 0; j < n; ++j) {
    if (0 < j) {
      characters
	= (0 == j%columns)? fprintf(stream, "\n"): fprintf(stream, " ");
      if (0 <= characters) { total += characters; } else { return -1; }
      }
    characters = nml_ifprintf(
      stream, pv[j], width, precision, flags);
    if (0 <= characters) { total += characters; } else { return -1; }
    }
  characters = fprintf(stream, "\n");
  if (0 <= characters) { total += characters; } else { return -1; }
  return total;
  }

/* Constructor								*/

#ifdef EXT_MEM_HANDLING
nml_ivector*	(nml_iv_new)(nml_extent n) {
   nml_idarray * p ;
   nml_extent * pointerCopy ;

   p = ( nml_idarray * ) nml_calloc(nml_ioffsetLength +
                                       n * sizeof(nml_iscalar) , 1 ) ;
  if (NULL == p) {
    nml_message("In function nml_iv_new(nml_extent):\n"
      "unable to allocate memory.");
       }
  else {
      pointerCopy = (nml_extent * ) p ;
      *pointerCopy = n ;
       } 
   return( (nml_ivector*) &p->start  ) ;
  }
#else  /*   EXT_MEM_HANDLING   */
nml_ivector*	(nml_iv_new)(nml_extent n) {
     nml_iscalar*	p
       = (nml_iscalar*)nml_calloc(1+n, sizeof(nml_iscalar));
  if (NULL == p) {
    nml_message("In function nml_iv_new(nml_extent):\n"
      "unable to allocate memory.");
    }
  else {
       *((nml_extent*)p)	= n; ++p;
    }
  return p;
  }
#endif  /*   EXT_MEM_HANDLING   */

nml_ivector*	(nml_iv_clone)(const nml_ivector* pv) {
  nml_extent		n = nml_iv_extent(pv);
  return nml_iv_vcpy(nml_iv_new(n), pv);
  }

/* Destructor								*/
#ifdef EXT_MEM_HANDLING
void			(nml_iv_delete)(nml_ivector* pv) {
   nml_idarray * pointerToDarray ;
	if (NULL != pv ) {
   pointerToDarray = 
   ( nml_idarray * ) ( (char*) pv - nml_ioffsetLength ) ;
   nml_free( pointerToDarray ) ;
	}
  }
#else  /*   EXT_MEM_HANDLING   */
void			(nml_iv_delete)(nml_ivector* pv) {
  if (NULL != pv)
  nml_free(pv - 1);
  }
#endif  /*   EXT_MEM_HANDLING   */



/* Reconstructor							*/
nml_ivector*	(nml_iv_resize)(nml_ivector* *ppv, nml_extent n) {
  if (NULL != ppv) {
    (nml_iv_delete)(*ppv);
    *ppv = (nml_iv_new)(n);
    return *ppv;
    }
  else {
    nml_message(
"In function nml_iv_resize(nml_ivector**, nml_extent):\n"
"invalid pointer to nml_ivector*.");
    return NULL;
    }
  }

