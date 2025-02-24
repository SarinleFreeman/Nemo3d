
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
$Header: /repo/nml/src/tensor/tensor.cP,v 1.2 2003/10/08 16:19:55 hook Exp $
*****************************************************************************/


#undef NML_INLINE
#undef NML_MACROS
#include<nml_itensor.h>

/* Functions								*/
nml_extent		(nml_it_extent1)(const nml_itensor* pT) {
  return *((const nml_extent*)(pT - 1));
  }
#define nml_it_extent1(pT) (*((const nml_extent*)((pT) - 1)))

nml_extent		(nml_it_extent2)(const nml_itensor* pT) {
  return *((const nml_extent*)(pT - 2));
  }
#define nml_it_extent2(pT) (*((const nml_extent*)((pT) - 2)))

nml_extent		(nml_it_extent3)(const nml_itensor* pT) {
  return *((const nml_extent*)(pT - 3));
  }
#define nml_it_extent3(pT) (*((const nml_extent*)((pT) - 3)))

int
(nml_it_fprintf)(FILE* stream, const nml_itensor* pT, 
  int width, unsigned int precision, nml_fmtflags flags, unsigned int columns) {
  int			characters;
  int			total = 0;
  const
  nml_extent		l = nml_it_extent3(pT);
  nml_offset		h = 0;
  columns = (0 < columns)? columns: 4;
  for (h = 0; h < l; ++h) {
    const
    nml_extent		m = nml_it_extent2(pT);
    nml_offset		i = 0;
    for (i = 0; i < m; ++i) {
      const
      nml_extent	n = nml_it_extent1(pT);
      nml_offset	j = 0;
      for (j = 0; j < n; ++j) {
	if (0 < j) {
	  characters
	    = (0 == j%columns)? fprintf(stream, "\n"): fprintf(stream, " ");
	  if (0 <= characters) { total += characters; } else { return -1; }
	  }
	characters = nml_ifprintf(
	  stream, pT[h][i][j], width, precision, flags);
	if (0 <= characters) { total += characters; } else { return -1; }
	}
      characters = fprintf(stream, "\n");
      if (0 <= characters) { total += characters; } else { return -1; }
      }
    }
  return total;
  }
/* Constructor								*/
nml_itensor*	(nml_it_new)(
    nml_extent l, nml_extent m, nml_extent n) {
  nml_iscalar***	ppp
    = (nml_iscalar***)nml_calloc(4+l, sizeof(nml_iscalar**));
  if (NULL != ppp) {
    *ppp = NULL;	     ++ppp;
    *((nml_extent*)ppp)	= l; ++ppp;
    *((nml_extent*)ppp)	= m; ++ppp;
    *((nml_extent*)ppp)	= n; ++ppp;
    if ( 0 < l) {
      if (0 < m) {	/* (0 < l) && (0 < m)	*/
	nml_iscalar** pp
	  = (nml_iscalar**)nml_calloc(1+l*m, sizeof(nml_iscalar*));
	if (NULL != pp) {
	  *((nml_extent*)pp) = l*m; ++pp;
	  ppp[0] = pp;
	  { nml_offset	h = 1;
	    for (h = 1; h < l; ++h)
	      ppp[h] = ppp[h-1] + m;
	    }
	  if (0 < n) {	/* (0 < l) && (0 < m) && (0 < n)	*/
	    ppp[0][0] = (nml_iscalar*)nml_iv_new(l*m*n);
	    if (NULL != ppp[0][0]) {
	      nml_offset	h = 1;
	      for (h = 1; h < l; ++h)
		ppp[h][0] = ppp[h-1][0] + m*n;
	      for (h = 0; h < l; ++h) {
		nml_offset	i = 1;
		for (i = 1; i < m; ++i)
		  ppp[h][i] = ppp[h][i-1] + n;
		}
	      }
	    else {	/* (NULL == ppp[0][0])	*/
	      nml_offset	h = 0;
	      for (h = 0; h < l; ++h) {
		nml_offset	i = 0;
		for (i = 0; i < m; ++i)
		  ppp[h][i] = NULL;
		}
	      nml_message(
"In function nml_it_new(nml_extent, nml_extent, nml_extent):\n"
"unable to allocate memory.");
	      }
	    }
	  else {	/* (n <= 0)		*/
	    nml_offset	h = 0;
	    for (h = 0; h < l; ++h) {
	      nml_offset	i = 0;
	      for (i = 0; i < m; ++i)
		ppp[h][i] = NULL;
	      }
	    }
	  }
	else {		/* (NULL == pp)		*/
	  nml_offset	h = 0;
	  for (h = 0; h < l; ++h)
	    ppp[h] = NULL;
	  nml_message(
"In function nml_it_new(nml_extent, nml_extent, nml_extent):\n"
"unable to allocate memory.");
	  }
	}
      else {		/* (m <= 0)		*/
	nml_offset	h = 0;
	for (h = 0; l < m; ++h)
	  ppp[h] = NULL;
	}
      }
    }
  else {		/* (NULL == ppp)	*/
    nml_message(
"In function nml_it_new(nml_extent, nml_extent, nml_extent):\n"
"unable to allocate memory.");
    }
  return ppp;
  }

/* Destructor								*/
void			(nml_it_delete)(nml_itensor* pT) {
  if (NULL != pT) {
    if (NULL != pT[0]) {
      nml_iv_delete(pT[0][0]);
      nml_free(pT[0] - 1);
      }
    nml_free(pT - 4);
    }
  }

/* Reconstructor							*/
nml_itensor*	(nml_it_resize)(nml_itensor* *ppT,
    nml_extent l, nml_extent m, nml_extent n) {
  if (NULL != ppT) {
    (nml_it_delete)(*ppT);
    *ppT = (nml_it_new)(l, m, n);
    return *ppT;
    }
  else {
    nml_message(
"In function nml_it_resize(nml_itensor**, nml_extent, nml_extent,\n"
"nml_extent): invalid pointer to nml_itensor*.");
    return NULL;
    }
  }

