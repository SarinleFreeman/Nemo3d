
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
$Header: /repo/nml/src/matrix/matrix.cP,v 1.2 2003/10/08 16:19:34 hook Exp $
*****************************************************************************/

#undef NML_INLINE
#undef NML_MACROS
#include<nml_smatrix.h>

/* Functions								*/
nml_extent		(nml_sm_extent1)(const nml_smatrix* pM) {
  return *((const nml_extent*)(pM - 1));
  }
#define nml_sm_extent1(pM) (*((const nml_extent*)((pM) - 1)))

nml_extent		(nml_sm_extent2)(const nml_smatrix* pM) {
  return *((const nml_extent*)(pM - 2));
  }
#define nml_sm_extent2(pM) (*((const nml_extent*)((pM) - 2)))

int
(nml_sm_fprintf)(FILE* stream, const nml_smatrix* pM, 
  int width, unsigned int precision, nml_fmtflags flags, unsigned int columns) {
  int		characters;
  int		total = 0;
  const
  nml_extent	m = nml_sm_extent2(pM);
  nml_offset	i = 0;
  columns = (0 < columns)? columns: 4;
  for (i = 0; i < m; ++i) {
    const
    nml_extent	n = nml_sm_extent1(pM);
    nml_offset	j = 0;
    for (j = 0; j < n; ++j) {
      if (0 < j) {
	characters
	  = (0 == j%columns)? fprintf(stream, "\n"): fprintf(stream, " ");
	if (0 <= characters) { total += characters; } else { return -1; }
	}
      characters = nml_sfprintf(
	stream, pM[i][j], width, precision, flags);
      if (0 <= characters) { total += characters; } else { return -1; }
      }
    characters = fprintf(stream, "\n");
    if (0 <= characters) { total += characters; } else { return -1; }
    }
  return total;
  }

/* Constructor								*/
nml_smatrix*	(nml_sm_new)(nml_extent m, nml_extent n) {
  nml_sscalar**	pp
    = (nml_sscalar**)nml_calloc(3+m, sizeof(nml_sscalar*));
  if (NULL == pp) {
    nml_message("In function nml_sm_new(nml_extent, nml_extent):\n"
      "unable to allocate memory.");
    }
  else {		/* (NULL != pp)		*/
    *pp = NULL;		     ++pp;
    *((nml_extent*)pp)	= m; ++pp;
    *((nml_extent*)pp)	= n; ++pp;
    if (0 < m) {
      if (0 < n) {
	pp[0] = (nml_sscalar*)nml_sv_new(n*m);
	if (NULL == pp[0]) {
	  nml_offset	i = 1;
          for (i = 1; i < m; ++i)
	    pp[i] = NULL;
	  nml_message("In function nml_sm_new(nml_extent, nml_extent):\n"
	    "unable to allocate memory.");
	  }
	else {		/* (NULL != pp[0])	*/
	  nml_offset	i = 1;
	  for (i = 1; i < m; ++i)
	    pp[i] = pp[i-1] + n;
	  }
	}
      else {		/* (n <= 0)		*/
	nml_offset	i = 0;
        for (i = 0; i < m; ++i)
	  pp[i] = NULL;
	}
      }
    }
  return pp;
  }

/* Destructor								*/
void			(nml_sm_delete)(nml_smatrix* pM) {
  if (NULL != pM) {
    nml_sv_delete(pM[0]);
    nml_free(pM - 3);
    }
  }

/* Reconstructor							*/
nml_smatrix*	(nml_sm_resize)(nml_smatrix* *ppM,
    nml_extent m, nml_extent n) {
  if (NULL != ppM) {
    (nml_sm_delete)(*ppM);
    *ppM = (nml_sm_new)(m, n);
    return *ppM;
    }
  else {
    nml_message(
"In function nml_sm_resize(nml_smatrix**, nml_extent, nml_extent):\n"
"invalid pointer to nml_smatrix*.");
    return NULL;
    }
  }

