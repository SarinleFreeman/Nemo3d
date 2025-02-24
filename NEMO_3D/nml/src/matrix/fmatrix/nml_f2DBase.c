
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
$Header: /repo/nml/src/matrix/2DBase.cP,v 1.2 2003/10/08 16:19:30 hook Exp $
*****************************************************************************/

#undef NML_INLINE
#undef NML_MACROS
#include<nml_f2DBase.h>

/* Functions								*/
nml_index		(nml_f2D_lower2)(const nml_f2DBase* pM) {
  return *((const nml_index*)(pM - 1));
  }
#define nml_f2D_lower2(pM) (*((const nml_index*)((pM) - 1)))

nml_index		(nml_f2D_upper2)(const nml_f2DBase* pM) {
  return *((const nml_index*)(pM - 2));
  }
#define nml_f2D_upper2(pM) (*((const nml_index*)((pM) - 2)))

nml_extent		(nml_f2D_extent2)(const nml_f2DBase* pM) {
  return 1+nml_f2D_upper2(pM)-nml_f2D_lower2(pM);
  }
#define nml_f2D_extent2(pM) (1+nml_f2D_upper2(pM)-nml_f2D_lower2(pM))

nml_index		(nml_f2D_lower1)(const nml_f2DBase* pM) {
  return *((const nml_index*)(pM - 3));
  }
#define nml_f2D_lower1(pM) (*((const nml_index*)((pM) - 3)))

nml_index		(nml_f2D_upper1)(const nml_f2DBase* pM) {
  return *((const nml_index*)(pM - 4));
  }
#define nml_f2D_upper1(pM) (*((const nml_index*)((pM) - 4)))

nml_extent		(nml_f2D_extent1)(const nml_f2DBase* pM) {
  return 1+nml_f2D_upper1(pM)-nml_f2D_lower1(pM);
  }
#define nml_f2D_extent1(pM) (1+nml_f2D_upper1(pM)-nml_f2D_lower1(pM))

int
(nml_f2D_fprintf)(FILE* stream, const nml_f2DBase* pM, 
  int width, unsigned int precision, nml_fmtflags flags, unsigned int columns) {
  int		characters;
  int		total = 0;
  const
  nml_extent	m = nml_f2D_extent2(pM);
  nml_offset	i = 0;
  columns = (0 < columns)? columns: 4;
  for (i = 0; i < m; ++i) {
    const
    nml_extent	n = nml_f2D_extent1(pM);
    nml_offset	j = 0;
    for (j = 0; j < n; ++j) {
      nml_index k = (nml_index)j + nml_f2D_lower1(pM);
      if (0 < j) {
	characters
	  = (0 == j%columns)? fprintf(stream, "\n"): fprintf(stream, " ");
	if (0 <= characters) { total += characters; } else { return -1; }
	}
      characters = nml_ffprintf(
	stream, pM[i][k], width, precision, flags);
      if (0 <= characters) { total += characters; } else { return -1; }
      }
    characters = fprintf(stream, "\n");
    if (0 <= characters) { total += characters; } else { return -1; }
    }
  return total;
  }

/* Constructor								*/
nml_f2DBase*	(nml_f2D_new)(
    nml_index lb2, nml_index ub2, nml_index lb1, nml_index ub1) {
  nml_fscalar**	pp = NULL;
  nml_index		m = 1 + ub2 - lb2;
  nml_index		n = 1 + ub1 - lb1;
  if (m < 0) {
    nml_message(
"In function nml_f2D_new(nml_index, nml_index, nml_index, nml_index):\n"
"extent2 is negative.");
    }
  else
  if (n < 0) {
    nml_message(
"In function nml_f2D_new(nml_index, nml_index, nml_index, nml_index):\n"
"extent1 is negative.");
    }
  else {		/* (0 <= m) && (0 <= n)	*/
    pp = (nml_fscalar**)nml_calloc(4+m, sizeof(nml_fscalar*));
    if (NULL != pp) {
      *((nml_index*)pp) = ub1; ++pp;
      *((nml_index*)pp) = lb1; ++pp;
      *((nml_index*)pp) = ub2; ++pp;
      *((nml_index*)pp) = lb2; ++pp;
      if (0 < m) {
	if (0 < n) {
	  nml_fscalar*	p
	    = (nml_fscalar*)nml_calloc(m*n, sizeof(nml_fscalar));
	  if (NULL != p) {
	    pp[0] = p - lb1;
	    { nml_index	i = 1;
	      for (i = 1; i < m; ++i)
		pp[i] = pp[i-1] + n;
	      }
	    }
	  else {	/* (NULL == p)		*/
	    int	i = 0;
            for (i = 0; i < m; ++i)
	      pp[i] = NULL;
	    nml_message(
"In function nml_f2D_new(nml_index, nml_index, nml_index, nml_index):\n"
"unable to allocate memory.");
	    }
	  }
	else {		/* (n <= 0)		*/
	  int	i = 0;
          for (i = 0; i < m; ++i)
	    pp[i] = NULL;
	  }
	}
      }
    else {		/* (NULL == pp)		*/
      nml_message(
"In function nml_f2D_new(nml_index, nml_index, nml_index, nml_index):\n"
"unable to allocate memory.");
      }
    }
  return pp;
  }

/* Destructor								*/
void			(nml_f2D_delete)(nml_f2DBase* pM) {
  if (NULL != pM) {
    if (NULL != pM[0]) {
      nml_index   lb1 = nml_f2D_lower1(pM);
      nml_free(pM[0] + lb1);
      }
    nml_free(pM - 4);
    }
  }

/* Reconstructor							*/
nml_f2DBase*	(nml_f2D_resize)(nml_f2DBase* *ppM,
    nml_index lb2, nml_index ub2, nml_index lb1, nml_index ub1) {
  if (NULL != ppM) {
    (nml_f2D_delete)(*ppM);
    *ppM = (nml_f2D_new)(lb2, ub2, lb1, ub1);
    return *ppM;
    }
  else {
    nml_message(
"In function nml_f2D_resize(nml_f2DBase**, nml_index, nml_index,\n"
"  nml_index, nml_index): invalid pointer to nml_f2DBase*.");
    return NULL;
    }
  }

