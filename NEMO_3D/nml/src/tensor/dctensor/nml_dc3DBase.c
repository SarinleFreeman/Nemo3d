
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
$Header: /repo/nml/src/tensor/3DBase.cP,v 1.2 2003/10/08 16:19:53 hook Exp $
*****************************************************************************/


#undef NML_INLINE
#undef NML_MACROS
#include<nml_dc3DBase.h>

/* Functions								*/
nml_index		(nml_dc3D_lower3)(const nml_dc3DBase* pT) {
  return *((const nml_index*)(pT - 1));
  }
#define nml_dc3D_lower3(pT) (*((const nml_index*)((pT) - 1)))

nml_index		(nml_dc3D_upper3)(const nml_dc3DBase* pT) {
  return *((const nml_index*)(pT - 2));
  }
#define nml_dc3D_upper3(pT) (*((const nml_index*)((pT) - 2)))

nml_extent		(nml_dc3D_extent3)(const nml_dc3DBase* pT) {
  return 1+nml_dc3D_upper3(pT)-nml_dc3D_lower3(pT);
  }
#define nml_dc3D_extent3(pT) (1+nml_dc3D_upper3(pT)-nml_dc3D_lower3(pT))

nml_index		(nml_dc3D_lower2)(const nml_dc3DBase* pT) {
  return *((const nml_index*)(pT - 3));
  }
#define nml_dc3D_lower2(pT) (*((const nml_index*)((pT) - 3)))

nml_index		(nml_dc3D_upper2)(const nml_dc3DBase* pT) {
  return *((const nml_index*)(pT - 4));
  }
#define nml_dc3D_upper2(pT) (*((const nml_index*)((pT) - 4)))

nml_extent		(nml_dc3D_extent2)(const nml_dc3DBase* pT) {
  return 1+nml_dc3D_upper2(pT)-nml_dc3D_lower2(pT);
  }
#define nml_dc3D_extent2(pT) (1+nml_dc3D_upper2(pT)-nml_dc3D_lower2(pT))

nml_index		(nml_dc3D_lower1)(const nml_dc3DBase* pT) {
  return *((const nml_index*)(pT - 5));
  }
#define nml_dc3D_lower1(pT) (*((const nml_index*)((pT) - 5)))

nml_index		(nml_dc3D_upper1)(const nml_dc3DBase* pT) {
  return *((const nml_index*)(pT - 6));
  }
#define nml_dc3D_upper1(pT) (*((const nml_index*)((pT) - 6)))

nml_extent		(nml_dc3D_extent1)(const nml_dc3DBase* pT) {
  return 1+nml_dc3D_upper1(pT)-nml_dc3D_lower1(pT);
  }
#define nml_dc3D_extent1(pT) (1+nml_dc3D_upper1(pT)-nml_dc3D_lower1(pT))
int
(nml_dc3D_fprintf)(FILE* stream, const nml_dc3DBase* pT, 
  int width, unsigned int precision, nml_fmtflags flags, unsigned int columns) {
  int		characters;
  int		total = 0;
  const
  nml_extent	l = nml_dc3D_extent3(pT);
  nml_offset	h = 0;
  columns = (0 < columns)? columns: 4;
  for (h = 0; h < l; ++h) {
    const
    nml_extent	m = nml_dc3D_extent2(pT);
    nml_offset	i = 0;
    for (i = 0; i < m; ++i) {
      nml_index k2 = (nml_index)i + nml_dc3D_lower2(pT);
      const
      nml_extent	n = nml_dc3D_extent1(pT);
      nml_offset	j = 0;
      for (j = 0; j < n; ++j) {
	nml_index k1 = (nml_index)j + nml_dc3D_lower1(pT);
	if (0 < j) {
	  characters
	    = (0 == j%columns)? fprintf(stream, "\n"): fprintf(stream, " ");
	  if (0 <= characters) { total += characters; } else { return -1; }
	  }
	characters = nml_dcfprintf(
	  stream, pT[h][k2][k1], width, precision, flags);
	if (0 <= characters) { total += characters; } else { return -1; }
	}
      characters = fprintf(stream, "\n");
      if (0 <= characters) { total += characters; } else { return -1; }
      }
    }
  return total;
  }

/* Constructor								*/
nml_dc3DBase*	(nml_dc3D_new)(nml_index lb3, nml_index ub3,
    nml_index lb2, nml_index ub2, nml_index lb1, nml_index ub1) {
  nml_dcscalar***	ppp = NULL;
  nml_index		l = 1 + ub3 - lb3;
  nml_index		m = 1 + ub2 - lb2;
  nml_index		n = 1 + ub1 - lb1;
  if (l < 0) {
    nml_message(
"In function nml_dc3D_new(nml_index, nml_index, nml_index, nml_index,\n"
"  nml_index, nml_index): extent3 is negative.");
    }
  else
  if (m < 0) {
    nml_message(
"In function nml_dc3D_new(nml_index, nml_index, nml_index, nml_index,\n"
"  nml_index, nml_index): extent2 is negative.");
    }
  else
  if (n < 0) {
    nml_message(
"In function nml_dc3D_new(nml_index, nml_index, nml_index, nml_index,\n"
"  nml_index, nml_index): extent1 is negative.");
    }
  else {		/* (0 <= l) && (0 <= m) && (0 <= n)	*/
    ppp = (nml_dcscalar***)nml_calloc(6 + l, sizeof(nml_dcscalar**));
    if (NULL != ppp) {
      *((nml_index*)ppp) = ub1; ++ppp;
      *((nml_index*)ppp) = lb1; ++ppp;
      *((nml_index*)ppp) = ub2; ++ppp;
      *((nml_index*)ppp) = lb2; ++ppp;
      *((nml_index*)ppp) = ub3; ++ppp;
      *((nml_index*)ppp) = lb3; ++ppp;
      if (0 < l) {
	if (0 < m) {	/* (0 < l) && (0 < m)	*/
	  nml_dcscalar** pp
	    = (nml_dcscalar**)nml_calloc(l*m, sizeof(nml_dcscalar*));
	  if (NULL != pp) {
	    ppp[0] = pp - lb2;
	    { nml_index	h = 1;
	      for (h = 1; h < l; ++h)
		ppp[h] = ppp[h-1] + m;
	      }
	    if (0 < n) {
	      nml_dcscalar*	p
		= (nml_dcscalar*)nml_calloc(l*m*n, sizeof(nml_dcscalar));
	      if (NULL != p) {
		ppp[0][lb2] = p - lb1;
		{ nml_index	h = 1;
		  for (h = 1; h < l; ++h)
		    ppp[h][lb2] = ppp[h-1][lb2] + m*n;
		  for (h = 0; h < l; ++h) {
		    nml_index	i = lb2 + 1;
		    for (i = lb2 + 1; i <= ub2; ++i)
		      ppp[h][i] = ppp[h][i-1] + n;
		    }
		  }
		}

	      else {	/* (NULL == p)		*/
		nml_index	h = 0;
		for (h = 0; h < l; ++h) {
		  nml_index	i = lb2;
		  for (i = lb2; i <= ub2; ++i)
		    ppp[h][i] = NULL;
		  }
		nml_message(
"In function nml_dc3D_new(nml_index, nml_index, nml_index, nml_index,\n"
"  nml_index, nml_index): unable to allocate memory.");
		}
	      }
	    else {	/* (n <= 0)		*/
	      nml_index	h = 0;
	      for (h = 0; h < l; ++h) {
		nml_index	i = lb2;
		for (i = lb2; i <= ub2; ++i)
		  ppp[h][i] = NULL;
		}
	      }
	    }
	  else {	/* (NULL == pp)		*/
	    nml_index	h = 0;
	    for (h = 0; h < l; ++h)
	      ppp[h] = NULL;
	    nml_message(
"In function nml_dc3D_new(nml_index, nml_index, nml_index, nml_index,\n"
"  nml_index, nml_index): unable to allocate memory.");
	    }
	  }
	else {		/* (m <= 0)		*/
	  nml_index	h = 0;
	  for (h = 0; h < l; ++h)
	    ppp[h] = NULL;
	  }
	}
      }
    else {		/* (NULL == ppp)	*/
      nml_message(
"In function nml_dc3D_new(nml_index, nml_index, nml_index, nml_index,\n"
"  nml_index, nml_index): unable to allocate memory.");
      }
    }
  return ppp;
  }

/* Destructor								*/
void			(nml_dc3D_delete)(nml_dc3DBase* pT) {
  if (NULL != pT) {
    if (NULL != pT[0]) {
      nml_index	lb2 = nml_dc3D_lower2(pT);
      if (NULL != pT[0][lb2]) {
	nml_index	lb1 = nml_dc3D_lower1(pT);
	nml_free(pT[0][lb2] + lb1);
	}
      nml_free(pT[0] + lb2);
      }
    nml_free(pT - 6);
    }
  }

/* Reconstructor							*/
nml_dc3DBase*	(nml_dc3D_resize)(nml_dc3DBase* *ppT,
    nml_index lb3, nml_index ub3, nml_index lb2, nml_index ub2,
    nml_index lb1, nml_index ub1) {
  if (NULL != ppT) {
    (nml_dc3D_delete)(*ppT);
    *ppT = (nml_dc3D_new)(lb3, ub3, lb2, ub2, lb1, ub1);
    return *ppT;
    }
  else {
    nml_message(
"In function nml_dc3D_resize(nml_dc3DBase**,\n"
"  nml_index, nml_index, nml_index, nml_index, nml_index, nml_index):\n"
"invalid pointer to nml_dc3DBase*.");
    return NULL;
    }
  }

