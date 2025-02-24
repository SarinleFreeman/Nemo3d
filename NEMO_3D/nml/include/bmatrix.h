#ifndef	BMATRIX_H
#define	BMATRIX_H 1

/*****************************************************************************
Copyright (C) 2002-2004 California Institute of Technology (Caltech)

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

Written by Marek J. Korkusinski

*****************************************************************************/


#include "nml_bmatrix.h"
#include "matrix.h"

typedef nml_bmatrix* bmatrix;

#ifdef	NML_INLINE
inline static
bmatrix	Bmatrix(int n, int m) {
  return nml_bm_new(n, m);
  }
inline static
void	rm_bmatrix(bmatrix* ppM) {
  if (NULL != ppM) {
    nml_bm_delete(*ppM);
    *ppM = NULL;
    }
  }
inline static
bmatrix	check_bmatrix(bmatrix pM, int n, int m) {
  if (NULL == pM)
    pM = nml_bm_new(n, m);
  else {		/* (NULL != pM)	*/
    if ((nml_bm_extent2(pM) != n) || (nml_bm_extent1(pM) != m))
      pM = nml_bm_resize(&pM, n, m);
    }
  return pM;
  }
#else /*NML_INLINE	*/
/* Constructor								*/
bmatrix	Bmatrix(int, int);
/* Destructor								*/
void	rm_bmatrix(bmatrix*);
/* Reconstructor							*/
bmatrix	check_bmatrix(bmatrix, int, int);
bmatrix	check_bmatrix_min(bmatrix, int, int);
bmatrix	check_clear_bmatrix(bmatrix, int);
bmatrix	clear_bmatrix(bmatrix, int);
bmatrix	cp_bmatrix(bmatrix, const bmatrix);
bmatrix	cp_bmatrix2to1(const bmatrix, const bmatrix);
bmatrix	cp_bmatrix3to1(const bmatrix, const bmatrix, const bmatrix);
bmatrix	cp_bmatrix4to1(const bmatrix, const bmatrix, const bmatrix, const bmatrix);
#endif/*NML_INLINE	*/

#ifdef	NML_MACROS
#define	Bmatrix(n, m) nml_bm_new((n), (m))
#define rm_bmatrix(ppM) \
do { \
  if (NULL != (ppM)) { \
    nml_bm_delete(*(ppM)); \
    *(ppM) = NULL; \
    } \
  } while(0)
#define check_bmatrix(pM, n, m) \
((NULL == (pM))? nml_bm_new((n), (m)): \
((nml_bm_extent2(pM) == (n)) && (nml_bm_extent1(pM) == (m)))? (pM): \
nml_bm_resize(&(pM), (n), (m)))

#endif/*NML_MACROS	*/

#endif/*BMATRIX_H	*/

