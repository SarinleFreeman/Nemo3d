#ifndef	BVECTOR_H
#define	BVECTOR_H 1

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


#include "vector.h"
#include "nml_bvector.h"

typedef nml_bvector* bvectr;

#ifdef	NML_INLINE
inline static
bvectr	Bvectr(int n) {
  return nml_bv_new(n);
  }
inline static
void	rm_bvectr(bvectr* ppv) {
  if (NULL != ppv) {
    nml_bv_delete(*ppv);
    *ppv = NULL;
    }
  }
inline static
bvectr	check_bvectr(bvectr pv, int n) {
  if (NULL == pv)
    pv = nml_bv_new(n);
  else			/* (NULL != pv)	*/
    if (nml_bv_extent(pv) != (nml_extent)n)
      pv = nml_bv_resize(&pv, n);
  return pv;
  }
#else /*NML_INLINE	*/
/* Constructor								*/
bvectr	Bvectr(int);
/* Destructor								*/
void	rm_bvectr(bvectr*);
/* Reconstructor							*/
bvectr	check_bvectr(bvectr, int);
bvectr	check_bvectr_min(bvectr, int);
bvectr	check_clear_bvectr(bvectr, int);
bvectr	clear_bvectr(bvectr, int);
bvectr	cp_bvectr(bvectr, const bvectr);
bvectr	cp_bvectr2to1(const bvectr, const bvectr);
bvectr	cp_bvectr3to1(const bvectr, const bvectr, const bvectr);
bvectr	cp_bvectr4to1(const bvectr, const bvectr, const bvectr, const bvectr);
#endif/*NML_INLINE	*/
bvectr	cp_bvectr(bvectr, const bvectr);
bvectr	cp_bvectr2to1(const bvectr, const bvectr);
bvectr	cp_bvectr3to1(const bvectr, const bvectr, const bvectr);
bvectr	cp_bvectr4to1(const bvectr, const bvectr, const bvectr, const bvectr);

#ifdef	NML_MACROS
#define	Bvectr(n) nml_bv_new(n)
#define rm_bvectr(ppv) \
do { \
  if (NULL != (ppv)) { \
    nml_bv_delete(*(ppv)); \
    *(ppv) = NULL; \
    } \
  } while(0)
#define check_bvectr(pv, n) \
((NULL == (pv))? nml_bv_new(n): \
(nml_bv_extent(pv) == (n))? (pv): nml_bv_resize(&(pv), (n)))
#endif/*NML_MACROS	*/

#endif/*BVECTOR_H	*/

