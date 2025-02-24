#ifndef	FVECTOR_H
#define	FVECTOR_H 1

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
#include "nml_fvector.h"

typedef nml_fvector* fvectr;

#ifdef	NML_INLINE
inline static
fvectr	Fvectr(int n) {
  return nml_fv_new(n);
  }
inline static
void	rm_fvectr(fvectr* ppv) {
  if (NULL != ppv) {
    nml_fv_delete(*ppv);
    *ppv = NULL;
    }
  }
inline static
fvectr	check_fvectr(fvectr pv, int n) {
  if (NULL == pv)
    pv = nml_fv_new(n);
  else			/* (NULL != pv)	*/
    if (nml_fv_extent(pv) != (nml_extent)n)
      pv = nml_fv_resize(&pv, n);
  return pv;
  }
#else /*NML_INLINE	*/
/* Constructor								*/
fvectr	Fvectr(int);
/* Destructor								*/
void	rm_fvectr(fvectr*);
/* Reconstructor							*/
fvectr	check_fvectr(fvectr, int);
fvectr	check_fvectr_min(fvectr, int);
fvectr	check_clear_fvectr(fvectr, int);
fvectr	clear_fvectr(fvectr, int);
fvectr	cp_fvectr(fvectr, const fvectr);
fvectr	cp_fvectr2to1(const fvectr, const fvectr);
fvectr	cp_fvectr3to1(const fvectr, const fvectr, const fvectr);
fvectr	cp_fvectr4to1(const fvectr, const fvectr, const fvectr, const fvectr);
#endif/*NML_INLINE	*/
fvectr	cp_fvectr(fvectr, const fvectr);
fvectr	cp_fvectr2to1(const fvectr, const fvectr);
fvectr	cp_fvectr3to1(const fvectr, const fvectr, const fvectr);
fvectr	cp_fvectr4to1(const fvectr, const fvectr, const fvectr, const fvectr);

#ifdef	NML_MACROS
#define	Fvectr(n) nml_fv_new(n)
#define rm_fvectr(ppv) \
do { \
  if (NULL != (ppv)) { \
    nml_fv_delete(*(ppv)); \
    *(ppv) = NULL; \
    } \
  } while(0)
#define check_fvectr(pv, n) \
((NULL == (pv))? nml_fv_new(n): \
(nml_fv_extent(pv) == (n))? (pv): nml_fv_resize(&(pv), (n)))
#endif/*NML_MACROS	*/

#endif/*FVECTOR_H	*/

