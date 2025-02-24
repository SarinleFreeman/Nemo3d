#ifndef	IVECTOR_H
#define	IVECTOR_H 1

/*****************************************************************************
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
$Header: /repo/nml/include/ivector.h,v 1.3 2003/10/08 16:19:18 hook Exp $
*****************************************************************************/


#include "vector.h"
#include "nml_ivector.h"

typedef nml_ivector* ivectr;

#ifdef	NML_INLINE
inline static
ivectr	Ivectr(int n) {
  return nml_iv_new(n);
  }
inline static
void	rm_ivectr(ivectr* ppv) {
  if (NULL != ppv) {
    nml_iv_delete(*ppv);
    *ppv = NULL;
    }
  }
inline static
ivectr	check_ivectr(ivectr pv, int n) {
  if (NULL == pv)
    pv = nml_iv_new(n);
  else			/* (NULL != pv)	*/
    if (nml_iv_extent(pv) != (nml_extent)n)
      pv = nml_iv_resize(&pv, n);
  return pv;
  }
#else /*NML_INLINE	*/
/* Constructor								*/
ivectr	Ivectr(int);
/* Destructor								*/
void	rm_ivectr(ivectr*);
/* Reconstructor							*/
ivectr	check_ivectr(ivectr, int);
ivectr	check_ivectr_min(ivectr, int);
ivectr	check_clear_ivectr(ivectr, int);
ivectr	clear_ivectr(ivectr, int);
ivectr	cp_ivectr(ivectr, const ivectr);
ivectr	cp_ivectr2to1(const ivectr, const ivectr);
ivectr	cp_ivectr3to1(const ivectr, const ivectr, const ivectr);
ivectr	cp_ivectr4to1(const ivectr, const ivectr, const ivectr, const ivectr);
#endif/*NML_INLINE	*/

#ifdef	NML_MACROS
#define	Ivectr(n) nml_iv_new(n)
#define rm_ivectr(ppv) \
do { \
  if (NULL != (ppv)) { \
    nml_iv_delete(*(ppv)); \
    *(ppv) = NULL; \
    } \
  } while(0)
#define check_ivectr(pv, n) \
((NULL == (pv))? nml_iv_new(n): \
(nml_iv_extent(pv) == (n))? (pv): nml_iv_resize(&(pv), (n)))
#endif/*NML_MACROS	*/

#endif/*IVECTOR_H	*/

