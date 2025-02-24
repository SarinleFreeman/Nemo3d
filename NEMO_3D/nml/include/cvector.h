#ifndef	CVECTOR_H
#define	CVECTOR_H 1

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
$Header: /repo/nml/include/cvector.h,v 1.3 2003/10/08 16:19:13 hook Exp $
*****************************************************************************/

#include "nml_dcvector.h"

typedef nml_dcvector* cvectr;

#ifdef	NML_INLINE
inline static
cvectr	Cvectr(int n) {
  return nml_dcv_new(n);
  }
inline static
void	rm_cvectr(cvectr* ppv) {
  if (NULL != ppv) {
    nml_dcv_delete(*ppv);
    *ppv = NULL;
    }
  }
inline static
rvectr	check_cvectr(cvectr pv, int n) {
  if (NULL == pv)
    pv = nml_dcv_new(n);
  else {		/* (NULL != pv)	*/
    if (nml_dcv_extent(pv) != n)
      pv = nml_dcv_resize(&pv, n)
  return pv;
  }
#else /*NML_INLINE	*/
/* Constructor								*/
cvectr	Cvectr(int);
/* Destructor								*/
void	rm_cvectr(cvectr*);
/* Reconstructor							*/
cvectr	check_cvectr(cvectr, int);
cvectr	check_cvectr_min(cvectr, int);
cvectr	check_clear_cvectr(cvectr, int);
cvectr	check_clear(cvectr, int);
cvectr	clear_cvectr(cvectr, int);
cvectr	cp_cvectr(cvectr, const cvectr);
cvectr	cp_cvectr2to1(const cvectr, const cvectr);
cvectr	cp_cvectr3to1(const cvectr, const cvectr, const cvectr);
cvectr	cp_cvectr4to1(const cvectr, const cvectr, const cvectr, const cvectr);
#endif/*NML_INLINE	*/

#ifdef	NML_MACROS
#define	Cvectr(n) nml_dcv_new(n)
#define rm_cvectr(ppv) \
do { \
  if (NULL != (ppv)) { \
    nml_dcv_delete(*(ppv)); \
    *(ppv) = NULL; \
    } \
  } while(0)
#define check_cvectr(pv, n) \
((NULL == (pv))? nml_dcv_new(n):  \
(nml_dcv_extent(pv) == (n))? (pv): nml_dcv_resize(&(pv), (n)))
#endif/*NML_MACROS	*/

#endif/*CVECTOR_H	*/

