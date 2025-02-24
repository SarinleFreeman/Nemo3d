#ifndef	RVECTOR_H
#define	RVECTOR_H 1

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
$Header: /repo/nml/include/rvector.h,v 1.4 2004/12/08 18:40:04 dseaman Exp $
*****************************************************************************/


#include "vector.h"
#include "nml_dvector.h"

typedef nml_dvector* rvectr;

#ifdef	NML_INLINE
inline static
rvectr	Rvectr(int n) {
  return nml_dv_new(n);
  }
inline static
void	rm_rvectr(rvectr* ppv) {
  if (NULL != ppv) {
    nml_dv_delete(*ppv);
    *ppv = NULL;
    }
  }
inline static
rvectr	check_rvectr(rvectr pv, int n) {
  if (NULL == pv)
    pv = nml_dv_new(n);
  else			/* (NULL != pv)	*/
    if (nml_dv_extent(pv) != (nml_extent)n)
      pv = nml_dv_resize(&pv, n);
  return pv;
  }
#else /*NML_INLINE	*/
/* Constructor								*/
rvectr	Rvectr(int);
/* Destructor								*/
void	rm_rvectr(rvectr*);
/* Reconstructor							*/
rvectr	check_rvectr(rvectr, int);
rvectr	check_rvectr_min(rvectr, int);
rvectr	check_clear_rvectr(rvectr, int);
rvectr	clear_rvectr(rvectr, int);
#ifdef _AIX
rvectr	cp_rvectr(rvectr, rvectr);
#else
rvectr	cp_rvectr(rvectr, const rvectr);
#endif
rvectr	cp_rvectr2to1(const rvectr, const rvectr);
rvectr	cp_rvectr3to1(const rvectr, const rvectr, const rvectr);
rvectr	cp_rvectr4to1(const rvectr, const rvectr, const rvectr, const rvectr);
#endif/*NML_INLINE	*/
rvectr	cp_rvectr(rvectr, const rvectr);
rvectr	cp_rvectr2to1(const rvectr, const rvectr);
rvectr	cp_rvectr3to1(const rvectr, const rvectr, const rvectr);
rvectr	cp_rvectr4to1(const rvectr, const rvectr, const rvectr, const rvectr);

#ifdef	NML_MACROS
#define	Rvectr(n) nml_dv_new(n)
#define rm_rvectr(ppv) \
do { \
  if (NULL != (ppv)) { \
    nml_dv_delete(*(ppv)); \
    *(ppv) = NULL; \
    } \
  } while(0)
#define check_rvectr(pv, n) \
((NULL == (pv))? nml_dv_new(n): \
(nml_dv_extent(pv) == (n))? (pv): nml_dv_resize(&(pv), (n)))
#endif/*NML_MACROS	*/

#endif/*RVECTOR_H	*/

