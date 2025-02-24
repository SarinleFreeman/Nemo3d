#ifndef	SVECTOR_H
#define	SVECTOR_H 1

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
#include "nml_svector.h"

typedef nml_svector* svectr;

#ifdef	NML_INLINE
inline static
svectr	Svectr(int n) {
  return nml_sv_new(n);
  }
inline static
void	rm_svectr(svectr* ppv) {
  if (NULL != ppv) {
    nml_sv_delete(*ppv);
    *ppv = NULL;
    }
  }
inline static
svectr	check_svectr(svectr pv, int n) {
  if (NULL == pv)
    pv = nml_sv_new(n);
  else			/* (NULL != pv)	*/
    if (nml_sv_extent(pv) != (nml_extent)n)
      pv = nml_sv_resize(&pv, n);
  return pv;
  }
#else /*NML_INLINE	*/
/* Constructor								*/
svectr	Svectr(int);
/* Destructor								*/
void	rm_svectr(svectr*);
/* Reconstructor							*/
svectr	check_svectr(svectr, int);
svectr	check_svectr_min(svectr, int);
svectr	check_clear_svectr(svectr, int);
svectr	clear_svectr(svectr, int);
svectr	cp_svectr(svectr, const svectr);
svectr	cp_svectr2to1(const svectr, const svectr);
svectr	cp_svectr3to1(const svectr, const svectr, const svectr);
svectr	cp_svectr4to1(const svectr, const svectr, const svectr, const svectr);
#endif/*NML_INLINE	*/
svectr	cp_svectr(svectr, const svectr);
svectr	cp_svectr2to1(const svectr, const svectr);
svectr	cp_svectr3to1(const svectr, const svectr, const svectr);
svectr	cp_svectr4to1(const svectr, const svectr, const svectr, const svectr);

#ifdef	NML_MACROS
#define	Svectr(n) nml_sv_new(n)
#define rm_svectr(ppv) \
do { \
  if (NULL != (ppv)) { \
    nml_sv_delete(*(ppv)); \
    *(ppv) = NULL; \
    } \
  } while(0)
#define check_svectr(pv, n) \
((NULL == (pv))? nml_sv_new(n): \
(nml_sv_extent(pv) == (n))? (pv): nml_sv_resize(&(pv), (n)))
#endif/*NML_MACROS	*/

#endif/*SVECTOR_H	*/

