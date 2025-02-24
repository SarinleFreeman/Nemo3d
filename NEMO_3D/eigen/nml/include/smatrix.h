#ifndef	SMATRIX_H
#define	SMATRIX_H 1

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


#include "nml_smatrix.h"
#include "matrix.h"

typedef nml_smatrix* smatrix;

#ifdef	NML_INLINE
inline static
smatrix	Smatrix(int n, int m) {
  return nml_sm_new(n, m);
  }
inline static
void	rm_smatrix(smatrix* ppM) {
  if (NULL != ppM) {
    nml_sm_delete(*ppM);
    *ppM = NULL;
    }
  }
inline static
smatrix	check_smatrix(smatrix pM, int n, int m) {
  if (NULL == pM)
    pM = nml_sm_new(n, m);
  else {		/* (NULL != pM)	*/
    if ((nml_sm_extent2(pM) != n) || (nml_sm_extent1(pM) != m))
      pM = nml_sm_resize(&pM, n, m);
    }
  return pM;
  }
#else /*NML_INLINE	*/
/* Constructor								*/
smatrix	Smatrix(int, int);
/* Destructor								*/
void	rm_smatrix(smatrix*);
/* Reconstructor							*/
smatrix	check_smatrix(smatrix, int, int);
smatrix	check_smatrix_min(smatrix, int, int);
smatrix	check_clear_smatrix(smatrix, int);
smatrix	clear_smatrix(smatrix, int);
smatrix	cp_smatrix(smatrix, const smatrix);
smatrix	cp_smatrix2to1(const smatrix, const smatrix);
smatrix	cp_smatrix3to1(const smatrix, const smatrix, const smatrix);
smatrix	cp_smatrix4to1(const smatrix, const smatrix, const smatrix, const smatrix);
#endif/*NML_INLINE	*/

#ifdef	NML_MACROS
#define	Smatrix(n, m) nml_sm_new((n), (m))
#define rm_smatrix(ppM) \
do { \
  if (NULL != (ppM)) { \
    nml_sm_delete(*(ppM)); \
    *(ppM) = NULL; \
    } \
  } while(0)
#define check_smatrix(pM, n, m) \
((NULL == (pM))? nml_sm_new((n), (m)): \
((nml_sm_extent2(pM) == (n)) && (nml_sm_extent1(pM) == (m)))? (pM): \
nml_sm_resize(&(pM), (n), (m)))

#endif/*NML_MACROS	*/

#endif/*SMATRIX_H	*/

