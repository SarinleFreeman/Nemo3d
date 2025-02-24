#ifndef	FMATRIX_H
#define	FMATRIX_H 1

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


#include "nml_fmatrix.h"
#include "matrix.h"

typedef nml_fmatrix* fmatrix;

#ifdef	NML_INLINE
inline static
fmatrix	Fmatrix(int n, int m) {
  return nml_fm_new(n, m);
  }
inline static
void	rm_fmatrix(fmatrix* ppM) {
  if (NULL != ppM) {
    nml_fm_delete(*ppM);
    *ppM = NULL;
    }
  }
inline static
fmatrix	check_fmatrix(fmatrix pM, int n, int m) {
  if (NULL == pM)
    pM = nml_fm_new(n, m);
  else {		/* (NULL != pM)	*/
    if ((nml_fm_extent2(pM) != n) || (nml_fm_extent1(pM) != m))
      pM = nml_fm_resize(&pM, n, m);
    }
  return pM;
  }
#else /*NML_INLINE	*/
/* Constructor								*/
fmatrix	Fmatrix(int, int);
/* Destructor								*/
void	rm_fmatrix(fmatrix*);
/* Reconstructor							*/
fmatrix	check_fmatrix(fmatrix, int, int);
fmatrix	check_fmatrix_min(fmatrix, int, int);
fmatrix	check_clear_fmatrix(fmatrix, int);
fmatrix	clear_fmatrix(fmatrix, int);
fmatrix	cp_fmatrix(fmatrix, const fmatrix);
fmatrix	cp_fmatrix2to1(const fmatrix, const fmatrix);
fmatrix	cp_fmatrix3to1(const fmatrix, const fmatrix, const fmatrix);
fmatrix	cp_fmatrix4to1(const fmatrix, const fmatrix, const fmatrix, const fmatrix);
#endif/*NML_INLINE	*/

#ifdef	NML_MACROS
#define	Fmatrix(n, m) nml_fm_new((n), (m))
#define rm_fmatrix(ppM) \
do { \
  if (NULL != (ppM)) { \
    nml_fm_delete(*(ppM)); \
    *(ppM) = NULL; \
    } \
  } while(0)
#define check_fmatrix(pM, n, m) \
((NULL == (pM))? nml_fm_new((n), (m)): \
((nml_fm_extent2(pM) == (n)) && (nml_fm_extent1(pM) == (m)))? (pM): \
nml_fm_resize(&(pM), (n), (m)))

#endif/*NML_MACROS	*/

#endif/*FMATRIX_H	*/

