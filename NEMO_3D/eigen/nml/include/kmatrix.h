#ifndef	KMATRIX_H
#define	KMATRIX_H 1

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


#include "nml_fcmatrix.h"
#include "matrix.h"

typedef nml_fcmatrix* kmatrix;

#ifdef	NML_INLINE
inline static
kmatrix	Kmatrix(int n, int m) {
  return nml_fcm_new(n, m);
  }
inline static
void	rm_kmatrix(kmatrix* ppM) {
  if (NULL != ppM) {
    nml_fcm_delete(*ppM);
    *ppM = NULL;
    }
  }
inline static
kmatrix	check_kmatrix(kmatrix pM, int n, int m) {
  if (NULL == pM)
    pM = nml_fcm_new(n, m);
  else {		/* (NULL != pM)	*/
    if ((nml_fcm_extent2(pM) != n) || (nml_fcm_extent1(pM) != m))
      pM = nml_fcm_resize(&pM, n, m);
    }
  return pM;
  }
#else /*NML_INLINE	*/
/* Constructor								*/
kmatrix	Kmatrix(int, int);
/* Destructor								*/
void	rm_kmatrix(kmatrix*);
/* Reconstructor							*/
kmatrix	check_kmatrix(kmatrix, int, int);
kmatrix	check_kmatrix_min(kmatrix, int, int);
kmatrix	check_clear_kmatrix(kmatrix, int);
kmatrix	clear_kmatrix(kmatrix, int);
kmatrix	cp_kmatrix(kmatrix, const kmatrix);
kmatrix	cp_kmatrix2to1(const kmatrix, const kmatrix);
kmatrix	cp_kmatrix3to1(const kmatrix, const kmatrix, const kmatrix);
kmatrix	cp_kmatrix4to1(const kmatrix, const kmatrix, const kmatrix, const kmatrix);
#endif/*NML_INLINE	*/

#ifdef	NML_MACROS
#define	Kmatrix(n, m) nml_fcm_new((n), (m))
#define rm_kmatrix(ppM) \
do { \
  if (NULL != (ppM)) { \
    nml_fcm_delete(*(ppM)); \
    *(ppM) = NULL; \
    } \
  } while(0)
#define check_kmatrix(pM, n, m) \
((NULL == (pM))? nml_fcm_new((n), (m)): \
((nml_fcm_extent2(pM) == (n)) && (nml_fcm_extent1(pM) == (m)))? (pM): \
nml_fcm_resize(&(pM), (n), (m)))

#endif/*NML_MACROS	*/

#endif/*KMATRIX_H	*/

