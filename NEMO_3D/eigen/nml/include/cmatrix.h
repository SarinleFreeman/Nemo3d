#ifndef	CMATRIX_H
#define	CMATRIX_H 1

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
$Header: /repo/nml/include/cmatrix.h,v 1.4 2004/12/08 18:40:04 dseaman Exp $
*****************************************************************************/


#include "nml_dcmatrix.h"
#include "matrix.h"

typedef nml_dcmatrix* cmatrix;

#ifdef	NML_INLINE
inline static
cmatrix	Cmatrix(int n, int m) {
  return nml_dcm_new(n, m);
  }
inline static
void	rm_cmatrix(cmatrix* ppM) {
  if (NULL != ppM) {
    nml_dcm_delete(*ppM);
    *ppM = NULL;
    }
  }
inline static
cmatrix	check_cmatrix(cmatrix pM, int n, int m) {
  if (NULL == pM)
    pM = nml_dcm_new(n, m);
  else {		/* (NULL != pM)	*/
    if ((nml_dcm_extent2(pM) != n) || (nml_dcm_extent1(pM) != m))
      pM = nml_dcm_resize(&pM, n, m);
    }
  return pM;
  }
#else /*NML_INLINE	*/
/* Constructor								*/
cmatrix	Cmatrix(int, int);
/* Destructor								*/
void	rm_cmatrix(cmatrix*);
/* Reconstructor							*/
cmatrix	check_cmatrix(cmatrix, int, int);
cmatrix	check_cmatrix_min(cmatrix, int, int);
cmatrix	check_clear_cmatrix(cmatrix, int);
cmatrix	clear_cmatrix(cmatrix, int);
#ifdef _AIX
cmatrix	cp_cmatrix(cmatrix, cmatrix);
#else
cmatrix	cp_cmatrix(cmatrix, const cmatrix);
#endif
cmatrix	cp_cmatrix2to1(const cmatrix, const cmatrix);
cmatrix	cp_cmatrix3to1(const cmatrix, const cmatrix, const cmatrix);
cmatrix	cp_cmatrix4to1(const cmatrix, const cmatrix, const cmatrix, const cmatrix);
#endif/*NML_INLINE	*/

#ifdef	NML_MACROS
#define	Cmatrix(n, m) nml_dcm_new((n), (m))
#define rm_cmatrix(ppM) \
do { \
  if (NULL != (ppM)) { \
    nml_dcm_delete(*(ppM)); \
    *(ppM) = NULL; \
    } \
  } while(0)
#define check_cmatrix(pM, n, m) \
((NULL == (pM))? nml_dcm_new((n), (m)): \
((nml_dcm_extent2(pM) == (n)) && (nml_dcm_extent1(pM) == (m)))? (pM): \
nml_dcm_resize(&(pM), (n), (m)))

#endif/*NML_MACROS	*/

#endif/*CMATRIX_H	*/

