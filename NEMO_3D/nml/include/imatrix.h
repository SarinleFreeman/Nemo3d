#ifndef	IMATRIX_H
#define	IMATRIX_H 1

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
$Header: /repo/nml/include/imatrix.h,v 1.3 2003/10/08 16:19:16 hook Exp $
*****************************************************************************/


#include "nml_imatrix.h"
#include "matrix.h"

typedef nml_imatrix* imatrix;

#ifdef	NML_INLINE
inline static
imatrix	Imatrix(int n, int m) {
  return nml_im_new(n, m);
  }
inline static
void	rm_imatrix(imatrix* ppM) {
  if (NULL != ppM) {
    nml_im_delete(*ppM);
    *ppM = NULL;
    }
  }
inline static
imatrix	check_imatrix(imatrix pM, int n, int m) {
  if (NULL == pM)
    pM = nml_im_new(n, m);
  else {		/* (NULL != pM)	*/
    if ((nml_im_extent2(pM) != n) || (nml_im_extent1(pM) != m))
      pM = nml_im_resize(&pM, n, m);
    }
  return pM;
  }
#else /*NML_INLINE	*/
/* Constructor								*/
imatrix	Imatrix(int, int);
/* Destructor								*/
void	rm_imatrix(imatrix*);
/* Reconstructor							*/
imatrix	check_imatrix(imatrix, int, int);
imatrix	check_imatrix_min(imatrix, int, int);
imatrix	check_clear_imatrix(imatrix, int);
imatrix	clear_imatrix(imatrix, int);
imatrix	cp_imatrix(imatrix, const imatrix);
imatrix	cp_imatrix2to1(const imatrix, const imatrix);
imatrix	cp_imatrix3to1(const imatrix, const imatrix, const imatrix);
imatrix	cp_imatrix4to1(const imatrix, const imatrix, const imatrix, const imatrix);
#endif/*NML_INLINE	*/

#ifdef	NML_MACROS
#define	Imatrix(n, m) nml_im_new((n), (m))
#define rm_imatrix(ppM) \
do { \
  if (NULL != (ppM)) { \
    nml_im_delete(*(ppM)); \
    *(ppM) = NULL; \
    } \
  } while(0)
#define check_imatrix(pM, n, m) \
((NULL == (pM))? nml_im_new((n), (m)): \
((nml_im_extent2(pM) == (n)) && (nml_im_extent1(pM) == (m)))? (pM): \
nml_im_resize(&(pM), (n), (m)))
#endif/*NML_MACROS	*/

#endif/*IMATRIX_H	*/

