#ifndef	RMATRIX_H
#define	RMATRIX_H 1

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
$Header: /repo/nml/include/rmatrix.h,v 1.3 2003/10/08 16:19:22 hook Exp $
*****************************************************************************/


#include "nml_dmatrix.h"
#include "matrix.h"

typedef nml_dmatrix* rmatrix;

#ifdef	NML_INLINE
inline static
rmatrix	Rmatrix(int n, int m) {
  return nml_dm_new(n, m);
  }
inline static
void	rm_rmatrix(rmatrix* ppM) {
  if (NULL != ppM) {
    nml_dm_delete(*ppM);
    *ppM = NULL;
    }
  }
inline static
rmatrix	check_rmatrix(rmatrix pM, int n, int m) {
  if (NULL == pM)
    pM = nml_dm_new(n, m);
  else {		/* (NULL != pM)	*/
    if ((nml_dm_extent2(pM) != n) || (nml_dm_extent1(pM) != m))
      pM = nml_dm_resize(&pM, n, m);
    }
  return pM;
  }
#else /*NML_INLINE	*/
/* Constructor								*/
rmatrix	Rmatrix(int, int);
/* Destructor								*/
void	rm_rmatrix(rmatrix*);
/* Reconstructor							*/
rmatrix	check_rmatrix(rmatrix, int, int);
rmatrix	check_rmatrix_min(rmatrix, int, int);
rmatrix	check_clear_rmatrix(rmatrix, int);
rmatrix	clear_rmatrix(rmatrix, int);
rmatrix	cp_rmatrix(rmatrix, const rmatrix);
rmatrix	cp_rmatrix2to1(const rmatrix, const rmatrix);
rmatrix	cp_rmatrix3to1(const rmatrix, const rmatrix, const rmatrix);
rmatrix	cp_rmatrix4to1(const rmatrix, const rmatrix, const rmatrix, const rmatrix);
#endif/*NML_INLINE	*/

#ifdef	NML_MACROS
#define	Rmatrix(n, m) nml_dm_new((n), (m))
#define rm_rmatrix(ppM) \
do { \
  if (NULL != (ppM)) { \
    nml_dm_delete(*(ppM)); \
    *(ppM) = NULL; \
    } \
  } while(0)
#define check_rmatrix(pM, n, m) \
((NULL == (pM))? nml_dm_new((n), (m)): \
((nml_dm_extent2(pM) == (n)) && (nml_dm_extent1(pM) == (m)))? (pM): \
nml_dm_resize(&(pM), (n), (m)))

#endif/*NML_MACROS	*/

#endif/*RMATRIX_H	*/

