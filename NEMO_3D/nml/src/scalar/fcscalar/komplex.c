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


#include "nml_fcscalar.h"
#undef NML_INLINE
#undef NML_MACROS
#include "komplex.h"

/* extern const COMPLEX czero={0.0,0.0};				*/
/* extern const COMPLEX cone ={1.0,0.0};				*/
/* global definition of a number that is used often			*/

/* Constructor								*/
nml_fcscalar	Komplex(nml_fscalar r, nml_fscalar i) {
  return nml_fcmplx(r, i);
  }
#define Komplex(r, i) nml_fcmplx((r), (i))

nml_fcscalar kmul(nml_fcscalar x1, nml_fcscalar x2) {
  nml_fcscalar result;
  NML_FCMUL(result, x1, x2);  
  return result;
  }

nml_fcscalar kadd(nml_fcscalar x1, nml_fcscalar x2) {
  nml_fcscalar result;
  NML_FCADD(result, x1, x2);  
  return result;
  }

nml_fcscalar ksub(nml_fcscalar x1, nml_fcscalar x2) {
  nml_fcscalar result;
  NML_FCSUB(result, x1, x2);  
  return result;
  }

nml_fcscalar kxdiv(nml_fcscalar x1, nml_fcscalar x2) {
  nml_fcscalar result;
  NML_FCDIV(result, x1, x2);  
  return result;
  }

nml_fcscalar kconj(nml_fcscalar x1) {
  nml_fcscalar result = x1;
  result.r = +x1.r;
  result.i = -x1.i;
  return result;
  }

nml_fscalar k_abs(nml_fcscalar x) {
  nml_fscalar result = sqrt(x.r*x.r + x.i*x.i);
  return result;
  }
#define k_abs(x) sqrt((x).r*(x).r + (x).i*(x).i)

nml_fcscalar kxsqrt(nml_fcscalar x) {
  nml_fcscalar result = kzero;
  if (0.0 != x.i) {
    nml_fscalar norm = k_abs(x);
    nml_fscalar sgn = (x.i < 0.0)? -1.0: +1.0;
    result = Komplex(sqrt(0.5*(norm + x.r)), sgn*sqrt(0.5*(norm - x.r)));
    }
  else {			/* (0.0 == x)				*/
    if (0.0 < x.r){
      result = Komplex(sqrt(+x.r), 0.0);
      }
    else {
      result = Komplex(0.0, sqrt(-x.r));
      }
    }
  return result;
  }	/* nml_fcscalar cxsqrt	*/

nml_fcscalar RKmul(nml_fscalar r, nml_fcscalar x) {
  nml_fcscalar result;
  result.r = r*x.r;
  result.i = r*x.i;
  return result;
  }
 
