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
$Header: /repo/nml/src/scalar/dcscalar/complex.c,v 1.3 2003/10/08 16:19:50 hook Exp $
*****************************************************************************/


#include "nml_dcscalar.h"
#undef NML_INLINE
#undef NML_MACROS
#include "complex.h"

/* extern const COMPLEX czero={0.0,0.0};				*/
/* extern const COMPLEX cone ={1.0,0.0};				*/
/* global definition of a number that is used often			*/

/* Constructor								*/
nml_dcscalar	Complex(nml_dscalar r, nml_dscalar i) {
  return nml_dcmplx(r, i);
  }
#define Complex(r, i) nml_dcmplx((r), (i))

nml_dcscalar cmul(nml_dcscalar x1, nml_dcscalar x2) {
  nml_dcscalar result;
  NML_DCMUL(result, x1, x2);  
  return result;
  }

nml_dcscalar cadd(nml_dcscalar x1, nml_dcscalar x2) {
  nml_dcscalar result;
  NML_DCADD(result, x1, x2);  
  return result;
  }

nml_dcscalar csub(nml_dcscalar x1, nml_dcscalar x2) {
  nml_dcscalar result;
  NML_DCSUB(result, x1, x2);  
  return result;
  }

nml_dcscalar cxdiv(nml_dcscalar x1, nml_dcscalar x2) {
  nml_dcscalar result;
  NML_DCDIV(result, x1, x2);  
  return result;
  }

nml_dcscalar cconj(nml_dcscalar x1) {
  nml_dcscalar result = x1;
  result.r = +x1.r;
  result.i = -x1.i;
  return result;
  }

nml_dscalar c_abs(nml_dcscalar x) {
  nml_dscalar result = sqrt(x.r*x.r + x.i*x.i);
  return result;
  }
#define c_abs(x) sqrt((x).r*(x).r + (x).i*(x).i)

nml_dcscalar cxsqrt(nml_dcscalar x) {
  nml_dcscalar result = czero;
  if (0.0 != x.i) {
    nml_dscalar norm = c_abs(x);
    nml_dscalar sgn = (x.i < 0.0)? -1.0: +1.0;
    result = Complex(sqrt(0.5*(norm + x.r)), sgn*sqrt(0.5*(norm - x.r)));
    }
  else {			/* (0.0 == x)				*/
    if (0.0 < x.r){
      result = Complex(sqrt(+x.r), 0.0);
      }
    else {
      result = Complex(0.0, sqrt(-x.r));
      }
    }
  return result;
  }	/* nml_dcscalar cxsqrt	*/

nml_dcscalar RCmul(nml_dscalar r, nml_dcscalar x) {
  nml_dcscalar result;
  result.r = r*x.r;
  result.i = r*x.i;
  return result;
  }
 
