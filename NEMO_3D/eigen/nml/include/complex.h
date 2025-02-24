#ifndef	COMPLEX_H
#define	COMPLEX_H 1

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
$Header: /repo/nml/include/complex.h,v 1.3 2003/10/08 16:19:11 hook Exp $
*****************************************************************************/

#include "nml_dcscalar.h"


/* For compatibility with previous complex.h	*/
/* typedef COMPLEX complex;			*/
typedef nml_dcscalar complex;
/* extern const nml_dscalar czero;		*/
/* extern const nml_dscalar cone;		*/
const nml_dcscalar czero = {0.0, 0.0};
const nml_dcscalar cone  = {1.0, 0.0};

nml_dcscalar Complex(nml_dscalar r, nml_dscalar i );
nml_dcscalar cmul(nml_dcscalar x1, nml_dcscalar x2);
nml_dcscalar cadd(nml_dcscalar x1, nml_dcscalar x2);
nml_dcscalar csub(nml_dcscalar x1, nml_dcscalar x2);
nml_dcscalar cxdiv(nml_dcscalar x1, nml_dcscalar x2);
nml_dcscalar cconj(nml_dcscalar x);
nml_dscalar c_abs(nml_dcscalar x);
nml_dcscalar cxsqrt (nml_dcscalar x);
nml_dcscalar RCmul(nml_dscalar r, nml_dcscalar x);

#endif/*COMPLEX_H	*/

