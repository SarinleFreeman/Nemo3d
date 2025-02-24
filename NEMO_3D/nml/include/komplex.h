#ifndef	KOMPLEX_H
#define	KOMPLEX_H 1

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


/* typedef COMPLEX complex;			*/
typedef nml_fcscalar komplex;
/* extern const nml_dscalar czero;		*/
/* extern const nml_dscalar cone;		*/
const nml_fcscalar kzero = {0.0, 0.0};
const nml_fcscalar kone  = {1.0, 0.0};

nml_fcscalar Komplex(nml_fscalar r, nml_fscalar i );
nml_fcscalar kmul(nml_fcscalar x1, nml_fcscalar x2);
nml_fcscalar kadd(nml_fcscalar x1, nml_fcscalar x2);
nml_fcscalar ksub(nml_fcscalar x1, nml_fcscalar x2);
nml_fcscalar kxdiv(nml_fcscalar x1, nml_fcscalar x2);
nml_fcscalar kconj(nml_fcscalar x);
nml_fscalar k_abs(nml_fcscalar x);
nml_fcscalar kxsqrt (nml_fcscalar x);
nml_fcscalar RKmul(nml_fscalar r, nml_fcscalar x);

#endif/*KOMPLEX_H	*/

