#ifndef _nml_bscalar_h
#define _nml_bscalar_h 1

/*****************************************************************************
The NEMO Math Library.
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
$Header: /repo/nml/src/scalar/scalar.hP,v 1.4 2004/10/22 13:37:23 marek Exp $
*****************************************************************************/



/* Include global type definitions.	*/
#include<nml_global.h>

typedef nml_char	nml_bscalar;

#define NML_bZERO ((nml_bscalar)0)
#define NML_bEQ(x, y) ((x) == (y))
#define NML_bNE(x, y) ((x) != (y))
#define NML_bMUL(t, x, y) ((t) = (x)*(y))
#define NML_bDIV(t, x, y) ((t) = (x)/(y))
#define NML_bADD(t, x, y) ((t) = (x) + (y))
#define NML_bSUB(t, x, y) ((t) = (x) - (y))
#define NML_bSQR(x) ((x)*(x))
#define NML_bABS(x)  abs(x)

int
(nml_bfprintf)(FILE *stream, nml_bscalar x,
    int width, unsigned int precision, nml_fmtflags flags);

#endif/*_nml_bscalar_h */

