
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
$Header: /repo/nml/src/scalar/scalar.cP,v 1.2 2003/10/08 16:19:47 hook Exp $
*****************************************************************************/


#undef NML_INLINE
#undef NML_MACROS
/* Include float complex scalar type definitions.	*/
#include<nml_fcscalar.h>

nml_fcomplex
(nml_fcmplx)(nml_float r, nml_float i) {
  nml_fcomplex	c;
  c.r = r; c.i = i;
  return c; }

nml_float
(nml_fcreal)(nml_fcomplex c) {
  return c.r; }
#define nml_fcreal_(c) ((c).r)

nml_float
(nml_fcimag)(nml_fcomplex c) {
  return c.i; }
#define nml_fcimag_(c) ((c).i)

int
(nml_fcfprintf)(FILE *stream, nml_fcomplex c,
    int width, unsigned int precision, nml_fmtflags flags) {
  int  total = 0;
  int
  characters = fprintf(stream, "(");
  if (0 <= characters) { total += characters; } else { return -1; }
  characters =
    nml_ffprintf(stream, nml_fcreal_(c), width, precision, flags);
  if (0 <= characters) { total += characters; } else { return -1; }
  characters = fprintf(stream, ",");
  if (0 <= characters) { total += characters; } else { return -1; }
  characters =
    nml_ffprintf(stream, nml_fcimag_(c), width, precision, flags);
  if (0 <= characters) { total += characters; } else { return -1; }
  characters = fprintf(stream, ")");
  if (0 <= characters) { total += characters; } else { return -1; }
  return total;
  }

