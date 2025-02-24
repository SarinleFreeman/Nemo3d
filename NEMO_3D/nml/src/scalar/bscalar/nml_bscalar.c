
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
/* Include char  scalar type definitions.	*/
#include<nml_bscalar.h>

int
(nml_bfprintf)(FILE *stream, nml_bscalar x,
    int width, unsigned int precision, nml_fmtflags flags) {
  char format[64];
  int   characters = sprintf(format, "%%%d.%ue", width, precision);
  if (0 <= characters) {
    return fprintf(stream, format, (double)x); } else { return -1; }
  }

