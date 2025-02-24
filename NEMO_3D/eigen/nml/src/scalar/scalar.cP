dnl This is a class prototype input file for .../<ns>/bin/genclass.
pushdef(`index', `ifelse($#, 0, ``index'', ``index($*)'')')dnl
pushdef(`shift', `ifelse($#, 0, ``shift'', ``shift($*)'')')dnl
pushdef(`format', `ifelse($#, 0, ``format'', ``format($*)'')')dnl

/*****************************************************************************
The <namespace>.
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


#undef <NS>INLINE
#undef <NS>MACROS
/* Include <precision> <system> scalar type definitions.	*/
#include<<ns><p><s>scalar.h>
ifelse(_system_, complex, `
<ns><p>complex
(<ns><p><s>mplx)(<ns><precision> r, <ns><precision> i) {
  <ns><p>complex	c;
  c.r = r; c.i = i;
  return c; }

<ns><precision>
(<ns><p><s>real)(<ns><p>complex c) {
  return c.r; }
#define <ns><p><s>real_(c) ((c).r)

<ns><precision>
(<ns><p><s>imag)(<ns><p>complex c) {
  return c.i; }
#define <ns><p><s>imag_(c) ((c).i)

int
(<ns><p><s>fprintf)(FILE *stream, <ns><p>complex c,
    int width, unsigned int precision, <ns>fmtflags flags) {
  int  total = 0;
  int
  characters = fprintf(stream, "(");
  if (0 <= characters) { total += characters; } else { return -1; }
  characters =
    <ns><p>fprintf(stream, <ns><p><s>real_(c), width, precision, flags);
  if (0 <= characters) { total += characters; } else { return -1; }
  characters = fprintf(stream, ",");
  if (0 <= characters) { total += characters; } else { return -1; }
  characters =
    <ns><p>fprintf(stream, <ns><p><s>imag_(c), width, precision, flags);
  if (0 <= characters) { total += characters; } else { return -1; }
  characters = fprintf(stream, ")");
  if (0 <= characters) { total += characters; } else { return -1; }
  return total;
  }
', `
int
(<ns><p><s>fprintf)(FILE *stream, <ns><p><s>scalar x,
    int width, unsigned int precision, <ns>fmtflags flags) {
  char format[64];
  int   characters = sprintf(format, "%%%d.%ue", width, precision);
  if (0 <= characters) {
    return fprintf(stream, format, (double)x); } else { return -1; }
  }
')dnl

