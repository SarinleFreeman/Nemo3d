#ifndef _nml_fcscalar_h
#define _nml_fcscalar_h 1

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



/* Include real scalar type definitions.	*/
#include<nml_fscalar.h>

typedef struct {
  /* public: */
  nml_float	r, i;
  } nml_fcomplex;

typedef nml_fcomplex	nml_fcscalar;


typedef nml_fcscalar KOMPLEX; 




#ifdef NML_INLINE
inline static
nml_fcomplex
(nml_fcmplx)(nml_float r, nml_float i) {
  nml_fcomplex	c;
  c.r = r; c.i = i;
  return c; }
inline static
nml_fscalar
(nml_fcreal)(nml_fcomplex c) {
  return c.r; }
inline static
nml_fscalar
(nml_fcimag)(nml_fcomplex c) {
  return c.i; }
#else /*NML_INLINE	*/
nml_fcomplex
(nml_fcmplx)(nml_float r, nml_float i);
nml_float
(nml_fcreal)(nml_fcomplex c);
nml_float
(nml_fcimag)(nml_fcomplex c);
#endif/*NML_INLINE	*/

#ifdef	NML_MACROS
#define nml_fcreal(c) ((c).r)
#define nml_fcimag(c) ((c).i)
#endif/*NML_MACROS	*/

#define NML_FCZERO nml_fcmplx((nml_fscalar)0, (nml_fscalar)0)
#define NML_FCEQ(x, y) (((x.r) == (y.r)) && ((x.i) == (y.i)))
#define NML_FCNE(x, y) (((x.r) != (y.r)) || ((x.i) != (y.i)))
#define NML_FCSET(c, r, i) do { (c).r = (r); (c).i = (i); }
#define NML_FCMUL(t, x, y) \
do { \
  (t).r = (x).r*(y).r - (x).i*(y).i; \
  (t).i = (x).r*(y).i + (x).i*(y).r; \
  } while (0)
#define NML_FCDIV(t, x, y) \
do { \
  (t).i = ((y).r*(y).r + (y).i*(y).i); \
  (t).r = ((x).r*(y).r + (x).i*(y).i)/(t).i; \
  (t).i = ((x).i*(y).r - (x).r*(y).i)/(t).i; \
  } while (0)
#define NML_FCADD(t, x, y) \
do { \
  (t).r = (x).r + (y).r; \
  (t).i = (x).i + (y).i; \
  } while (0)
#define NML_FCSUB(t, x, y) \
do { \
  (t).r = (x).r - (y).r; \
  (t).i = (x).i - (y).i; \
  } while (0)
#define NML_FCSQR(x) ((x).r*(x).r + (x).i*(x).i)
#define NML_FCABS(x) sqrt(NML_FCSQR(x))

int
(nml_fcfprintf)(FILE *stream, nml_fcscalar x,
    int width, unsigned int precision, nml_fmtflags flags);

#endif/*_nml_fcscalar_h */

