#ifndef _nml_dcscalar_h
#define _nml_dcscalar_h 1

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
#include<nml_dscalar.h>

typedef struct {
  /* public: */
  nml_double	r, i;
  } nml_dcomplex;

typedef nml_dcomplex	nml_dcscalar;



typedef nml_dcscalar COMPLEX; 



#ifdef NML_INLINE
inline static
nml_dcomplex
(nml_dcmplx)(nml_double r, nml_double i) {
  nml_dcomplex	c;
  c.r = r; c.i = i;
  return c; }
inline static
nml_dscalar
(nml_dcreal)(nml_dcomplex c) {
  return c.r; }
inline static
nml_dscalar
(nml_dcimag)(nml_dcomplex c) {
  return c.i; }
#else /*NML_INLINE	*/
nml_dcomplex
(nml_dcmplx)(nml_double r, nml_double i);
nml_double
(nml_dcreal)(nml_dcomplex c);
nml_double
(nml_dcimag)(nml_dcomplex c);
#endif/*NML_INLINE	*/

#ifdef	NML_MACROS
#define nml_dcreal(c) ((c).r)
#define nml_dcimag(c) ((c).i)
#endif/*NML_MACROS	*/

#define NML_DCZERO nml_dcmplx((nml_dscalar)0, (nml_dscalar)0)
#define NML_DCEQ(x, y) (((x.r) == (y.r)) && ((x.i) == (y.i)))
#define NML_DCNE(x, y) (((x.r) != (y.r)) || ((x.i) != (y.i)))
#define NML_DCSET(c, r, i) do { (c).r = (r); (c).i = (i); }
#define NML_DCMUL(t, x, y) \
do { \
  (t).r = (x).r*(y).r - (x).i*(y).i; \
  (t).i = (x).r*(y).i + (x).i*(y).r; \
  } while (0)
#define NML_DCDIV(t, x, y) \
do { \
  (t).i = ((y).r*(y).r + (y).i*(y).i); \
  (t).r = ((x).r*(y).r + (x).i*(y).i)/(t).i; \
  (t).i = ((x).i*(y).r - (x).r*(y).i)/(t).i; \
  } while (0)
#define NML_DCADD(t, x, y) \
do { \
  (t).r = (x).r + (y).r; \
  (t).i = (x).i + (y).i; \
  } while (0)
#define NML_DCSUB(t, x, y) \
do { \
  (t).r = (x).r - (y).r; \
  (t).i = (x).i - (y).i; \
  } while (0)
#define NML_DCSQR(x) ((x).r*(x).r + (x).i*(x).i)
#define NML_DCABS(x) sqrt(NML_DCSQR(x))

int
(nml_dcfprintf)(FILE *stream, nml_dcscalar x,
    int width, unsigned int precision, nml_fmtflags flags);

#endif/*_nml_dcscalar_h */

