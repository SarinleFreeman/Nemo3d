/*****************************************************************************
The Jet Propulsion Laboratory (JPL) NanoElectronicMOdeling-3D package.
Copyright (C) 2002 California Institute of Technology (Caltech)

This application is free software, which you can redistribute and/or modify
under the terms of the GNU Lesser General Public License as published by the
Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this library; see the file COPYING. If not, write to the
Free Software Foundation, Inc.,
59 Temple Place, Suite 330,
Boston, MA  02111-1307  USA

For additional information, please contact
  Gerhard Klimeck (gekco@jpl.nasa.gov)
  Fabiano Oyafuso (fabiano@jpl.nasa.gov)

Written by:  Chris Bowen
             Gerhard Klimeck
             Fabiano Oyafuso
             Seungwon Lee
             Olga Lazarenkova
             Hook Hua

This product includes software developed by the Apache Software Foundation
(http://www.apache.org/).

*****************************************************************************
$Header: /repo/nemo3d/src/math/math_vdim.h,v 1.2 2003/10/08 16:18:04 hook Exp $
*****************************************************************************/

#ifndef MATH_VDIM_H__
#define MATH_VDIM_H__

#include <vector.h>

/*T_HEADER_SHORT{vdim1}
  Returns the first dimension of a matrix.
  Macro that returns the address of a pointer reduced by {\em one} pointer step
  size. This macro uses in turn the macro \verb|*int_ptr|
  REFERENCE{int_ptr,vector.h}. T*/
#define vdim1(p) *int_ptr((p - 1))

/*T_HEADER_SHORT{vdim2}
  Returns second dimension of a tensor.
  Macro that returns the address of a pointer reduced by {\em two} pointer step
  sizes. This macro uses in turn the macro \verb|*int_ptr|
  REFERENCE{int_ptr,vector.h}. T*/
#define vdim2(p) *int_ptr((p - 2))

/*T_HEADER_SHORT{vdim3}
  Returns the third dimension of a tensor.
  Macro that returns the address of a pointer reduced by {\em three} pointer
  step sizes. This macro uses in turn the macro \verb|*int_ptr|
  REFERENCE{int_ptr,vector.h}. T*/
#define vdim3(p) *int_ptr((p - 3))

/*T_HEADER_SHORT{vdim4}
  Returns the third dimension of a tensor.
  Macro that returns the address of a pointer reduced by {\em three} pointer
  step sizes. This macro uses in turn the macro \verb|*int_ptr|
  REFERENCE{int_ptr,vector.h}. T*/
#define vdim4(p) *int_ptr((p - 4))

#endif
