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
$Header: /repo/nemo3d/src/math/tensor.h,v 1.2 2003/10/08 16:18:21 hook Exp $
*****************************************************************************/


#ifndef TENSOR_H
#define TENSOR_H 1


#include <vector.h>



const int NR_END = 0;


/*
 * Macro returns the lower bound of a tensor's first dimension.
 * It accesses the address of the pointer minus one.
 */
#define d_1l(p) *int_ptr((p-1))

/*
 * Macro returns the upper bound of a tensor's first dimension.
 * It accesses the address of the pointer minus two.
 */
#define d_1h(p) *int_ptr((p-2))



/*
 * Macro returns the lower bound of a tensor's first dimension.
 * It accesses the address of the pointer minus three.
 */
#define d_2l(p) *int_ptr((p-3))

/*
 * Macro returns the upper bound of a tensor's first dimension.
 * It accesses the address of the pointer minus four.
 */
#define d_2h(p) *int_ptr((p-4))



/*
 * Macro returns the lower bound of a tensor's first dimension.
 * It accesses the address of the pointer minus five.
 */
#define d_3l(p) *int_ptr((p-5))

/*
 * Macro returns the upper bound of a tensor's first dimension.
 * It accesses the address of the pointer minus six.
 */
#define d_3h(p) *int_ptr((p-6))


#endif

