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
$Header: /repo/nemo3d/src/math/r3tensor.h,v 1.3 2004/10/27 19:43:10 gekco Exp $
*****************************************************************************/

#ifndef R3TENSOR_H
#define R3TENSOR_H 1

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "nml_global.h"
#include "realtype.h"
#include "tensor.h"
#include "util_memory.h"
#include "util_messages.h"

typedef real ***r3tensor;

r3tensor R3tensor(int dimRowLow, int dimRowHigh, int dimColLow, int dimColHigh,
                  int dimDepthLow, int dimDepthHigh);
void rm_r3tensor(r3tensor *ptrR3t, int dimRowLow);
r3tensor check_r3tensor(r3tensor r3t, int d1l, int d1h, int d2l, int d2h,
                        int d3l, int d3h);

#endif /* R3TENSOR_H */
