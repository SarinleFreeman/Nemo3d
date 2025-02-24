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
$Header: /repo/nemo3d/src/math/c3tensor.h,v 1.4 2004/10/27 19:43:10 gekco Exp $
*****************************************************************************/

#ifndef C3TENSOR_H
#define C3TENSOR_H 1

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "util_messages.h"
#include "nml_global.h"
#include "tensor.h"
#include "complex.h"

typedef complex ***c3tensor;

c3tensor C3tensor(int dimRowLow, int dimRowHigh, int dimColLow, int dimColHigh, int dimDepthLow, int dimDepthHigh);
void rm_c3tensor(c3tensor *ptrC3t, int dimRowLow);
c3tensor check_c3tensor(c3tensor c3t, int d1l, int d1h, int d2l, int d2h, int d3l, int d3h);

#endif /* C3TENSOR_H */

