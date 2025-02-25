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
$Header: /repo/nemo3d/src/math/i3matrix.h,v 1.3 2004/10/27 19:43:10 gekco Exp $
*****************************************************************************/

#ifndef I3MATRIX_H
#define I3MATRIX_H 1

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <vector.h>

#include "math_vdim.h"
#include "nml_global.h"
#include "realtype.h"
#include "util_memory.h"

typedef int ***i3matrix;

i3matrix I3matrix(int n1, int n2, int n3);

void rm_i3matrix(i3matrix *t_ptr);

i3matrix check_i3matrix(i3matrix, int n1, int n2, int n3);

#endif
