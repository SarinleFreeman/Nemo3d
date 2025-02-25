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
$Header: /repo/nemo3d/src/math/rmatrixrvector.h,v 1.3 2004/10/27 19:43:10 gekco
Exp $
*****************************************************************************/

#ifndef Rmatrixrvector_H
#define Rmatrixrvector_H 1

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "io_utils.h"
#include "math_vdim.h"
#include "nml_global.h"
#include "realtype.h"
#include <rmatrix.h>
#include <rvector.h>

typedef rmatrix *rmatrixrvectr;

rmatrixrvectr Rmatrixrvectr(int d1, int d2);
void rm_rmatrixrvectr(rmatrixrvectr *t_ptr);
void rm_rmatrixrvectr_shell_only(rmatrixrvectr *t_ptr);
int add_elem_rmatrixrvectr(rvectr data_vector, int pos1, int pos2,
                           rmatrixrvectr hold_matrix);

#endif
