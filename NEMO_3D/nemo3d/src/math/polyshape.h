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
$Header: /repo/nemo3d/src/math/polyshape.h,v 1.2 2003/10/08 16:18:06 hook Exp $ 
*****************************************************************************/


#ifndef POLYSHAPE_H
#define POLYSHAPE_H

#include "util_memory.h"
#include "geopointvectorvector.h"
#include "realtype.h"

#include <rvector.h>
#include <ivector.h>
#include <rmatrix.h>
#include <imatrix.h>


typedef struct POLYSHAPESTRUCT {
    imatrix rv;
    imatrix rf;
    ivectr nvert;
    rmatrix vert;
} *polyshape;

void rm_polyshape( polyshape *p_ptr );
polyshape Polyshape_Empty(void);

polyshape Polyshape_Rhombohedron(real delta, real ah);
polyshape Polyshape_Cube(real ax, real ay, real az, 
			 real delta_x, real delta_y, real delta_z);

geopointvectrvectr Polyshape_Cube_gp(real ax, real ay, real az, 
				     real delta_x, real delta_y, real delta_z);

#endif

