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
$Header: /repo/nemo3d/src/math/ptinply3.h,v 1.3 2004/10/18 23:39:23 gekco Exp $ 
*****************************************************************************/

#ifndef PTINPLY3_H
#define PTINPLY3_H 1

#include <math.h>
#include <stdlib.h>

#include "geopoint.h"
#include "geopointvectorvector.h"
#include "constants_nemo.h"
#include "realtype.h"

#include <rmatrix.h>
#include <ivector.h>

#define GeoZeroVec(v) ((v).x = (v).y = (v).z = 0.0)
#define GeoMultVec(a,b,c) do {(c).x = a*(b).x; (c).y = a*(b).y;	(c).z = a*(b).z; } while (0)
#define Geo_Vet(a,b,c) do {(c).x = (b).x-(a).x; (c).y = (b).y-(a).y; (c).z = (b).z-(a).z;} while (0)

#define Geo_Add(a,b,c) {(c).x = (b).x+(a).x; (c).y = (b).y+(a).y; (c).z = (b).z+(a).z;}

real GeoDotProd ( geopoint *vec0, geopoint *vec1 );
void GeoCrossProd ( geopoint *in0, geopoint *in1, geopoint *out );
real GeoTripleProd ( geopoint *vec0, geopoint *vec1, geopoint *vec2 );
real GeoVecLen ( geopoint *vec );
int GeoPolyNormal ( int	n_verts, geopoint *verts, geopoint *n );
real geo_solid_angle ( int n_vert, geopoint *verts, geopoint *p );
int ptinpolyhedron( rmatrix v, ivectr nv, real x, real y, real z );
int ptinpolyhedron_gp( geopointvectrvectr gpvv, geopoint p );

#endif

