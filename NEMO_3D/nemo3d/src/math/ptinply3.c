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
$Header: /repo/nemo3d/src/math/ptinply3.c,v 1.4 2004/12/02 12:09:56 marek Exp $
*****************************************************************************/

#include "ptinply3.h"

#ifndef max
#define max(a, b) ((a) > (b) ? (a) : (b))
#endif

#ifndef min
#define min(a, b) ((a) < (b) ? (a) : (b))
#endif

real GeoDotProd(geopoint *vec0, geopoint *vec1) {
  return (vec0->x * vec1->x + vec0->y * vec1->y + vec0->z * vec1->z);
}

void GeoCrossProd(geopoint *in0, geopoint *in1, geopoint *out) {
  out->x = (in0->y * in1->z) - (in0->z * in1->y);
  out->y = (in0->z * in1->x) - (in0->x * in1->z);
  out->z = (in0->x * in1->y) - (in0->y * in1->x);
}

real GeoTripleProd(geopoint *vec0, geopoint *vec1, geopoint *vec2) {
  geopoint tmp;

  GeoCrossProd(vec0, vec1, &tmp);
  return (GeoDotProd(&tmp, vec2));
}

real GeoVecLen(geopoint *vec) { return sqrt(GeoDotProd(vec, vec)); }

int GeoPolyNormal(int n_verts, geopoint *verts, geopoint *n) {
  int i;
  float n_size;
  geopoint v0, v1, p;

  GeoZeroVec(*n);
  Geo_Vet(verts[0], verts[1], v0);
  for (i = 2; i < n_verts; i++) {
    Geo_Vet(verts[0], verts[i], v1);
    GeoCrossProd(&v0, &v1, &p);
    n->x += p.x;
    n->y += p.y;
    n->z += p.z;
    v0 = v1;
  }

  n_size = GeoVecLen(n);
  if (n_size > 0.0) {
    GeoMultVec(1 / n_size, *n, *n);
    return 1;
  } else
    return 0;
}

/*=========================  geo_solid_angle  =========================*/
/*
  Calculates the solid angle given by the spherical projection of
  a 3D plane polygon
*/

real geo_solid_angle(int n_vert,      /* number of vertices */
                     geopoint *verts, /* vertex coordinates list */
                     geopoint *p)     /* point to be tested */
{
  int i;
  real area = 0.0, ang, s, l1, l2;
  geopoint p1, p2, r1, a, b, n1, n2;
  geopoint plane;

  if (n_vert < 3)
    return 0.0;

  GeoPolyNormal(n_vert, verts, &plane);

  /*
     WARNING: at this point, a practical implementation should check
     whether p is too close to the polygon plane. If it is, then
     there are two possibilities:
       a) if the projection of p onto the plane is outside the
          polygon, then area zero should be returned;
       b) otherwise, p is on the polyhedron boundary.
  */

  p2 = verts[n_vert - 1]; /* last vertex */
  p1 = verts[0];          /* first vertex */
  Geo_Vet(p1, p2, a);     /* a = p2 - p1 */

  for (i = 0; i < n_vert; i++) {
    Geo_Vet(*p, p1, r1);
    p2 = verts[(i + 1) % n_vert];
    Geo_Vet(p1, p2, b);
    GeoCrossProd(&a, &r1, &n1);
    GeoCrossProd(&r1, &b, &n2);

    l1 = GeoVecLen(&n1);
    l2 = GeoVecLen(&n2);
    s = GeoDotProd(&n1, &n2) / (l1 * l2);
    ang = acos(max(-1.0, min(1.0, s)));
    s = GeoTripleProd(&b, &a, &plane);
    area += s > 0.0 ? Pi - ang : Pi + ang;

    GeoMultVec(-1.0, b, a);
    p1 = p2;
  }

  area -= Pi * (n_vert - 2);

  return (GeoDotProd(&plane, &r1) > 0.0) ? -area : area;
}

/* Function returns 1 if the point (x,y,z) is inside the polyhedron
   defined by the list of verticies in rmatrix v[i][j].  Otherwise
   zero is returned.  The first index of rmatrix labels vertices of
   the polyhedron.  The second index labels the  x,y and z coordinates
   (j=0 is x, j=1 is y, j=2 is z) of each vertex. */
int ptinpolyhedron(rmatrix v, ivectr nv, real x, real y, real z) {
  int i, j, cnt, inout, n = nml_iv_extent(nv);
  real area;
  geopoint verts[100], p;

  p.x = x;
  p.y = y;
  p.z = z;

  cnt = 0;
  area = 0.0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < nv[i]; j++) {
      verts[j].x = v[cnt][0];
      verts[j].y = v[cnt][1];
      verts[j].z = v[cnt][2];
      cnt++;
    }

    area += geo_solid_angle(nv[i], verts, &p);
  }

  if ((area > 2 * Pi) || (area < -2 * Pi))
    inout = 1;
  else
    inout = 0;

  return inout;
}

/* Function returns 1 if the point (x,y,z) is inside the polyhedron
   defined by the list of verticies in rmatrix v[i][j].  Otherwise
   zero is returned.  The first index of rmatrix labels vertices of
   the polyhedron.  The second index labels the  x,y and z coordinates
   (j=0 is x, j=1 is y, j=2 is z) of each vertex. */
int ptinpolyhedron_gp(geopointvectrvectr gpvv, geopoint p) {
  int i, inout, n = vdim(gpvv);
  real area;

  area = 0.0;
  for (i = 0; i < n; i++) {
    area += geo_solid_angle(vdim(gpvv[i]), gpvv[i], &p);
  }

  if ((area > 2 * Pi) || (area < -2 * Pi))
    inout = 1;
  else
    inout = 0;

  return inout;
}
