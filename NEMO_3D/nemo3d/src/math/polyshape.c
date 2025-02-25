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
$Header: /repo/nemo3d/src/math/polyshape.c,v 1.3 2004/10/27 19:43:10 gekco Exp $
*****************************************************************************/

#include "polyshape.h"

polyshape Polyshape_Empty(void) {

  polyshape p = (polyshape)nml_calloc(1, sizeof(polyshape));
  return p;
}

void rm_polyshape(polyshape *p_ptr) {

  polyshape p;
  if (!p_ptr || !(p = *p_ptr))
    return;

  rm_imatrix(&p->rv);
  rm_imatrix(&p->rf);
  rm_ivectr(&p->nvert);
  rm_rmatrix(&p->vert);

  nml_free(p);
  *p_ptr = NULL;

  return;
}

polyshape Polyshape_Rhombohedron(real delta, real ah) {
  imatrix rv, rf;
  ivectr nvert;
  rmatrix vert;
  int i, j, k;

  polyshape p = Polyshape_Empty();

  p->rv = rv = Imatrix(8, 3);
  p->rf = rf = Imatrix(6, 4);
  p->nvert = nvert = Ivectr(6);
  p->vert = vert = Rmatrix(24, 3);

  /* Store the integer vertices of a rhombohedron */
  rv[0][0] = 2;
  rv[0][1] = 2;
  rv[0][2] = 2;

  rv[1][0] = 2;
  rv[1][1] = 1;
  rv[1][2] = 1;

  rv[2][0] = 1;
  rv[2][1] = 1;
  rv[2][2] = 2;

  rv[3][0] = 1;
  rv[3][1] = 2;
  rv[3][2] = 1;

  rv[4][0] = 1;
  rv[4][1] = 1;
  rv[4][2] = 0;

  rv[5][0] = 1;
  rv[5][1] = 0;
  rv[5][2] = 1;

  rv[6][0] = 0;
  rv[6][1] = 1;
  rv[6][2] = 1;

  rv[7][0] = 0;
  rv[7][1] = 0;
  rv[7][2] = 0;

  /* Store the 6 faces of the rhombohedron in terms of the vertices.  The
     face vertices are ordered such that the cross product of any two face
     edge vectors is directed outside of the rhombohedron.  This ordering must
     be consistent in order for the function ptinpolyhedron() to work properly.
   */
  rf[0][0] = 0;
  rf[0][1] = 2;
  rf[0][2] = 5;
  rf[0][3] = 1;
  rf[1][0] = 0;
  rf[1][1] = 3;
  rf[1][2] = 6;
  rf[1][3] = 2;
  rf[2][0] = 0;
  rf[2][1] = 1;
  rf[2][2] = 4;
  rf[2][3] = 3;
  rf[3][0] = 7;
  rf[3][1] = 4;
  rf[3][2] = 1;
  rf[3][3] = 5;
  rf[4][0] = 7;
  rf[4][1] = 5;
  rf[4][2] = 2;
  rf[4][3] = 6;
  rf[5][0] = 7;
  rf[5][1] = 6;
  rf[5][2] = 3;
  rf[5][3] = 4;

  for (i = 0; i < 6; i++) {
    nvert[i] = 4;
    for (j = 0; j < nvert[i]; j++)
      for (k = 0; k < 3; k++)
        vert[i * 4 + j][k] = ah * rv[rf[i][j]][k] + delta;
  }

  return p;
}

polyshape Polyshape_Cube(real ax, real ay, real az, real delta_x, real delta_y,
                         real delta_z) {
  imatrix rv, rf;
  ivectr nvert;
  rmatrix vert;
  int i, j;

  polyshape p = Polyshape_Empty();

  p->rv = rv = Imatrix(8, 3);
  p->rf = rf = Imatrix(6, 4);
  p->nvert = nvert = Ivectr(6);
  p->vert = vert = Rmatrix(24, 3);

  /* Store the integer vertices of a cube */
  rv[0][0] = 1;
  rv[0][1] = 1;
  rv[0][2] = 1;

  rv[1][0] = 1;
  rv[1][1] = 1;
  rv[1][2] = 0;

  rv[2][0] = 0;
  rv[2][1] = 1;
  rv[2][2] = 0;

  rv[3][0] = 0;
  rv[3][1] = 1;
  rv[3][2] = 1;

  rv[4][0] = 1;
  rv[4][1] = 0;
  rv[4][2] = 1;

  rv[5][0] = 1;
  rv[5][1] = 0;
  rv[5][2] = 0;

  rv[6][0] = 0;
  rv[6][1] = 0;
  rv[6][2] = 0;

  rv[7][0] = 0;
  rv[7][1] = 0;
  rv[7][2] = 1;

  /* Store the 6 faces of the cube in terms of the vertices.  The
     face vertices are ordered such that the cross product of any two face
     edge vectors is directed outside of the rhombohedron.  This ordering must
     be consistent in order for the function ptinpolyhedron() to work properly.
   */
  rf[0][0] = 0;
  rf[0][1] = 1;
  rf[0][2] = 2;
  rf[0][3] = 3;

  rf[1][0] = 0;
  rf[1][1] = 4;
  rf[1][2] = 5;
  rf[1][3] = 1;

  rf[2][0] = 4;
  rf[2][1] = 7;
  rf[2][2] = 6;
  rf[2][3] = 5;

  rf[3][0] = 3;
  rf[3][1] = 2;
  rf[3][2] = 6;
  rf[3][3] = 7;

  rf[4][0] = 0;
  rf[4][1] = 3;
  rf[4][2] = 7;
  rf[4][3] = 4;

  rf[5][0] = 1;
  rf[5][1] = 5;
  rf[5][2] = 6;
  rf[5][3] = 2;

  for (i = 0; i < 6; i++) {
    nvert[i] = 4;
    for (j = 0; j < nvert[i]; j++) {
      vert[i * 4 + j][0] = ax * rv[rf[i][j]][0] + delta_x;
      vert[i * 4 + j][1] = ay * rv[rf[i][j]][1] + delta_y;
      vert[i * 4 + j][2] = az * rv[rf[i][j]][2] + delta_z;
    }
  }

  return p;
}

geopointvectrvectr Polyshape_Cube_gp(real ax, real ay, real az, real delta_x,
                                     real delta_y, real delta_z) {
  int i, j;
  geopointvectrvectr gpvv = Geopointvectrvectr(6);
  geopointvectr gpv_temp = NULL;

  imatrix rv = Imatrix(8, 3);
  imatrix rf = Imatrix(6, 4);

  /* Store the integer vertices of a cube */
  rv[0][0] = 1;
  rv[0][1] = 1;
  rv[0][2] = 1;

  rv[1][0] = 1;
  rv[1][1] = 1;
  rv[1][2] = 0;

  rv[2][0] = 0;
  rv[2][1] = 1;
  rv[2][2] = 0;

  rv[3][0] = 0;
  rv[3][1] = 1;
  rv[3][2] = 1;

  rv[4][0] = 1;
  rv[4][1] = 0;
  rv[4][2] = 1;

  rv[5][0] = 1;
  rv[5][1] = 0;
  rv[5][2] = 0;

  rv[6][0] = 0;
  rv[6][1] = 0;
  rv[6][2] = 0;

  rv[7][0] = 0;
  rv[7][1] = 0;
  rv[7][2] = 1;

  /* Store the 6 faces of the cube in terms of the vertices.  The
     face vertices are ordered such that the cross product of any two face
     edge vectors is directed outside of the rhombohedron.  This ordering must
     be consistent in order for the function ptinpolyhedron() to work properly.
   */
  rf[0][0] = 0;
  rf[0][1] = 1;
  rf[0][2] = 2;
  rf[0][3] = 3;

  rf[1][0] = 0;
  rf[1][1] = 4;
  rf[1][2] = 5;
  rf[1][3] = 1;

  rf[2][0] = 4;
  rf[2][1] = 7;
  rf[2][2] = 6;
  rf[2][3] = 5;

  rf[3][0] = 3;
  rf[3][1] = 2;
  rf[3][2] = 6;
  rf[3][3] = 7;

  rf[4][0] = 0;
  rf[4][1] = 3;
  rf[4][2] = 7;
  rf[4][3] = 4;

  rf[5][0] = 1;
  rf[5][1] = 5;
  rf[5][2] = 6;
  rf[5][3] = 2;

  for (i = 0; i < 6; i++) {
    gpv_temp = Geopointvectr(4);
    for (j = 0; j < 4; j++) {
      gpv_temp[j].x = ax * rv[rf[i][j]][0] + delta_x;
      gpv_temp[j].y = ay * rv[rf[i][j]][1] + delta_y;
      gpv_temp[j].z = az * rv[rf[i][j]][2] + delta_z;
    }
    add_elem_geopointvectrvectr(gpv_temp, i, gpvv);
  }

  rm_imatrix(&rf);
  rm_imatrix(&rv);

  return gpvv;
}
