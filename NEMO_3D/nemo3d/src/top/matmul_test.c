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
$Header: /repo/nemo3d/src/top/matmul_test.c,v 1.4 2005/05/27 17:34:33 marek Exp
$
*****************************************************************************/

#include "../base/matmul.c"

int main(int argc, char *argv[]) {

  int dim = 20, i, j, offset = 0, n;
  cmatrix a = Cmatrix(dim, dim);
  cvectr x = Cvectr(dim);
  cvectr y = Cvectr(dim);

  for (n = 0; n < 100000; n++) {
    for (i = 0; i < dim; i++) {
      x[i].r = i;
      x[i].i = i + 100;
      for (j = 0; j < dim; j++) {
        a[i][j].r = (i + j) / dim;
        a[i][j].i = (i - j) / dim;
      }
    }

    for (i = 0; i < dim; i++) {
      y[i] = czero;
    }

    cmatmul_spds_20(y, a, x, offset, offset, dim);
  }

  for (i = 0; i < dim; i++) {
    printf("%d (%g,%g)\n", i, y[i].r, y[i].i);
  }

  return 0;
  ;
}
