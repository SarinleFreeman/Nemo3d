/*****************************************************************************
Copyright (C) 2002 California Institute of Technology (Caltech)

This file is part of
The NanoElectronic MOdeling (NEMO) Math Library.

This library is free software which you can redistribute and/or modify
under the terms of the GNU Library General Public License
as published by the Free Software Foundation;
either version 2, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this library; see the file COPYING. If not, write to the
Free Software Foundation, Inc.,
59 Temple Place, Suite 330,
Boston, MA  02111-1307  USA

Written by E. Robert Tisdale

*****************************************************************************
$Header: /repo/nml/examples/example.c,v 1.2 2003/10/08 16:19:07 hook Exp $
*****************************************************************************/

#include<r3tensor.h>

int
main(int argc, char* argv[]) {
  const
  nml_index	lb3 = -2;
  const
  nml_index	ub3 = +3;
  const
  nml_index	lb2 = +2;
  const
  nml_index	ub2 = +9;
  const
  nml_index	lb1 = +3;
  const
  nml_index	ub1 = +7;
  r3tensor	T = R3tensor(lb3, ub3, lb2, ub2, lb1, ub1);

  nml_index	h = lb3;

  for (h = lb3; h <= ub3; ++h) {
    nml_index	i = lb2;
    for (i = lb2; i <= ub2; ++i) {
      nml_index	j = lb1;
      for (j = lb1; j <= ub1; ++j) {
	T[h][i][j] = j + i*10 + h*100;
	}
      }
    }

  fprintf(stdout, "T:\n");
  nml_d3D_fprintf(stdout, T + lb3, 10, 3, 0, 1 + ub1 - lb1);

  rm_r3tensor(&T, lb3);

  return 0;
  }
