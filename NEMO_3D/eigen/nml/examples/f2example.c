/*****************************************************************************
Copyright (C) 2002-2004 California Institute of Technology (Caltech)

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

Written by Marek J. Korkusinski

*****************************************************************************/

#include<f2tensor.h>

int
main(int argc, char* argv[]) {
  const
  nml_index	lb2 = 2;
  const
  nml_index	ub2 = 9;
  const
  nml_index	lb1 = 3;
  const
  nml_index	ub1 = 7;

  nml_memory_report("Before matrix allocation.") ;

  f2tensor	M = F2tensor(lb2, ub2, lb1, ub1);

  nml_memory_report("After matrix allocation.") ;

  nml_index	i = lb2;
  for (i = lb2; i <= ub2; ++i) {
    nml_index	j = lb1;
    for (j = lb1; j <= ub1; ++j)
      M[i][j] = j + i*10;
    }

  fprintf(stdout, "M:\n");
  nml_f2D_fprintf(stdout, M + lb2, 8, 1, 0, 1 + ub1 - lb1);

  rm_f2tensor(&M, lb2);

  nml_memory_report("After matrix deallocation.") ;

  return 0;
  }
