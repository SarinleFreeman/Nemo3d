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
$Header: /repo/nml/examples/rmexample.c,v 1.2 2003/10/08 16:19:09 hook Exp $
*****************************************************************************/

#include<rmatrix.h>

int main(int argc, char* argv[]) {
    const nml_extent m = 4;
    const nml_extent n = 7;
    rmatrix M = Rmatrix(m, n);

    nml_offset i, j;
    for (i = 0; i < m; ++i) {
        for (j = 0; j < n; ++j) {
        M[i][j] = j + i*10;
        }
    }

    fprintf(stdout, "M:\n");
    nml_dm_fprintf(stdout, M, 10, 3, 0, n);

    rm_rmatrix(&M);

    return 0;
}
