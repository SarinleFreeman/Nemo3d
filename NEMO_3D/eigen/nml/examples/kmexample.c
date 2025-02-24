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

#include<kmatrix.h>

int main(int argc, char* argv[]) {
    const
    nml_extent	 m = 7;
    const
    nml_extent	 n = 4;

    nml_memory_report("before") ;

    kmatrix	 M = Kmatrix(m, n);

    nml_memory_report("after") ;

    nml_offset	 i, j;
    for (i = 0; i < m; ++i) {
        for (j = 0; j < n; ++j) {
        M[i][j] = nml_fcmplx((j + i*10), (i + j*10));
        }
    }

    fprintf(stdout, "M:\n");
    nml_fcm_fprintf(stdout, M, 8, 1, 0, n);

    rm_kmatrix(&M);

    nml_memory_report("the end") ;

    return 0;
}
