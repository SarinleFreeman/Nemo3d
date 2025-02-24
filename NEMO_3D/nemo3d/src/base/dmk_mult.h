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
$Header: /repo/nemo3d/src/base/dmk_mult.h,v 1.1 2004/02/27 04:09:05 lol Exp $
*****************************************************************************/

#ifndef DMK_MULT_H
#define DMK_MULT_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "realtype.h"
#include "MPI_Timing.h"


#include "qd_struct.h"
#include "dmk_mult.h"


void av(qd_struct d,
        int num_states,//number of states on the processor =3*N_atoms_on_proc
        cvectr v,//initial vector
        cvectr w//the rezult of vector-matrix multiplication w<==DM*v
        );
void av_mpi(qd_struct d,
            cvectr v,
            cvectr w);
     //Calls av to compute (1) residential w<==DM*v,
     //                    (2) wR, wL to be sent to right/left neighbors
     //Communicates with neighbors (1) send wR, wL
     //                            (2) receive vR, vL
     //Computes the final value of w
void qsortEph(cvectr Eph,//[nconv==iparam[4]] vector of converged eigenvalues
                    int beg, int end//lower and upper bounds of the arrays
                     );
     /*Modification of quick_sort to
      sort obtained phonon eigenvalues
      in ascending order of its real part.
     */
int find_min_ivec_mpi(qd_struct d,ivectr iv);


#endif /* DMK_H */
