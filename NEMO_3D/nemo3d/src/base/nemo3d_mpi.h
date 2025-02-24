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
$Header: /repo/nemo3d/src/base/nemo3d_mpi.h,v 1.4 2003/10/08 16:16:54 hook Exp $
*****************************************************************************/

#ifndef NEMO3D_MPI_H
#define NEMO3d_MPI_H 1


/*
 * for the embedding of NEMO 3D we need to have just the function definitions
 * at the top level defined.  Not the whole include path should be followed
 * down in order to make the embedding of the overall nemo3d library into
 * python simpler.
 */
#ifndef NEMO3D_DB

#include "eigsys3d.h"
#include "eigsys3d_par.h"
#include "math3d.h"
#include "ham3d.h"
#include "opt3d.h"
#include "mat3d.h"
#include "ham3d_final.h"
#include "util_command_BX.h"
#include "qd_struct.h"

#include "i_mat.h"
#include "io_utils.h"

#include "MPI_Timing.h"

#endif  /* NEMO3D_DB */

void nemo3d_mpi_init(int argc, char *argv[]);
void nemo3d_mpi_finalize(void);


#endif /* NEMO3D_MPI_H */
