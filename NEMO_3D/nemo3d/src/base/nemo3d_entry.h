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
$Header: /repo/nemo3d/src/base/nemo3d_entry.h,v 1.6 2003/10/08 16:16:52 hook Exp $
*****************************************************************************/

#ifndef NEMO3d_ENTRY_H
#define NEMO3d_ENTRY_H 1

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

#include "nemo3d_mpi.h"

#include "i_mat.h"
#include "io_utils.h"

#ifdef MPI3d
#include "run3d_mpi.h"
#endif

#include "MPI_Timing.h"
#include "top_struct.H"
#include "phon_struct.h"
#include "phon_Keating.h"

#ifdef Embedded_Postprocess
#include "Postprocess.h"
#endif

#endif /* NEMO3D_DB */

int nemo3d(const char* inputFile, const char *inFileName);

int nemo3d_wrapper(const char* str,  char *filename);

#ifdef __cplusplus
extern "C" {
#endif
int nemo3d_wrapper_c(const char* str,  char *filename);
#ifdef __cplusplus
}
#endif


#endif /* NEMO3d_ENTRY_H */
