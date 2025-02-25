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
$Header: /repo/nemo3d/src/io/i_mat.h,v 1.5 2004/02/16 03:36:56 gekco Exp $
*****************************************************************************/

#ifndef I_MAT_H
#define I_MAT_H

#include "mat_def.h"
#include "mb_ham_spds_micro.h"
#include "qd_struct.h"
#include "run3d_mpi.h"
#include "rvector_lib.h"
/*#include "nemo_mb_models.h" */
#include "io_utils.h"
#include "util_messages.h"
#include "util_stringops.h"

int str_to_atomID(const char *);
char *atomID_to_str(int);
void set_sMatList(qd_struct);

#endif /* I_MAT_H */
