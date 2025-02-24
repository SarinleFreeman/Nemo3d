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
$Header: /repo/nemo3d/src/base/mat3d.h,v 1.3 2004/02/16 03:36:56 gekco Exp $ 
*****************************************************************************/

#ifndef MAT3D_H
#define MAT3D_H 1

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "io_utils.h"
#include "util_messages.h"

#include "realtype.h"
#include <rvector.h>
#include "rvector_lib.h"
#include "rmatrixrvector.h"

#include "constants_nemo.h"
#include "qd_struct.h"
#include "run3d_mpi.h"
#include "mat_def.h"
#include "i_mat.h"
#include "mb_ham_spds_micro.h"
   /* #include "mat_tokens.h" */

#ifdef OK
int unit_cell(I_Dev_CrystalStruct structure, int site, int anion, int cation, int cation1, real ratio1, int cation2, real ratio2);
#endif

void strnmat_init(qd_struct d, const vector<double>& mb_strain , int cation, int anion);
void param_init( qd_struct d );
int mat_param_init ( qd_struct d );

int trivial_mat_init( qd_struct d );

#endif







