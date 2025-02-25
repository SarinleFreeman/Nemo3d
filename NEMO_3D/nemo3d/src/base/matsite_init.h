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
$Header: /repo/nemo3d/src/base/matsite_init.h,v 1.4 2004/10/27 19:43:10 gekco
Exp $
*****************************************************************************/

#ifndef MATSITE_INIT_H
#define MATSITE_INIT_H

#include "Boolean.h"
#include "ivectorvector.h"
#include "nml_global.h"
#include "realtype.h"
#include "util_stringops.h"
#include <rvector.h>

#include "mat_def.h"
#include "qd_struct.h"
#include "random.h"
#include "run3d_mpi.h"

#include "Cell.h"

#include "i_3dout_files.h"
#include "i_mat.h"

typedef struct MatSite_Init_struct {
  int imat;
  ivectr cation, anion;
  real cation_weight, anion_weight;
  ivectrvectr cellatom;
  rvectr randweight;
  int random_alloy;
  ivectr species_count;
} *matsite_init_struct;

matsite_init_struct make_matsite_init_struct(qd_struct d, const Cell &cell,
                                             int matid3d);
void rm_matsite_init_struct(matsite_init_struct *mat_ptr);
void print_species_count(qd_struct d, matsite_init_struct mat);

#endif
