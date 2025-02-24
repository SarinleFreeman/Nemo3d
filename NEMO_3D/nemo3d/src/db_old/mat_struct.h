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
$Header: /repo/nemo3d/src/db_old/mat_struct.h,v 1.1 2004/02/16 03:36:56 gekco Exp $
*****************************************************************************/

#ifndef MAT_STRUCT_H
#define MAT_STRUCT_H 1
#ifdef DATABASE_CODE_FRAGMENT_TO_BE_FIXED_LATER

#include <r2tensor.h>
#include <rvector.h>
#include <ivector.h>
#include <cmatrix.h>
#include "realtype.h"


const int MB_Param_Max  = 38;
const int B2_Param_Max  =  3;
const int MSB_Param_Max =  6;


typedef struct SEMI_MATERIAL {

    char   *material;
    char   *band_model;
    char   *comp1_mat;
    char   *comp2_mat;

    real    xval, yval;
    real    tkel;
    real    Ev_offset;
    real    Ec_offset;

    int     matid3d;
    char   *cation, *anion;
    ivectr  cation_id, anion_id;

    real   *mb_p;  
    int     mb_new_3d;
    cmatrix mb_d, mb_s;
    rvectr  mb_param;
    int     mb_bdim, mb_abdim, mb_spin_orbit, mb_band_inc, mb_nbr, mb_param_n, mb_spin_degeneracy, mb_sigma_singular;
    char   *mb_model;
    real    mb_unit_cell;
    real    unstrnd_cubic_cell_length;
    real    strain_alpha, strain_beta;
    real    mbn_p[MB_Param_Max];
    real    mb2_p[B2_Param_Max]; 
    rvectr  mb_pv;
    rvectr  mb_strain;

    char   *msb_label[MSB_Param_Max];

    int    *Abort_result;

} *semiMaterial;

#endif /* #ifdef DATABASE_CODE_FRAGMENT_TO_BE_FIXED_LATER  */

#endif /* MAT_STRUCT_H */


