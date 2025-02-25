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
$Header: /repo/nemo3d/src/db_old/nemo_mb_models.h,v 1.1 2004/02/16 03:36:56
gekco Exp $
*****************************************************************************/

#ifndef NEMO_MB_MODELS_H
#define NEMO_MB_MODELS_H 1
#ifdef DATABASE_CODE_FRAGMENT_TO_BE_FIXED_LATER

#ifndef LIMIT_THEORY
#include "Boolean.h"
#include "mat_struct.h"
#include "mat_tokens.h"
#include "util_messages.h"
#include "util_stringops.h"

const int SideBands = 1;

Boolean isIndependentSingleBand(char *);

int setMBModelParam(
    char *band_model,
    int *param_n,            /* Number of associated multiband parameters. */
    char **model, int *bdim, /* Block dimension of the bulk Hamiltonian. */
    int *band_inc,           /* Number of incoming bands. */
    int *nbr,                /* Number of neighbors. */
    int *abdim,              /* Block dimension for the RGF algorithm. */
    int *spin_orbit,         /* spin orbit coupling */
    int *spin_degeneracy,
    int *sigma_singular /* singular tight binding boundary condition. */
);

semiMaterial adjustMBParamPointer(semiMaterial sm);

semiMaterial initMBSemiMaterial(semiMaterial sm);

#endif /* LIMIT_THEORY */
#endif /* DATABASE_CODE_FRAGMENT_TO_BE_FIXED_LATER  */

#endif /* NEMO_MB_MODELS_H */
