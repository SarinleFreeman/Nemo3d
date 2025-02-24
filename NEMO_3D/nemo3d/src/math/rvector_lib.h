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
$Header: /repo/nemo3d/src/math/rvector_lib.h,v 1.3 2004/10/27 19:43:10 gekco Exp $
*****************************************************************************/

#ifndef RVECTOR_LIB_H
#define RVECTOR_LIB_H 1

#include <rvector.h>
#include "realtype.h"
#include "nml_global.h"
#include "util_messages.h"


/* Used by NEMO3D/nemo/math/rmatrix_lib.[ch]'s luDecompositionBackSubstitution() function */
typedef rvectr *roperator;
rvectr cp_rvectr (rvectr rvTarget, rvectr rvSource);
real* cp_rvectr2rarray(real *raTarget, rvectr rvSource);



/* Used by
 *   "NEMO3D/nemo/base/mat3d.c"
 *   "NEMO3D/nemo/io/i_mat."
 */
void cp_rarray2rarray(real* raSource, real* raTarget, int arraySize);


#endif /* RVECTOR_LIB_H */
