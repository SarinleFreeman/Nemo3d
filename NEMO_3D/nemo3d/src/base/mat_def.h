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
$Header: /repo/nemo3d/src/base/mat_def.h,v 1.4 2006/10/17 14:10:14 hoonryu Exp $
*****************************************************************************/

#ifndef MAT_DEF_H
#define MAT_DEF_H

#include "io_utils.h"
#include "util_stringops.h"

/*! \file mat_def.h
 *  \brief definition of available materials.
 *  defines available materials in NEMO3D.
 */

#define EL_None 0
/*!< 0 = No material */
#define EL_B 1
/*!< 1 = Boron */
#define EL_C 2
/*!< 2 = Carbon */
#define EL_N 3
/*!< 3 = Nitride */
#define EL_O 4
/*!< 4 = Oxygen */
#define EL_F 5
/*!< 5 = Fluorine */
#define EL_Al 6
/*!< 6 = Aluminium */
#define EL_Si 7
/*!< 7 = Silicon */
#define EL_P 8
/*!< 8 = Phosporous */
#define EL_Zn 9
/*!< 9 = Zinc */
#define EL_Ga 10
/*!< 10 = Galium */
#define EL_Ge 11
/*!< 11 = Germanium */
#define EL_As 12
/*!< 12 = Arsenic */
#define EL_Cd 13
/*!< 13 = Cadmium */
#define EL_In 14
/*!< 14 = Indium */
#define EL_Sn 15
/*!< 15 = Stannum */
#define EL_Sb 16
/*!< 16 = Stibium */
#define EL_H 17
/*!< 17 = Hydrogen */
#define EL_S 18
/*!< 18 = Sulphur */
#define EL_Se 19
/*!< 19 = Selenium */
#define EL_EXPLICIT 20
/*!< 20 = Not quire sure about this */
#define EL_USERDEF 5
/*!< 5 = User-defined Material */
#define EL_TOTAL EL_EXPLICIT + EL_USERDEF
#endif
