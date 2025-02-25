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
$Header: /repo/nemo3d/src/base/qd_struct_global.c,v 1.2 2003/10/08 16:17:11 hook
Exp $
*****************************************************************************/

#include "qd_struct_global.h"

/*T_HEADER_SHORT{I_N3D_global}
  Basic input structure needed in the input deck handling.  T*/
qd_struct QD_global = NULL;

/*T_HEADER_SHORT{Get_QD_Global}
  Return the pointer to the overall simulation structure
  of the type \verb|qd_struct|. T*/
qd_struct Get_QD_Global(void) { return QD_global; }

/*T_HEADER_SHORT{Set_QD_Global}
  Set the pointer to the overall simulation structure
  of the type \verb|qd_struct| to the given input structure.  T*/
void Set_QD_Global(qd_struct QD) {
  QD_global = QD;

  return;
}
