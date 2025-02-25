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
$Header: /repo/nemo3d/src/base/MaterialHandle.h,v 1.4 2003/10/08 16:16:19 hook
Exp $
*****************************************************************************/

#ifndef MATERIAL_HNDL_H
#define MATERIAL_HNDL_H

#include "Material_struct.H"
#include "realtype.h"

class MaterialHandle {
public:
  Material_struct *hndl;
  int matid3d;
  vector<int> anion_id, cation_id;
  string anion, cation;
  real cation_weight, anion_weight;

  MaterialHandle(Material_struct *h = NULL) { hndl = h; }

  friend ostream &operator<<(ostream &o, const MaterialHandle &mh);
};

#endif
