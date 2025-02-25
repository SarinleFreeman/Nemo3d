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
$Header: /repo/nemo3d/src/postprocessing/TimeReversal.cpp,v 1.1 2004/08/26
21:28:07 swlee Exp $
*****************************************************************************/

#include "TimeReversal.h"
#include "SimpleParser.h"

void TimeReversal::get_data(complex *coef) { _cf = coef; }

void TimeReversal::apply_time_reversal_operator(char *evec_file) {

  for (int g = g_min; g < g_max; g++) {
    for (int orbital = 0; orbital < o_max; orbital++) {
      cf_tr(g, 1, orbital).r = cf(g, 0, orbital).r;
      cf_tr(g, 1, orbital).i = -cf(g, 0, orbital).i;
      cf_tr(g, 0, orbital).r = -cf(g, 1, orbital).r;
      cf_tr(g, 0, orbital).i = cf(g, 1, orbital).i;
    }
  }

  char filename[123];
  sprintf(filename, "%s_tr", evec_file);
  cout_master << "Writing time-reversal wave function to " << filename << endl;
  writeSiteInfo(d, "basis", filename, "{real wf.r, real wf.i}", _cf_tr);
}
