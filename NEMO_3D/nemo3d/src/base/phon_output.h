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
$Header: /repo/nemo3d/src/base/phon_output.h,v 1.5 2004/02/16 03:36:56 gekco Exp $ 
*****************************************************************************/


/* some of print-outs for phon_struct*/

#ifndef PHON_OUTPUT_H
#define PHON_OUTPUT_H

#include <stdlib.h>
#include <nml_dcmatrix.h>
#include <nml_dmatrix.h>
#include <nml_imatrix.h>

#include <i2tensor.h>
#include <i3tensor.h>
#include <r3tensor.h>
#include <rvector.h>
#include <cvector.h>
#include <rmatrix.h>
#include <cmatrix.h>
#include <imatrix.h>
#include "realtype.h"
#include "ivectori3matrix.h"
   /* #include "mat_struct.h" */
#include "ivectorvector.h"
#include "cmatrixmatrix.h"
#include "rmatrixrvector.h"

#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
//#include <fstream.h>

#include "qd_struct.h"
#include "phon_parpack.h"
#include "phon_struct.h"

void print_phonon_dispersion(ostream& output_file,qd_struct d,
                             real V_over_V0//Delete me!!!!!!
                             );
void print_phonon_dispersion(string file_name,qd_struct d,real V_over_V0/*Delete me!!!*/);
//Prints phonon dispersion in file_name
void DM_print_matlab(qd_struct d,
                     int proc_to_print);
//Prints DM stored on the processor for Matlab
void Map_print(qd_struct d);//Prints map for DM
void nbr_on_proc_print(qd_struct d,
                       int cell0,int cellF);//print neighbor indexes
#endif
