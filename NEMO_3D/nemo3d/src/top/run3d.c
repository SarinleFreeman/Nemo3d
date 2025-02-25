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
$Header: /repo/nemo3d/src/top/run3d.c,v 1.22 2008/02/07 15:47:45 stevenclark Exp
$
*****************************************************************************/

#include "run3d.h"

//! This is the main function. The execution of NEMO3D starts here.
/*! The first step the program takes is checking its own name.
 *  If the name is "fmtdat", then the function fmtdat() is called.
 *  If it is "postprocess", then the function nemo3d_postprocessing()
 *  is called. Finally, if the name is "nemo3d", then the following
 *  sequence of activities takes place:
 * -# initialize MPI;
 * -# copy the contents of the entire input file to a string;
 * -# call the function nemo3d() which does all the calculations;
 * -# finalize MPI.
 */
int main(int argc, char *argv[]) {
  for (int i = 0; i < argc; i++) {
    if (strstr(argv[i], "fmtdat")) {
      fmtdat(argc - i, argv + i);
      nml_memory_report("End of fmtdat");
      exit(0);
    }
#ifdef Embedded_Postprocess
    else if (strstr(argv[i], "postprocess")) {
      nemo3d_postprocessing(argc, argv);
      nml_memory_report("End of postprocessing");
      exit(0);
    }
#endif
    else if (strstr(argv[i], "nemo3d")) {
      nemo3d_mpi_init(argc, argv);

      setup_global_streams();

      // copy file to a string
      fstream fs(CmdOptions.input, ios::in);
      string input;
      getline(fs, input, '\0');
      fs.close();

      nemo3d(CmdOptions.input, CmdOptions.input);
      nml_memory_report("End of nemo3d");

      //         nemo3d_wrapper(input.c_str(), CmdOptions.input);

      nemo3d_mpi_finalize();
      return 0;
    }
  }
  return 0;
}
