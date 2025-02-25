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
$Header: /repo/nemo3d/src/top/serv3d.c,v 1.5 2005/05/27 17:34:33 marek Exp $
*****************************************************************************/

#include "serv3d.h"

int Server3D(qd_struct d) {
  int run_loop = 1;
  char *line = NULL;

  while (run_loop) {
    line = get_line(stdin);
    if (!strlcmp(line, "STOPserver")) {
      break;
    } else if (!strlcmp(line, "DOsim_io")) {
      d->I_N3D->sptr = IO_Tcl_Server(d->I_N3D->sptr, d->I_N3D_sptr_Store);
    }
    /*    else if (!strlcmp(line,"DOstruct_io")){
          d->I_N3D->sptr = IO_Tcl_Server(d->I_Dev->sptr, d->I_Dev_sptr_Store);
          }
    */
    else if (!strlcmp(line, "SAVEinput")) {
      i_save_shell(d, d->inputfile, 1, &d->Abort_result);
    } else {
      printf("Received an unknown command.\n"
             "Known Commands:\n"
             "STOPserver, DOsim_io, SAVEinput,\n");
    }
  }

  return 0;
}

int main(int argc, char *argv[]) {
  char *inputfile;
  int Abort_result = 0;
  int i, num;
  real dx, min, max, *target;
  qd_struct d = NULL;
  StructDescript_ptr sptr_store = NULL;

  /*T Process any command-line arguments. T*/
  /*T_FUNCTION{argv_process,argv_B.c} T*/
  argv_process(args, argc, argv);

  if (!argv_check())
    exit(1);

  /*T Create memory for the central simulation structure. T*/
  d = new QD_struct;
  /*T Set the global pointer to the central simulation structure. T*/
  Set_QD_Global(d);

  /* Assign the input filename from the command line options. T*/
  d->inputfile = copy_str(CmdOptions.input);

  dx = .5;

  i_read(d->inputfile, d, &Abort_result);

#ifdef OLDSERVER
  d->I_N3D->sptr = IO_Tcl_Server(d->I_N3D->sptr, d->I_N3D_sptr_Store);

  i_save_shell(d, d->inputfile, 1, &Abort_result);
#endif

  Server3D(d);

  delete d;

  return 0;
}
