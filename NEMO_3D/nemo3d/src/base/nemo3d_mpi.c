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
$Header: /repo/nemo3d/src/base/nemo3d_mpi.c,v 1.4 2003/10/08 16:16:53 hook Exp $ 
*****************************************************************************/

#include "nemo3d_mpi.h"

#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <fstream>


void nemo3d_mpi_init(int argc, char *argv[]) 
{
   int argc0=argc;

#if (defined MPI3d && !defined FAKE_MPI)
   int _is_mpi_initialized=0;
   // Initialize the MPI message passing interface for parallel execution.
   MPI_Initialized(&_is_mpi_initialized);
   if (_is_mpi_initialized) 
      printf("One process already initialized!\n"), fflush(stdout);
   printf("x"), fflush(stdout);
   MPI_Init(&argc, &argv);
#else /* not (defined MPI3d && !defined FAKE_MPI) */
   printf("Initializing Fake MPI.\n");
#endif /* (defined MPI3d && !defined FAKE_MPI) */
   MPI_Comm_size(MPI_COMM_WORLD, &mpi_n3d_numprocs);
   MPI_Comm_rank(MPI_COMM_WORLD, &mpi_n3d_id);

   if (argc0==argc && mpi_n3d_numprocs==1) {
      char *env=NULL;
      env = getenv("PWD");
      chdir(env);
   }

   printf("\nProcess [%i]/%i:  activated\n",mpi_n3d_id+1,mpi_n3d_numprocs);
   fflush(stdout);
   
   // Process any command-line arguments.
   argv_process(args, argc, argv);
   
   if (!argv_check())
      exit(1);
   
#define LAUNCH_DEBUGGER
#ifdef LAUNCH_DEBUGGER
   if (CmdOptions.debuglaunch) {
      
#ifdef Linux
      {
         char leg[100];
         char DISP[100];
         /* gethostname(HOST,100); */
         /* sprintf(DISP,"volcanoes"); */
         sprintf(DISP,CmdOptions.debuglaunch);
         sprintf(leg,"ddd %s %d -display %s:0 &",argv[0],(int) getpid(),DISP);
         printf(leg);
         fflush(stdout);
         /* sprintf(leg,"ddd  %d --auto --debugger /usr/bin/gdb &",(int) getpid()); */
         system(leg);
#ifndef NO_SYSTEM_CALL
         system("sleep 10");
#else /* NO_SYSTEM_CALL is defined */
         {
            int debug=0;
            
            for (;;){
               if (debug)
                  break;
            }
         }
#endif /* NO_SYSTEM_CALL */
      }
      
#endif /* Linux */
#ifdef IRIX64
      {
         char leg[100];
         sprintf(leg,"cvd -pid %d &",(int) getpid());
         system(leg);
         system("sleep 20");
      }
#endif /* IRIX64 */
   }
   
#endif /* LAUNCH_DEBUGGER */
   
   
#ifdef MPI_TIMING
   mpiTiming.Initialize();
   MPI_TIC(mpiTiming.total_0);
   MPI_TIC(mpiTiming.setup_0);
#endif /* MPI_TIMING */
}


void nemo3d_mpi_finalize() 
{
#ifdef MPI3d
   fflush(stdout);
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
#endif /* MPI3d */
}

 
