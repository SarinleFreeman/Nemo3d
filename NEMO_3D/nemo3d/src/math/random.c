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
$Header: /repo/nemo3d/src/math/random.c,v 1.5 2004/06/18 02:49:57 swlee Exp $
*****************************************************************************/

#include "random.h"
#include "run3d_mpi.h"
#include <iostream>
extern "C" {
#include <sys/time.h>
#include <sys/times.h>
#include <sys/types.h>
#include <unistd.h>
}

static int seed_set = 0;
int set_random_seed(int seed) {
  if (seed_set)
    return 0;

  srand(seed);
  seed_set = 1;

  return 0;
}

int reset_random_seed(int seed) {
  srand(seed);
  seed_set = 1;

  return 0;
}

int reset_random_seed_with_time(void) {
  int seed;
  if (mpi_n3d_id == 0) {
    struct timeval tvbuf;  // Values from call to gettimeofday
    struct timezone tzbuf; // Timezone
    gettimeofday(&tvbuf, &tzbuf);
    unsigned long last_secs = tvbuf.tv_sec;
    seed = last_secs % 16081;
  }
  MPI_Bcast(&seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
  srand(seed);
  return 0;
}
