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
$Header: /repo/nemo3d/src/io/io_utils.c,v 1.13 2003/11/11 19:23:45 hook Exp $
*****************************************************************************/
#define SET_GLOBAL_STREAMS

#include "io_utils.h"

using std::cout;
using std::ios;

void die(char *s, ...) {
  char dum[1000];
  va_list argptr;

  va_start(argptr, s);
  vsprintf(dum, s, argptr);
  va_end(argptr);

  if (dum[0] != '\0')
    fprintf(stderr, "ERROR:  %s\n", dum);
  fflush(stderr);
  MPI_Finalize();
  exit(1);
}

void masterPrint(char *s, ...) {
  if (mpi_n3d_id != mpi_n3d_masterid)
    return;

  char dum[1000];
  va_list argptr;

  va_start(argptr, s);
  vsprintf(dum, s, argptr);
  va_end(argptr);

  fprintf(stdout, "%s\n", dum);
  fflush(stdout);
}

string format(const char *fmt, ...) {
  char dum[256];

  va_list args;
  va_start(args, fmt);
  vsnprintf(dum, 256, fmt, args);
  va_end(args);

  return string(dum);
}

void setup_global_streams() {

  if (mpi_n3d_id) {
    cout_master.clear(ios::badbit);
    cout_slave.rdbuf(cout.rdbuf());
  } else {
    cout_slave.clear(ios::badbit);
    cout_master.rdbuf(cout.rdbuf());

    // filebuf* fbp = new filebuf(1);
    // filebuf fb;
    // fb.open ("test.txt",ios::out);
    // ostream os(&fb);
    // os << "Test sentence\n";
    // fb.close();
  }
}

// std::ostream& operator<<(std::ostream& o, const string_debug& s)
//{
//    o << "DEBUG:  " << string(s);
//    return o;
// }
//
//
// std::ostream& operator<<(std::ostream& o, const DEATH& s)
//{
//    o << "ERROR:  " << string(s) << endl;
//    MPI_Finalize();
//    exit(1);
//
//    return o;
// }
//
// template<typename T>
// ostream_debug operator<<(ostream_debug& o, const T& t)
//{
//    o << "DEBUG:  " << o << endl;
// }

// The definition below establishes an unmangled C function compiled with C++
// such that this function can be called from a C-linked executable
#ifdef __cplusplus
extern "C" {
#endif
void setup_global_streams_c(void) { return; }
#ifdef __cplusplus
}
#endif
