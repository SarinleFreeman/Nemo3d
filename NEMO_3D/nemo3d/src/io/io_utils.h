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
$Header: /repo/nemo3d/src/io/io_utils.h,v 1.12 2003/10/31 22:37:17 hook Exp $
*****************************************************************************/

#ifndef _IO_UTILS_H
#define _IO_UTILS_H

/*
 * for the embedding of NEMO 3D we need to have just the function definitions
 * at the top level defined.  Not the whole include path should be followed
 * down in order to make the embedding of the overall nemo3d library into
 * python simpler.
 */
#ifndef NEMO3D_DB

#include <iostream>
#include <iomanip>
#include <string>

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "run3d_mpi.h"

using std::string;

#ifdef PRINT_DEBUG_COM
#define print_com2(A,B) "printf(A,B), fflush(stdout);"
#define print_com3(A,B,C) "printf(A,B,C), fflush(stdout);"
#else
#define print_com2(A,B) ;
#define print_com3(A,B,C) ;
#endif


// global variable
#ifdef SET_GLOBAL_STREAMS
#define SCOPE_STREAM
#else
#define SCOPE_STREAM extern
#endif

class streambuf_N3D: public std::streambuf {
};

/* NOTE for Intel C++ compilers:
 * class ostream_N3D: private virtual streambuf_N3D, public std::ostream {
 * Though it will compile fine, the inheritance of a VIRTUAL specified class
 * will cause "Segmentation fault (core dumped)" during RUNTIME.
 * Removing the VIRTUAL will solve the problem.
 */
class ostream_N3D: private streambuf_N3D, public std::ostream {
    public:
        ostream_N3D() : std::ostream(this) {
        }
};

SCOPE_STREAM ostream_N3D cout_master, cout_slave;

//class ostream_debug: ostream_N3D
//{
//   std::string
// public:
//   ostream_debug() : ostream_N3D()  {    }
//
//   friend template<typename T>
//      ostream_debug operator<<(ostream_debug& o, const T& t);
//};
;//
//
//class string_debug : public std::string
//{
// public:
//   string_debug(const std::string& s) : std::string(s)  {  }
//   friend std::ostream& operator<<(std::ostream& o, const string_debug& s);
//};
;//
//
//class DEATH : public std::string
//{
// public:
//   DEATH(const std::string& s) : std::string(s)  {  }
//   friend std::ostream& operator<<(std::ostream& o, const DEATH& s);
//};


void die(char*, ...);
void masterPrint(char*, ...);

string format(const char* fmt, ...);

void setup_global_streams();
double* getBinFileData(char* filename, int prec, int *rdim, int *cdim);

#endif  /* NEMO3D_DB */

// The definition below establishes an unmangled C function compiled with C++
// such that this function can be called from a C-linked executable
#ifdef __cplusplus
extern "C" {
#endif
  void setup_global_streams_c(void);
#ifdef __cplusplus
}
#endif


#endif /* _IO_UTILS_H */
