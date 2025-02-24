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
$Header: /repo/nemo3d/src/util/system.h,v 1.3 2005/04/15 19:09:36 dseaman Exp $
*****************************************************************************/


#ifndef SYSTEM_H
#define SYSTEM_H 1

#include <stdio.h>


#define LOGFILE stderr 
#define CUTOFF_PREC 1.0e-8


#if (__STDC__ | WINNT | __VMS | VMS | Linux | NEW_XRT | MACOSX)
#include <float.h>
#define XRT3D_HUGE_VAL     DBL_MAX
#else
#include <values.h>
#define XRT3D_HUGE_VAL     MAXDOUBLE
#endif /* (__STDC__ | WINNT | __VMS | VMS | Linux | NEW_XRT) */


#endif /* SYSTEM_H */
