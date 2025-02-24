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
$Header: /repo/nemo3d/src/math/rmatrix_lib.h,v 1.2 2003/10/08 16:18:16 hook Exp $
*****************************************************************************/


#ifndef RMATRIX_LIB_H
#define RMATRIX_LIB_H 1

#include <rmatrix.h>
#include <rvector.h>
#include "rvector_lib.h" /* references typedef roperator and vector copy functions */

#define Swap(g, h)  { temp=(g); (g)=(h); (h)=temp; }

#ifdef __cplusplus
extern "C" {
#endif

#ifdef FORTRAN_UNDERSCORE
#define dgetrf            dgetrf_
#define dgetrs            dgetrs_
#endif

#ifdef FORTRAN_UNDERSCORE2
#define dgetrf            dgetrf__
#define dgetrs            dgetrs__
#endif

#ifdef FORTRAN_ALLCAPS
#define dgetrf            DGETRF
#define dgetrs            DGETRS
#endif

    void dgetrf(int*, int*, nml_dscalar*, int*, int**, int*);
    void dgetrs(char*, int*, int*, nml_dscalar*, int*, int**, rvectr, int*, int*);

#ifdef __cplusplus
}
#endif


roperator Roperator(int n);
rvectr luDecompositionBackSubstitution(roperator a0, roperator a, rvectr x, rvectr b);

#endif /* RMATRIX_LIB_H  */
