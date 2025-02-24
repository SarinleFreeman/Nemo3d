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
$Header: /repo/nemo3d/src/math/rvector_lib.c,v 1.5 2004/12/07 19:47:59 marek Exp $
*****************************************************************************/

#include "rvector_lib.h"



/**
 * Copies rvector "rvSource" to the rvector "raTarget".
 * If rvTarget is NULL, a new Rvector is allocated with the same dimension as rvSource.
 * If rvTarget exists and is different dimension from rvSource, memory for rvTarget
 * is deallocated, and a new Rvector is allocated with the same dimensions as rvSource.
 */
rvectr cp_rvectr (rvectr rvTarget, rvectr rvSource) {
    int i;

    if (!rvSource) {
        n3d_warning("In cp_rvectr():: Error - Source rvector 'rvSource' is NULL.\n");
        return rvSource;
    }

    if (!rvTarget) {
      rvTarget = Rvectr( nml_dv_extent(rvSource) );
    }
    else { /* rvTarget exists */
        /* Now check if different dimensions */
        if ( nml_dv_extent(rvTarget) != nml_dv_extent(rvSource) ) {
            rm_rvectr( &rvTarget );
                rvTarget = Rvectr( nml_dv_extent(rvSource) );
        }
    }

    for (i=0; i < nml_dv_extent(rvSource); i++) {
        rvTarget[i] = rvSource[i];
    }
    return rvTarget;
}



/**
 * Copies rvector "rvSource"  to the real array "raTarget".
 */
real* cp_rvectr2rarray(real *raTarget, rvectr rvSource) {
    int i;

    if (!rvSource){
        n3d_warning("In cp_rvectr2rarray():: Error - Source rvector 'rvSource' is NULL.\n");
        return raTarget;
    }

    if (!raTarget) {
        n3d_warning("In cp_rvectr2rarray():: Error - Target real array 'raTarget' is NULL.\n");
        return raTarget;
    }

    for (i=0; i < nml_dv_extent(rvSource); i++) {
        raTarget[i] = rvSource[i];
    }
    return raTarget;
}



/*
 * Copies real array "raSource" to the real array "raTarget".
 * It copies arraySize number of elements from raSource to raTarget.
 */
void cp_rarray2rarray(real* raSource, real* raTarget, int arraySize) {
    int i;

    if (!raSource){
        n3d_warning("In cp_rarray2rarray():: Error - Source real array 'raSource' is NULL.\n");
        return;
    }

    if (!raTarget) {
        n3d_warning("In cp_rarray2rarray():: Error - Target real array 'raTarget' is NULL.\n");
        return;
    }

    for (i=0; i<arraySize; i++){
        raTarget[i] = raSource[i];
    }
    return;
}


