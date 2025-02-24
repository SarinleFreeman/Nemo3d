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
$Header: /repo/nemo3d/src/math/cvector_lib.c,v 1.6 2004/12/07 19:47:59 marek Exp $
*****************************************************************************/

#include "cvector_lib.h"



/*
 * Set elements to 0.0.
 */
cvectr cvectrClear(cvectr cv) {  
    int i;

    for (i= 0; i < nml_dcv_extent(cv); i++) {
        cv[i].r = 0.0;
        cv[i].i = 0.0;
    }

    return cv;
}




/*
 * Copies cvector "cvSource" to the cvector "cvTarget".
 * It copies arraySize number of elements from cvSource to cvTarget.
 * If cvTarget is NULL, a new Cvector is allocated with the same dimension as cvSource.
 * If cvTarget exists and is different dimension from cvSource, memory for cvTarget
 * is deallocated, and a new Cvector is allocated with the same dimensions as cvSource.
 */
cvectr cp_cvectr_dim(cvectr cvTarget, cvectr cvSource, int arraySize) {
    int i;

    if (!cvSource) {
        n3d_warning("In cp_cvectr_dim():: Error - Source cvector 'cvSource' is NULL.\n");
        return cvSource;
    }

    if (!cvTarget) {
        cvTarget = Cvectr( nml_dcv_extent(cvSource) );
    }
    else { /* cvTarget exists */
        if ( nml_dcv_extent(cvSource) != nml_dcv_extent(cvTarget) ) {
            rm_cvectr( &cvTarget );
            cvTarget = Cvectr( nml_dcv_extent(cvSource) );
        }
    }

    if ( nml_dcv_extent(cvSource) < arraySize ) {
        n3d_warning("In cp_cvectr_dim():: specified arraySize to copy is larger than dimension of source cvector 'cvSource'.\n");
        arraySize = nml_dcv_extent(cvSource);
    }

    for (i = 0; i < arraySize; i++) {
	cvTarget[i] = cvSource[i];
    }
    return cvTarget;
}


/*
 * Sets the complex vector "cRandVector" with random numbers.
 */
int crandvectr(cvectr cRandVector) {
    int i;

    for ( i=0; i < nml_dcv_extent(cRandVector); i++ ) {
        cRandVector[i] = Complex( rand()/RAND_MAX, rand()/RAND_MAX );
    }

    return 0;
   
}


/*
 * Return the complex dot product of two complex vectors.
 */
complex cvectrmul(cvectr x, cvectr y) {
    int i;
    complex cz = czero;
    real a, b, c, d;

    for ( i=0; i < nml_dcv_extent(x) ; i++ ) {

        a = x[i].r;
        b = x[i].i;

        c = y[i].r;
        d = y[i].i;

        cz.r += a*c - b*d;
        cz.i += b*c + a*d;
    }

    return cz;
}




/*
 * Copies complex array "caSource" to the complex array "caTarget".
 * It copies arraySize number of elements from caSource to caTarget.
 */
void cp_carray2carray (complex* caSource, complex* caTarget, int arraySize) {
    int i;

    if (!caSource){
        n3d_warning("In cp_carray2carray():: Error - Source complex array 'caSource' is NULL.\n");
        return;
    }

    if (!caTarget) {
        n3d_warning("In cp_carray2carray():: Error - Target complex array 'caTarget' is NULL.\n");
        return;
    }

    for (i=0; i<arraySize; i++){
        caTarget[i] = caSource[i];
    }
    return;
}

