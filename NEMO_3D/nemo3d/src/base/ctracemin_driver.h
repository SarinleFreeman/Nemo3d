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
*****************************************************************************/

#ifndef CTRACEMIN_DRIVER_H
#define CTRACEMIN_DRIVER_H

//#define SMALLER_THEN_ALL_EIG -14 ////-40  //-14 //-60
//#define LARGER_THEN_ALL_EIG   40 //// 40  // 39 // 60

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus 

#ifdef FORTRAN_UNDERSCORE
#define zcopy                        zcopy_
#define zdotc                        zdotc_
#define zlarnv                       zlarnv_
#define dznrm2                       dznrm2_
#define zgtracemin_wrapper           zgtracemin_wrapper_
#define map_chebyshev_shift_wrapper  map_chebyshev_shift_wrapper_
#define pcextract_eigenpairs_wrapper pcextract_eigenpairs_wrapper_
#endif

#ifdef FORTRAN_UNDERSCORE2
#define zcopy                        zcopy__
#define zdotc                        zdotc__
#define zlarnv                       zlarnv__
#define dznrm2                       dznrm2__
#define zgtracemin_wrapper           zgtracemin_wrapper__
#define map_chebyshev_shift_wrapper  map_chebyshev_shift_wrapper__
#define pcextract_eigenpairs_wrapper pcextract_eigenpairs_wrapper__
#endif

#ifdef FORTRAN_ALLCAPS
#define zcopy                        ZCOPY
#define zdotc                        ZDOTC
#define zlarnv                       ZLARNV
#define dznrm2                       DZNRM2
#define zgtracemin_wrapper           ZGTRACEMIN_WRAPPER             
#define map_chebyshev_shift_wrapper  MAP_CHEBYSHEV_SHIFT_WRAPPER
#define pcextract_eigenpairs_wrapper PCEXTRACT_EIGENPAIRS_WRAPPER
#endif

//LAPACK Functions
void zcopy (int *, complex *, int *, complex *, int *);                  //subroutine used to copy vectors
void zlarnv(int *, int *, int *, complex *);                             //subroutine used to generate random vector 
real dznrm2(int *, complex *, int *);                                    //functions used to compute norm of two vectors
void zdotc (complex *, int *, complex *, int *, complex *, int *);       //function used to compute complex conjugate dot product of two vectors

//TRACEMIN WRAPPER FUNCTIONS
void map_chebyshev_shift_wrapper(int *,int *,int *,complex *,complex *,complex *,complex *,complex *,int *,real *,real *,real *,int *,int *);
void zgtracemin_wrapper(int *,int *,complex *,real *,complex *, complex *,complex *,complex *,complex *,complex *,real *,int *,complex *,real *,int *,int *,int *,int *,int *,int *,int *);
void pcextract_eigenpairs_wrapper(complex *,real *,int *,int *,complex *,complex *,complex *,real *,complex *,int *,int *,int *);

#ifdef __cplusplus
}
#endif // __cplusplus 


int ctracemin_driver( qd_struct d );

#endif // CTRACEMIN_DRIVER_H 
