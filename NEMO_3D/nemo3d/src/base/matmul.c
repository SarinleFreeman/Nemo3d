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
$Header: /repo/nemo3d/src/base/matmul.c,v 1.9 2006/12/18 21:40:04 gekco Exp $ 
*****************************************************************************/


#include "matmul.h"


/*
#ifdef Linux
#define CMATMUL_EXPLICIT_INLINE_20
#endif
*/

#ifndef CMATMUL_EXPLICIT_INLINE_20
#define CMATMUL_POINTER_REDUCE_20
#endif

#ifndef CMATMUL_EXPLICIT_INLINE_10
#define CMATMUL_POINTER_REDUCE_10
#endif


int cmatmul_spds( cvectr y, cmatrix h, cvectr x, int isy, int isx, int hdim )
{
    if (hdim==20){
	cmatmul_spds_20( y, h, x, isy, isx, hdim );
	return 0;
    } else if (hdim==10){
	cmatmul_spds_10( y, h, x, isy, isx, hdim );
	return 0;
    } else {
	{
#if 0
	    int i, j;
	    for ( i=0; i < hdim; i++ )
		for ( j=0; j < hdim; j++ )
		    {
			y[i+isy].r += h[i][j].r*x[j+isx].r - h[i][j].i*x[j+isx].i;
			y[i+isy].i += h[i][j].r*x[j+isx].i + h[i][j].i*x[j+isx].r;
		    }
#else
	    int i, j;
	    register int ip,jp;
	    cvectr h_i;
	    register double y_ip_r,y_ip_i;
	    register real h_i_j_r, h_i_j_i;
	    for ( i=0; i < hdim; i++ ){
		ip= i+isy;
		h_i = h[i];
		y_ip_r = y[ip].r;
		y_ip_i = y[ip].i;
		for ( j=0; j < hdim; j++ ) {
		    jp=j+isx;
		    h_i_j_r = h_i[j].r, h_i_j_i = h_i[j].i;
		    y_ip_r += h_i_j_r*x[jp].r - h_i_j_i*x[jp].i;
		    y_ip_i += h_i_j_r*x[jp].i + h_i_j_i*x[jp].r;
		}
		y[ip].r = y_ip_r;
		y[ip].i = y_ip_i;
	    }
#endif
	}
    }
    return 0;
}

/* same as cmatmul_spds but perform an adjoint operation on the matrix h first. */
int cmatmul_spds_hc( cvectr y, cmatrix h, cvectr x, int isy, int isx, int hdim )
{
#if 0
    int i, j;
    for ( j=0; j < hdim; j++ )
	for ( i=0; i < hdim; i++ )
	    {
		y[i+isy].r += h[j][i].r*x[j+isx].r + h[j][i].i*x[j+isx].i;
		y[i+isy].i += h[j][i].r*x[j+isx].i - h[j][i].i*x[j+isx].r;
	    }
#else
    int i, j;
    register int jp;
    cvectr h_j;
    register double x_jp_r,x_jp_i;
    register real h_j_i_r, h_j_i_i;
    for ( j=0; j < hdim; j++ ){
	jp=j+isx;
	x_jp_r=x[jp].r, x_jp_i=x[jp].i;
	h_j=h[j];
	for ( i=0; i < hdim; i++ ){
	    h_j_i_r = h_j[i].r, h_j_i_i = h_j[i].i;
	    y[i+isy].r += h_j_i_r*x_jp_r + h_j_i_i*x_jp_i;
	    y[i+isy].i += h_j_i_r*x_jp_i - h_j_i_i*x_jp_r;
	}
    }
#endif
    return 0;
}

int cmatmul_spds_hc_20( cvectr y, cmatrix h, cvectr x, int isy, int isx, int hdim )
{
#if 0
    int i, j;
    register int jp;
    cvectr h_j;
    register double x_jp_r,x_jp_i;
    register real h_j_i_r, h_j_i_i;
    for ( j=0; j < 20; j++ ){
	jp=j+isx;
	x_jp_r=x[jp].r, x_jp_i=x[jp].i;
	h_j=h[j];
	for ( i=0; i < 20; i++ ){
	    h_j_i_r = h_j[i].r, h_j_i_i = h_j[i].i;
	    y[i+isy].r += h_j_i_r*x_jp_r + h_j_i_i*x_jp_i;
	    y[i+isy].i += h_j_i_r*x_jp_i - h_j_i_i*x_jp_r;
	}
    }
#else
    int i, j;
    register int jp;
    cvectr h_j;
    register double x_jp_r,x_jp_i;
    register real h_j_i_r, h_j_i_i;
    static real y_r[20],y_i[20];
    for ( i=0; i < 20; i++ ){
	y_r[i] = y[i+isy].r;
	y_i[i] = y[i+isy].i;
    }
    for ( j=0; j < 20; j++ ){
	jp=j+isx;
	x_jp_r=x[jp].r, x_jp_i=x[jp].i;
	h_j=h[j];
	for ( i=0; i < 20; i++ ){
	    h_j_i_r = h_j[i].r, h_j_i_i = h_j[i].i;
	    y_r[i] += h_j_i_r*x_jp_r + h_j_i_i*x_jp_i;
	    y_i[i] += h_j_i_r*x_jp_i - h_j_i_i*x_jp_r;
	}
    }
    for ( i=0; i < 20; i++ ){
	y[i+isy].r = y_r[i];
	y[i+isy].i = y_i[i];
    }
#endif
	    
    return 0;
}

int cmatmul_spds_hc_10( cvectr y, cmatrix h, cvectr x, int isy, int isx, int hdim )
{
    int i, j;
    register int jp;
    cvectr h_j;
    register double x_jp_r,x_jp_i;
    register real h_j_i_r, h_j_i_i;
    for ( j=0; j < 10; j++ ){
	jp=j+isx;
	x_jp_r=x[jp].r, x_jp_i=x[jp].i;
	h_j=h[j];
	for ( i=0; i < 10; i++ ){
	    h_j_i_r = h_j[i].r, h_j_i_i = h_j[i].i;
	    y[i+isy].r += h_j_i_r*x_jp_r + h_j_i_i*x_jp_i;
	    y[i+isy].i += h_j_i_r*x_jp_i - h_j_i_i*x_jp_r;
	}
    }
	    
    return 0;
}


int cmatmul_spds_20( cvectr y, cmatrix h, cvectr x, int isy, int isx, int hdim )
{
#ifdef CMATMUL_POINTER_REDUCE_20
    register int i, j;
    register int ip,jp;
    cvectr h_i;
    register double y_ip_r,y_ip_i;
    register double h_i_j_r, h_i_j_i;
    static double xjp_r[20], xjp_i[20];
    for ( j=0; j < 20; j++ ) {
	jp=j+isx;
	xjp_r[j]=x[jp].r;
	xjp_i[j]=x[jp].i;
    }
    for ( i=0; i < 20; i++ ){
	ip= i+isy;
	h_i = h[i];
	y_ip_r = y[ip].r;
	y_ip_i = y[ip].i;
	for ( j=0; j < 20; j++ ) {
	    h_i_j_r = h_i[j].r, h_i_j_i = h_i[j].i;
	    y_ip_r += h_i_j_r*xjp_r[j] - h_i_j_i*xjp_i[j];
	    y_ip_i += h_i_j_r*xjp_i[j] + h_i_j_i*xjp_r[j];
	}
	y[ip].r = y_ip_r;
	y[ip].i = y_ip_i;
    }
    return 0;
#endif

#ifdef CMATMUL_EXPLICIT_INLINE_20
    register int i, j;
    register int ip,jp;
    cvectr h_i;
    register double y_ip_r,y_ip_i;
    register double h_i_j_r, h_i_j_i;
    register double xjp_r_00, xjp_i_00;
    register double xjp_r_01, xjp_i_01;
    register double xjp_r_02, xjp_i_02;
    register double xjp_r_03, xjp_i_03;
    register double xjp_r_04, xjp_i_04;
    register double xjp_r_05, xjp_i_05;
    register double xjp_r_06, xjp_i_06;
    register double xjp_r_07, xjp_i_07;
    register double xjp_r_08, xjp_i_08;
    register double xjp_r_09, xjp_i_09;
    register double xjp_r_10, xjp_i_10;
    register double xjp_r_11, xjp_i_11;
    register double xjp_r_12, xjp_i_12;
    register double xjp_r_13, xjp_i_13;
    register double xjp_r_14, xjp_i_14;
    register double xjp_r_15, xjp_i_15;
    register double xjp_r_16, xjp_i_16;
    register double xjp_r_17, xjp_i_17;
    register double xjp_r_18, xjp_i_18;
    register double xjp_r_19, xjp_i_19;

    jp=isx;
    xjp_r_00=x[jp].r;
    xjp_i_00=x[jp].i;
    
    jp++;;
    xjp_r_01=x[jp].r;
    xjp_i_01=x[jp].i;
    
    jp++;;
    xjp_r_02=x[jp].r;
    xjp_i_02=x[jp].i;
    
    jp++;;
    xjp_r_03=x[jp].r;
    xjp_i_03=x[jp].i;
    
    jp++;;
    xjp_r_04=x[jp].r;
    xjp_i_04=x[jp].i;
    
    jp++;;
    xjp_r_05=x[jp].r;
    xjp_i_05=x[jp].i;
    
    jp++;;
    xjp_r_06=x[jp].r;
    xjp_i_06=x[jp].i;
    
    jp++;;
    xjp_r_07=x[jp].r;
    xjp_i_07=x[jp].i;
    
    jp++;;
    xjp_r_08=x[jp].r;
    xjp_i_08=x[jp].i;
    
    jp++;;
    xjp_r_09=x[jp].r;
    xjp_i_09=x[jp].i;
    
    jp++;;
    xjp_r_10=x[jp].r;
    xjp_i_10=x[jp].i;
    
    jp++;;
    xjp_r_11=x[jp].r;
    xjp_i_11=x[jp].i;
    
    jp++;;
    xjp_r_12=x[jp].r;
    xjp_i_12=x[jp].i;
    
    jp++;;
    xjp_r_13=x[jp].r;
    xjp_i_13=x[jp].i;
    
    jp++;;
    xjp_r_14=x[jp].r;
    xjp_i_14=x[jp].i;
    
    jp++;;
    xjp_r_15=x[jp].r;
    xjp_i_15=x[jp].i;
    
    jp++;;
    xjp_r_16=x[jp].r;
    xjp_i_16=x[jp].i;
    
    jp++;;
    xjp_r_17=x[jp].r;
    xjp_i_17=x[jp].i;
    
    jp++;;
    xjp_r_18=x[jp].r;
    xjp_i_18=x[jp].i;
    
    jp++;;
    xjp_r_19=x[jp].r;
    xjp_i_19=x[jp].i;


    for ( i=0; i < 20; i++ ){
	ip= i+isy;
	h_i = h[i];
	y_ip_r = y[ip].r;
	y_ip_i = y[ip].i;

	j=0;
	h_i_j_r = h_i[j].r, h_i_j_i = h_i[j].i;
	y_ip_r += h_i_j_r*xjp_r_00 - h_i_j_i*xjp_i_00;
	y_ip_i += h_i_j_r*xjp_i_00 + h_i_j_i*xjp_r_00;

	j++;
	h_i_j_r = h_i[j].r, h_i_j_i = h_i[j].i;
	y_ip_r += h_i_j_r*xjp_r_01 - h_i_j_i*xjp_i_01;
	y_ip_i += h_i_j_r*xjp_i_01 + h_i_j_i*xjp_r_01;

	j++;
	h_i_j_r = h_i[j].r, h_i_j_i = h_i[j].i;
	y_ip_r += h_i_j_r*xjp_r_02 - h_i_j_i*xjp_i_02;
	y_ip_i += h_i_j_r*xjp_i_02 + h_i_j_i*xjp_r_02;

	j++;
	h_i_j_r = h_i[j].r, h_i_j_i = h_i[j].i;
	y_ip_r += h_i_j_r*xjp_r_03 - h_i_j_i*xjp_i_03;
	y_ip_i += h_i_j_r*xjp_i_03 + h_i_j_i*xjp_r_03;

	j++;
	h_i_j_r = h_i[j].r, h_i_j_i = h_i[j].i;
	y_ip_r += h_i_j_r*xjp_r_04 - h_i_j_i*xjp_i_04;
	y_ip_i += h_i_j_r*xjp_i_04 + h_i_j_i*xjp_r_04;

	j++;
	h_i_j_r = h_i[j].r, h_i_j_i = h_i[j].i;
	y_ip_r += h_i_j_r*xjp_r_05 - h_i_j_i*xjp_i_05;
	y_ip_i += h_i_j_r*xjp_i_05 + h_i_j_i*xjp_r_05;

	j++;
	h_i_j_r = h_i[j].r, h_i_j_i = h_i[j].i;
	y_ip_r += h_i_j_r*xjp_r_06 - h_i_j_i*xjp_i_06;
	y_ip_i += h_i_j_r*xjp_i_06 + h_i_j_i*xjp_r_06;

	j++;
	h_i_j_r = h_i[j].r, h_i_j_i = h_i[j].i;
	y_ip_r += h_i_j_r*xjp_r_07 - h_i_j_i*xjp_i_07;
	y_ip_i += h_i_j_r*xjp_i_07 + h_i_j_i*xjp_r_07;

	j++;
	h_i_j_r = h_i[j].r, h_i_j_i = h_i[j].i;
	y_ip_r += h_i_j_r*xjp_r_08 - h_i_j_i*xjp_i_08;
	y_ip_i += h_i_j_r*xjp_i_08 + h_i_j_i*xjp_r_08;

	j++;
	h_i_j_r = h_i[j].r, h_i_j_i = h_i[j].i;
	y_ip_r += h_i_j_r*xjp_r_09 - h_i_j_i*xjp_i_09;
	y_ip_i += h_i_j_r*xjp_i_09 + h_i_j_i*xjp_r_09;

	j++;
	h_i_j_r = h_i[j].r, h_i_j_i = h_i[j].i;
	y_ip_r += h_i_j_r*xjp_r_10 - h_i_j_i*xjp_i_10;
	y_ip_i += h_i_j_r*xjp_i_10 + h_i_j_i*xjp_r_10;

	j++;
	h_i_j_r = h_i[j].r, h_i_j_i = h_i[j].i;
	y_ip_r += h_i_j_r*xjp_r_11 - h_i_j_i*xjp_i_11;
	y_ip_i += h_i_j_r*xjp_i_11 + h_i_j_i*xjp_r_11;

	j++;
	h_i_j_r = h_i[j].r, h_i_j_i = h_i[j].i;
	y_ip_r += h_i_j_r*xjp_r_12 - h_i_j_i*xjp_i_12;
	y_ip_i += h_i_j_r*xjp_i_12 + h_i_j_i*xjp_r_12;

	j++;
	h_i_j_r = h_i[j].r, h_i_j_i = h_i[j].i;
	y_ip_r += h_i_j_r*xjp_r_13 - h_i_j_i*xjp_i_13;
	y_ip_i += h_i_j_r*xjp_i_13 + h_i_j_i*xjp_r_13;

	j++;
	h_i_j_r = h_i[j].r, h_i_j_i = h_i[j].i;
	y_ip_r += h_i_j_r*xjp_r_14 - h_i_j_i*xjp_i_14;
	y_ip_i += h_i_j_r*xjp_i_14 + h_i_j_i*xjp_r_14;

	j++;
	h_i_j_r = h_i[j].r, h_i_j_i = h_i[j].i;
	y_ip_r += h_i_j_r*xjp_r_15 - h_i_j_i*xjp_i_15;
	y_ip_i += h_i_j_r*xjp_i_15 + h_i_j_i*xjp_r_15;

	j++;
	h_i_j_r = h_i[j].r, h_i_j_i = h_i[j].i;
	y_ip_r += h_i_j_r*xjp_r_16 - h_i_j_i*xjp_i_16;
	y_ip_i += h_i_j_r*xjp_i_16 + h_i_j_i*xjp_r_16;

	j++;
	h_i_j_r = h_i[j].r, h_i_j_i = h_i[j].i;
	y_ip_r += h_i_j_r*xjp_r_17 - h_i_j_i*xjp_i_17;
	y_ip_i += h_i_j_r*xjp_i_17 + h_i_j_i*xjp_r_17;

	j++;
	h_i_j_r = h_i[j].r, h_i_j_i = h_i[j].i;
	y_ip_r += h_i_j_r*xjp_r_18 - h_i_j_i*xjp_i_18;
	y_ip_i += h_i_j_r*xjp_i_18 + h_i_j_i*xjp_r_18;

	j++;
	h_i_j_r = h_i[j].r, h_i_j_i = h_i[j].i;
	y_ip_r += h_i_j_r*xjp_r_19 - h_i_j_i*xjp_i_19;
	y_ip_i += h_i_j_r*xjp_i_19 + h_i_j_i*xjp_r_19;




	y[ip].r = y_ip_r;
	y[ip].i = y_ip_i;
    }
    return 0;
#endif


#if 0
    /* Use a pointer to the matrix element h_i_p with incremental 
       Pointers.  This does not seem to work very well on the 
       Intel Pentium III. */
    register int i, j;
    register int ip,jp;
    cvectr h_i;
    register complex *h_i_p;
    register double y_ip_r,y_ip_i;
    register double h_i_j_r, h_i_j_i;
    static double xjp_r[20], xjp_i[20];
    for ( j=0; j < 20; j++ ) {
	jp=j+isx;
	xjp_r[j]=x[jp].r;
	xjp_i[j]=x[jp].i;
    }
    for ( i=0; i < 20; i++ ){
	ip= i+isy;
	h_i = h[i];
	y_ip_r = y[ip].r;
	y_ip_i = y[ip].i;
	h_i_p = h_i;
	for ( j=0; j < 20; j++ ) {
	    
	    h_i_j_r = (*h_i_p).r, h_i_j_i = (*h_i_p).i;
	    y_ip_r += h_i_j_r*xjp_r[j] - h_i_j_i*xjp_i[j];
	    y_ip_i += h_i_j_r*xjp_i[j] + h_i_j_i*xjp_r[j];
	    h_i_p++;
	}
	y[ip].r = y_ip_r;
	y[ip].i = y_ip_i;
    }
    return 0;
#endif
#if 0
    /* use a temporary complex register....
       need to do pointer arithmatic on the inner loop....*/
    int i, j;
    register int ip,jp;
    cvectr h_i;
    register complex y_ip;
    register double h_i_j_r, h_i_j_i;
        static double xjp_r[20], xjp_i[20];
        for ( j=0; j < 20; j++ ) {
            jp=j+isx;
            xjp_r[j]=x[jp].r;
            xjp_i[j]=x[jp].i;
        }
        for ( i=0; i < 20; i++ ){
            ip= i+isy;
            h_i = h[i];
            y_ip = y[ip];
            for ( j=0; j < 20; j++ ) {
                    h_i_j_r = h_i[j].r, h_i_j_i = h_i[j].i;
                    y_ip.r += h_i_j_r*xjp_r[j] - h_i_j_i*xjp_i[j];
                    y_ip.i += h_i_j_r*xjp_i[j] + h_i_j_i*xjp_r[j];
            }
            y[ip] = y_ip;
        }
        return 0;
#endif
}

int cmatmul_spds_10( cvectr y, cmatrix h, cvectr x, int isy, int isx, int hdim )
{
#ifdef CMATMUL_POINTER_REDUCE_10
    register int i, j;
    register int ip,jp;
    cvectr h_i;
    register double y_ip_r,y_ip_i;
    register double h_i_j_r, h_i_j_i;
    static double xjp_r[10], xjp_i[10];
    for ( j=0; j < 10; j++ ) {
	jp=j+isx;
	xjp_r[j]=x[jp].r;
	xjp_i[j]=x[jp].i;
    }
    for ( i=0; i < 10; i++ ){
	ip= i+isy;
	h_i = h[i];
	y_ip_r = y[ip].r;
	y_ip_i = y[ip].i;
	for ( j=0; j < 10; j++ ) {
	    h_i_j_r = h_i[j].r, h_i_j_i = h_i[j].i;
	    y_ip_r += h_i_j_r*xjp_r[j] - h_i_j_i*xjp_i[j];
	    y_ip_i += h_i_j_r*xjp_i[j] + h_i_j_i*xjp_r[j];
	}
	y[ip].r = y_ip_r;
	y[ip].i = y_ip_i;
    }
    return 0;
#endif

#ifdef CMATMUL_EXPLICIT_INLINE_10
    register int i, j;
    register int ip,jp;
    cvectr h_i;
    register double y_ip_r,y_ip_i;
    register double h_i_j_r, h_i_j_i;
    register double xjp_r_00, xjp_i_00;
    register double xjp_r_01, xjp_i_01;
    register double xjp_r_02, xjp_i_02;
    register double xjp_r_03, xjp_i_03;
    register double xjp_r_04, xjp_i_04;
    register double xjp_r_05, xjp_i_05;
    register double xjp_r_06, xjp_i_06;
    register double xjp_r_07, xjp_i_07;
    register double xjp_r_08, xjp_i_08;
    register double xjp_r_09, xjp_i_09;

    jp=isx;
    xjp_r_00=x[jp].r;
    xjp_i_00=x[jp].i;
    
    jp++;;
    xjp_r_01=x[jp].r;
    xjp_i_01=x[jp].i;
    
    jp++;;
    xjp_r_02=x[jp].r;
    xjp_i_02=x[jp].i;
    
    jp++;;
    xjp_r_03=x[jp].r;
    xjp_i_03=x[jp].i;
    
    jp++;;
    xjp_r_04=x[jp].r;
    xjp_i_04=x[jp].i;
    
    jp++;;
    xjp_r_05=x[jp].r;
    xjp_i_05=x[jp].i;
    
    jp++;;
    xjp_r_06=x[jp].r;
    xjp_i_06=x[jp].i;
    
    jp++;;
    xjp_r_07=x[jp].r;
    xjp_i_07=x[jp].i;
    
    jp++;;
    xjp_r_08=x[jp].r;
    xjp_i_08=x[jp].i;
    
    jp++;;
    xjp_r_09=x[jp].r;
    xjp_i_09=x[jp].i;
    



    for ( i=0; i < 10; i++ ){
	ip= i+isy;
	h_i = h[i];
	y_ip_r = y[ip].r;
	y_ip_i = y[ip].i;

	j=0;
	h_i_j_r = h_i[j].r, h_i_j_i = h_i[j].i;
	y_ip_r += h_i_j_r*xjp_r_00 - h_i_j_i*xjp_i_00;
	y_ip_i += h_i_j_r*xjp_i_00 + h_i_j_i*xjp_r_00;

	j++;
	h_i_j_r = h_i[j].r, h_i_j_i = h_i[j].i;
	y_ip_r += h_i_j_r*xjp_r_01 - h_i_j_i*xjp_i_01;
	y_ip_i += h_i_j_r*xjp_i_01 + h_i_j_i*xjp_r_01;

	j++;
	h_i_j_r = h_i[j].r, h_i_j_i = h_i[j].i;
	y_ip_r += h_i_j_r*xjp_r_02 - h_i_j_i*xjp_i_02;
	y_ip_i += h_i_j_r*xjp_i_02 + h_i_j_i*xjp_r_02;

	j++;
	h_i_j_r = h_i[j].r, h_i_j_i = h_i[j].i;
	y_ip_r += h_i_j_r*xjp_r_03 - h_i_j_i*xjp_i_03;
	y_ip_i += h_i_j_r*xjp_i_03 + h_i_j_i*xjp_r_03;

	j++;
	h_i_j_r = h_i[j].r, h_i_j_i = h_i[j].i;
	y_ip_r += h_i_j_r*xjp_r_04 - h_i_j_i*xjp_i_04;
	y_ip_i += h_i_j_r*xjp_i_04 + h_i_j_i*xjp_r_04;

	j++;
	h_i_j_r = h_i[j].r, h_i_j_i = h_i[j].i;
	y_ip_r += h_i_j_r*xjp_r_05 - h_i_j_i*xjp_i_05;
	y_ip_i += h_i_j_r*xjp_i_05 + h_i_j_i*xjp_r_05;

	j++;
	h_i_j_r = h_i[j].r, h_i_j_i = h_i[j].i;
	y_ip_r += h_i_j_r*xjp_r_06 - h_i_j_i*xjp_i_06;
	y_ip_i += h_i_j_r*xjp_i_06 + h_i_j_i*xjp_r_06;

	j++;
	h_i_j_r = h_i[j].r, h_i_j_i = h_i[j].i;
	y_ip_r += h_i_j_r*xjp_r_07 - h_i_j_i*xjp_i_07;
	y_ip_i += h_i_j_r*xjp_i_07 + h_i_j_i*xjp_r_07;

	j++;
	h_i_j_r = h_i[j].r, h_i_j_i = h_i[j].i;
	y_ip_r += h_i_j_r*xjp_r_08 - h_i_j_i*xjp_i_08;
	y_ip_i += h_i_j_r*xjp_i_08 + h_i_j_i*xjp_r_08;

	j++;
	h_i_j_r = h_i[j].r, h_i_j_i = h_i[j].i;
	y_ip_r += h_i_j_r*xjp_r_09 - h_i_j_i*xjp_i_09;
	y_ip_i += h_i_j_r*xjp_i_09 + h_i_j_i*xjp_r_09;



	y[ip].r = y_ip_r;
	y[ip].i = y_ip_i;
    }
    return 0;
#endif



}

int cmatmul_spds_sc( cvectr y, cmatrix h, cvectr x, int isy, int isx, int hdim, complex scalar )
{
    int i, j;

#ifdef _OLD_eikHx
    register int ip,jp;
    register complex a;
    cvectr h_i;

    for ( i=0; i < hdim; i++ ){
	ip= i+isy;
	h_i = h[i];
        for ( j=0; j < hdim; j++ ) {
	    a.r = scalar.r*h_i[j].r - scalar.i*h_i[j].i;
	    a.i = scalar.r*h_i[j].i + scalar.i*h_i[j].r;  
	    jp=j+isx;
	    y[ip].r += a.r*x[jp].r - a.i*x[jp].i;
	    y[ip].i += a.r*x[jp].i + a.i*x[jp].r;
	} 
    }
#endif

    for ( i=0; i < hdim; i++ ){
       register complex Hx = Complex(0.0,0.0);
       cvectr h_i = h[i];
       complex* y_i = y + i + isy;
       for ( j=0; j < hdim; j++ ) {
          complex* x_j = x + j + isx;
          complex* H_ij = h_i + j;
          Hx.r += H_ij->r*x_j->r - H_ij->i*x_j->i;
          Hx.i += H_ij->r*x_j->i + H_ij->i*x_j->r;
       }
       y_i->r += scalar.r*Hx.r - scalar.i*Hx.i;
       y_i->i += scalar.r*Hx.i + scalar.i*Hx.r;
    }
    return(0);
}





#ifdef CMATMUL_SPDS_1_MACRO 
#else /* CMATMUL_SPDS_1_MACRO */
int cmatmul_spds_1( cvectr y, cmatrix h, cvectr x, int isy, int isx, int hdim )
{
    y[isy].r += h[0][0].r*x[isx].r - h[0][0].i*x[isx].i;
    y[isy].i += h[0][0].r*x[isx].i + h[0][0].i*x[isx].r;

    return 0;

}
/* same as cmatmul_spds but perform an adjoint operation on the matrix h first. */
int cmatmul_spds_hc_1( cvectr y, cmatrix h, cvectr x, int isy, int isx, int hdim )
{
    y[isy].r += h[0][0].r*x[isx].r + h[0][0].i*x[isx].i;
    y[isy].i += h[0][0].r*x[isx].i - h[0][0].i*x[isx].r;
    return 0;
    
}
int cmatmul_sparse_spds_1( cvectr y, double h, cvectr x, int isy, int isx )
{
    y[isy].r += h*x[isx].r ;
    y[isy].i += h*x[isx].i ;

    return 0;

}
/* same as cmatmul_spds but perform an adjoint operation on the matrix h first. */
int cmatmul_sparse_spds_hc_1( cvectr y, double h, cvectr x, int isy, int isx )
{
    y[isy].r += h*x[isx].r;
    y[isy].i += h*x[isx].i;
    return 0;
    
}
#endif /* CMATMUL_SPDS_1_MACRO */

