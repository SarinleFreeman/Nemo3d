#ifndef _nml_global_h
#define _nml_global_h 1

/*****************************************************************************
The NEMO Math Library.
Copyright (C) 2002-2004 California Institute of Technology (Caltech)

This file is part of
The NanoElectronic MOdeling (NEMO) Math Library.

This library is free software which you can redistribute and/or modify
under the terms of the GNU Library General Public License
as published by the Free Software Foundation;
either version 2, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this library; see the file COPYING. If not, write to the
Free Software Foundation, Inc.,
59 Temple Place, Suite 330,
Boston, MA  02111-1307  USA

Written by E. Robert Tisdale
Modified by Marek J. Korkusinski, October 2004

*****************************************************************************/

#include<math.h>
#include<float.h>
#include<limits.h>
#include<stdio.h>
#include<stdlib.h>

#define LOGFILE stderr

#ifndef LONG_8Byte
#define LONG_8Byte long long
#ifdef mem_64BIT
#define LONG_8Byte long
#endif
#endif

/* global type definitions	*/

typedef   signed int	ptrdiff_nml_t;

typedef    size_t	nml_offset;
typedef    size_t	nml_extent;
typedef ptrdiff_nml_t	nml_stride;
typedef ptrdiff_nml_t	nml_index;

typedef	     float	nml_single;
typedef      double	nml_double;
typedef	     signed int	nml_integer;
typedef	   unsigned int	nml_logical;
typedef     short       nml_short ;
typedef      unsigned char        nml_char ;
typedef     float       nml_float ;



void * nml_malloc( nml_extent ) ;
void * nml_calloc( nml_extent , nml_extent ) ; 
void   nml_free( void * ) ; 
void   nml_memory_report( char *tag) ;


void
nml_message(const char* m);

#define nml_message(m) fprintf(stderr, (m))

typedef unsigned int nml_fmtflags;

#define nml_skipws	= (1 <<  0) /* skip whitespace on input		   */
				    /* field adjustment:		   */
#define nml_left	= (1 <<  1) /* pad after value			   */
#define nml_right	= (1 <<  2) /* pad before value			   */
#define nml_internal	= (1 <<  3) /* pad between sign and value	   */

#define nml_boolalpha	= (1 <<  4) /* represent true and false symbolically*/
				    /* integer base:			   */
#define nml_dec	  	= (1 <<  5) /* base 10 output (decimal)		   */
#define nml_hex	  	= (1 <<  6) /* base 16 output (hexadecimal)	   */
#define nml_oct	  	= (1 <<  7) /* base  8 output (octal)		   */

#define nml_scientific  = (1 <<  8) /* floating-point notation d.ddddddEdd */
#define nml_fixed	= (1 <<  9) /* fixed-point    notation dddd.dd	   */

#define nml_showbase	= (1 << 10) /* prefix oct by 0 and hex by 0x on output*/
#define nml_showpoint	= (1 << 11) /* print trailing zeros		   */
#define nml_showpos	= (1 << 12) /* explicit '+' for positive ints	   */
#define nml_uppercase	= (1 << 13) /* 'E', 'X' rather than 'e', 'x'	   */
				    /* flags related to:		   */
#define nml_adjustfield = (1 << 14) /* adjustment			   */
#define nml_basefield	= (1 << 15) /* integer base			   */
#define nml_floatfield	= (1 << 16) /* floating-point output		   */

#endif /* _nml_global_h */
