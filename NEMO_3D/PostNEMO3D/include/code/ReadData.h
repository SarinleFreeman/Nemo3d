/*****************************************************************************
The Jet Propulsion Laboratory (JPL) NanoElectronicMOdeling-3D
PostProcessing package.
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
  Seungwon Lee (Seungwon.Lee@jpl.nasa.gov)

Written by:  Seungwon Lee
*****************************************************************************/


#ifndef _read_data_h
#define _read_data_h 

#include <stdio.h>
#include <stdlib.h>
#include <stream.h>
#include <complex>
using namespace std;

void die(char* s, ...);
void parseRecordFormat(char*, int*, int*, int*, int*, int*, char*);
FILE* readHeader(const char*, char*, int*, int*, char*, int*, char*);
void read_natoms_norbitals(const char*, int*, int*, const int*);
void read_wf(const char*, int, int* eList, complex<double>**, int, int, int, int);
void read_write_wf_10(const char*, int, int* eList, int);

void read_data_for_postprocessing(const char*, int, const int*, int, const int*,
                                  complex<double>**, complex<double>**, int**, int**,
				  double**, double**, double**,  int, int, int, int);
void read_data_for_postprocessing_new(const char*, int, const int*, int, const int*,
                                  int**, int**, double**, double**, double**);

#endif
