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
$Header: /repo/nemo3d/src/io/fileformat.h,v 1.11 2006/11/30 22:12:52 nkharche
Exp $
*****************************************************************************/

#ifndef FILEFORMAT_H
#define FILEFORMAT_H

/*
 * for the embedding of NEMO 3D we need to have just the function definitions
 * at the top level defined.  Not the whole include path should be followed
 * down in order to make the embedding of the overall nemo3d library into
 * python simpler.
 */

#ifndef NEMO3D_DB

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "complex.h"
#include "io_utils.h"
#include "realtype.h"

#define MAX_VARS 10

#define FORMATLONG "%8i"
#define FORMATORIG "%5i"
#define FORMATOUT FORMATLONG

#define PDBFORMATATOM                                                          \
  "%5i" /* Used to write pdb file in rasmol readable format */
#define PDBFORMATCONCT                                                         \
  "%4i" /* Used to write pdb file in rasmol readable format */

/* For the structure output Wei Qiao (Purdue) decided that he does not need/
 *  want the values "  1.00  0.00" in the output file.   */
#define COORDINATEWRITEOLD "  1.00  0.00"
#define COORDINATEWRITENEW ""
#define COORDINATEWRITE COORDINATEWRITENEW

/*if COORDINATEWRITEOLD is chosen, COORDINATEREADOLD should be chosen either */
#define COORDINATEREADOLD "%f%f"
#define COORDINATEREADNEW ""
#define COORDINATEREAD COORDINATEREADNEW

using std::cout;
using std::endl;
using std::string;

void parseRecordFormat(const char *label, int *Nint, int *Nreal, int *Nfloat,
                       int *Ncol_var, int *sz_var, int *Nvar, char *fieldType);
FILE *readHeader(const char *filename, const char *recFmt, int *Ndim, int *dim,
                 char *dim_lbl[], int *Nuser, char *usr_comment[]);
void writeHeader(const char *file, const char *recFmt, int *dim,
                 char *dim_lbl[], char *usr_comment[]);
void displayHeader(const char *filename);

void fmtdat(int, char *[]);
void reformatBinaryData(char *filename, const char *flags, int prec);
void convertBinToAscii(char *filename, int prec, int _inclHdr);
void convertWfToPdb(char *filename, char *filename2);
void convertAsciiToBin(char *filename);
double *getBinFileData(const char *filename, int prec, int *rdim, int *cdim);
void writeDXfile(char *filename);
void read_oneLine(FILE *fp, char *xbuffer, char *ybuffer, char *zbuffer,
                  char *databuffer, int *x, int *y, int *z, double *value);
void writeDXHeader_For_WF_Profile(FILE *fp, int cellxmax, int cellymax,
                                  int cellzmax, int evN, double eLvalue,
                                  double lattice_x, double lattice_y,
                                  double lattice_z);
void writeDXHeader_For_Tension_Profile(FILE *fp, int cellxmax, int cellymax,
                                       int cellzmax, double lattice_x,
                                       double lattice_y, double lattice_z,
                                       char *tension);
void writeDXTail(FILE *fp);

#endif /* NEMO3D_DB */

// The definition below establishes an unmangled C function compiled with C++
// such that this function can be called from a C-linked executable
#ifdef __cplusplus
extern "C" {
#endif
double *getBinFileData_c(const char *filename, int prec, int *rdim, int *cdim);
#ifdef __cplusplus
}
#endif

#endif /* FILEFORMAT_H */
