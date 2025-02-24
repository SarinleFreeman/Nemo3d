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
$Header: /repo/nemo3d/src/io/i_3dout_files.h,v 1.23 2008/09/24 02:09:19 hoonryu Exp $ 
*****************************************************************************/

#ifndef I_3DOUT_FILES_H
#define I_3DOUT_FILES_H 1


#include "io_utils.h"
#include "Boolean.h"
#include "qd_struct.h"
#include "mat_def.h"
#include "run3d_mpi.h"
#include "i_mat.h"
#include "MPI_Timing.h"
#include "fileformat.h"
#include "nml_global.h"
#include "fvector.h"
#include "i3tensor.h"
#include "b3tensor.h"

#if (!defined MPI3d && defined FAKE_MPI)
#include "mpi_fake.h"

#endif /* (defined MPI3d && !defined FAKE_MPI) */

#define bzero(a,b)  memset(a,0,b)


string replaceSuffix(const char* ch, const char* sfx, int indx=-1);

void* malloc_strippedVector(qd_struct d, int sz_datstr, void* vIn);
char *getfileext(char *fname);


void writeHeader_siteInfo(qd_struct d, const char*, const char* filename, 
                          const char* label);
void writeSiteInfo(qd_struct d, const char*, const char* filename, 
                   const char* label, void* dat0, ...);
void readSiteInfo(qd_struct d, const char*, const char* filename, void* dat0, ...);
int initStartVect(char*, qd_struct d, cvectr p);
void readWaveFunction(char* filename, qd_struct d, cvectr p);

void writeAtomInfo(qd_struct d);
void writeShapeInfo(qd_struct d);
void writeNbrIndx(qd_struct d);
void writePhaseInfoForNEMO1D(qd_struct d);

void writePos(qd_struct d, const char *type);
void writeEigenvector(qd_struct d, complex* wf, int indx, const char* suffix="");
void writeMomentumMatrix(qd_struct d, complex* wf, int indx, const char* suffix="");
void writeTransitionRate(qd_struct d, real *wf, int indx, const char* suffix="");
void writeOverLap(qd_struct d, complex *wf, int indx, const char* suffix="");

void writePsiSqr(cmatrix wfc, int n, qd_struct d, int wf_indx);
void printPsisq(qd_struct d, cvectr p, int n_offset);

void writeEigenvalues(qd_struct d, char mode);

void multicastDisplacement( qd_struct d);
void readDisplacement( qd_struct d, rvectr p);


int print_bond_lengths( qd_struct d );

void write_nbr( qd_struct d, int first_or_second );


void printHam ( qd_struct d, char* fmt );
int print_hmat_struct ( qd_struct d );

int print_Psi4_IPR( qd_struct d , rvectr wf4, int n, cvectr eigval);
void hin_atom_parse(qd_struct d, geo_struct* di);

/* obsolete stuff OR used only by serial code */

int print_isosurf_pdb( qd_struct d, rmatrix data, real isoval, int eigvallabel );
int print_struct_pdb( qd_struct d );
int print_struct_gmv( qd_struct d );

int scatter_print_all(qd_struct d, rmatrix psisq, FILE *f_out,int isx, int iex);
int scatter_print_unitcell(qd_struct d, rmatrix psisq, FILE *f_out,int isx, int iex);
void print_psisq_gnu_pdb(cmatrix wfc,int n,qd_struct d,int n_offset);
void print_psisq_gnu_pdb_k(cmatrix,int,qd_struct,int);
int print_psisqr_gnuplot( qd_struct d, rmatrix psisqr, int k );

void write_displacement( qd_struct d);
void read_displacement( qd_struct d, rvectr p);

int write_cubic_struct_wavefunctions( qd_struct d );
void writeOneVector(qd_struct d, complex* x, complex* y);

#endif /* I_3DOUT_FILES_H */

