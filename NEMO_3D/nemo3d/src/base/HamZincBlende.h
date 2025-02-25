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
$Header: /repo/nemo3d/src/base/HamZincBlende.h,v 1.12 2004/10/29 02:54:34 gekco
Exp $
*****************************************************************************/

#ifndef HAMZINCBLENDE_H
#define HAMZINCBLENDE_H

#include "complex.h"
#include "cvector.h"
#include "io_utils.h"
#include "ivector.h"
#include "mb_ham_spds_micro2.h"
#include "rvector.h"
#include <cmatrix.h>

#ifdef HamMemorySave
#include "fvector.h"
#include "komplex.h"
#include "kvector.h"
#define ham_mem_rvectr fvectr
#define ham_mem_Rvectr Fvectr
#define ham_mem_rm_rvectr rm_fvectr
#define ham_mem_real_ptr float *
#define ham_mem_real float

#define ham_mem_cvectr kvectr
#define ham_mem_Cvectr Kvectr
#define ham_mem_rm_cvectr rm_kvectr
#define ham_mem_complex_ptr komplex *
#define ham_mem_complex komplex

#define ham_mem_double (float)

#else
#define ham_mem_rvectr rvectr
#define ham_mem_Rvectr Rvectr
#define ham_mem_rm_rvectr rm_rvectr
#define ham_mem_real_ptr real *
#define ham_mem_real real

#define ham_mem_cvectr cvectr
#define ham_mem_Cvectr Cvectr
#define ham_mem_rm_cvectr rm_cvectr
#define ham_mem_complex_ptr complex *
#define ham_mem_complex complex

#endif

class HamZincBlende {
  void print_20_sp3d5ss_spin();
  void print_10_sp3ss_spin();
  void print_10_sp3d5ss_nospin();

public:
  int N_Basis;
  BM_Type BandModel;
  int Natom;
  int Natom_surf;
  bool MagneticFieldOn;

  complex HBxy; // off-diagonal element of H_{atom-atom} due to magnetic field
                // in x or y direction. (Bx-iBy)*bohr_magneton

  ham_mem_rvectr Hdd; // diagonal elements of H_{atom-atom}
  ham_mem_cvectr Hdu; // off-diagonal elements of H_{atom-atom}
  ham_mem_rvectr Hds; // off-diagonal shifts of H_{atom-atom} for surface atoms
  ham_mem_rvectr Hu;  // matrix of off-diag elements of H on this processor
  ham_mem_rvectr
      Ho; // matrix of off-diag elements of H on neighboring processor
  ham_mem_cvectr HphaseIn; // vector-potential phase of off-diag element of H on
                           // this processor
  ham_mem_cvectr HphaseOut; // vector-potential phase of off-diag element of H
                            // on neighboring processor

  ivectr indxHU_to_row;
  ivectr indxHU_to_col;
  ivectr indxHO_to_row;
  ivectr indxHO_to_col;
  ivectr indxHU_nbr;
  ivectr indxHO_nbr;

  HamZincBlende() {
    Hdd = 0;
    Hdu = 0;
    Hds = 0;
    Hu = 0;
    Ho = 0;
    HphaseIn = 0;
    HphaseOut = 0;
    indxHU_to_row = indxHU_to_col = indxHO_to_row = indxHO_to_col = indxHU_nbr =
        indxHO_nbr = 0;
  }
  void Initialize(BM_Type BandModel, int N_Basis, int Natom, int Natom_surf,
                  int Nout, int Nin, bool MagneticFieldOn, complex HBxy);
  void setDiag(const cmatrix hd, int Zatom, int I, int J, int Zsurf);
  void printStorageInfo();
  void print() {
    switch (BandModel) {
    case BM_10_sp3ss_spin:
      print_10_sp3ss_spin();
      break;
    case BM_20_sp3d5ss_spin:
      print_20_sp3d5ss_spin();
      break;
    case BM_10_sp3d5ss_nospin:
      print_10_sp3d5ss_nospin();
      break;
    default:
      die("Unidentified Bandstrucutre model in print\n");
    }
  }
  void Deallocate();
  bool isInitialized() { return this->Hdd != NULL; } // a kludge
};

#endif
