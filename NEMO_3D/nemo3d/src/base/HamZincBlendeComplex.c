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
$Header: /repo/nemo3d/src/base/HamZincBlendeComplex.c,v 1.13 2006/11/30 22:08:30
nkharche Exp $
*****************************************************************************/

#include "HamZincBlendeComplex.h"
#include "io_utils.h"
#include "realtype.h"

void HamZincBlendeComplex::Initialize(
    BM_Type BandModel, int N_Basis, int Natom, int Natom_surf, int Nout,
    int Nin, bool MagneticFieldOn, /* Magnetic field added to
                                      HamZincBlendeComplex class */
    complex HBxy) {
  nml_memory_report("Before Allocation of Hamiltonian in Hzb.Initialize - "
                    "HamZincBlendeComplex::Initialize");
  this->N_Basis = N_Basis;
  this->BandModel = BandModel;
  this->Natom = Natom;
  this->Natom_surf = Natom_surf;
  this->MagneticFieldOn =
      MagneticFieldOn; /* Magnetic field added to HamZincBlendeComplex class */
  this->HBxy = HBxy;
  double memory_estimate = 0;

  switch (BandModel) {
  case BM_20_sp3d5ss_spin:
    memory_estimate =
        (Natom * N_Basis + Natom_surf * 5) * sizeof(ham_mem_real) +
        Natom * (6 + 2 * N_Basis * N_Basis) * sizeof(ham_mem_complex) +
        Nin * 3 * sizeof(int) + Nout * 3 * sizeof(int);
    /* Additional storage required when Magnetic field is added */
    if (MagneticFieldOn)
      memory_estimate += Natom * 4 * 2 * sizeof(ham_mem_complex);
    /* End: Additional storage required */
    printf("mpi_id=%d memory estimate Hamiltonian storage: %g MB Natom=%d "
           "Natom_surface=%d\n",
           mpi_n3d_id, memory_estimate / 1.0e6, Natom, Natom_surf);
    fflush(stdout);
    this->Hdd = ham_mem_Rvectr(Natom * N_Basis);
    this->Hdu = ham_mem_Cvectr(Natom * 6);
    this->Hds = ham_mem_Rvectr(Natom_surf * 5);
    this->Hu = ham_mem_Cvectr(Natom * N_Basis * N_Basis);
    this->Ho = ham_mem_Cvectr(Natom * N_Basis * N_Basis);
    /* Allocate memory for HphaseIn/Out */
    if (MagneticFieldOn) {
      this->HphaseIn = ham_mem_Cvectr(/*Natom*4*/ Nin);
      this->HphaseOut = ham_mem_Cvectr(/*Natom*4*/ Nout);
    }
    /* End: Allocate memory for HphaseIn/Out */
    break;
  case BM_10_sp3ss_spin:
    memory_estimate =
        (Natom * N_Basis + Natom_surf * 5) * sizeof(ham_mem_real) +
        Natom * (6 + 2 * N_Basis * N_Basis) * sizeof(ham_mem_complex) +
        Nin * 3 * sizeof(int) + Nout * 3 * sizeof(int);
    /* Additional storage required when Magnetic field is added */
    if (MagneticFieldOn)
      memory_estimate += Natom * 4 * 2 * sizeof(complex);
    /* End: Additional storage required */
    printf("mpi_id=%d memory estimate Hamiltonian storage: %g MB Natom=%d "
           "Natom_surface=%d\n",
           mpi_n3d_id, memory_estimate / 1.0e6, Natom, Natom_surf);
    fflush(stdout);
    this->Hdd = ham_mem_Rvectr(Natom * N_Basis);
    this->Hdu = ham_mem_Cvectr(Natom * 6);
    this->Hds = ham_mem_Rvectr(Natom_surf * 5);
    this->Hu = ham_mem_Cvectr(Natom * N_Basis * N_Basis);
    this->Ho = ham_mem_Cvectr(Natom * N_Basis * N_Basis);
    /* Allocate memory for HphaseIn/Out */
    if (MagneticFieldOn) {
      this->HphaseIn = ham_mem_Cvectr(/*Natom*4*/ Nin);
      this->HphaseOut = ham_mem_Cvectr(/*Natom*4*/ Nout);
    }
    /* End: Allocate memory for HphaseIn/Out */
    break;
  case BM_10_sp3d5ss_nospin:
    memory_estimate =
        (Natom * N_Basis + Natom_surf * 6) * sizeof(ham_mem_real) +
        Natom * (2 * N_Basis * N_Basis * 4) * sizeof(ham_mem_complex) +
        Nin * 3 * sizeof(int) + Nout * 3 * sizeof(int);
    /* Additional storage required when Magnetic field is added */
    if (MagneticFieldOn)
      memory_estimate += Natom * 4 * 2 * sizeof(complex);
    /* End: Additional storage required */
    printf("mpi_id=%d memory estimate Hamiltonian storage: %g MB Natom=%d "
           "Natom_surface=%d\n",
           mpi_n3d_id, memory_estimate / 1.0e6, Natom, Natom_surf);
    fflush(stdout);
    this->Hdd = ham_mem_Rvectr(Natom * N_Basis);
    this->Hds = ham_mem_Rvectr(Natom_surf * 6);
    this->Hu = ham_mem_Cvectr(Natom * N_Basis * N_Basis * 4);
    this->Ho = ham_mem_Cvectr(Natom * N_Basis * N_Basis * 4);
    /* Allocate memory for HphaseIn/Out */
    if (MagneticFieldOn) {
      this->HphaseIn = ham_mem_Cvectr(/*Natom*4*/ Nin);
      this->HphaseOut = ham_mem_Cvectr(/*Natom*4*/ Nout);
    }
    /* End: Allocate memory for HphaseIn/Out */
    break;
  case BM_1_s_nospin:
    die("BM_1_s_nospin is not yet implemented "
        "HamZincBlendeComplex::Initialize");
    break;
  default:
    die("Unsupported bandstructure model specified.  Choose one from "
        "Bands_1_s_nospin, Bands_10_sp3ss_spin, Bands_20_sp3d5ss_spin, "
        "Bands_10_sp3d5ss_nospin\n");
    break;
  }

  this->indxHU_nbr = Ivectr(Nin);
  this->indxHU_to_row = Ivectr(Nin);
  this->indxHU_to_col = Ivectr(Nin);
  this->indxHO_nbr = Ivectr(Nout);
  this->indxHO_to_row = Ivectr(Nout);
  this->indxHO_to_col = Ivectr(Nout);

  nml_memory_report("After Allocation of Hamiltonian in Hzb.Initialize - "
                    "HamZincBlendeComplex::Initialize");
}

void HamZincBlendeComplex::setDiag(const cmatrix hd, int atom, int I, int J,
                                   int Zsurf) {
  // diagonal part of diagonal block
  for (int ii = 0; ii < this->N_Basis; ii++)
    this->Hdd[J + ii] = hd[ii][ii].r;

  // p-p interaction of diagonal part of diagonal block
  int offset = 6 * atom;

  switch (this->BandModel) {
  case BM_20_sp3d5ss_spin:
    this->Hdu[offset + 0].r = hd[2][3].r;
    this->Hdu[offset + 0].i = hd[2][3].i;
    this->Hdu[offset + 1].r = hd[12][13].r;
    this->Hdu[offset + 1].i = hd[12][13].i;
    this->Hdu[offset + 2].r = hd[4][12].r;
    this->Hdu[offset + 2].i = hd[4][12].i;
    this->Hdu[offset + 3].r = hd[4][13].r;
    this->Hdu[offset + 3].i = hd[4][13].i;
    this->Hdu[offset + 4].r = hd[2][14].r;
    this->Hdu[offset + 4].i = hd[2][14].i;
    this->Hdu[offset + 5].r = hd[3][14].r;
    this->Hdu[offset + 5].i = hd[3][14].i;
    break;
  case BM_10_sp3ss_spin:
    this->Hdu[offset + 0].r = hd[2][3].r;
    this->Hdu[offset + 0].i = hd[2][3].i;
    this->Hdu[offset + 1].r = hd[7][8].r;
    this->Hdu[offset + 1].i = hd[7][8].i;
    this->Hdu[offset + 2].r = hd[4][7].r;
    this->Hdu[offset + 2].i = hd[4][7].i;
    this->Hdu[offset + 3].r = hd[4][8].r;
    this->Hdu[offset + 3].i = hd[4][8].i;
    this->Hdu[offset + 4].r = hd[2][9].r;
    this->Hdu[offset + 4].i = hd[2][9].i;
    this->Hdu[offset + 5].r = hd[3][9].r;
    this->Hdu[offset + 5].i = hd[3][9].i;
    break;
  case BM_10_sp3d5ss_nospin:
    break; // no spin-orbit coupling
  case BM_1_s_nospin:
    die("BM_1_s_nospin is not yet implemented HamZincBlendeComplex::setDiag");
    break;
  default:
    die("Unsupported bandstructure model specified.  Choose one from "
        "Bands_1_s_nospin, Bands_10_sp3ss_spin, Bands_20_sp3d5ss_spin, "
        "Bands_10_sp3d5ss_nospin\n");
    break;
  }

  // surface-atom dangling-bond shift in the diagonal block
  if (Zsurf >= 0) {
    if (this->BandModel == BM_20_sp3d5ss_spin ||
        this->BandModel == BM_10_sp3ss_spin) {
      int offset = 5 * Zsurf;
      this->Hds[offset + 0] = hd[1][2].r;
      this->Hds[offset + 1] = hd[1][3].r;
      this->Hds[offset + 2] = hd[1][4].r;
      this->Hds[offset + 3] = hd[2][4].r;
      this->Hds[offset + 4] = hd[3][4].r;
    } else if (this->BandModel == BM_10_sp3d5ss_nospin) {
      int offset = 6 * Zsurf;
      this->Hds[offset + 0] = hd[1][2].r;
      this->Hds[offset + 1] = hd[1][3].r;
      this->Hds[offset + 2] = hd[1][4].r;
      this->Hds[offset + 3] = hd[2][3].r;
      this->Hds[offset + 4] = hd[2][4].r;
      this->Hds[offset + 5] = hd[3][4].r;
    }
  }
}

void HamZincBlendeComplex::Deallocate() {
  if (this->Hu) {
    ham_mem_rm_cvectr(&this->Hu);
    this->Hu = NULL;
  }
  if (this->Ho) {
    ham_mem_rm_cvectr(&this->Ho);
    this->Ho = NULL;
  }
  if (this->Hdd) {
    ham_mem_rm_rvectr(&this->Hdd);
    this->Hdd = NULL;
  }
  if (this->Hdu) {
    ham_mem_rm_cvectr(&this->Hdu);
    this->Hdu = NULL;
  }
  if (this->Hds) {
    ham_mem_rm_rvectr(&this->Hds);
    this->Hds = NULL;
  }
  /* Deallocate memory for HphaseIn/Out */
  if (this->HphaseIn) {
    ham_mem_rm_cvectr(&this->HphaseIn);
  }
  if (this->HphaseOut) {
    ham_mem_rm_cvectr(&this->HphaseOut);
  }
  /* Deallocate memory for HphaseIn/Out */
  if (this->indxHU_nbr) {
    rm_ivectr(&this->indxHU_nbr);
    this->indxHU_nbr = NULL;
  }
  if (this->indxHU_to_row) {
    rm_ivectr(&this->indxHU_to_row);
    this->indxHU_to_row = NULL;
  }
  if (this->indxHU_to_col) {
    rm_ivectr(&this->indxHU_to_col);
    this->indxHU_to_col = NULL;
  }
  if (this->indxHO_nbr) {
    rm_ivectr(&this->indxHO_nbr);
    this->indxHO_nbr = NULL;
  }
  if (this->indxHO_to_row) {
    rm_ivectr(&this->indxHO_to_row);
    this->indxHO_to_row = NULL;
  }
  if (this->indxHO_to_col) {
    rm_ivectr(&this->indxHO_to_col);
    this->indxHO_to_col = NULL;
  }
}

void HamZincBlendeComplex::printStorageInfo() {
  if (mpi_n3d_id)
    return; // only master prints this info

  real storage_main = 0.0;

  switch (this->BandModel) {
  case BM_20_sp3d5ss_spin:
    storage_main = sizeof(real) * Natom * (100 * 16 + 32) / 1048576.0;
    break;
  case BM_10_sp3ss_spin:
    storage_main = sizeof(real) * Natom * (25 * 16 + 22) / 1048576.0;
    break;
  case BM_10_sp3d5ss_nospin:
    storage_main = sizeof(real) * Natom * (100 * 16 + 16) / 1048576.0;
    break;
  case BM_1_s_nospin:
    die("Band_1_s_nospin is not yet implemented "
        "HamZincBlendeComplex::printStorageInfo");
    break;
  default:
    die("Unsupported bandstructure model specified.  Choose one from "
        "Bands_1_s_nospin, Bands_10_sp3ss_spin, Bands_20_sp3d5ss_spin, "
        "Bands_10_sp3d5ss_nospin\n");
    break;
  }

  printf("atom count:  total=%d  surface=%d\n", Natom, Natom_surf);
  /* printf("Estimated Hamiltonian Storage %.1fMB (main)\n", storage_main); */
}

void HamZincBlendeComplex::print_20_sp3d5ss_spin() {
  /* print the Hamiltonian to various files */
  int Nb_2 = this->N_Basis / 2;
  int atom, i, j, k;
  FILE *fp;

  fp = fopen("Hdd", "w");
  for (atom = 0; atom < this->Natom; atom++) {
    for (i = 0; i < this->N_Basis; i++) {
      int indx = i + atom * this->N_Basis;
      fprintf(fp, "%d %d %20.12e\n", indx, indx, this->Hdd[indx]);
    }
  }
  fclose(fp);

  fp = fopen("Hdu", "w");
  for (atom = 0; atom < this->Natom; atom++) {
    int indx = 6 * atom;
    int offset = atom * this->N_Basis;

    fprintf(fp, "%d %d %20.12e %20.12e\n", offset + 2, offset + 3,
            this->Hdu[indx].r, this->Hdu[indx].i);
    fprintf(fp, "%d %d %20.12e %20.12e\n", offset + 3, offset + 2,
            this->Hdu[indx].r, -this->Hdu[indx].i);
    indx++;
    fprintf(fp, "%d %d %20.12e %20.12e\n", offset + 12, offset + 13,
            this->Hdu[indx].r, this->Hdu[indx].i);
    fprintf(fp, "%d %d %20.12e %20.12e\n", offset + 13, offset + 12,
            this->Hdu[indx].r, -this->Hdu[indx].i);
    indx++;
    fprintf(fp, "%d %d %20.12e %20.12e\n", offset + 4, offset + 12,
            this->Hdu[indx].r, this->Hdu[indx].i);
    fprintf(fp, "%d %d %20.12e %20.12e\n", offset + 12, offset + 4,
            this->Hdu[indx].r, -this->Hdu[indx].i);
    indx++;
    fprintf(fp, "%d %d %20.12e %20.12e\n", offset + 4, offset + 13,
            this->Hdu[indx].r, this->Hdu[indx].i);
    fprintf(fp, "%d %d %20.12e %20.12e\n", offset + 13, offset + 4,
            this->Hdu[indx].r, -this->Hdu[indx].i);
    indx++;
    fprintf(fp, "%d %d %20.12e %20.12e\n", offset + 2, offset + 14,
            this->Hdu[indx].r, this->Hdu[indx].i);
    fprintf(fp, "%d %d %20.12e %20.12e\n", offset + 14, offset + 2,
            this->Hdu[indx].r, -this->Hdu[indx].i);
    indx++;
    fprintf(fp, "%d %d %20.12e %20.12e\n", offset + 3, offset + 14,
            this->Hdu[indx].r, this->Hdu[indx].i);
    fprintf(fp, "%d %d %20.12e %20.12e\n", offset + 14, offset + 3,
            this->Hdu[indx].r, -this->Hdu[indx].i);
    indx++;
  }
  fclose(fp);

  fp = fopen("Ho", "w");
  for (k = 0; k < nml_iv_extent(this->indxHO_to_row); k++) {
    for (i = 0; i < Nb_2; i++) {
      for (j = 0; j < Nb_2; j++) {
        int indx = Nb_2 * (Nb_2 * k + i) + j;
        fprintf(fp, "%d %d %20.12e %20.12e\n", this->indxHO_to_row[k] + i,
                this->indxHO_to_col[k] + j, this->Ho[indx].r, this->Ho[indx].i);
        fprintf(fp, "%d %d %20.12e %20.12e\n",
                this->indxHO_to_row[k] + i + Nb_2,
                this->indxHO_to_col[k] + j + Nb_2, this->Ho[indx].r,
                this->Ho[indx].i);
      }
    }
  }
  fclose(fp);

  fp = fopen("Hu", "w");
  for (k = 0; k < nml_iv_extent(this->indxHU_to_row); k++) {
    for (i = 0; i < Nb_2; i++) {
      for (j = 0; j < Nb_2; j++) {
        int indx = Nb_2 * (Nb_2 * k + i) + j;
        int I = this->indxHU_to_row[k] + i;
        int J = this->indxHU_to_col[k] + j;
        fprintf(fp, "%d %d %20.12e %20.12e\n", I, J, this->Hu[indx].r,
                this->Hu[indx].i);
        fprintf(fp, "%d %d %20.12e %20.12e\n", J, I, this->Hu[indx].r,
                this->Hu[indx].i);
        fprintf(fp, "%d %d %20.12e %20.12e\n", I + Nb_2, J + Nb_2,
                this->Hu[indx].r, this->Hu[indx].i);
        fprintf(fp, "%d %d %20.12e %20.12e\n", J + Nb_2, I + Nb_2,
                this->Hu[indx].r, this->Hu[indx].i);
      }
    }
  }
  fclose(fp);
}

void HamZincBlendeComplex::print_10_sp3d5ss_nospin() {
  /* print the Hamiltonian to various files */
  const int Nb_2 = this->N_Basis;
  int atom, i, j, k;
  FILE *fp;

  fp = fopen("Hdd", "w");
  for (atom = 0; atom < this->Natom; atom++) {
    for (i = 0; i < this->N_Basis; i++) {
      int indx = i + atom * this->N_Basis;
      fprintf(fp, "%d %d %20.12e\n", indx, indx, this->Hdd[indx]);
    }
  }
  fclose(fp);

  fp = fopen("Ho", "w");
  for (k = 0; k < nml_iv_extent(this->indxHO_to_row); k++) {
    for (i = 0; i < Nb_2; i++) {
      for (j = 0; j < Nb_2; j++) {
        int indx = Nb_2 * (Nb_2 * k + i) + j;
        fprintf(fp, "%d %d %20.12e %20.12e\n", this->indxHO_to_row[k] + i,
                this->indxHO_to_col[k] + j, this->Ho[indx].r, this->Ho[indx].i);
      }
    }
  }
  fclose(fp);

  fp = fopen("Hu", "w");
  for (k = 0; k < nml_iv_extent(this->indxHU_to_row); k++) {
    for (i = 0; i < Nb_2; i++) {
      for (j = 0; j < Nb_2; j++) {
        int indx = Nb_2 * (Nb_2 * k + i) + j;
        int I = this->indxHU_to_row[k] + i;
        int J = this->indxHU_to_col[k] + j;
        fprintf(fp, "%d %d %20.12e %20.12e\n", I, J, this->Hu[indx].r,
                this->Hu[indx].i);
        fprintf(fp, "%d %d %20.12e %20.12e\n", J, I, this->Hu[indx].r,
                this->Hu[indx].i);
      }
    }
  }
  fclose(fp);
}

void HamZincBlendeComplex::print_10_sp3ss_spin() {
  /* print the Hamiltonian to various files */
  const int Nb_2 = this->N_Basis / 2;
  int atom, i, j, k;
  char fil[100];
  FILE *fp;

  sprintf(fil, "Hdd_%d", mpi_n3d_id);
  fp = fopen(fil, "w");
  for (atom = 0; atom < this->Natom; atom++) {
    for (i = 0; i < this->N_Basis; i++) {
      int indx = i + atom * this->N_Basis;
      fprintf(fp, "%d %d %20.12e\n", indx, indx, this->Hdd[indx]);
    }
  }
  fclose(fp);

  sprintf(fil, "Hdu_%d", mpi_n3d_id);
  fp = fopen(fil, "w");
  for (atom = 0; atom < this->Natom; atom++) {
    int indx = 6 * atom;
    int offset = atom * this->N_Basis;

    fprintf(fp, "%d %d %20.12e %20.12e\n", offset + 2, offset + 3,
            this->Hdu[indx].r, this->Hdu[indx].i);
    fprintf(fp, "%d %d %20.12e %20.12e\n", offset + 3, offset + 2,
            this->Hdu[indx].r, -this->Hdu[indx].i);
    indx++;
    fprintf(fp, "%d %d %20.12e %20.12e\n", offset + 7, offset + 8,
            this->Hdu[indx].r, this->Hdu[indx].i);
    fprintf(fp, "%d %d %20.12e %20.12e\n", offset + 8, offset + 7,
            this->Hdu[indx].r, -this->Hdu[indx].i);
    indx++;
    fprintf(fp, "%d %d %20.12e %20.12e\n", offset + 4, offset + 7,
            this->Hdu[indx].r, this->Hdu[indx].i);
    fprintf(fp, "%d %d %20.12e %20.12e\n", offset + 7, offset + 4,
            this->Hdu[indx].r, -this->Hdu[indx].i);
    indx++;
    fprintf(fp, "%d %d %20.12e %20.12e\n", offset + 4, offset + 8,
            this->Hdu[indx].r, this->Hdu[indx].i);
    fprintf(fp, "%d %d %20.12e %20.12e\n", offset + 8, offset + 4,
            this->Hdu[indx].r, -this->Hdu[indx].i);
    indx++;
    fprintf(fp, "%d %d %20.12e %20.12e\n", offset + 2, offset + 9,
            this->Hdu[indx].r, this->Hdu[indx].i);
    fprintf(fp, "%d %d %20.12e %20.12e\n", offset + 9, offset + 2,
            this->Hdu[indx].r, -this->Hdu[indx].i);
    indx++;
    fprintf(fp, "%d %d %20.12e %20.12e\n", offset + 3, offset + 9,
            this->Hdu[indx].r, this->Hdu[indx].i);
    fprintf(fp, "%d %d %20.12e %20.12e\n", offset + 9, offset + 3,
            this->Hdu[indx].r, -this->Hdu[indx].i);
    indx++;
  }
  fclose(fp);

  sprintf(fil, "Ho_%d", mpi_n3d_id);
  fp = fopen(fil, "w");
  for (k = 0; k < nml_iv_extent(this->indxHO_to_row); k++) {
    for (i = 0; i < Nb_2; i++) {
      for (j = 0; j < Nb_2; j++) {
        int indx = Nb_2 * (Nb_2 * k + i) + j;
        fprintf(fp, "%d %d %20.12e %20.12e\n", this->indxHO_to_row[k] + i,
                this->indxHO_to_col[k] + j, this->Ho[indx].r, this->Ho[indx].i);
        fprintf(fp, "%d %d %20.12e %20.12e\n",
                this->indxHO_to_row[k] + i + Nb_2,
                this->indxHO_to_col[k] + j + Nb_2, this->Ho[indx].r,
                this->Ho[indx].i);
      }
    }
  }
  fclose(fp);

  sprintf(fil, "Hu_%d", mpi_n3d_id);
  fp = fopen(fil, "w");
  for (k = 0; k < nml_iv_extent(this->indxHU_to_row); k++) {
    for (i = 0; i < Nb_2; i++) {
      for (j = 0; j < Nb_2; j++) {
        int indx = Nb_2 * (Nb_2 * k + i) + j;
        int I = this->indxHU_to_row[k] + i;
        int J = this->indxHU_to_col[k] + j;
        fprintf(fp, "%d %d %20.12e %20.12e\n", I, J, this->Hu[indx].r,
                this->Hu[indx].i);
        fprintf(fp, "%d %d %20.12e %20.12e\n", J, I, this->Hu[indx].r,
                this->Hu[indx].i);
        fprintf(fp, "%d %d %20.12e %20.12e\n", I + Nb_2, J + Nb_2,
                this->Hu[indx].r, this->Hu[indx].i);
        fprintf(fp, "%d %d %20.12e %20.12e\n", J + Nb_2, I + Nb_2,
                this->Hu[indx].r, this->Hu[indx].i);
      }
    }
  }
  fclose(fp);
}
