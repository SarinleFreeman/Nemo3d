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
$Header: /repo/nemo3d/src/ham/mb_ham_spds.c,v 1.7 2003/10/08 16:17:21 hook Exp $
*****************************************************************************/

#include "mb_ham_spds.h"

#define max(A, B) ((A) > (B) ? (A) : (B))
#define min(A, B) ((A) < (B) ? (A) : (B))

int Make_spd_s_hamiltonian(real *param, real *mb_strain, real Ev_offset,
                           real kxa, real kya, cmatrix d, cmatrix s) {
  /*  Here we just copy the spds* paramters into the matrix d.  Valence
      band offset is added to the onsite energies.  Matrix elements
      are generated during the folding routine.  The spds     */

  int i;
  int nparam = sp3d5s_nparam;
  int nstrain = sp3d5s_nstrain;

  for (i = 0; i < 8; i++)
    d[0][i].r = param[i] + Ev_offset;
  d[0][4].r = d[0][5].r = 0;
  for (i = 8; i < nparam; i++)
    d[0][i].r = param[i];
  if (mb_strain) {
    for (i = 0; i < nstrain; i++)
      d[1][i].r = mb_strain[i];
  }

  return (0);
}

int Make_spds_s_hamiltonian(real *param, real *mb_strain, real Ev_offset,
                            real kxa, real kya, cmatrix d, cmatrix s) {
  /*  Here we just copy the spds* paramters into the matrix d.  Valence
      band offset is added to the onsite energies.  Matrix elements
      are generated during the folding routine.  The spds     */

  int i;
  int nparam = sp3d5s_nparam;
  int nstrain = sp3d5s_nstrain;
  real StrainShift = 0.0;
  int r = 0, c = 0;
#ifdef STRAIN_SHIFT_ENABLE_1D
  if (mb_strain) {
    if (sp3d5s_Ebound == sp3d5s_Ebound_undefine) {
      sp3d5s_Ebound = mb_strain[eta_cat_shift];
      printf("FIX:  currently use cation shift\n");
    }
    if (sp3d5s_Ebound != sp3d5s_Ebound_default) {
      printf("Only support a default value of %f as an overall\n"
             "energy shift for the strain corrections at this time.\n"
             "Adjust the strain parameters accordingly!!!\n",
             sp3d5s_Ebound_default);
      exit(-1);
    }
    StrainShift = sp3d5s_Ebound;
  } else {
    sp3d5s_Ebound = sp3d5s_Ebound_default;
    StrainShift = sp3d5s_Ebound;
  }
#endif

#define NEW_param_distr
#ifdef NEW_param_distr
  /* Chris decided to copy the parameters into the d matrices to conform to the
     old nemo standard of not carrying the paramters around after the
     construction of the d's and s's. So now we have to be careful in reusing
     this for bandstructure models whose basis is smaller than the parameter
     length. assume we have at least the size of sp3d5s_nparam_copy_size
  */
  for (i = 0; i < 8; i++) {
    d[r][c].r = param[i] + Ev_offset - StrainShift;
    c++;
  }
  for (i = 8; i < nparam; i++) {
    if (c == sp3d5s_nparam_copy_size) {
      c = 0;
      r++;
    }
    d[r][c].r = param[i];
    c++;
  }
  if (mb_strain) {
    r = 0;
    c = 0;
    for (i = 0; i < nstrain; i++) {
      if (c == sp3d5s_nparam_copy_size) {
        c = 0;
        r++;
      }
      d[r][c].i = mb_strain[i];
      c++;
    }
  }

#else
  for (i = 0; i < 8; i++)
    d[0][i].r = param[i] + Ev_offset - StrainShift;
  for (i = 8; i < nparam; i++)
    d[0][i].r = param[i];
  if (mb_strain) {
    for (i = 0; i < nstrain; i++)
      d[1][i].r = mb_strain[i];
  }
#endif

  return (0);
}

int Make_spdd_s_hamiltonian(real *param, real *mb_strain, real Ev_offset,
                            real kxa, real kya, cmatrix d, cmatrix s) {
  /*  Here we just copy the spds* paramters into the matrix d.  Valence
      band offset is added to the onsite energies.  Matrix elements
      are generated during the folding routine.  The spds     */

  int i;
  int nparam = sp3d5s_nparam;
  int nstrain = sp3d5s_nstrain;

  printf("Make_spdd_s_hamiltonian\n  Not fixed yet  \n");
  exit(-1);
  for (i = 0; i < 8; i++)
    d[0][i].r = param[i] + Ev_offset;
  for (i = 8; i < nparam; i++)
    d[0][i].r = param[i];
  if (mb_strain) {
    for (i = 0; i < nstrain; i++)
      d[1][i].r = mb_strain[i];
  }

  return (0);
}

/*T_HEADER_SHORT{fold_t_k_spd_nn_s}
  Fold transverse wavevectors into the bulk nearest neighbor
  $spds^*$ Hamiltonian. T*/
void fold_t_k_spd_s(real kxa, real kya, cmatrix d, cmatrix s, cmatrix d_fold,
                    cmatrix s_fold) {
  /*  The spds* Hamiltonian is generated and stored in d_fold and s_fold.
      the parameters are stored in the first row of d. */

  int i, n;
  real a, dx, dy, dz;
  ivectr basis, spin;
  rvectr k, param, mb_strain = NULL;
  real unstrnd_cubic_cell_length = 0.0;
  rmatrix nnv;
  rmatrix nnv0;

  n = rdim(d);
  nnv = Rmatrix(4, 3);
  nnv0 = Rmatrix(4, 3);
  k = Rvectr(3);
  param = Rvectr(sp3d5s_nparam);
  basis = Ivectr(36);
  spin = Ivectr(36);

  for (i = 0; i < 9; i++) {
    spin[i] = 1;
    spin[i + 18] = 1;
    basis[i] = basis[i + 9] = i + 1;
    basis[i + 18] = basis[i + 27] = i + 11;
  }

  /* Initialize parameter vector */
  for (i = 0; i < sp3d5s_nparam; i++)
    param[i] = d[0][i].r;

  if (d[1][2].r != 0.0) {
    {
      int nstrain = sp3d5s_nstrain;
      mb_strain = Rvectr(nstrain);
      for (i = 0; i < nstrain; i++)
        mb_strain[i] = d[1][i].r;
    }
  }

  /* Initialize crystal displacements and lattice constant */

  a = param[31];
  dx = param[32];
  dy = param[33];
  dz = param[34];
  unstrnd_cubic_cell_length = a;

  /* Initialize k-vector */
  k[0] = kxa * 2 * Pi / a;
  k[1] = kya * 2 * Pi / a;

  /* Initialize nearest-neighbor vectors */
  nnv0[0][0] = a / 4;
  nnv0[0][1] = a / 4;
  nnv0[0][2] = a / 4;

  nnv0[1][0] = a / 4;
  nnv0[1][1] = -a / 4;
  nnv0[1][2] = -a / 4;

  nnv0[2][0] = -a / 4;
  nnv0[2][1] = a / 4;
  nnv0[2][2] = -a / 4;

  nnv0[3][0] = -a / 4;
  nnv0[3][1] = -a / 4;
  nnv0[3][2] = a / 4;

  /* Initialize nearest-neighbor vectors */
  nnv[0][0] = dx * a / 4;
  nnv[0][1] = dy * a / 4;
  nnv[0][2] = dz * a / 4;

  nnv[1][0] = dx * a / 4;
  nnv[1][1] = -dy * a / 4;
  nnv[1][2] = -dz * a / 4;

  nnv[2][0] = -dx * a / 4;
  nnv[2][1] = dy * a / 4;
  nnv[2][2] = -dz * a / 4;

  nnv[3][0] = -dx * a / 4;
  nnv[3][1] = -dy * a / 4;
  nnv[3][2] = dz * a / 4;

  /* Generate hamiltonian and store in d_fold and s_fold */
  tbh_const_1d_nemo(d_fold, s_fold, nnv, nnv0, basis, spin, param, mb_strain,
                    unstrnd_cubic_cell_length, k);
  /* tbh_const_1d_nemo_old( d_fold, s_fold, nnv, param, k ); */

  rm_ivectr(&basis);
  rm_ivectr(&spin);
  rm_rvectr(&k);
  rm_rvectr(&param);
  rm_rvectr(&mb_strain);
  rm_rmatrix(&nnv);
  rm_rmatrix(&nnv0);
}

/*T_HEADER_SHORT{fold_t_k_sp3s_s_new}
  Fold transverse wavevectors into the bulk nearest neighbor
  $spds^*$ Hamiltonian. T*/
void fold_t_k_sp3s_s_new(real kxa, real kya, cmatrix d, cmatrix s,
                         cmatrix d_fold, cmatrix s_fold) {
  /*  The spds* Hamiltonian is generated and stored in d_fold and s_fold.
      the parameters are stored in the first row of d. */

  int i, n;
  real a, dx, dy, dz;
  ivectr basis, spin;
  rvectr k, param, mb_strain = NULL;
  real unstrnd_cubic_cell_length = 0.0;
  rmatrix nnv;
  rmatrix nnv0;

  n = rdim(d);
  nnv = Rmatrix(4, 3);
  nnv0 = Rmatrix(4, 3);
  k = Rvectr(3);
  param = Rvectr(sp3d5s_nparam);
  basis = Ivectr(20);
  spin = Ivectr(20);

  /* This really loops over all the orbitals that are included in the model
     here. in mb_ham_spds_micro we define the list of available orbitals STa Sa
     Pxa Pya etc.... Here we will take out only the orbitals STa=0 Sa=1 Pxa=2
     Pya=3 Pza=4  and STc=10 Sc=11 Pxc=12 Pyc=13 Pzc=14  and The order is as
     follows in the blocks: anion up-spin     0- 4 anion down-spin   5- 9 cation
     up-spin   10-14 cation down spin 15-19
   */
  for (i = 0; i < 5; i++) {
    spin[i] = 1;
    spin[i + 10] = 1;
    basis[i] = basis[i + 5] = i;
    basis[i + 10] = basis[i + 15] = i + 10;
  }

#ifdef NEW_param_distr
  {
    int r = 0, c = 0, len = sp3d5s_nparam_copy_size;
    i = 0;
    while (i < sp3d5s_nparam) {
      for (c = 0; c < len; (c++, i++)) {
        param[i] = d[r][c].r;
      }
      r++;
      len = min(sp3d5s_nparam - r * sp3d5s_nparam_copy_size,
                sp3d5s_nparam_copy_size);
    }
    if (d[0][2].i != 0.0 || d[0][5].i != 0.0) {
      int nstrain = sp3d5s_nstrain;
      len = sp3d5s_nparam_copy_size;
      r = 0, c = 0;

      mb_strain = Rvectr(nstrain);
      i = 0;
      while (i < nstrain) {
        for (c = 0; c < len; (c++, i++)) {
          mb_strain[i] = d[r][c].i;
        }
        r++;
        len =
            min(nstrain - r * sp3d5s_nparam_copy_size, sp3d5s_nparam_copy_size);
      }
    }
  }

#else
  for (i = 0; i < sp3d5s_nparam; i++)
    param[i] = d[0][i].r;

  if (d[1][2].r != 0.0) {
    {
      int nstrain = sp3d5s_nstrain;
      mb_strain = Rvectr(nstrain);
      for (i = 0; i < nstrain; i++)
        mb_strain[i] = d[1][i].r;
    }
  }
#endif

  /* Initialize crystal displacements and lattice constant */
  a = param[31];
  dx = param[32];
  dy = param[33];
  dz = param[34];
  unstrnd_cubic_cell_length = a;

  /* Initialize k-vector */
  k[0] = kxa * 2 * Pi / a;
  k[1] = kya * 2 * Pi / a;

  /* Initialize nearest-neighbor vectors */
  nnv0[0][0] = a / 4;
  nnv0[0][1] = a / 4;
  nnv0[0][2] = a / 4;

  nnv0[1][0] = a / 4;
  nnv0[1][1] = -a / 4;
  nnv0[1][2] = -a / 4;

  nnv0[2][0] = -a / 4;
  nnv0[2][1] = a / 4;
  nnv0[2][2] = -a / 4;

  nnv0[3][0] = -a / 4;
  nnv0[3][1] = -a / 4;
  nnv0[3][2] = a / 4;

  /* Initialize nearest-neighbor vectors */
  nnv[0][0] = dx * a / 4;
  nnv[0][1] = dy * a / 4;
  nnv[0][2] = dz * a / 4;

  nnv[1][0] = dx * a / 4;
  nnv[1][1] = -dy * a / 4;
  nnv[1][2] = -dz * a / 4;

  nnv[2][0] = -dx * a / 4;
  nnv[2][1] = dy * a / 4;
  nnv[2][2] = -dz * a / 4;

  nnv[3][0] = -dx * a / 4;
  nnv[3][1] = -dy * a / 4;
  nnv[3][2] = dz * a / 4;

  /* Generate hamiltonian and store in d_fold and s_fold */
  tbh_const_1d_nemo(d_fold, s_fold, nnv, nnv0, basis, spin, param, mb_strain,
                    unstrnd_cubic_cell_length, k);
  /* tbh_const_1d_nemo_old( d_fold, s_fold, nnv, param, k ); */

  rm_ivectr(&basis);
  rm_ivectr(&spin);
  rm_rvectr(&k);
  rm_rvectr(&param);
  rm_rvectr(&mb_strain);
  rm_rmatrix(&nnv);
  rm_rmatrix(&nnv0);
}

/*T_HEADER_SHORT{fold_t_k_spds_nn_s}
  Fold transverse wavevectors into the bulk nearest neighbor
  $spds^*$ Hamiltonian. T*/
void fold_t_k_spds_s(real kxa, real kya, cmatrix d, cmatrix s, cmatrix d_fold,
                     cmatrix s_fold) {
  /*  The spds* Hamiltonian is generated and stored in d_fold and s_fold.
      the parameters are stored in the first row of d. */

  int i, n;
  real a, dx, dy, dz;
  ivectr basis, spin;
  rvectr k, param, mb_strain = NULL;
  real unstrnd_cubic_cell_length = 0.0;
  rmatrix nnv, nnv0;

  n = rdim(d);
  nnv = Rmatrix(4, 3);
  nnv0 = Rmatrix(4, 3);
  k = Rvectr(3);
  param = Rvectr(sp3d5s_nparam);
  basis = Ivectr(40);
  spin = Ivectr(40);

  /* This really loops over all the orbitals that are included in the model
  here. in mb_ham_spds_micro we define the list of available orbitals STa Sa Pxa
  Pya etc.... Here we take all the 20 orbitals over s, p, d, and s* The order is
  as follows in the blocks: anion up-spin     0- 9 anion down-spin  10-19 cation
  up-spin   20-29 cation down spin 30-39
  */
  for (i = 0; i < 10; i++) {
    spin[i] = 1;
    spin[i + 20] = 1;
    basis[i] = basis[i + 10] = i;
    basis[i + 20] = basis[i + 30] = i + 10;
  }

  /* Initialize parameter vector */
#ifdef NEW_param_distr
  {
    int r = 0, c = 0, len = sp3d5s_nparam_copy_size;
    i = 0;
    while (i < sp3d5s_nparam) {
      for (c = 0; c < len; (c++, i++)) {
        param[i] = d[r][c].r;
      }
      r++;
      len = min(sp3d5s_nparam - r * sp3d5s_nparam_copy_size,
                sp3d5s_nparam_copy_size);
    }
    if (d[0][2].i != 0.0 || d[0][5].i != 0.0) {
      int nstrain = sp3d5s_nstrain;
      len = sp3d5s_nparam_copy_size;
      r = 0, c = 0;

      mb_strain = Rvectr(nstrain);
      i = 0;
      while (i < nstrain) {
        for (c = 0; c < len; (c++, i++)) {
          mb_strain[i] = d[r][c].i;
        }
        r++;
        len =
            min(nstrain - r * sp3d5s_nparam_copy_size, sp3d5s_nparam_copy_size);
      }
    }
  }

#else
  for (i = 0; i < sp3d5s_nparam; i++)
    param[i] = d[0][i].r;

  if (d[1][2].r != 0.0) {
    {
      int nstrain = sp3d5s_nstrain;
      mb_strain = Rvectr(nstrain);
      for (i = 0; i < nstrain; i++)
        mb_strain[i] = d[1][i].r;
    }
  }
#endif

  /* Initialize crystal displacements and lattice constant */
  a = param[31];
  dx = param[32];
  dy = param[33];
  dz = param[34];
  unstrnd_cubic_cell_length = a;

  /* Initialize k-vector */
  k[0] = kxa * 2 * Pi / a;
  k[1] = kya * 2 * Pi / a;

  /* Initialize unperturbed nearest-neighbor vectors */
  nnv0[0][0] = a / 4;
  nnv0[0][1] = a / 4;
  nnv0[0][2] = a / 4;

  nnv0[1][0] = a / 4;
  nnv0[1][1] = -a / 4;
  nnv0[1][2] = -a / 4;

  nnv0[2][0] = -a / 4;
  nnv0[2][1] = a / 4;
  nnv0[2][2] = -a / 4;

  nnv0[3][0] = -a / 4;
  nnv0[3][1] = -a / 4;
  nnv0[3][2] = a / 4;

  /* Initialize nearest-neighbor vectors */
  nnv[0][0] = dx * a / 4;
  nnv[0][1] = dy * a / 4;
  nnv[0][2] = dz * a / 4;

  nnv[1][0] = dx * a / 4;
  nnv[1][1] = -dy * a / 4;
  nnv[1][2] = -dz * a / 4;

  nnv[2][0] = -dx * a / 4;
  nnv[2][1] = dy * a / 4;
  nnv[2][2] = -dz * a / 4;

  nnv[3][0] = -dx * a / 4;
  nnv[3][1] = -dy * a / 4;
  nnv[3][2] = dz * a / 4;

  /* Generate hamiltonian and store in d_fold and s_fold */
  tbh_const_1d_nemo(d_fold, s_fold, nnv, nnv0, basis, spin, param, mb_strain,
                    unstrnd_cubic_cell_length, k);
  /* tbh_const_1d_nemo_old( d_fold, s_fold, nnv, param, k ); */

  rm_ivectr(&basis);
  rm_ivectr(&spin);
  rm_rvectr(&k);
  rm_rvectr(&param);
  rm_rvectr(&mb_strain);
  rm_rmatrix(&nnv);
  rm_rmatrix(&nnv0);
}

/*T_HEADER_SHORT{fold_t_k_spdd_nn_s}
  Fold transverse wavevectors into the bulk nearest neighbor
  $spds^*$ Hamiltonian. T*/
void fold_t_k_spdd_s(real kxa, real kya, cmatrix d, cmatrix s, cmatrix d_fold,
                     cmatrix s_fold) {
  /*  The spds* Hamiltonian is generated and stored in d_fold and s_fold.
      the parameters are stored in the first row of d. */

  int i, n;
  real a, dx, dy, dz;
  ivectr basis, spin;
  rvectr k, param, mb_strain = NULL;
  real unstrnd_cubic_cell_length = 0.0;
  rmatrix nnv, nnv0;

  printf("fold_t_k_spdd_s\nNot implemented yet\n");
  exit(-1);
  n = rdim(d);
  nnv = Rmatrix(4, 3);
  nnv0 = Rmatrix(4, 3);
  k = Rvectr(3);
  param = Rvectr(sp3d5s_nparam);
  basis = Ivectr(40);
  spin = Ivectr(40);

  for (i = 0; i < 10; i++) {
    spin[i] = 1;
    spin[i + 20] = 1;
    basis[i] = basis[i + 10] = i;
    basis[i + 20] = basis[i + 30] = i + 10;
  }

  /* Initialize parameter vector */
  for (i = 0; i < sp3d5s_nparam; i++)
    param[i] = d[0][i].r;

  if (d[1][2].r != 0.0) {
    {
      int nstrain = sp3d5s_nstrain;
      mb_strain = Rvectr(nstrain);
      for (i = 0; i < nstrain; i++)
        mb_strain[i] = d[1][i].r;
    }
  }

  /* Initialize crystal displacements and lattice constant */
  a = param[31];
  dx = param[32];
  dy = param[33];
  dz = param[34];
  unstrnd_cubic_cell_length = a;

  /* Initialize k-vector */
  k[0] = kxa * 2 * Pi / a;
  k[1] = kya * 2 * Pi / a;

  /* Initialize nearest-neighbor vectors */
  nnv0[0][0] = a / 4;
  nnv0[0][1] = a / 4;
  nnv0[0][2] = a / 4;

  nnv0[1][0] = a / 4;
  nnv0[1][1] = -a / 4;
  nnv0[1][2] = -a / 4;

  nnv0[2][0] = -a / 4;
  nnv0[2][1] = a / 4;
  nnv0[2][2] = -a / 4;

  nnv0[3][0] = -a / 4;
  nnv0[3][1] = -a / 4;
  nnv0[3][2] = a / 4;

  /* Initialize nearest-neighbor vectors */
  nnv[0][0] = dx * a / 4;
  nnv[0][1] = dy * a / 4;
  nnv[0][2] = dz * a / 4;

  nnv[1][0] = dx * a / 4;
  nnv[1][1] = -dy * a / 4;
  nnv[1][2] = -dz * a / 4;

  nnv[2][0] = -dx * a / 4;
  nnv[2][1] = dy * a / 4;
  nnv[2][2] = -dz * a / 4;

  nnv[3][0] = -dx * a / 4;
  nnv[3][1] = -dy * a / 4;
  nnv[3][2] = dz * a / 4;

  /* Generate hamiltonian and store in d_fold and s_fold */
  tbh_const_1d_nemo(d_fold, s_fold, nnv, nnv0, basis, spin, param, mb_strain,
                    unstrnd_cubic_cell_length, k);
  /* tbh_const_1d_nemo_old( d_fold, s_fold, nnv, param, k ); */

  rm_ivectr(&basis);
  rm_ivectr(&spin);
  rm_rvectr(&k);
  rm_rvectr(&mb_strain);
  rm_rvectr(&param);
  rm_rmatrix(&nnv);
  rm_rmatrix(&nnv0);
}

/* Returns the band energy.  Input vector contains everything necessary to
construct the spds hamiltonian.  It is structured as follows:
p[0] - p[2]       Contains kx, ky, and kz in units of 2Pi/a.
p[3]              Band index.
p[4] - p[38]      spds parameters.

This function is designed to be called by powell() to find band minima. */
real mb_energy_spds(rvectr p, void *data1, void *data2, void *data3) {
  int i, band, idum = 0;
  real a, dx, dy, dz, E;
  ivectr basis, spin;
  rvectr param, mb_strain, k, eigv;
  real unstrnd_cubic_cell_length;
  real *bond;
  rmatrix nnv, nnv0;
  cmatrix h;

  param = (rvectr)(data1);
  mb_strain = (rvectr)(data2);
  bond = (real *)(data3);
  unstrnd_cubic_cell_length = *bond;

  k = Rvectr(3);
  nnv = Rmatrix(4, 3);
  nnv0 = Rmatrix(4, 3);
  basis = Ivectr(40);
  spin = Ivectr(40);
  h = Coperator(40);
  eigv = Rvectr(40);

  for (i = 0; i < 10; i++) {
    spin[i] = 1;
    spin[i + 20] = 1;
    basis[i] = basis[i + 10] = i;
    basis[i + 20] = basis[i + 30] = i + 10;
  }

  a = param[31];
  dx = param[32];
  dy = param[33];
  dz = param[34];
  band = (int)(param[35]);

  for (i = 0; i < 3; i++)
    k[i] = p[i] * 2 * Pi / a;

  nnv0[0][0] = a / 4;
  nnv0[0][1] = a / 4;
  nnv0[0][2] = a / 4;

  nnv0[1][0] = a / 4;
  nnv0[1][1] = -a / 4;
  nnv0[1][2] = -a / 4;

  nnv0[2][0] = -a / 4;
  nnv0[2][1] = a / 4;
  nnv0[2][2] = -a / 4;

  nnv0[3][0] = -a / 4;
  nnv0[3][1] = -a / 4;
  nnv0[3][2] = a / 4;

  nnv[0][0] = dx * a / 4;
  nnv[0][1] = dy * a / 4;
  nnv[0][2] = dz * a / 4;

  nnv[1][0] = dx * a / 4;
  nnv[1][1] = -dy * a / 4;
  nnv[1][2] = -dz * a / 4;

  nnv[2][0] = -dx * a / 4;
  nnv[2][1] = dy * a / 4;
  nnv[2][2] = -dz * a / 4;

  nnv[3][0] = -dx * a / 4;
  nnv[3][1] = -dy * a / 4;
  nnv[3][2] = dz * a / 4;

  tbh_const_bulk(h, nnv, nnv0, basis, spin, param, mb_strain,
                 unstrnd_cubic_cell_length, k);

  eigval_h_full(eigv, h, &idum);

  E = eigv[band];

  printf("\n%1.15e %1.15e %1.15e  %1.15e", k[0], k[1], k[2], E);

  rm_cmatrix(&h);
  rm_rvectr(&eigv);
  rm_rmatrix(&nnv);
  rm_rmatrix(&nnv0);
  rm_ivectr(&basis);
  rm_ivectr(&spin);
  rm_rvectr(&k);

  return E;
}

/* Returns effective masses. */
real mb_newton_spds(rvectr p, rvectr param, rvectr mb_strain,
                    real unstrnd_cubic_cell_length) {
  int i, band, it;
  real a, dx, dy, dz, E, dkmax, fmax;
  ivectr basis, spin;
  rvectr k, eigv, f, dk;
  rmatrix nnv, nnv0, jac, jacd;
  cmatrix h;

  jac = Roperator(3);
  jacd = Roperator(3);
  f = Rvectr(3);
  k = Rvectr(3);
  dk = Rvectr(3);
  nnv = Rmatrix(4, 3);
  nnv0 = Rmatrix(4, 3);
  basis = Ivectr(40);
  spin = Ivectr(40);
  h = Coperator(40);
  eigv = Rvectr(40);

  for (i = 0; i < 10; i++) {
    spin[i] = 1;
    spin[i + 20] = 1;
    basis[i] = basis[i + 10] = i;
    basis[i + 20] = basis[i + 30] = i + 10;
  }
  /* for ( i=0; i < 20; i++ )
     basis[i] = i; */

  a = param[31];
  dx = param[32];
  dy = param[33];
  dz = param[34];
  band = (int)(param[35]);

  for (i = 0; i < 3; i++) {
    k[i] = p[i] * 2 * Pi / a;
  }

  nnv0[0][0] = a / 4;
  nnv0[0][1] = a / 4;
  nnv0[0][2] = a / 4;

  nnv0[1][0] = a / 4;
  nnv0[1][1] = -a / 4;
  nnv0[1][2] = -a / 4;

  nnv0[2][0] = -a / 4;
  nnv0[2][1] = a / 4;
  nnv0[2][2] = -a / 4;

  nnv0[3][0] = -a / 4;
  nnv0[3][1] = -a / 4;
  nnv0[3][2] = a / 4;

  nnv[0][0] = dx * a / 4;
  nnv[0][1] = dy * a / 4;
  nnv[0][2] = dz * a / 4;

  nnv[1][0] = dx * a / 4;
  nnv[1][1] = -dy * a / 4;
  nnv[1][2] = -dz * a / 4;

  nnv[2][0] = -dx * a / 4;
  nnv[2][1] = dy * a / 4;
  nnv[2][2] = -dz * a / 4;

  nnv[3][0] = -dx * a / 4;
  nnv[3][1] = -dy * a / 4;
  nnv[3][2] = dz * a / 4;

  it = 0;
  do {
    it++;
    E = Newton_micro(h, nnv, nnv0, basis, spin, param, mb_strain,
                     unstrnd_cubic_cell_length, k, band, f, jac);
    luDecompositionBackSubstitution(jac, jacd, dk, f);
    dkmax = fmax = -1e100;
    for (i = 0; i < 3; i++) {
      if (fabs(dk[i]) > dkmax)
        dkmax = fabs(dk[i]);
      if (fabs(f[i]) > fmax)
        fmax = fabs(f[i]);
      k[i] += dk[i];
    }
    printf("\niter=%i dkmax = %e  k = %e %e %e  E = %e f= %e", it, dkmax, k[0],
           k[1], k[2], E, fmax);
  } while (dkmax > 1e-8);

  for (i = 0; i < 3; i++)
    p[i] = k[i] * a / 2 / Pi;

  rm_cmatrix(&h);
  rm_rmatrix(&nnv);
  rm_rmatrix(&nnv0);
  rm_ivectr(&basis);
  rm_ivectr(&spin);
  rm_rvectr(&k);
  rm_rvectr(&dk);
  rm_rmatrix(&jac);
  rm_rmatrix(&jacd);
  rm_rvectr(&f);

  return (0);
}

/* Returns effective masses. */
real mb_mstar_spds(rvectr p, rvectr param, rvectr mb_strain,
                   real unstrnd_cubic_cell_length, int xdir, int ydir, int zdir,
                   rmatrix mtens) {
  int i, band;
  real a, dx, dy, dz;
  ivectr basis, spin;
  rvectr k;
  rmatrix nnv, nnv0;
  cmatrix h;

  k = Rvectr(3);
  nnv = Rmatrix(4, 3);
  nnv0 = Rmatrix(4, 3);
  basis = Ivectr(40);
  spin = Ivectr(40);
  h = Coperator(40);

  for (i = 0; i < 10; i++) {
    spin[i] = 1;
    spin[i + 20] = 1;
    basis[i] = basis[i + 10] = i;
    basis[i + 20] = basis[i + 30] = i + 10;
  }

  a = param[31];
  dx = param[32];
  dy = param[33];
  dz = param[34];
  band = (int)(param[35]);

  for (i = 0; i < 3; i++)
    k[i] = p[i] * 2 * Pi / a;

  nnv0[0][0] = a / 4;
  nnv0[0][1] = a / 4;
  nnv0[0][2] = a / 4;

  nnv0[1][0] = a / 4;
  nnv0[1][1] = -a / 4;
  nnv0[1][2] = -a / 4;

  nnv0[2][0] = -a / 4;
  nnv0[2][1] = a / 4;
  nnv0[2][2] = -a / 4;

  nnv0[3][0] = -a / 4;
  nnv0[3][1] = -a / 4;
  nnv0[3][2] = a / 4;

  nnv[0][0] = dx * a / 4;
  nnv[0][1] = dy * a / 4;
  nnv[0][2] = dz * a / 4;

  nnv[1][0] = dx * a / 4;
  nnv[1][1] = -dy * a / 4;
  nnv[1][2] = -dz * a / 4;

  nnv[2][0] = -dx * a / 4;
  nnv[2][1] = dy * a / 4;
  nnv[2][2] = -dz * a / 4;

  nnv[3][0] = -dx * a / 4;
  nnv[3][1] = -dy * a / 4;
  nnv[3][2] = dz * a / 4;

  mtensor(h, nnv, nnv0, basis, spin, param, mb_strain,
          unstrnd_cubic_cell_length, k, band, mtens);

  rm_cmatrix(&h);
  rm_rmatrix(&nnv);
  rm_rmatrix(&nnv0);
  rm_ivectr(&basis);
  rm_ivectr(&spin);
  rm_rvectr(&k);

  return (0);
}

/* Returns effective masses for 111 direction. */
real mb_mstar_spds_111(rvectr p, rvectr param, rvectr mb_strain,
                       real unstrnd_cubic_cell_length, rmatrix mtens) {
  int i, band, idum = 0;
  real a, dx, dy, dz, kdelta, kdivt, kdivl, E0, Ep1, Em1;
  ivectr basis, spin;
  rvectr k, eigv;
  rmatrix nnv, nnv0;
  cmatrix h;

  k = Rvectr(3);
  nnv = Rmatrix(4, 3);
  nnv0 = Rmatrix(4, 3);
  basis = Ivectr(40);
  spin = Ivectr(40);
  h = Coperator(40);
  eigv = Rvectr(40);

  kdelta = 1e-4;
  kdivl = 3 * kdelta * kdelta;
  kdivt = 2 * kdelta * kdelta;

  /* for ( i=0; i < 20; i++ )
     basis[i] = i; */
  for (i = 0; i < 10; i++) {
    spin[i] = 1;
    spin[i + 20] = 1;
    basis[i] = basis[i + 10] = i;
    basis[i + 20] = basis[i + 30] = i + 10;
  }

  a = param[31];
  dx = param[32];
  dy = param[33];
  dz = param[34];
  band = (int)(param[35]);

  for (i = 0; i < 3; i++)
    k[i] = p[i] * 2 * Pi / a;

  nnv0[0][0] = a / 4;
  nnv0[0][1] = a / 4;
  nnv0[0][2] = a / 4;

  nnv0[1][0] = a / 4;
  nnv0[1][1] = -a / 4;
  nnv0[1][2] = -a / 4;

  nnv0[2][0] = -a / 4;
  nnv0[2][1] = a / 4;
  nnv0[2][2] = -a / 4;

  nnv0[3][0] = -a / 4;
  nnv0[3][1] = -a / 4;
  nnv0[3][2] = a / 4;

  nnv[0][0] = dx * a / 4;
  nnv[0][1] = dy * a / 4;
  nnv[0][2] = dz * a / 4;

  nnv[1][0] = dx * a / 4;
  nnv[1][1] = -dy * a / 4;
  nnv[1][2] = -dz * a / 4;

  nnv[2][0] = -dx * a / 4;
  nnv[2][1] = dy * a / 4;
  nnv[2][2] = -dz * a / 4;

  nnv[3][0] = -dx * a / 4;
  nnv[3][1] = -dy * a / 4;
  nnv[3][2] = dz * a / 4;

  /* Offset kx */
  k[0] += 0.005 * 2 * Pi / a / sqrt(2.0);
  k[1] -= 0.005 * 2 * Pi / a / sqrt(2.0);
  tbh_const_bulk(h, nnv, nnv0, basis, spin, param, mb_strain,
                 unstrnd_cubic_cell_length, k);
  eigval_h_full(eigv, h, &idum);
  E0 = eigv[band];
  k[0] -= kdelta;
  k[1] -= kdelta;
  k[2] -= kdelta;
  tbh_const_bulk(h, nnv, nnv0, basis, spin, param, mb_strain,
                 unstrnd_cubic_cell_length, k);
  eigval_h_full(eigv, h, &idum);
  Em1 = eigv[band];
  k[0] += 2 * kdelta;
  k[1] += 2 * kdelta;
  k[2] += 2 * kdelta;
  tbh_const_bulk(h, nnv, nnv0, basis, spin, param, mb_strain,
                 unstrnd_cubic_cell_length, k);
  eigval_h_full(eigv, h, &idum);
  Ep1 = eigv[band];
  mtens[0][0] = 2 * HBAR_MASSFACTOR / ((Em1 - 2 * E0 + Ep1) / kdivl);

  k[0] -= kdelta;
  k[1] -= kdelta;
  k[2] -= kdelta;

  k[0] -= kdelta;
  k[1] += kdelta;
  tbh_const_bulk(h, nnv, nnv0, basis, spin, param, mb_strain,
                 unstrnd_cubic_cell_length, k);
  eigval_h_full(eigv, h, &idum);
  Em1 = eigv[band];
  k[0] += 2 * kdelta;
  k[1] -= 2 * kdelta;
  tbh_const_bulk(h, nnv, nnv0, basis, spin, param, mb_strain,
                 unstrnd_cubic_cell_length, k);
  eigval_h_full(eigv, h, &idum);
  Ep1 = eigv[band];
  mtens[1][1] = mtens[2][2] =
      2 * HBAR_MASSFACTOR / ((Em1 - 2 * E0 + Ep1) / kdivt);

  rm_cmatrix(&h);
  rm_rmatrix(&nnv);
  rm_rmatrix(&nnv0);
  rm_ivectr(&basis);
  rm_ivectr(&spin);
  rm_rvectr(&k);
  rm_rvectr(&eigv);

  return (0);
}
