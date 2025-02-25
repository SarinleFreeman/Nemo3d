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
$Header: /repo/nemo3d/src/base/ham_local.c,v 1.24 2008/07/23 16:49:14 nkharche
Exp $
*****************************************************************************/

#include "ham_local.h"

#ifndef max
#define max(a, b) ((a) > (b) ? (a) : (b))
#endif

#ifndef min
#define min(a, b) ((a) < (b) ? (a) : (b))
#endif

#ifdef TEST
#include "Atom.h"
void dummy_func {
  Neighborhood cat(Neighborhood::zincblende_cation);
  Neighborhood an(Neighborhood::zincblende_anion);
  cout << cat.NbrPos_a(1, 0) << endl;
  cout << an.NbrPos_a(1, 0) << endl;
}
#endif

void local_Hoff_zincBlende(qd_struct d, cmatrix h, int l, int m,
                           int map_mm_to_n[8][8], int map_mb_to_m[8][4],
                           real a_lattice[3], int block, int nbr5, int nbr6,
                           int nbr7, real phase5, real phase6, real phase7) {
  rvectr nnv = Rvectr(3);
  rvectr nnv0 = Rvectr(3);

  real r5 = cos(Pi * phase5);
  real i5 = sin(Pi * phase5);

  real r6 = cos(Pi * phase6);
  real i6 = sin(Pi * phase6);

  real r7 = cos(Pi * phase7);
  real i7 = sin(Pi * phase7);

  int i = d->geo.cell__ijk[l][0];
  int j = d->geo.cell__ijk[l][1];
  int k = d->geo.cell__ijk[l][2];

  int m_nbr = map_mb_to_m[m][block];
  int n = map_mm_to_n[m][m_nbr];
  int n_nbr = map_mm_to_n[m_nbr][m];

  int i_nbr = i + d->geo.NbrCell(m, n, 0);
  int j_nbr = j + d->geo.NbrCell(m, n, 1);
  int k_nbr = k + d->geo.NbrCell(m, n, 2);

  int l_nbr = d->geo.ijk__cell[i_nbr][j_nbr][k_nbr];

#ifdef ELIMINATE_SSMAP
  if (l_nbr < 0 && (l_nbr + d->geo.N_Cell + 1) >= 0) {
    l_nbr = l_nbr + d->geo.N_Cell + 1;
  }
  if (d->geo.is_tilted)
    l_nbr = d->geo.get_cindx_tilted(l, m, n);

#else  /* ELIMINATE_SSMAP */
  if (l_nbr < 0 && d->geo.ssmap[l][m][n] != -1) {
    l_nbr = d->geo.ssmap[l][m][n];
  }
#endif /* ELIMINATE_SSMAP */

  /* correct cation (l,m) diagonal energies (first zero hd) */
  for (int r = 0; r < d->NBasisStates; r++) {
    for (int c = 0; c < d->NBasisStates; c++) {
      d->hd[r][c].r = d->hd[r][c].i = 0.0;
    }
  }
  H_DiagStrainCorr_OffDiag_Wrapper(d, l, m, (int)d->geo.AtomType[l][m], n,
                                   l_nbr, m_nbr, nnv, nnv0, a_lattice);
  for (int ii = 4; ii < 8; ii++) {
    for (int r = 0; r < d->NBasisStates; r++) {
      for (int c = 0; c < d->NBasisStates; c++) {
        h[r + ii * d->NBasisStates][c + ii * d->NBasisStates].r +=
            d->hd[r][c].r;
        h[r + ii * d->NBasisStates][c + ii * d->NBasisStates].i +=
            d->hd[r][c].i;
      }
    }
  }

  /* compute cation-anion interaction (l,m)-(l_nbr,m_nbr) */
  for (int r = 0; r < d->NBasisStates; r++) {
    for (int c = 0; c < d->NBasisStates; c++) {
      int cp, rp;
      real hr = d->ho[r][c].r;
      real hi = d->ho[r][c].i;

      cp = 4 * d->NBasisStates + c;
      rp = block * d->NBasisStates + r;
      h[rp][cp].r = hr;
      h[rp][cp].i = hi;
      h[cp][rp].r = h[rp][cp].r;
      h[cp][rp].i = -h[rp][cp].i;

      cp = 5 * d->NBasisStates + c;
      rp = nbr5 * d->NBasisStates + r;
      h[rp][cp].r = hr * r5 - hi * i5;
      h[rp][cp].i = hi * r5 + hr * i5;
      h[cp][rp].r = h[rp][cp].r;
      h[cp][rp].i = -h[rp][cp].i;

      cp = 6 * d->NBasisStates + c;
      rp = nbr6 * d->NBasisStates + r;
      h[rp][cp].r = hr * r6 - hi * i6;
      h[rp][cp].i = hi * r6 + hr * i6;
      h[cp][rp].r = h[rp][cp].r;
      h[cp][rp].i = -h[rp][cp].i;

      cp = 7 * d->NBasisStates + c;
      rp = nbr7 * d->NBasisStates + r;
      h[rp][cp].r = hr * r7 - hi * i7;
      h[rp][cp].i = hi * r7 + hr * i7;
      h[cp][rp].r = h[rp][cp].r;
      h[cp][rp].i = -h[rp][cp].i;
    }
  }

  /* correct anion (l_nbr, m_nbr) diagonal energies  (first zero hd) */
  for (int r = 0; r < d->NBasisStates; r++) {
    for (int c = 0; c < d->NBasisStates; c++) {
      d->hd[r][c].r = d->hd[r][c].i = 0.0;
    }
  }
  H_DiagStrainCorr_OffDiag_Wrapper(d, l_nbr, m_nbr,
                                   (int)d->geo.AtomType[l_nbr][m_nbr], n_nbr, l,
                                   m, nnv, nnv0, a_lattice);
  for (int ii = 0; ii < 4; ii++) {
    for (int r = 0; r < d->NBasisStates; r++) {
      for (int c = 0; c < d->NBasisStates; c++) {
        h[r + ii * d->NBasisStates][c + ii * d->NBasisStates].r +=
            d->hd[r][c].r;
        h[r + ii * d->NBasisStates][c + ii * d->NBasisStates].i +=
            d->hd[r][c].i;
      }
    }
  }

  rm_rvectr(&nnv);
  rm_rvectr(&nnv0);
}

/* need to assume S^3 symmetry for now */
static void compute_local_hamiltonian_at_cation(cmatrix h, qd_struct d, int l,
                                                int m, real a_lattice[3],
                                                int map_mm_to_n[8][8],
                                                int map_mb_to_m[8][4]) {
  int AtomType_this = (int)d->geo.AtomType[l][m];

  int i = d->geo.cell__ijk[l][0];
  int j = d->geo.cell__ijk[l][1];
  int k = d->geo.cell__ijk[l][2];

  int sz = 8 * d->NBasisStates;

  /* zero h */
  for (int r = 0; r < sz; r++) {
    for (int c = 0; c < sz; c++) {
      h[r][c].r = h[r][c].i = 0.0;
    }
  }

  /* there is only one cation in the primitive tetrahedral cell.
     The cubic cell contains three other copies of the same
     cation.  Thus, the diagonal blocks h(4,4), h(5,5),
     h(6,6), h(7,7) must be identical. */
  H_Diag_Wrapper(d, l, m, i, j, k, AtomType_this, d->BasisCation);
  for (int block = 4; block < 8; block++) {
    int offset = block * d->NBasisStates;
    for (int r = 0; r < d->NBasisStates; r++) {
      for (int c = 0; c < d->NBasisStates; c++) {
        complex *h_r1_c1 = &h[r + offset][c + offset];
        h_r1_c1->r = d->hd[r][c].r;
        h_r1_c1->i = d->hd[r][c].i;
      }
    }
  }

  /* fill anion diagonal blocks:  h(0,0), h(1,1), h(2,2), h(3,3)
   */
  for (int block = 0; block < 4; block++) {
    int offset = block * d->NBasisStates;

    int m_nbr = map_mb_to_m[m][block];
    int n = map_mm_to_n[m][m_nbr];

    int i_nbr = i + d->geo.NbrCell(m, n, 0);
    int j_nbr = j + d->geo.NbrCell(m, n, 1);
    int k_nbr = k + d->geo.NbrCell(m, n, 2);

    /* Get the unit cell and atomic index for this neighbor */
    int l_nbr = d->geo.ijk__cell[i_nbr][j_nbr][k_nbr];

#ifdef ELIMINATE_SSMAP
    if (l_nbr < 0 && (l_nbr + d->geo.N_Cell + 1) >= 0) {
      l_nbr = l_nbr + d->geo.N_Cell + 1;
    }
    if (d->geo.is_tilted)
      l_nbr = d->geo.get_cindx_tilted(l, m, n);
#else  /* ELIMINATE_SSMAP */
    if (l_nbr < 0 && d->geo.ssmap[l][m][n] != -1) {
      l_nbr = d->geo.ssmap[l][m][n];
    }
#endif /* ELIMINATE_SSMAP */

    int AtomType_nbr = (int)d->geo.AtomType[l_nbr][m_nbr];

    /* reassign {ijk}_nbr since l_nbr may have changed */
    /*
    i_nbr = d->geo.cell__ijk[l_nbr][0];
    j_nbr = d->geo.cell__ijk[l_nbr][1];
    k_nbr = d->geo.cell__ijk[l_nbr][2];
    */

    /* for local bandstruct version, all cation neighbors of
       this anion are of the same species within this virtual
       primitive cell.  Also, for now, assume periodic BC so
       that we don't have to worry about shifting energies at
       dangling bonds. */
    for (int p = 0; p < sp3d5s_nparam; p++)
      d->param[p] = d->parmat[AtomType_nbr][AtomType_this][p];

    h3d_diag(d->hd, d->BasisAnion, d->Basis_Spin, d->param, d->NBasisStates);

    for (int r = 0; r < d->NBasisStates; r++) {
      for (int c = 0; c < d->NBasisStates; c++) {
        h[r + offset][c + offset].r = d->hd[r][c].r;
        h[r + offset][c + offset].i = d->hd[r][c].i;
      }
    }

#ifdef STRAIN_SHIFT_ENABLE
    if (d->opt.ExecParam.ElCalc.include_strain_in_Hamiltonian) {
      real StrainShift = d->strnmat[AtomType_this][AtomType_nbr][eta_an_shift];
      for (int r = 0; r < d->NBasisStates; r++) {
        int diag = r + block * d->NBasisStates;
        h[diag][diag].r += StrainShift;
      }
    }
#endif
  }

  /* do off-diagonal elements h(*,1) and h(1,*) coupling
     cation (l,m) with its neighoring anions */
  local_Hoff_zincBlende(d, h, l, m, map_mm_to_n, map_mb_to_m, a_lattice, 0, 2,
                        1, 3, d->kxL + d->kyL, d->kxL + d->kzL,
                        d->kyL + d->kzL);
  local_Hoff_zincBlende(d, h, l, m, map_mm_to_n, map_mb_to_m, a_lattice, 1, 3,
                        0, 2, d->kyL, 0.0, d->kyL);
  local_Hoff_zincBlende(d, h, l, m, map_mm_to_n, map_mb_to_m, a_lattice, 2, 0,
                        3, 1, 0.0, d->kzL, d->kzL);
  local_Hoff_zincBlende(d, h, l, m, map_mm_to_n, map_mb_to_m, a_lattice, 3, 1,
                        2, 0, d->kxL, d->kxL, 0.0);
}

/* for now assume zinc-blende structure with cubic basis */
/* also, for now, just compute E(k=0) */
void local_bandstruct_zincBlende_cubic(qd_struct d) {
  MPI_TIME_INIT(Tstart);

  int Ncat = 0;
  int isx = d->cell_s[mpi_n3d_id];
  int iex = isx + d->cell_ln[mpi_n3d_id];
  int report_n = 0;
  int atom_n = 0;
  int atom_total = 0;

  int sz_ham = d->NBasisStates * 8;
  cmatrix h = Coperator(sz_ham);
  cmatrix h_vca = Coperator(sz_ham);
  cvectr E = Cvectr(sz_ham);
  cvectr E_vca = Cvectr(sz_ham);

  int map_mm_to_n[8][8];
  int map_mb_to_m[8][4]; /* map (cation m, block label of nbr) to m of nbr */

  char filetag[100];
  FILE *fpLocal = NULL;

  real a_lattice[3];

  a_lattice[0] = d->geo.lattice_x;
  a_lattice[1] = d->geo.lattice_y;
  a_lattice[2] = d->geo.lattice_z;

  /* need map from nbr indx (n) to atom indx (m) */
  for (int m = 0; m < d->geo.AtomsPerCellMax(); m++) {
    for (int n = 0; n < d->geo.Neighbors(m); n++) {
      int m_nbr = d->geo.NbrCell(m, n, 3);
      map_mm_to_n[m][m_nbr] = n;
    }
  }

  for (int m = 0; m < 4; m++) {
    map_mb_to_m[m][0] = 4 + m;
    map_mb_to_m[m][1] = 4 + (7 - m) % 4;
    map_mb_to_m[m][2] = 4 + (2 + m) % 4;
    map_mb_to_m[m][3] = 4 + (5 - m) % 4;
  }
  for (int m = 4; m < 8; m++) {
    map_mb_to_m[m][0] = (m - 2) % 4;
    map_mb_to_m[m][1] = (9 - m) % 4;
    map_mb_to_m[m][2] = m % 4;
    map_mb_to_m[m][3] = (7 - m) % 4;
  }

  /* zero h_vca */
  for (int r = 0; r < sz_ham; r++) {
    for (int c = 0; c < sz_ham; c++) {
      h_vca[r][c].r = h_vca[r][c].i = 0.0;
    }
  }

  if (d->_bandstruct_on_local_domain) {
    char *filename = n3d_strdup_n(d->inputfile);
    sprintf(filetag, "nd_Ek_proc_%d_kx%.3f_ky_%.3f_kz_%.3f_r", mpi_n3d_id,
            d->kxL, d->kyL, d->kzL);
    n3d_FileTypeSet(&filename, filetag, TRUE);

    if ((fpLocal = fopen(filename, "w")) == NULL)
      die("Could not open the file \"%s\"", filename);
    str_free(filename);
  }

  atom_total = (iex - isx) * d->geo.AtomsPerCellMax() / 2;

  report_n = max(1, atom_total / 100000);
  atom_n = 0;

  // Loop over all the cations in each unit cell belonging to this processor
  for (int l = isx; l < iex; l++) {
    for (int m = 0; m < d->geo.AtomsPerCellMax() / 2; m++) {

      atom_n++;
      if (!(atom_n % report_n))
        printf("#%d:  %d / %d\n", mpi_n3d_id, atom_n, atom_total);

      double x = d->geo.lattice_x *
                     (d->geo.cell__ijk[l][0] + d->geo.PositionInCell(m, 0)) +
                 d->geo.l_m_2_atomdsp[l][m][0];
      double y = d->geo.lattice_y *
                     (d->geo.cell__ijk[l][1] + d->geo.PositionInCell(m, 1)) +
                 d->geo.l_m_2_atomdsp[l][m][1];
      double z = d->geo.lattice_z *
                     (d->geo.cell__ijk[l][2] + d->geo.PositionInCell(m, 2)) +
                 d->geo.l_m_2_atomdsp[l][m][2];

      /*  Check if atom is not in simulation domain, skip */
      if (!d->geo.AtomType[l][m])
        continue;

      /* only look at those cations within this box */
      if (x < d->opt.ExecParam.Subdomain.xMin ||
          x > d->opt.ExecParam.Subdomain.xMax)
        continue;
      if (y < d->opt.ExecParam.Subdomain.yMin ||
          y > d->opt.ExecParam.Subdomain.yMax)
        continue;
      if (z < d->opt.ExecParam.Subdomain.zMin ||
          z > d->opt.ExecParam.Subdomain.zMax)
        continue;

      Ncat++;

      // we are lazy; toss out if we or neighbor are on boundary
      bool flag = true;
      {
        int i = d->geo.cell__ijk[l][0];
        int j = d->geo.cell__ijk[l][1];
        int k = d->geo.cell__ijk[l][2];

        if (d->geo.isOnSurface(l, m))
          flag = false;
        else {
          for (int n = 0; n < d->geo.Neighbors(m); n++) {
            int id = d->geo.NbrCell(m, n, 0);
            int jd = d->geo.NbrCell(m, n, 1);
            int kd = d->geo.NbrCell(m, n, 2);
            if (d->geo.ijk__cell[i + id][j + jd][k + kd] < 0) {
              flag = false;
            }
          }
        }
      }

      if (flag) {
        compute_local_hamiltonian_at_cation(h, d, l, m, a_lattice, map_mm_to_n,
                                            map_mb_to_m);
      }

      if (d->opt.ExecParam.ElCalc.ElOut.LocalHam) {
        FILE *fpHam;
        char *filename = n3d_strdup_n(d->inputfile);

        n3d_FileTypeSet(&filename, "nd_ham_k_r", TRUE);
        fpHam = fopen(filename, "a");

        for (int r = 0; r < sz_ham; r++) {
          for (int c = 0; c < sz_ham; c++) {
            if (h[r][c].r != 0.0 || h[r][c].i != 0.0) {
              fprintf(fpHam, "%d %d %e %e\n", r, c, h[r][c].r, h[r][c].i);
            }
          }
        }
        fclose(fpHam);
        str_free(filename);
      }

      if (d->_bandstruct_on_local_domain && fpLocal) {

        ceigval_ns_full(E, h);

        for (int i = 0; i < sz_ham; i++) {
          fprintf(fpLocal, "%e %e %e %e %e %e %e   %d\n", d->kxL, d->kyL,
                  d->kzL, x, y, z, E[i].r, (int)d->geo.AtomType[l][m]);
        }
      }

      if (d->_bandstruct_on_full_domain) {
        for (int r = 0; r < sz_ham; r++) {
          for (int c = 0; c < sz_ham; c++) {
            h_vca[r][c].r += h[r][c].r;
            h_vca[r][c].i += h[r][c].i;
          }
        }
      }
    }
  }

  if (d->_bandstruct_on_local_domain && fpLocal)
    fclose(fpLocal);

  /* compute vca */
  if (d->_bandstruct_on_full_domain) {
    int Ncat_tot;
    cmatrix h_vca_tot = Coperator(sz_ham);

#if (defined MPI3d && !defined FAKE_MPI)
    MPI_Allreduce(&h_vca[0][0], &h_vca_tot[0][0], 2 * sz_ham * sz_ham,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&Ncat, &Ncat_tot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
    h_vca_tot[0][0] = h_vca[0][0];
    Ncat_tot = Ncat;
#endif

    if (mpi_n3d_id == mpi_n3d_masterid) {
      FILE *fp_vca;
      char *filename = n3d_strdup_n(d->inputfile);
      ceigval_ns_full(E_vca, h_vca_tot);

      sprintf(filetag, "nd_Ek_vca_kx=%.3f_ky=%.3f_kz=%.3f", d->kxL, d->kyL,
              d->kzL);
      n3d_FileTypeSet(&filename, filetag, TRUE);

      if ((fp_vca = fopen(filename, "w")) == NULL) {
        die("Could not open the file \"%s\"", filename);
      }
      str_free(filename);

      for (int i = 0; i < sz_ham; i++)
        fprintf(fp_vca, "%e %e %e %e\n", d->kxL, d->kyL, d->kzL,
                E_vca[i].r / Ncat_tot);
      fclose(fp_vca);
    }
    rm_cmatrix(&h_vca_tot);
  }

  rm_cmatrix(&h);
  rm_cmatrix(&h_vca);
  rm_cvectr(&E);
  rm_cvectr(&E_vca);
}

void compute_local_hamiltonian_at_cell(cmatrix h, qd_struct d, int l) {
  for (int r = 0; r < 20; r++) {
    for (int c = 0; c < 20; c++) {
      d->hd[r][c].r = d->hd[r][c].i = 0.0;
    }
  }

  /* These are real space indexes for the unit cell */
  int i = d->geo.cell__ijk[l][0];
  int j = d->geo.cell__ijk[l][1];
  int k = d->geo.cell__ijk[l][2];

  /* Loop through all the atoms in the unit cell */
  for (int m = 0; m < d->geo.AtomsPerCellMax(); m++) {

    for (int r = 0; r < 20; r++) {
      for (int c = 0; c < 20; c++) {
        d->ho[r][c].r = d->ho[r][c].i = 0.0;
      }
    }

    int I = 20 * m;
    int AtomType_l_m = (int)d->geo.AtomType[l][m];
    /*  Check if atom is not in simulation domain, skip */
    if (!AtomType_l_m)
      continue;

    /* Hamiltonian diagonal (before strain correction).
       If anion send anion basis vector, otherwise use
       cation vector . */

    //==========================================================
    //  was H_diag_Wrapper
    //==========================================================

    ivectr Basis = (d->geo.isCation(l, m) ? d->BasisCation : d->BasisAnion);

    // I don't think this averaging is necessary for the
    // onsite terms  -- F.O.
    /* We need a weighted average of the orbital energies */
    for (int p = 0; p < sp3d5s_nparam; p++)
      d->param[p] = 0;

    real strainShift = 0.0;
    int nbrcnt = 0;
    for (int n = 0; n < d->geo.Neighbors(m); n++) {
      int id = d->geo.NbrCell(m, n, 0);
      int jd = d->geo.NbrCell(m, n, 1);
      int kd = d->geo.NbrCell(m, n, 2);
      int aindx = d->geo.NbrCell(m, n, 3);
#ifdef USE_ACTUAL_CELL
      int cindx = d->geo.ijk__cell[i + id][j + jd][k + kd];

#ifdef ELIMINATE_SSMAP
      if (cindx < 0 && (cindx + d->geo.N_Cell + 1) >= 0) {
        cindx = cindx + d->geo.N_Cell + 1;
      }
      if (d->geo.is_tilted)
        cindx = d->geo.get_cindx_tilted(l, m, n);
#else  /* ELIMINATE_SSMAP */
      if (cindx < 0 && d->geo.ssmap[l][m][n] != -1) {
        cindx = d->geo.ssmap[l][m][n];
      }
#endif /* ELIMINATE_SSMAP */

#else
      int cindx = l;
#endif

      if (cindx >= 0) {
        /* The new data storage can crash on the dereferencing of the casted
         * data */
        /* It appears to crash on the new gcc compilers                   */
        /* const int& AtomType_nbr = (int) d->geo.AtomType[cindx][aindx]; */
        int AtomType_nbr = (int)d->geo.AtomType[cindx][aindx];
        if (!AtomType_nbr)
          continue;

        int cation, anion;
        if (d->geo.isCation(l, m)) {
          cation = AtomType_l_m;
          anion = AtomType_nbr;
#ifdef STRAIN_SHIFT_ENABLE
          if (d->opt.ExecParam.ElCalc.include_strain_in_Hamiltonian)
            strainShift += d->strnmat[cation][anion][eta_cat_shift];
#endif
        } else {
          cation = AtomType_nbr;
          anion = AtomType_l_m;
#ifdef STRAIN_SHIFT_ENABLE
          if (d->opt.ExecParam.ElCalc.include_strain_in_Hamiltonian)
            strainShift += d->strnmat[cation][anion][eta_an_shift];
#endif
        }

        rvectr param_local = d->parmat[cation][anion];
        for (int p = 0; p < sp3d5s_nparam; p++)
          d->param[p] += param_local[p];
        nbrcnt++;
      }
    }
    strainShift /= nbrcnt;

    for (int p = 0; p < sp3d5s_nparam; p++)
      d->param[p] /= nbrcnt;

    for (int p = 0; p < sp3d5s_nparam; p++) {
      for (int n = 0; n < d->geo.Neighbors(m); n++) {
        int id = d->geo.NbrCell(m, n, 0);
        int jd = d->geo.NbrCell(m, n, 1);
        int kd = d->geo.NbrCell(m, n, 2);
        int aindx = d->geo.NbrCell(m, n, 3);
#ifdef USE_ACTUAL_CELL
        int cindx = d->geo.ijk__cell[i + id][j + jd][k + kd];

#ifdef ELIMINATE_SSMAP
        if (cindx < 0 && (cindx + d->geo.N_Cell + 1) >= 0) {
          cindx = cindx + d->geo.N_Cell + 1;
        }
        if (d->geo.is_tilted)
          cindx = d->geo.get_cindx_tilted(l, m, n);

#else  /* ELIMINATE_SSMAP */
        if (cindx < 0 && d->geo.ssmap[l][m][n] != -1) {
          cindx = d->geo.ssmap[l][m][n];
        }
#endif /* ELIMINATE_SSMAP */

#else
        int cindx = l;
#endif
        if (cindx >= 0) {
          /* The new data storage can crash on the dereferencing of the casted
           * data */
          /* It appears to crash on the new gcc compilers                   */
          /* const int& AtomType_nbr = (int) d->geo.AtomType[cindx][aindx]; */
          int AtomType_nbr = (int)d->geo.AtomType[cindx][aindx];
          if (!AtomType_nbr)
            continue;

          int cation, anion;
          if (d->geo.isCation(l, m)) {
            cation = AtomType_l_m;
            anion = AtomType_nbr;
          } else {
            cation = AtomType_nbr;
            anion = AtomType_l_m;
          }
        }
      }
    }

    h3d_diag(d->hd, Basis, d->Basis_Spin, d->param, d->NBasisStates);

    if (d->geo.isCellWurtzite()) {
      double pz_addition_energy = -0.04;
      if (pz_addition_energy != 0.0)
        h3d_diag_pz_addition_for_wurtzite_crystal_splitting(d,
                                                            pz_addition_energy);
    }

    if (d->geo.isOnSurface(l, m)) {
      if (d->opt.ExecParam.Phys.SurfacePassivation ==
          Phys_struct::PassivateOrbital)
        h3d_diag_surface_atom_shift(d);
      if (d->opt.ExecParam.Phys.SurfacePassivation ==
          Phys_struct::PassivateDanglingBond)
        h3d_diag_dangling_bond_shift(d, l, m);
    }

#ifdef STRAIN_SHIFT_ENABLE
    for (int ii = 0; ii < d->NBasisStates; ii++) {
      d->hd[ii][ii].r += strainShift;
    }
#endif

    //==========================================================
    //  end of H_diag_Wrapper
    //==========================================================

    /* Hamiltonian off-diagonal:  Here we loop through all
       this atom's nearest neighbors, get the coupling
       hamiltonians */
    for (int n = 0; n < d->geo.Neighbors(m); n++) {

      int J = 20 * d->geo.NbrCell(m, n, 3);
      complex eik;
      real nnv[3];
      real nnv0[3];

      // Get relative position of unit cell where this nbr lives.
      int i_nbr = i + d->geo.NbrCell(m, n, 0);
      int j_nbr = j + d->geo.NbrCell(m, n, 1);
      int k_nbr = k + d->geo.NbrCell(m, n, 2);

      // Get the unit cell and atomic index for this neighbor
      int aindx = d->geo.NbrCell(m, n, 3);
      int cindx = d->geo.ijk__cell[i_nbr][j_nbr][k_nbr];
      if (cindx != l) {
        real phase =
            d->geo.phasePeriodic(d->kxL, d->kyL, d->kzL, i_nbr, j_nbr, k_nbr);

        eik.r = cos(phase);
        eik.i = sin(phase);
        cindx = l;
      } else {
        eik.r = 1.0;
        eik.i = 0.0;
      }

      // if nbr's cell is not in simulation domain, skip it
      if (!d->geo.AtomType[cindx][aindx])
        continue;

      // ==================================================
      // was H_DiagStrainCorr_OffDiag_Wrapper
      // ==================================================

      /* The new data storage can crash on the dereferencing of the casted data
       */
      /* It appears to crash on the new gcc compilers                   */
      /* const int& AtomType_nbr = (int) d->geo.AtomType[cindx][aindx];*/
      int AtomType_nbr = (int)d->geo.AtomType[cindx][aindx];

      real a_lattice_unstrnd =
          d->parmat[AtomType_l_m][AtomType_nbr][pV_unstr_latt];

      d->geo.getRelativePosVec(nnv0, nnv, a_lattice_unstrnd, cindx, aindx, l, m,
                               n);

      ivectr basis_this, basis_nbr;
      rvectr mb_strain = NULL, param;
      if (d->geo.isCation(l, m)) {
        basis_this = d->BasisCation;
        basis_nbr = d->BasisAnion;
        param = d->parmat[AtomType_l_m][AtomType_nbr];
        if (d->opt.ExecParam.ElCalc.include_strain_in_Hamiltonian)
          mb_strain = d->strnmat[AtomType_l_m][AtomType_nbr];

        for (int p = 0; p < 3; p++) {
          nnv[p] *= -1;
          nnv0[p] *= -1;
        }
      } else {
        basis_this = d->BasisAnion;
        basis_nbr = d->BasisCation;
        param = d->parmat[AtomType_nbr][AtomType_l_m];
        if (d->opt.ExecParam.ElCalc.include_strain_in_Hamiltonian)
          mb_strain = d->strnmat[AtomType_nbr][AtomType_l_m];
      }

      /* Generate off-diagonal Hamiltonian */
      /* returns correct d->ho (original contents destroyed),
         d->hd (new) = d->hd (old) + update */

      h3d_offdiag_diag_corr(d->ho, d->hd, nnv, nnv0, basis_this, basis_nbr,
                            d->Basis_Spin, param, mb_strain, a_lattice_unstrnd,
                            d->NBasisStates);

      // ==================================================
      // end of H_DiagStrainCorr_OffDiag_Wrapper
      // ==================================================

      /* Matrix-vector multiply (same proc, full store) */
      for (int r = 0; r < 20; r++) {
        for (int c = 0; c < 20; c++) {
          h[J + r][I + c].r = d->ho[r][c].r;
          h[J + r][I + c].i = d->ho[r][c].i;
        }
      }
    }

    for (int r = 0; r < 20; r++) {
      for (int c = 0; c < 20; c++) {
        h[I + r][I + c].r = d->hd[r][c].r;
        h[I + r][I + c].i = d->hd[r][c].i;
      }
    }
  }
}

void local_bandstruct_cubic(qd_struct d) {
  MPI_TIME_INIT(Tstart);

  int Ncat = 0;
  int isx = d->cell_s[mpi_n3d_id];
  int iex = isx + d->cell_ln[mpi_n3d_id];
  int atom_total = 0;

  int sz_ham = d->NBasisStates * 8;
  cmatrix h = Coperator(sz_ham);
  cmatrix h_vca = Coperator(sz_ham);
  cvectr E = Cvectr(sz_ham);
  cvectr E_vca = Cvectr(sz_ham);

  int map_mm_to_n[8][8];
  int map_mb_to_m[8][4]; /* map (cation m, block label of nbr) to m of nbr */

  char filetag[100];
  FILE *fpLocal = NULL;

  real a_lattice[3];

  a_lattice[0] = d->geo.lattice_x;
  a_lattice[1] = d->geo.lattice_y;
  a_lattice[2] = d->geo.lattice_z;

  /* need map from nbr indx (n) to atom indx (m) */
  for (int m = 0; m < d->geo.AtomsPerCellMax(); m++) {
    for (int n = 0; n < d->geo.Neighbors(m); n++) {
      int m_nbr = d->geo.NbrCell(m, n, 3);
      map_mm_to_n[m][m_nbr] = n;
    }
  }

  for (int m = 0; m < 4; m++) {
    map_mb_to_m[m][0] = 4 + m;
    map_mb_to_m[m][1] = 4 + (7 - m) % 4;
    map_mb_to_m[m][2] = 4 + (2 + m) % 4;
    map_mb_to_m[m][3] = 4 + (5 - m) % 4;
  }
  for (int m = 4; m < 8; m++) {
    map_mb_to_m[m][0] = (m - 2) % 4;
    map_mb_to_m[m][1] = (9 - m) % 4;
    map_mb_to_m[m][2] = m % 4;
    map_mb_to_m[m][3] = (7 - m) % 4;
  }

  /* zero h_vca */
  for (int r = 0; r < sz_ham; r++) {
    for (int c = 0; c < sz_ham; c++) {
      h_vca[r][c].r = h_vca[r][c].i = 0.0;
    }
  }

  if (d->_bandstruct_on_local_domain) {
    char *filename = n3d_strdup_n(d->inputfile);
    sprintf(filetag, "nd_Ek_proc_%d_kx%.3f_ky_%.3f_kz_%.3f_r", mpi_n3d_id,
            d->kxL, d->kyL, d->kzL);
    n3d_FileTypeSet(&filename, filetag, TRUE);

    if ((fpLocal = fopen(filename, "w")) == NULL)
      die("Could not open the file \"%s\"", filename);
    str_free(filename);
  }

  atom_total = (iex - isx) * d->geo.AtomsPerCellMax() / 2;

  // Loop over all the cells belonging to this processor
  for (int l = isx; l < iex; l++) {
    if (((l - isx) % 100) == 0) {
      printf("atom %d / %d\n", l - isx, iex - isx);
    }
    bool flag = true;
    // all atoms in cell must exist; else throw out whole cell
    for (int m = 0; m < d->geo.AtomsPerCellMax(); m++) {
      if (!d->geo.AtomType[l][m]) {
        flag = false;
        break;
      }

      double x = d->geo.lattice_x *
                     (d->geo.cell__ijk[l][0] + d->geo.PositionInCell(m, 0)) +
                 d->geo.l_m_2_atomdsp[l][m][0];
      double y = d->geo.lattice_y *
                     (d->geo.cell__ijk[l][1] + d->geo.PositionInCell(m, 1)) +
                 d->geo.l_m_2_atomdsp[l][m][1];
      double z = d->geo.lattice_z *
                     (d->geo.cell__ijk[l][2] + d->geo.PositionInCell(m, 2)) +
                 d->geo.l_m_2_atomdsp[l][m][2];

      /* all atoms in cell have to be inside subdomain box */
      if (x < d->opt.ExecParam.Subdomain.xMin ||
          x > d->opt.ExecParam.Subdomain.xMax ||
          y < d->opt.ExecParam.Subdomain.yMin ||
          y > d->opt.ExecParam.Subdomain.yMax ||
          z < d->opt.ExecParam.Subdomain.zMin ||
          z > d->opt.ExecParam.Subdomain.zMax) {
        flag = false;
        break;
      }

      int i = d->geo.cell__ijk[l][0];
      int j = d->geo.cell__ijk[l][1];
      int k = d->geo.cell__ijk[l][2];

      if (d->geo.isOnSurface(l, m)) {
        flag = false;
        break;
      } else {
        for (int n = 0; n < d->geo.Neighbors(m); n++) {
          int id = d->geo.NbrCell(m, n, 0);
          int jd = d->geo.NbrCell(m, n, 1);
          int kd = d->geo.NbrCell(m, n, 2);
          if (d->geo.ijk__cell[i + id][j + jd][k + kd] < 0) {
            flag = false;
          }
        }
      }
    }

    if (flag) {
      compute_local_hamiltonian_at_cell(h, d, l);
    } else {
      continue;
    }

    if (d->opt.ExecParam.ElCalc.ElOut.LocalHam) {
      FILE *fpHam;
      char *filename = n3d_strdup_n(d->inputfile);

      n3d_FileTypeSet(&filename, "nd_ham_k_r", TRUE);
      fpHam = fopen(filename, "a");

      for (int r = 0; r < sz_ham; r++) {
        for (int c = 0; c < sz_ham; c++) {
          if (h[r][c].r != 0.0 || h[r][c].i != 0.0) {
            fprintf(fpHam, "%d %d %e %e\n", r, c, h[r][c].r, h[r][c].i);
          }
        }
      }
      fclose(fpHam);
      str_free(filename);
    }

    if (d->_bandstruct_on_local_domain && fpLocal) {

      ceigval_ns_full(E, h);

      double x = d->geo.lattice_x * d->geo.cell__ijk[l][0];
      double y = d->geo.lattice_y * d->geo.cell__ijk[l][1];
      double z = d->geo.lattice_z * d->geo.cell__ijk[l][2];

      for (int i = 0; i < sz_ham; i++) {
        fprintf(fpLocal, "%e %e %e %e %e %e %e   0\n", d->kxL, d->kyL, d->kzL,
                x, y, z, E[i].r);
      }
    }

    if (d->_bandstruct_on_full_domain) {
      for (int r = 0; r < sz_ham; r++) {
        for (int c = 0; c < sz_ham; c++) {
          h_vca[r][c].r += h[r][c].r;
          h_vca[r][c].i += h[r][c].i;
        }
      }
    }
  }

  rm_cmatrix(&h);

  if (d->_bandstruct_on_local_domain && fpLocal)
    fclose(fpLocal);

  /* compute vca */
  if (d->_bandstruct_on_full_domain) {
    int Ncat_tot;
    cmatrix h_vca_tot = Coperator(sz_ham);

#if (defined MPI3d && !defined FAKE_MPI)
    MPI_Allreduce(&h_vca[0][0], &h_vca_tot[0][0], 2 * sz_ham * sz_ham,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&Ncat, &Ncat_tot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
    h_vca_tot[0][0] = h_vca[0][0];
    Ncat_tot = Ncat;
#endif

    if (mpi_n3d_id == mpi_n3d_masterid) {
      FILE *fp_vca;
      char *filename = n3d_strdup_n(d->inputfile);
      ceigval_ns_full(E_vca, h_vca_tot);

      sprintf(filetag, "nd_Ek_vca_kx=%.3f_ky=%.3f_kz=%.3f", d->kxL, d->kyL,
              d->kzL);
      n3d_FileTypeSet(&filename, filetag, TRUE);

      if ((fp_vca = fopen(filename, "w")) == NULL) {
        die("Could not open the file \"%s\"", filename);
      }
      str_free(filename);

      for (int i = 0; i < sz_ham; i++)
        fprintf(fp_vca, "%e %e %e %e\n", d->kxL, d->kyL, d->kzL,
                E_vca[i].r / Ncat_tot);
      fclose(fp_vca);
    }
    rm_cmatrix(&h_vca_tot);
  }

  rm_cmatrix(&h);
  rm_cmatrix(&h_vca);
  rm_cvectr(&E);
  rm_cvectr(&E_vca);
}
