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
$Header: /repo/nemo3d/src/base/Hmult_stored_old.c,v 1.1 2006/11/23 19:44:28
hoonryu Exp $
*****************************************************************************/

#include "Hmult_stored.h"

#ifdef INTEL_MODIFIED
#include <emmintrin.h>
#include <pmmintrin.h>
#endif

#ifdef HamMemorySave
#include "Hmult_stored__Hv_stored_20_Hzb_Vmem.c"
#endif

// #define MPI_TIMEmeasure
#undef MPI_TIMEmeasure

#ifdef MPI_TIMEmeasure
#include "nemo3d_mpi.h"
int counter = 0, max_counter = 10;
double T____diag = 0.0, T____mag = 0.0, T____Hdu = 0.0, T____Hds = 0.0,
       T____Hu = 0.0, T____Ho = 0.0;
double T_V0_diag = 0.0, T_V0_mag = 0.0, T_V0_Hdu = 0.0, T_V0_Hds = 0.0,
       T_V0_Hu = 0.0, T_V0_Ho = 0.0;
double T_V1_diag = 0.0, T_V1_mag = 0.0, T_V1_Hdu = 0.0, T_V1_Hds = 0.0,
       T_V1_Hu = 0.0, T_V1_Ho = 0.0;
#include "Hmult_stored__Hv_stored_20_Hzb_V0.c"
#include "Hmult_stored__Hv_stored_20_Hzb_V1.c"
#endif

/* Computes y=(H-Is)x where H is the Hamiltonian, s is input shift, and x is
   input vector.  In this function x is a subvector starting at isx and ending
   at iex.  isx and iex are atomic (not actual) indices.  On output y contains
   the [isx:isy] subvector.  The rest of y is contained in yc.
   It is arranged as follows yc[i*d-nb] = y[ycmap[i][0]*d->NBasisStates]. */
int Hmult_spds_col_stored(cvectr y, cvectr yc, qd_struct d, real shift,
                          cvectr x, int proc) {
  MPI_TIME_INIT(Tstart);

  int i;

#ifdef MPI3d
  for (i = 0; i < d->nvmap[proc] * d->NBasisStates; i++)
    yc[i].r = yc[i].i = 0.0;
#endif

#ifdef INTEL_MODIFIED
#include "Hmult_stored_SSE_1.c"
#else
  for (i = 0; i < d->seg_ln[proc]; i++) {
    y[i].r = y[i].i = 0.0;
  }
#endif

  MPI_TIC(Tstart);

  switch (d->BandModel) {
  case BM_20_sp3d5ss_spin:
    if (d->HzbComplex.isInitialized()) {
      Hv_stored_20_HzbComplex(y, yc, x, d);
      break;
    } else {
      Hv_stored_20_Hzb(y, yc, x, d);

#ifdef MPI_TIMEmeasure
      counter++;
      Hv_stored_20_Hzb_V0(y, yc, x, d);
      Hv_stored_20_Hzb_V1(y, yc, x, d);
      printf("CPU %3d Vector-time:  T____diag=%g T____mag=%g T____Hdu=%g "
             "T____Hds=%g T____Hu=%g T____Ho=%g\n",
             mpi_n3d_id, T____diag / counter, T____mag / counter,
             T____Hdu / counter, T____Hds / counter, T____Hu / counter,
             T____Ho / counter),
          fflush(stdout);
      printf("CPU %3d Vector-time:  T_V0_diag=%g T_V0_mag=%g T_V0_Hdu=%g "
             "T_V0_Hds=%g T_V0_Hu=%g T_V0_Ho=%g\n",
             mpi_n3d_id, T_V0_diag / counter, T_V0_mag / counter,
             T_V0_Hdu / counter, T_V0_Hds / counter, T_V0_Hu / counter,
             T_V0_Ho) /
          counter,
          fflush(stdout);
      printf("CPU %3d Vector-time:  T_V1_diag=%g T_V1_mag=%g T_V1_Hdu=%g "
             "T_V1_Hds=%g T_V1_Hu=%g T_V1_Ho=%g\n",
             mpi_n3d_id, T_V1_diag / counter, T_V1_mag / counter,
             T_V1_Hdu / counter, T_V1_Hds / counter, T_V1_Hu / counter,
             T_V1_Ho) /
          counter,
          fflush(stdout);

      if (counter >= max_counter) {
        nemo3d_mpi_finalize();
        exit(-1);
      }

#endif
      break;
    }
  case BM_10_sp3ss_spin:
    if (d->HzbComplex.isInitialized()) {
      Hv_stored_10_HzbComplex(y, yc, x, d);
      break;
    } else {
      Hv_stored_10_Hzb(y, yc, x, d);
      break;
    }
  case BM_10_sp3d5ss_nospin:
    if (d->HzbComplex.isInitialized()) {
      Hv_stored_10_sp3d5s_HzbComplex(y, yc, x, d);
      break;
    } else {
      Hv_stored_10_sp3d5s_Hzb(y, yc, x, d);
      break;
    }
  case BM_1_s_nospin:
    die("ERROR in Hmult_spds_col_stored: Bands_1_s_nospin is not yet "
        "implemented");
    break;
  default:
    die("ERROR in Hmult_spds_col_stored: not yet implemented "
        "Hmult_spds_col_stored");
    break;
  }

  MPI_TOC(mpiTiming.matmul_spds, Tstart);

  if (shift != 0.0) {
    for (i = 0; i < d->seg_ln[proc]; i++) {

      y[i].r -= shift * x[i].r;
      y[i].i -= shift * x[i].i;
    }
  }

  return d->nvmap[proc];
}

/* Matrix-vector multiplication customized for the
   spin-up and spin-down sp3d5s* model Bands_20_sp3d5ss_spin */
void Hv_stored_20_Hzb(cvectr y, cvectr yc, cvectr x, qd_struct d) {
  const int Nb_2 = 10;
  int atom, i, j, k, thislength;
  register int indx;
  bool has_periodicity = d->geo.hasPeriodicity();
  ham_mem_real_ptr H_IJ;
  ham_mem_complex_ptr Hc_IJ;
#ifdef MPI_TIMEmeasure
  MPI_TIME_INIT(Tstart);
  MPI_TIC(Tstart);
#endif

  /* y += Hdd * x */
#ifdef INTEL_MODIFIED
#include "Hmult_stored_SSE_2.c"
#else
  indx = 0;
  for (atom = 0; atom < d->Hzb.Natom; atom++) {
    for (i = 0; i < d->NBasisStates; i++) {
      y[indx].r += d->Hzb.Hdd[indx] * x[indx].r;
      y[indx].i += d->Hzb.Hdd[indx] * x[indx].i;
      indx++;
    }
  }

#endif /* INTEL_MODIFIED*/

#ifdef MPI_TIMEmeasure
  MPI_TOC(T____diag, Tstart);
  MPI_TIC(Tstart);
#endif

  /* y += HBxy * x */
  if (d->opt.ExecParam.Phys.MagneticFieldOn) {
    int nb = d->NBasisStates / 2;
    for (atom = 0; atom < d->Hzb.Natom; atom++) {
      int offset = atom * d->NBasisStates;
      for (int i = 0; i < nb; i++) {
        y[offset + i].r += d->Hzb.HBxy.r * x[offset + i + nb].r -
                           d->Hzb.HBxy.i * x[offset + i + nb].i;
        y[offset + i].i += d->Hzb.HBxy.r * x[offset + i + nb].i +
                           d->Hzb.HBxy.i * x[offset + i + nb].r;
        y[offset + i + nb].r +=
            d->Hzb.HBxy.r * x[offset + i].r + d->Hzb.HBxy.i * x[offset + i].i;
        y[offset + i + nb].i +=
            d->Hzb.HBxy.r * x[offset + i].i - d->Hzb.HBxy.i * x[offset + i].r;
      }
    }
  }

#ifdef MPI_TIMEmeasure
  MPI_TOC(T____mag, Tstart);
  MPI_TIC(Tstart);
#endif

  /* y += Hdu * x */
  Hc_IJ = &d->Hzb.Hdu[0];
  for (atom = 0; atom < d->Hzb.Natom; atom++) {
    int offset = atom * d->NBasisStates;
    int I, J;

    I = offset + 2;
    J = offset + 3;
    y[I].r += Hc_IJ->r * x[J].r - Hc_IJ->i * x[J].i;
    y[I].i += Hc_IJ->r * x[J].i + Hc_IJ->i * x[J].r;
    y[J].r += Hc_IJ->r * x[I].r + Hc_IJ->i * x[I].i;
    y[J].i += Hc_IJ->r * x[I].i - Hc_IJ->i * x[I].r;
    Hc_IJ++;

    I = offset + 12;
    J = offset + 13;
    y[I].r += Hc_IJ->r * x[J].r - Hc_IJ->i * x[J].i;
    y[I].i += Hc_IJ->r * x[J].i + Hc_IJ->i * x[J].r;
    y[J].r += Hc_IJ->r * x[I].r + Hc_IJ->i * x[I].i;
    y[J].i += Hc_IJ->r * x[I].i - Hc_IJ->i * x[I].r;
    Hc_IJ++;

    I = offset + 4;
    J = offset + 12;
    y[I].r += Hc_IJ->r * x[J].r - Hc_IJ->i * x[J].i;
    y[I].i += Hc_IJ->r * x[J].i + Hc_IJ->i * x[J].r;
    y[J].r += Hc_IJ->r * x[I].r + Hc_IJ->i * x[I].i;
    y[J].i += Hc_IJ->r * x[I].i - Hc_IJ->i * x[I].r;
    Hc_IJ++;

    I = offset + 4;
    J = offset + 13;
    y[I].r += Hc_IJ->r * x[J].r - Hc_IJ->i * x[J].i;
    y[I].i += Hc_IJ->r * x[J].i + Hc_IJ->i * x[J].r;
    y[J].r += Hc_IJ->r * x[I].r + Hc_IJ->i * x[I].i;
    y[J].i += Hc_IJ->r * x[I].i - Hc_IJ->i * x[I].r;
    Hc_IJ++;

    I = offset + 2;
    J = offset + 14;
    y[I].r += Hc_IJ->r * x[J].r - Hc_IJ->i * x[J].i;
    y[I].i += Hc_IJ->r * x[J].i + Hc_IJ->i * x[J].r;
    y[J].r += Hc_IJ->r * x[I].r + Hc_IJ->i * x[I].i;
    y[J].i += Hc_IJ->r * x[I].i - Hc_IJ->i * x[I].r;
    Hc_IJ++;

    I = offset + 3;
    J = offset + 14;
    y[I].r += Hc_IJ->r * x[J].r - Hc_IJ->i * x[J].i;
    y[I].i += Hc_IJ->r * x[J].i + Hc_IJ->i * x[J].r;
    y[J].r += Hc_IJ->r * x[I].r + Hc_IJ->i * x[I].i;
    y[J].i += Hc_IJ->r * x[I].i - Hc_IJ->i * x[I].r;
    Hc_IJ++;
  }

#ifdef MPI_TIMEmeasure
  MPI_TOC(T____Hdu, Tstart);
  MPI_TIC(Tstart);
#endif

  /* y += Hds * x */
  H_IJ = &d->Hzb.Hds[0];

  for (int Zsurf = 0; Zsurf < d->Hzb.Natom_surf; Zsurf++) {
    int I, J;
    int offset = d->surfaceAtoms[Zsurf] * d->NBasisStates;

    I = offset + 1;
    J = offset + 2;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;

    I = offset + 11;
    J = offset + 12;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;
    H_IJ++;

    I = offset + 1;
    J = offset + 3;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;

    I = offset + 11;
    J = offset + 13;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;
    H_IJ++;

    I = offset + 1;
    J = offset + 4;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;

    I = offset + 11;
    J = offset + 14;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;
    H_IJ++;

    I = offset + 4;
    J = offset + 2;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;

    I = offset + 14;
    J = offset + 12;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;
    H_IJ++;

    I = offset + 4;
    J = offset + 3;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;

    I = offset + 14;
    J = offset + 13;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;
    H_IJ++;
  }

#ifdef MPI_TIMEmeasure
  MPI_TOC(T____Hds, Tstart);
  MPI_TIC(Tstart);
#endif

  /* y += Hu * x */
  indx = 0;
#ifdef INTEL_MODIFIED
#include "Hmult_stored_SSE_3a.c"
#endif /* INTEL_MODIFIED*/

  thislength = nml_iv_extent(d->Hzb.indxHU_to_row);
  for (k = 0; k < thislength; k++) {

    real kdotl = 0.0;
    if (has_periodicity) {
      int aindx = d->Hzb.indxHU_to_col[k] / d->NBasisStates;
      kdotl = d->geo.phaseBloch(d->kxL, d->kyL, d->kzL, d->Atom_to_lm[aindx][0],
                                d->Atom_to_lm[aindx][1], d->Hzb.indxHU_nbr[k]);
    }

    if (has_periodicity && kdotl != 0.0) {
      real cs = cos(kdotl);
      real sn = sin(kdotl);
      //==================================================================================================
      for (i = 0; i < Nb_2; i++) {
        int row = d->Hzb.indxHU_to_row[k] + i;
        complex *x_I = x + row;
        complex *y_I = y + row;
        complex *x_I2 = x + row + Nb_2;
        complex *y_I2 = y + row + Nb_2;
        for (j = 0; j < Nb_2; j++) {
          int col = d->Hzb.indxHU_to_col[k] + j;
          complex *x_J = x + col;
          complex *y_J = y + col;
          complex *x_J2 = x + col + Nb_2;
          complex *y_J2 = y + col + Nb_2;
          register double Hu_indx_r, Hu_indx_i;
          if (d->opt.ExecParam.Phys.MagneticFieldOn) {
            Hu_indx_r = d->Hzb.Hu[indx] *
                        (cs * d->Hzb.HphaseIn[k].r - sn * d->Hzb.HphaseIn[k].i);
            Hu_indx_i = d->Hzb.Hu[indx] *
                        (cs * d->Hzb.HphaseIn[k].i + sn * d->Hzb.HphaseIn[k].r);
          } else {
            Hu_indx_r = d->Hzb.Hu[indx] * cs;
            Hu_indx_i = d->Hzb.Hu[indx] * sn;
          }

          /* do lower (upper?) triangular part */
          y_I->r += Hu_indx_r * x_J->r - Hu_indx_i * x_J->i;
          y_I->i += Hu_indx_r * x_J->i + Hu_indx_i * x_J->r;
          y_I2->r += Hu_indx_r * x_J2->r - Hu_indx_i * x_J2->i;
          y_I2->i += Hu_indx_r * x_J2->i + Hu_indx_i * x_J2->r;

          /* do upper (lower?) triangular part */
          y_J->r += Hu_indx_r * x_I->r + Hu_indx_i * x_I->i;
          y_J->i += Hu_indx_r * x_I->i - Hu_indx_i * x_I->r;
          y_J2->r += Hu_indx_r * x_I2->r + Hu_indx_i * x_I2->i;
          y_J2->i += Hu_indx_r * x_I2->i - Hu_indx_i * x_I2->r;

          indx++;
        }
      }
      //==================================================================================================

    } else {
      //==================================================================================================

#ifdef INTEL_MODIFIED
#include "Hmult_stored_SSE_3.c"
#else

      if (d->opt.ExecParam.Phys.MagneticFieldOn) {
        for (i = 0; i < Nb_2; i++) {
          int row = d->Hzb.indxHU_to_row[k] + i;
          complex *x_I = x + row;
          complex *y_I = y + row;
          complex *x_I2 = x + row + Nb_2;
          complex *y_I2 = y + row + Nb_2;
          for (j = 0; j < Nb_2; j++) {
            int col = d->Hzb.indxHU_to_col[k] + j;
            complex *x_J = x + col;
            complex *y_J = y + col;
            complex *x_J2 = x + col + Nb_2;
            complex *y_J2 = y + col + Nb_2;
            // inner if
            double Hu_indx_r = d->Hzb.Hu[indx] * d->Hzb.HphaseIn[k].r;
            double Hu_indx_i = d->Hzb.Hu[indx] * d->Hzb.HphaseIn[k].i;
            // do lower (upper?) triangular part
            y_I->r += Hu_indx_r * x_J->r - Hu_indx_i * x_J->i;
            y_I->i += Hu_indx_r * x_J->i + Hu_indx_i * x_J->r;
            y_I2->r += Hu_indx_r * x_J2->r - Hu_indx_i * x_J2->i;
            y_I2->i += Hu_indx_r * x_J2->i + Hu_indx_i * x_J2->r;
            // do upper (lower?) triangular part
            y_J->r += Hu_indx_r * x_I->r + Hu_indx_i * x_I->i;
            y_J->i += Hu_indx_r * x_I->i - Hu_indx_i * x_I->r;
            y_J2->r += Hu_indx_r * x_I2->r + Hu_indx_i * x_I2->i;
            y_J2->i += Hu_indx_r * x_I2->i - Hu_indx_i * x_I2->r;
            // end inner if
            indx++;
          }
        }
      } else {

        double tmp;

        complex *x_J;
        complex *y_J;
        complex *x_J2;
        complex *y_J2;

        for (i = 0; i < Nb_2; i++) {
          int row = d->Hzb.indxHU_to_row[k] + i;
          complex *x_I = x + row;
          complex *y_I = y + row;
          complex *x_I2 = x_I + Nb_2;
          complex *y_I2 = y_I + Nb_2;

          int col = d->Hzb.indxHU_to_col[k];
          for (j = 0; j < Nb_2; j++) {
            // int col = d->Hzb.indxHU_to_col[k] + j;

            x_J = x + col;
            y_J = y + col;
            x_J2 = x_J + Nb_2;
            y_J2 = y_J + Nb_2;

            tmp = (double)d->Hzb.Hu[indx];

            // inner if
            //  do lower (upper?) triangular part

            y_I->r += tmp * x_J->r;
            y_I->i += tmp * x_J->i;
            y_I2->r += tmp * x_J2->r;
            y_I2->i += tmp * x_J2->i;

            // do upper (lower?) triangular part
            y_J->r += tmp * x_I->r;
            y_J->i += tmp * x_I->i;
            y_J2->r += tmp * x_I2->r;
            y_J2->i += tmp * x_I2->i;
            // end inner if
            indx++;
            col++;
          }
        }
      }

#endif /* INTEL_MODIFIED*/
      //==================================================================================================
    }
  }

#ifdef MPI_TIMEmeasure
  MPI_TOC(T____Hu, Tstart);
  MPI_TIC(Tstart);
#endif
  /* yc += Ho * x */
  thislength = nml_iv_extent(d->Hzb.indxHO_to_row);
  for (k = 0; k < thislength; k++) {

    real kdotl = 0.0;
    if (has_periodicity) {
      int aindx = d->Hzb.indxHO_to_col[k] / d->NBasisStates;
      kdotl = d->geo.phaseBloch(d->kxL, d->kyL, d->kzL, d->Atom_to_lm[aindx][0],
                                d->Atom_to_lm[aindx][1], d->Hzb.indxHO_nbr[k]);
    }

    if (has_periodicity && kdotl != 0.0) {
      real cs = cos(kdotl);
      real sn = sin(kdotl);

      for (i = 0; i < Nb_2; i++) {
        int I = d->Hzb.indxHO_to_row[k] + i;
        for (j = 0; j < Nb_2; j++) {
          int indx = Nb_2 * (Nb_2 * k + i) + j;
          int J = d->Hzb.indxHO_to_col[k] + j;
          register double Ho_indx_r, Ho_indx_i;
          if (d->opt.ExecParam.Phys.MagneticFieldOn) {
            Ho_indx_r = d->Hzb.Ho[indx] * (cs * d->Hzb.HphaseOut[k].r -
                                           sn * d->Hzb.HphaseOut[k].i);
            Ho_indx_i = d->Hzb.Ho[indx] * (cs * d->Hzb.HphaseOut[k].i +
                                           sn * d->Hzb.HphaseOut[k].r);
          } else {
            Ho_indx_r = d->Hzb.Ho[indx] * cs;
            Ho_indx_i = d->Hzb.Ho[indx] * sn;
          }

          yc[I].r += Ho_indx_r * x[J].r - Ho_indx_i * x[J].i;
          yc[I].i += Ho_indx_r * x[J].i + Ho_indx_i * x[J].r;

          yc[I + Nb_2].r +=
              Ho_indx_r * x[J + Nb_2].r - Ho_indx_i * x[J + Nb_2].i;
          yc[I + Nb_2].i +=
              Ho_indx_r * x[J + Nb_2].i + Ho_indx_i * x[J + Nb_2].r;
        }
      }
    } else {
      for (i = 0; i < Nb_2; i++) {
        int I = d->Hzb.indxHO_to_row[k] + i;
        for (j = 0; j < Nb_2; j++) {
          int indx = Nb_2 * (Nb_2 * k + i) + j;
          int J = d->Hzb.indxHO_to_col[k] + j;
          if (d->opt.ExecParam.Phys.MagneticFieldOn) {
            yc[I].r += d->Hzb.Ho[indx] * (x[J].r * d->Hzb.HphaseOut[k].r -
                                          x[J].i * d->Hzb.HphaseOut[k].i);
            yc[I].i += d->Hzb.Ho[indx] * (x[J].i * d->Hzb.HphaseOut[k].r +
                                          x[J].r * d->Hzb.HphaseOut[k].i);

            yc[I + Nb_2].r +=
                d->Hzb.Ho[indx] * (x[J + Nb_2].r * d->Hzb.HphaseOut[k].r -
                                   x[J + Nb_2].i * d->Hzb.HphaseOut[k].i);
            yc[I + Nb_2].i +=
                d->Hzb.Ho[indx] * (x[J + Nb_2].r * d->Hzb.HphaseOut[k].i +
                                   x[J + Nb_2].i * d->Hzb.HphaseOut[k].r);
          } else {
            yc[I].r += d->Hzb.Ho[indx] * x[J].r;
            yc[I].i += d->Hzb.Ho[indx] * x[J].i;

            yc[I + Nb_2].r += d->Hzb.Ho[indx] * x[J + Nb_2].r;
            yc[I + Nb_2].i += d->Hzb.Ho[indx] * x[J + Nb_2].i;
          }
        }
      }
    }
  }
#ifdef MPI_TIMEmeasure
  MPI_TOC(T____Ho, Tstart);
  /* printf("Vector-time:  T____diag=%g T____mag=%g T____Hdu=%g T____Hds=%g
   * T____Hu=%g T____Ho=%g\n",T____diag, T____mag, T____Hdu, T____Hds, T____Hu,
   * T____Ho), fflush(stdout); */
#endif
}

/************************************
 ************************************
 */
/* Matrix-vector multiplication customized for the
   spin-less sp3d5s* model Bands_10_sp3d5ss_nospin */
void Hv_stored_10_sp3d5s_Hzb(cvectr y, cvectr yc, cvectr x, qd_struct d) {
  const int Nb = 10;
  int atom, indx, i, j, k, thislength;
  bool has_periodicity = d->geo.hasPeriodicity();
  ham_mem_real_ptr H_IJ;

  /* y += Hdd * x */
  indx = 0;
  for (atom = 0; atom < d->Hzb.Natom; atom++) {
    for (i = 0; i < d->NBasisStates; i++) {
      y[indx].r += d->Hzb.Hdd[indx] * x[indx].r;
      y[indx].i += d->Hzb.Hdd[indx] * x[indx].i;
      indx++;
    }
  }

  /* The spinless sp3d5s* model does not have spin
     So, the Zeeman energy cannot be added (y += HBxy * x)
  */

  /* The spinless sp3d5s* model does not have any
     spin-orbit coupling.
     The matrix Hdu therefore does not have to be utilized at all.
  */

  /* y += Hds * x */
  H_IJ = &d->Hzb.Hds[0];
  for (int Zsurf = 0; Zsurf < d->Hzb.Natom_surf; Zsurf++) {
    int I, J;
    int offset = d->surfaceAtoms[Zsurf] * d->NBasisStates;

    I = offset + 1;
    J = offset + 2;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;

    H_IJ++;

    I = offset + 1;
    J = offset + 3;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;

    H_IJ++;

    I = offset + 1;
    J = offset + 4;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;

    H_IJ++;

    I = offset + 2;
    J = offset + 3;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;

    H_IJ++;

    I = offset + 2;
    J = offset + 4;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;

    H_IJ++;

    I = offset + 3;
    J = offset + 4;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;

    H_IJ++;
  }

  /* y += Hu * x */
  indx = 0;
  thislength = nml_iv_extent(d->Hzb.indxHU_to_row);
  for (k = 0; k < thislength; k++) {

    real kdotl = 0.0;
    if (has_periodicity) {
      int aindx = d->Hzb.indxHU_to_col[k] / d->NBasisStates;
      kdotl = d->geo.phaseBloch(d->kxL, d->kyL, d->kzL, d->Atom_to_lm[aindx][0],
                                d->Atom_to_lm[aindx][1], d->Hzb.indxHU_nbr[k]);
    }

    if (has_periodicity && kdotl != 0.0) {
      real cs = cos(kdotl);
      real sn = sin(kdotl);

      for (i = 0; i < Nb; i++) {
        int row = d->Hzb.indxHU_to_row[k] + i;
        complex *x_I = x + row;
        complex *y_I = y + row;
        for (j = 0; j < Nb; j++) {
          int col = d->Hzb.indxHU_to_col[k] + j;
          complex *x_J = x + col;
          complex *y_J = y + col;
          register double Hu_indx_r, Hu_indx_i;
          if (d->opt.ExecParam.Phys.MagneticFieldOn) {
            Hu_indx_r = d->Hzb.Hu[indx] *
                        (cs * d->Hzb.HphaseIn[k].r - sn * d->Hzb.HphaseIn[k].i);
            Hu_indx_i = d->Hzb.Hu[indx] *
                        (cs * d->Hzb.HphaseIn[k].i + sn * d->Hzb.HphaseIn[k].r);
          } else {
            Hu_indx_r = d->Hzb.Hu[indx] * cs;
            Hu_indx_i = d->Hzb.Hu[indx] * sn;
          }

          /* do lower (upper?) triangular part */
          y_I->r += Hu_indx_r * x_J->r - Hu_indx_i * x_J->i;
          y_I->i += Hu_indx_r * x_J->i + Hu_indx_i * x_J->r;

          /* do upper (lower?) triangular part */
          y_J->r += Hu_indx_r * x_I->r + Hu_indx_i * x_I->i;
          y_J->i += Hu_indx_r * x_I->i - Hu_indx_i * x_I->r;

          indx++;
        }
      }
    } else {
      for (i = 0; i < Nb; i++) {
        int row = d->Hzb.indxHU_to_row[k] + i;
        complex *x_I = x + row;
        complex *y_I = y + row;
        for (j = 0; j < Nb; j++) {
          int col = d->Hzb.indxHU_to_col[k] + j;
          complex *x_J = x + col;
          complex *y_J = y + col;

          if (d->opt.ExecParam.Phys.MagneticFieldOn) {
            register double Hu_indx_r = d->Hzb.Hu[indx] * d->Hzb.HphaseIn[k].r;
            register double Hu_indx_i = d->Hzb.Hu[indx] * d->Hzb.HphaseIn[k].i;
            /* do lower (upper?) triangular part */
            y_I->r += Hu_indx_r * x_J->r - Hu_indx_i * x_J->i;
            y_I->i += Hu_indx_r * x_J->i + Hu_indx_i * x_J->r;

            /* do upper (lower?) triangular part */
            y_J->r += Hu_indx_r * x_I->r + Hu_indx_i * x_I->i;
            y_J->i += Hu_indx_r * x_I->i - Hu_indx_i * x_I->r;
          } else {
            /* do lower (upper?) triangular part */
            y_I->r += d->Hzb.Hu[indx] * x_J->r;
            y_I->i += d->Hzb.Hu[indx] * x_J->i;

            /* do upper (lower?) triangular part */
            y_J->r += d->Hzb.Hu[indx] * x_I->r;
            y_J->i += d->Hzb.Hu[indx] * x_I->i;
          }
          indx++;
        }
      }
    }
  }

  /* yc += Ho * x */
  thislength = nml_iv_extent(d->Hzb.indxHO_to_row);
  for (k = 0; k < thislength; k++) {

    real kdotl = 0.0;
    if (has_periodicity) {
      int aindx = d->Hzb.indxHO_to_col[k] / d->NBasisStates;
      kdotl = d->geo.phaseBloch(d->kxL, d->kyL, d->kzL, d->Atom_to_lm[aindx][0],
                                d->Atom_to_lm[aindx][1], d->Hzb.indxHO_nbr[k]);
    }

    if (has_periodicity && kdotl != 0.0) {
      real cs = cos(kdotl);
      real sn = sin(kdotl);

      for (i = 0; i < Nb; i++) {
        int I = d->Hzb.indxHO_to_row[k] + i;
        for (j = 0; j < Nb; j++) {
          int indx = Nb * (Nb * k + i) + j;
          int J = d->Hzb.indxHO_to_col[k] + j;
          register double Ho_indx_r, Ho_indx_i;
          if (d->opt.ExecParam.Phys.MagneticFieldOn) {
            Ho_indx_r = d->Hzb.Ho[indx] * (cs * d->Hzb.HphaseOut[k].r -
                                           sn * d->Hzb.HphaseOut[k].i);
            Ho_indx_i = d->Hzb.Ho[indx] * (cs * d->Hzb.HphaseOut[k].i +
                                           sn * d->Hzb.HphaseOut[k].r);
          } else {
            Ho_indx_r = d->Hzb.Ho[indx] * cs;
            Ho_indx_i = d->Hzb.Ho[indx] * sn;
          }

          yc[I].r += Ho_indx_r * x[J].r - Ho_indx_i * x[J].i;
          yc[I].i += Ho_indx_r * x[J].i + Ho_indx_i * x[J].r;
        }
      }
    } else {
      for (i = 0; i < Nb; i++) {
        int I = d->Hzb.indxHO_to_row[k] + i;
        for (j = 0; j < Nb; j++) {
          int indx = Nb * (Nb * k + i) + j;
          int J = d->Hzb.indxHO_to_col[k] + j;
          if (d->opt.ExecParam.Phys.MagneticFieldOn) {
            yc[I].r += d->Hzb.Ho[indx] * (x[J].r * d->Hzb.HphaseOut[k].r -
                                          x[J].i * d->Hzb.HphaseOut[k].i);
            yc[I].i += d->Hzb.Ho[indx] * (x[J].i * d->Hzb.HphaseOut[k].r +
                                          x[J].r * d->Hzb.HphaseOut[k].i);

          } else {
            yc[I].r += d->Hzb.Ho[indx] * x[J].r;
            yc[I].i += d->Hzb.Ho[indx] * x[J].i;
          }
        }
      }
    }
  }
}

void Hv_stored_10_Hzb(cvectr y, cvectr yc, cvectr x, qd_struct d) {
  const int Nb_2 = 5;
  int atom, indx, i, j, k, thislength;
  bool has_periodicity = d->geo.hasPeriodicity();
  ham_mem_real_ptr H_IJ;
  ham_mem_complex_ptr Hc_IJ;

  /* y += Hdd * x */
  indx = 0;
  for (atom = 0; atom < d->Hzb.Natom; atom++) {
    for (i = 0; i < d->NBasisStates; i++) {
      y[indx].r += d->Hzb.Hdd[indx] * x[indx].r;
      y[indx].i += d->Hzb.Hdd[indx] * x[indx].i;
      indx++;
    }
  }

  /* y += HBxy * x */
  if (d->opt.ExecParam.Phys.MagneticFieldOn) {
    int nb = d->NBasisStates / 2;
    for (atom = 0; atom < d->Hzb.Natom; atom++) {
      int offset = atom * d->NBasisStates;
      for (int i = 0; i < nb; i++) {
        y[offset + i].r += d->Hzb.HBxy.r * x[offset + i + nb].r -
                           d->Hzb.HBxy.i * x[offset + i + nb].i;
        y[offset + i].i += d->Hzb.HBxy.r * x[offset + i + nb].i +
                           d->Hzb.HBxy.i * x[offset + i + nb].r;
        y[offset + i + nb].r +=
            d->Hzb.HBxy.r * x[offset + i].r + d->Hzb.HBxy.i * x[offset + i].i;
        y[offset + i + nb].i +=
            d->Hzb.HBxy.r * x[offset + i].i - d->Hzb.HBxy.i * x[offset + i].r;
      }
    }
  }

  /* y += Hdu * x */
  Hc_IJ = &d->Hzb.Hdu[0];
  for (atom = 0; atom < d->Hzb.Natom; atom++) {
    int offset = atom * d->NBasisStates;
    int I, J;

    I = offset + 2;
    J = offset + 3;
    y[I].r += Hc_IJ->r * x[J].r - Hc_IJ->i * x[J].i;
    y[I].i += Hc_IJ->r * x[J].i + Hc_IJ->i * x[J].r;
    y[J].r += Hc_IJ->r * x[I].r + Hc_IJ->i * x[I].i;
    y[J].i += Hc_IJ->r * x[I].i - Hc_IJ->i * x[I].r;
    Hc_IJ++;

    I = offset + 7;
    J = offset + 8;
    y[I].r += Hc_IJ->r * x[J].r - Hc_IJ->i * x[J].i;
    y[I].i += Hc_IJ->r * x[J].i + Hc_IJ->i * x[J].r;
    y[J].r += Hc_IJ->r * x[I].r + Hc_IJ->i * x[I].i;
    y[J].i += Hc_IJ->r * x[I].i - Hc_IJ->i * x[I].r;
    Hc_IJ++;

    I = offset + 4;
    J = offset + 7;
    y[I].r += Hc_IJ->r * x[J].r - Hc_IJ->i * x[J].i;
    y[I].i += Hc_IJ->r * x[J].i + Hc_IJ->i * x[J].r;
    y[J].r += Hc_IJ->r * x[I].r + Hc_IJ->i * x[I].i;
    y[J].i += Hc_IJ->r * x[I].i - Hc_IJ->i * x[I].r;
    Hc_IJ++;

    I = offset + 4;
    J = offset + 8;
    y[I].r += Hc_IJ->r * x[J].r - Hc_IJ->i * x[J].i;
    y[I].i += Hc_IJ->r * x[J].i + Hc_IJ->i * x[J].r;
    y[J].r += Hc_IJ->r * x[I].r + Hc_IJ->i * x[I].i;
    y[J].i += Hc_IJ->r * x[I].i - Hc_IJ->i * x[I].r;
    Hc_IJ++;

    I = offset + 2;
    J = offset + 9;
    y[I].r += Hc_IJ->r * x[J].r - Hc_IJ->i * x[J].i;
    y[I].i += Hc_IJ->r * x[J].i + Hc_IJ->i * x[J].r;
    y[J].r += Hc_IJ->r * x[I].r + Hc_IJ->i * x[I].i;
    y[J].i += Hc_IJ->r * x[I].i - Hc_IJ->i * x[I].r;
    Hc_IJ++;

    I = offset + 3;
    J = offset + 9;
    y[I].r += Hc_IJ->r * x[J].r - Hc_IJ->i * x[J].i;
    y[I].i += Hc_IJ->r * x[J].i + Hc_IJ->i * x[J].r;
    y[J].r += Hc_IJ->r * x[I].r + Hc_IJ->i * x[I].i;
    y[J].i += Hc_IJ->r * x[I].i - Hc_IJ->i * x[I].r;
    Hc_IJ++;
  }

  /* y += Hds * x */
  H_IJ = &d->Hzb.Hds[0];
  for (int Zsurf = 0; Zsurf < d->Hzb.Natom_surf; Zsurf++) {
    int I, J;
    int offset = d->surfaceAtoms[Zsurf] * d->NBasisStates;

    I = offset + 1;
    J = offset + 2;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;

    I = offset + 6;
    J = offset + 7;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;
    H_IJ++;

    I = offset + 1;
    J = offset + 3;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;

    I = offset + 6;
    J = offset + 8;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;
    H_IJ++;

    I = offset + 1;
    J = offset + 4;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;

    I = offset + 6;
    J = offset + 9;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;
    H_IJ++;

    I = offset + 4;
    J = offset + 2;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;

    I = offset + 9;
    J = offset + 7;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;
    H_IJ++;

    I = offset + 4;
    J = offset + 3;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;

    I = offset + 9;
    J = offset + 8;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;
    H_IJ++;
  }

  /* y += Hu * x */
  indx = 0;
  thislength = nml_iv_extent(d->Hzb.indxHU_to_row);
  for (k = 0; k < thislength; k++) {

    real kdotl = 0.0;
    if (has_periodicity) {
      int aindx = d->Hzb.indxHU_to_col[k] / d->NBasisStates;
      kdotl = d->geo.phaseBloch(d->kxL, d->kyL, d->kzL, d->Atom_to_lm[aindx][0],
                                d->Atom_to_lm[aindx][1], d->Hzb.indxHU_nbr[k]);
    }

    if (has_periodicity && kdotl != 0.0) {
      real cs = cos(kdotl);
      real sn = sin(kdotl);

      for (i = 0; i < Nb_2; i++) {
        int row = d->Hzb.indxHU_to_row[k] + i;
        complex *x_I = x + row;
        complex *y_I = y + row;
        complex *x_I2 = x + row + Nb_2;
        complex *y_I2 = y + row + Nb_2;
        for (j = 0; j < Nb_2; j++) {
          int col = d->Hzb.indxHU_to_col[k] + j;
          complex *x_J = x + col;
          complex *y_J = y + col;
          complex *x_J2 = x + col + Nb_2;
          complex *y_J2 = y + col + Nb_2;
          register double Hu_indx_r, Hu_indx_i;
          if (d->opt.ExecParam.Phys.MagneticFieldOn) {
            Hu_indx_r = d->Hzb.Hu[indx] *
                        (cs * d->Hzb.HphaseIn[k].r - sn * d->Hzb.HphaseIn[k].i);
            Hu_indx_i = d->Hzb.Hu[indx] *
                        (cs * d->Hzb.HphaseIn[k].i + sn * d->Hzb.HphaseIn[k].r);
          } else {
            Hu_indx_r = d->Hzb.Hu[indx] * cs;
            Hu_indx_i = d->Hzb.Hu[indx] * sn;
          }

          /* do lower (upper?) triangular part */
          y_I->r += Hu_indx_r * x_J->r - Hu_indx_i * x_J->i;
          y_I->i += Hu_indx_r * x_J->i + Hu_indx_i * x_J->r;
          y_I2->r += Hu_indx_r * x_J2->r - Hu_indx_i * x_J2->i;
          y_I2->i += Hu_indx_r * x_J2->i + Hu_indx_i * x_J2->r;

          /* do upper (lower?) triangular part */
          y_J->r += Hu_indx_r * x_I->r + Hu_indx_i * x_I->i;
          y_J->i += Hu_indx_r * x_I->i - Hu_indx_i * x_I->r;
          y_J2->r += Hu_indx_r * x_I2->r + Hu_indx_i * x_I2->i;
          y_J2->i += Hu_indx_r * x_I2->i - Hu_indx_i * x_I2->r;

          indx++;
        }
      }
    } else {
      for (i = 0; i < Nb_2; i++) {
        int row = d->Hzb.indxHU_to_row[k] + i;
        complex *x_I = x + row;
        complex *y_I = y + row;
        complex *x_I2 = x + row + Nb_2;
        complex *y_I2 = y + row + Nb_2;
        for (j = 0; j < Nb_2; j++) {
          int col = d->Hzb.indxHU_to_col[k] + j;
          complex *x_J = x + col;
          complex *y_J = y + col;
          complex *x_J2 = x + col + Nb_2;
          complex *y_J2 = y + col + Nb_2;
          if (d->opt.ExecParam.Phys.MagneticFieldOn) {
            register double Hu_indx_r = d->Hzb.Hu[indx] * d->Hzb.HphaseIn[k].r;
            register double Hu_indx_i = d->Hzb.Hu[indx] * d->Hzb.HphaseIn[k].i;

            /* do lower (upper?) triangular part */
            y_I->r += Hu_indx_r * x_J->r - Hu_indx_i * x_J->i;
            y_I->i += Hu_indx_r * x_J->i + Hu_indx_i * x_J->r;
            y_I2->r += Hu_indx_r * x_J2->r - Hu_indx_i * x_J2->i;
            y_I2->i += Hu_indx_r * x_J2->i + Hu_indx_i * x_J2->r;

            /* do upper (lower?) triangular part */
            y_J->r += Hu_indx_r * x_I->r + Hu_indx_i * x_I->i;
            y_J->i += Hu_indx_r * x_I->i - Hu_indx_i * x_I->r;
            y_J2->r += Hu_indx_r * x_I2->r + Hu_indx_i * x_I2->i;
            y_J2->i += Hu_indx_r * x_I2->i - Hu_indx_i * x_I2->r;
          } else {
            /* do lower (upper?) triangular part */
            y_I->r += d->Hzb.Hu[indx] * x_J->r;
            y_I->i += d->Hzb.Hu[indx] * x_J->i;
            y_I2->r += d->Hzb.Hu[indx] * x_J2->r;
            y_I2->i += d->Hzb.Hu[indx] * x_J2->i;

            /* do upper (lower?) triangular part */
            y_J->r += d->Hzb.Hu[indx] * x_I->r;
            y_J->i += d->Hzb.Hu[indx] * x_I->i;
            y_J2->r += d->Hzb.Hu[indx] * x_I2->r;
            y_J2->i += d->Hzb.Hu[indx] * x_I2->i;
          }
          indx++;
        }
      }
    }
  }

  /* yc += Ho * x */
  thislength = nml_iv_extent(d->Hzb.indxHO_to_row);
  for (k = 0; k < thislength; k++) {

    real kdotl = 0.0;
    if (has_periodicity) {
      int aindx = d->Hzb.indxHO_to_col[k] / d->NBasisStates;
      kdotl = d->geo.phaseBloch(d->kxL, d->kyL, d->kzL, d->Atom_to_lm[aindx][0],
                                d->Atom_to_lm[aindx][1], d->Hzb.indxHO_nbr[k]);
    }

    if (has_periodicity && kdotl != 0.0) {
      real cs = cos(kdotl);
      real sn = sin(kdotl);

      for (i = 0; i < Nb_2; i++) {
        int I = d->Hzb.indxHO_to_row[k] + i;
        for (j = 0; j < Nb_2; j++) {
          int indx = Nb_2 * (Nb_2 * k + i) + j;
          int J = d->Hzb.indxHO_to_col[k] + j;
          register double Ho_indx_r, Ho_indx_i;
          if (d->opt.ExecParam.Phys.MagneticFieldOn) {
            Ho_indx_r = d->Hzb.Ho[indx] * (cs * d->Hzb.HphaseOut[k].r -
                                           sn * d->Hzb.HphaseOut[k].i);
            Ho_indx_i = d->Hzb.Ho[indx] * (cs * d->Hzb.HphaseOut[k].i +
                                           sn * d->Hzb.HphaseOut[k].r);
          } else {
            Ho_indx_r = d->Hzb.Ho[indx] * cs;
            Ho_indx_i = d->Hzb.Ho[indx] * sn;
          }

          yc[I].r += Ho_indx_r * x[J].r - Ho_indx_i * x[J].i;
          yc[I].i += Ho_indx_r * x[J].i + Ho_indx_i * x[J].r;

          yc[I + Nb_2].r +=
              Ho_indx_r * x[J + Nb_2].r - Ho_indx_i * x[J + Nb_2].i;
          yc[I + Nb_2].i +=
              Ho_indx_r * x[J + Nb_2].i + Ho_indx_i * x[J + Nb_2].r;
        }
      }
    } else {
      for (i = 0; i < Nb_2; i++) {
        int I = d->Hzb.indxHO_to_row[k] + i;
        for (j = 0; j < Nb_2; j++) {
          int indx = Nb_2 * (Nb_2 * k + i) + j;
          int J = d->Hzb.indxHO_to_col[k] + j;
          if (d->opt.ExecParam.Phys.MagneticFieldOn) {
            yc[I].r += d->Hzb.Ho[indx] * (x[J].r * d->Hzb.HphaseOut[k].r -
                                          x[J].i * d->Hzb.HphaseOut[k].i);
            yc[I].i += d->Hzb.Ho[indx] * (x[J].i * d->Hzb.HphaseOut[k].r +
                                          x[J].r * d->Hzb.HphaseOut[k].i);

            yc[I + Nb_2].r +=
                d->Hzb.Ho[indx] * (x[J + Nb_2].r * d->Hzb.HphaseOut[k].r -
                                   x[J + Nb_2].i * d->Hzb.HphaseOut[k].i);
            yc[I + Nb_2].i +=
                d->Hzb.Ho[indx] * (x[J + Nb_2].r * d->Hzb.HphaseOut[k].i +
                                   x[J + Nb_2].i * d->Hzb.HphaseOut[k].r);
          } else {
            yc[I].r += d->Hzb.Ho[indx] * x[J].r;
            yc[I].i += d->Hzb.Ho[indx] * x[J].i;

            yc[I + Nb_2].r += d->Hzb.Ho[indx] * x[J + Nb_2].r;
            yc[I + Nb_2].i += d->Hzb.Ho[indx] * x[J + Nb_2].i;
          }
        }
      }
    }
  }
}

void Hv_stored_20_HzbComplex(cvectr y, cvectr yc, cvectr x, qd_struct d) {
  const int Nb_2_20 = 10;
  int atom, indx, i, j, k, thislength;
  bool has_periodicity = d->geo.hasPeriodicity();
  ham_mem_real_ptr H_IJ;
  ham_mem_complex_ptr Hc_IJ;
#undef MPI_TIMEmeasure
#ifdef MPI_TIMEmeasure
  double T_diag = 0.0, T_mag = 0.0, T_Hdu = 0.0, T_Hds = 0.0, T_Hu = 0.0,
         T_Ho = 0.0;
  MPI_TIME_INIT(Tstart);
  MPI_TIC(Tstart);
#endif

  /* y += Hdd * x */
  indx = 0;
  for (atom = 0; atom < d->HzbComplex.Natom; atom++) {
    for (i = 0; i < d->NBasisStates; i++) {
      y[indx].r += d->HzbComplex.Hdd[indx] * x[indx].r;
      y[indx].i += d->HzbComplex.Hdd[indx] * x[indx].i;
      indx++;
    }
  }

#ifdef MPI_TIMEmeasure
  MPI_TOC(T_diag, Tstart);
  MPI_TIC(Tstart);
#endif
  /* y += HBxy * x */
  if (d->opt.ExecParam.Phys.MagneticFieldOn) {
    int nb = d->NBasisStates / 2;
    for (atom = 0; atom < d->HzbComplex.Natom; atom++) {
      int offset = atom * d->NBasisStates;
      for (int i = 0; i < nb; i++) {
        y[offset + i].r += d->HzbComplex.HBxy.r * x[offset + i + nb].r -
                           d->HzbComplex.HBxy.i * x[offset + i + nb].i;
        y[offset + i].i += d->HzbComplex.HBxy.r * x[offset + i + nb].i +
                           d->HzbComplex.HBxy.i * x[offset + i + nb].r;
        y[offset + i + nb].r += d->HzbComplex.HBxy.r * x[offset + i].r +
                                d->HzbComplex.HBxy.i * x[offset + i].i;
        y[offset + i + nb].i += d->HzbComplex.HBxy.r * x[offset + i].i -
                                d->HzbComplex.HBxy.i * x[offset + i].r;
      }
    }
  }

#ifdef MPI_TIMEmeasure
  MPI_TOC(T_mag, Tstart);
  MPI_TIC(Tstart);
#endif

  /* y += Hdu * x */
  Hc_IJ = &d->HzbComplex.Hdu[0];
  for (atom = 0; atom < d->HzbComplex.Natom; atom++) {
    int offset = atom * d->NBasisStates;
    int I, J;

    I = offset + 2;
    J = offset + 3;
    y[I].r += Hc_IJ->r * x[J].r - Hc_IJ->i * x[J].i;
    y[I].i += Hc_IJ->r * x[J].i + Hc_IJ->i * x[J].r;
    y[J].r += Hc_IJ->r * x[I].r + Hc_IJ->i * x[I].i;
    y[J].i += Hc_IJ->r * x[I].i - Hc_IJ->i * x[I].r;
    Hc_IJ++;

    I = offset + 12;
    J = offset + 13;
    y[I].r += Hc_IJ->r * x[J].r - Hc_IJ->i * x[J].i;
    y[I].i += Hc_IJ->r * x[J].i + Hc_IJ->i * x[J].r;
    y[J].r += Hc_IJ->r * x[I].r + Hc_IJ->i * x[I].i;
    y[J].i += Hc_IJ->r * x[I].i - Hc_IJ->i * x[I].r;
    Hc_IJ++;

    I = offset + 4;
    J = offset + 12;
    y[I].r += Hc_IJ->r * x[J].r - Hc_IJ->i * x[J].i;
    y[I].i += Hc_IJ->r * x[J].i + Hc_IJ->i * x[J].r;
    y[J].r += Hc_IJ->r * x[I].r + Hc_IJ->i * x[I].i;
    y[J].i += Hc_IJ->r * x[I].i - Hc_IJ->i * x[I].r;
    Hc_IJ++;

    I = offset + 4;
    J = offset + 13;
    y[I].r += Hc_IJ->r * x[J].r - Hc_IJ->i * x[J].i;
    y[I].i += Hc_IJ->r * x[J].i + Hc_IJ->i * x[J].r;
    y[J].r += Hc_IJ->r * x[I].r + Hc_IJ->i * x[I].i;
    y[J].i += Hc_IJ->r * x[I].i - Hc_IJ->i * x[I].r;
    Hc_IJ++;

    I = offset + 2;
    J = offset + 14;
    y[I].r += Hc_IJ->r * x[J].r - Hc_IJ->i * x[J].i;
    y[I].i += Hc_IJ->r * x[J].i + Hc_IJ->i * x[J].r;
    y[J].r += Hc_IJ->r * x[I].r + Hc_IJ->i * x[I].i;
    y[J].i += Hc_IJ->r * x[I].i - Hc_IJ->i * x[I].r;
    Hc_IJ++;

    I = offset + 3;
    J = offset + 14;
    y[I].r += Hc_IJ->r * x[J].r - Hc_IJ->i * x[J].i;
    y[I].i += Hc_IJ->r * x[J].i + Hc_IJ->i * x[J].r;
    y[J].r += Hc_IJ->r * x[I].r + Hc_IJ->i * x[I].i;
    y[J].i += Hc_IJ->r * x[I].i - Hc_IJ->i * x[I].r;
    Hc_IJ++;
  }

#ifdef MPI_TIMEmeasure
  MPI_TOC(T_Hdu, Tstart);
  MPI_TIC(Tstart);
#endif
  /* y += Hds * x */
  H_IJ = &d->HzbComplex.Hds[0];
  for (int Zsurf = 0; Zsurf < d->HzbComplex.Natom_surf; Zsurf++) {
    int I, J;
    int offset = d->surfaceAtoms[Zsurf] * d->NBasisStates;

    I = offset + 1;
    J = offset + 2;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;

    I = offset + 11;
    J = offset + 12;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;
    H_IJ++;

    I = offset + 1;
    J = offset + 3;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;

    I = offset + 11;
    J = offset + 13;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;
    H_IJ++;

    I = offset + 1;
    J = offset + 4;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;

    I = offset + 11;
    J = offset + 14;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;
    H_IJ++;

    I = offset + 4;
    J = offset + 2;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;

    I = offset + 14;
    J = offset + 12;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;
    H_IJ++;

    I = offset + 4;
    J = offset + 3;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;

    I = offset + 14;
    J = offset + 13;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;
    H_IJ++;
  }
#ifdef MPI_TIMEmeasure
  MPI_TOC(T_Hds, Tstart);
  MPI_TIC(Tstart);
#endif

  /* y += Hu * x */
  indx = 0;
  thislength = nml_iv_extent(d->HzbComplex.indxHU_to_row);
  for (k = 0; k < thislength; k++) {

    real kdotl = 0.0;
    if (has_periodicity) {
      int aindx = d->HzbComplex.indxHU_to_col[k] / d->NBasisStates;
      kdotl = d->geo.phaseBloch(d->kxL, d->kyL, d->kzL, d->Atom_to_lm[aindx][0],
                                d->Atom_to_lm[aindx][1],
                                d->HzbComplex.indxHU_nbr[k]);
    }

    if (has_periodicity && kdotl != 0.0) {
      real cs = cos(kdotl);
      real sn = sin(kdotl);

      for (i = 0; i < Nb_2_20; i++) {
        int row = d->HzbComplex.indxHU_to_row[k] + i;
        complex *x_I = x + row;
        complex *y_I = y + row;
        complex *x_I2 = x + row + Nb_2_20;
        complex *y_I2 = y + row + Nb_2_20;
        for (j = 0; j < Nb_2_20; j++) {
          int col = d->HzbComplex.indxHU_to_col[k] + j;
          complex *x_J = x + col;
          complex *y_J = y + col;
          complex *x_J2 = x + col + Nb_2_20;
          complex *y_J2 = y + col + Nb_2_20;
          register double Hu_indx_r =
              cs * d->HzbComplex.Hu[indx].r - sn * d->HzbComplex.Hu[indx].i;
          register double Hu_indx_i =
              sn * d->HzbComplex.Hu[indx].r + cs * d->HzbComplex.Hu[indx].i;

          /* do lower (upper?) triangular part */
          y_I->r += Hu_indx_r * x_J->r - Hu_indx_i * x_J->i;
          y_I->i += Hu_indx_r * x_J->i + Hu_indx_i * x_J->r;
          y_I2->r += Hu_indx_r * x_J2->r - Hu_indx_i * x_J2->i;
          y_I2->i += Hu_indx_r * x_J2->i + Hu_indx_i * x_J2->r;

          /* do upper (lower?) triangular part */
          y_J->r += Hu_indx_r * x_I->r + Hu_indx_i * x_I->i;
          y_J->i += Hu_indx_r * x_I->i - Hu_indx_i * x_I->r;
          y_J2->r += Hu_indx_r * x_I2->r + Hu_indx_i * x_I2->i;
          y_J2->i += Hu_indx_r * x_I2->i - Hu_indx_i * x_I2->r;

          indx++;
        }
      }
    } else {
      for (i = 0; i < Nb_2_20; i++) {
        int row = d->HzbComplex.indxHU_to_row[k] + i;
        complex *x_I = x + row;
        complex *y_I = y + row;
        complex *x_I2 = x + row + Nb_2_20;
        complex *y_I2 = y + row + Nb_2_20;
        for (j = 0; j < Nb_2_20; j++) {
          int col = d->HzbComplex.indxHU_to_col[k] + j;
          complex *x_J = x + col;
          complex *y_J = y + col;
          complex *x_J2 = x + col + Nb_2_20;
          complex *y_J2 = y + col + Nb_2_20;
          register double Hu_indx_r = d->HzbComplex.Hu[indx].r;
          register double Hu_indx_i = d->HzbComplex.Hu[indx].i;

          /* do lower (upper?) triangular part */
          y_I->r += Hu_indx_r * x_J->r - Hu_indx_i * x_J->i;
          y_I->i += Hu_indx_r * x_J->i + Hu_indx_i * x_J->r;
          y_I2->r += Hu_indx_r * x_J2->r - Hu_indx_i * x_J2->i;
          y_I2->i += Hu_indx_r * x_J2->i + Hu_indx_i * x_J2->r;

          /* do upper (lower?) triangular part */
          y_J->r += Hu_indx_r * x_I->r + Hu_indx_i * x_I->i;
          y_J->i += Hu_indx_r * x_I->i - Hu_indx_i * x_I->r;
          y_J2->r += Hu_indx_r * x_I2->r + Hu_indx_i * x_I2->i;
          y_J2->i += Hu_indx_r * x_I2->i - Hu_indx_i * x_I2->r;

          indx++;
        }
      }
    }
  }
#ifdef MPI_TIMEmeasure
  MPI_TOC(T_Hu, Tstart);
  MPI_TIC(Tstart);
#endif

  /* yc += Ho * x */
  thislength = nml_iv_extent(d->HzbComplex.indxHO_to_row);
  for (k = 0; k < thislength; k++) {

    real kdotl = 0.0;
    if (has_periodicity) {
      int aindx = d->HzbComplex.indxHO_to_col[k] / d->NBasisStates;
      kdotl = d->geo.phaseBloch(d->kxL, d->kyL, d->kzL, d->Atom_to_lm[aindx][0],
                                d->Atom_to_lm[aindx][1],
                                d->HzbComplex.indxHO_nbr[k]);
    }

    if (has_periodicity && kdotl != 0.0) {
      real cs = cos(kdotl);
      real sn = sin(kdotl);

      for (i = 0; i < Nb_2_20; i++) {
        for (j = 0; j < Nb_2_20; j++) {
          int indx = Nb_2_20 * (Nb_2_20 * k + i) + j;
          int I = d->HzbComplex.indxHO_to_row[k] + i;
          int J = d->HzbComplex.indxHO_to_col[k] + j;
          register double Ho_indx_r =
              cs * d->HzbComplex.Ho[indx].r - sn * d->HzbComplex.Ho[indx].i;
          register double Ho_indx_i =
              sn * d->HzbComplex.Ho[indx].r + cs * d->HzbComplex.Ho[indx].i;

          yc[I].r += Ho_indx_r * x[J].r - Ho_indx_i * x[J].i;
          yc[I].i += Ho_indx_r * x[J].i + Ho_indx_i * x[J].r;

          yc[I + Nb_2_20].r +=
              Ho_indx_r * x[J + Nb_2_20].r - Ho_indx_i * x[J + Nb_2_20].i;
          yc[I + Nb_2_20].i +=
              Ho_indx_r * x[J + Nb_2_20].i + Ho_indx_i * x[J + Nb_2_20].r;
        }
      }
    } else {
      for (i = 0; i < Nb_2_20; i++) {
        for (j = 0; j < Nb_2_20; j++) {
          int indx = Nb_2_20 * (Nb_2_20 * k + i) + j;
          int I = d->HzbComplex.indxHO_to_row[k] + i;
          int J = d->HzbComplex.indxHO_to_col[k] + j;

          yc[I].r += d->HzbComplex.Ho[indx].r * x[J].r -
                     d->HzbComplex.Ho[indx].i * x[J].i;
          yc[I].i += d->HzbComplex.Ho[indx].r * x[J].i +
                     d->HzbComplex.Ho[indx].i * x[J].r;

          yc[I + Nb_2_20].r += d->HzbComplex.Ho[indx].r * x[J + Nb_2_20].r -
                               d->HzbComplex.Ho[indx].i * x[J + Nb_2_20].i;
          yc[I + Nb_2_20].i += d->HzbComplex.Ho[indx].r * x[J + Nb_2_20].i +
                               d->HzbComplex.Ho[indx].i * x[J + Nb_2_20].r;
        }
      }
    }
  }
#ifdef MPI_TIMEmeasure
  MPI_TOC(T_Ho, Tstart);
  printf("Vector-time:  T_diag=%g T_mag=%g T_Hdu=%g T_Hds=%g T_Hu=%g T_Ho=%g\n",
         T_diag, T_mag, T_Hdu, T_Hds, T_Hu, T_Ho),
      fflush(stdout);
#endif
}

void Hv_stored_10_sp3d5s_HzbComplex(cvectr y, cvectr yc, cvectr x,
                                    qd_struct d) {
  const int Nb_2_20 = 10;
  int atom, indx, i, j, k, thislength;
  bool has_periodicity = d->geo.hasPeriodicity();
  ham_mem_real_ptr H_IJ;

  /* y += Hdd * x */
  indx = 0;
  for (atom = 0; atom < d->HzbComplex.Natom; atom++) {
    for (i = 0; i < d->NBasisStates; i++) {
      y[indx].r += d->HzbComplex.Hdd[indx] * x[indx].r;
      y[indx].i += d->HzbComplex.Hdd[indx] * x[indx].i;
      indx++;
    }
  }

  /* In the spin-less system there is no spin-up and spin-down coupling */
  /* y += HBxy * x */

  /* In the spin-less system we do not have any spin-orbit coupling
     Do not need to consider the matrix Hdu at all. */
  /* y += Hdu * x */

  /* y += Hds * x */
  H_IJ = &d->HzbComplex.Hds[0];
  for (int Zsurf = 0; Zsurf < d->HzbComplex.Natom_surf; Zsurf++) {
    int I, J;
    int offset = d->surfaceAtoms[Zsurf] * d->NBasisStates;

    I = offset + 1;
    J = offset + 2;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;

    H_IJ++;

    I = offset + 1;
    J = offset + 3;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;

    H_IJ++;

    I = offset + 1;
    J = offset + 4;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;

    H_IJ++;

    I = offset + 2;
    J = offset + 3;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;

    H_IJ++;

    I = offset + 2;
    J = offset + 4;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;

    H_IJ++;

    I = offset + 3;
    J = offset + 4;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;

    H_IJ++;
  }

  /* y += Hu * x */
  indx = 0;
  thislength = nml_iv_extent(d->HzbComplex.indxHU_to_row);
  for (k = 0; k < thislength; k++) {

    real kdotl = 0.0;
    if (has_periodicity) {
      int aindx = d->HzbComplex.indxHU_to_col[k] / d->NBasisStates;
      kdotl = d->geo.phaseBloch(d->kxL, d->kyL, d->kzL, d->Atom_to_lm[aindx][0],
                                d->Atom_to_lm[aindx][1],
                                d->HzbComplex.indxHU_nbr[k]);
    }

    if (has_periodicity && kdotl != 0.0) {
      real cs = cos(kdotl);
      real sn = sin(kdotl);

      for (i = 0; i < Nb_2_20; i++) {
        int row = d->HzbComplex.indxHU_to_row[k] + i;
        complex *x_I = x + row;
        complex *y_I = y + row;
        for (j = 0; j < Nb_2_20; j++) {
          int col = d->HzbComplex.indxHU_to_col[k] + j;
          complex *x_J = x + col;
          complex *y_J = y + col;
          register double Hu_indx_r =
              cs * d->HzbComplex.Hu[indx].r - sn * d->HzbComplex.Hu[indx].i;
          register double Hu_indx_i =
              sn * d->HzbComplex.Hu[indx].r + cs * d->HzbComplex.Hu[indx].i;

          /* do lower (upper?) triangular part */
          y_I->r += Hu_indx_r * x_J->r - Hu_indx_i * x_J->i;
          y_I->i += Hu_indx_r * x_J->i + Hu_indx_i * x_J->r;

          /* do upper (lower?) triangular part */
          y_J->r += Hu_indx_r * x_I->r + Hu_indx_i * x_I->i;
          y_J->i += Hu_indx_r * x_I->i - Hu_indx_i * x_I->r;

          indx++;
        }
      }
    } else {
      for (i = 0; i < Nb_2_20; i++) {
        int row = d->HzbComplex.indxHU_to_row[k] + i;
        complex *x_I = x + row;
        complex *y_I = y + row;
        for (j = 0; j < Nb_2_20; j++) {
          int col = d->HzbComplex.indxHU_to_col[k] + j;
          complex *x_J = x + col;
          complex *y_J = y + col;
          register double Hu_indx_r = d->HzbComplex.Hu[indx].r;
          register double Hu_indx_i = d->HzbComplex.Hu[indx].i;

          /* do lower (upper?) triangular part */
          y_I->r += Hu_indx_r * x_J->r - Hu_indx_i * x_J->i;
          y_I->i += Hu_indx_r * x_J->i + Hu_indx_i * x_J->r;

          /* do upper (lower?) triangular part */
          y_J->r += Hu_indx_r * x_I->r + Hu_indx_i * x_I->i;
          y_J->i += Hu_indx_r * x_I->i - Hu_indx_i * x_I->r;

          indx++;
        }
      }
    }
  }

  /* yc += Ho * x */
  thislength = nml_iv_extent(d->HzbComplex.indxHO_to_row);
  for (k = 0; k < thislength; k++) {

    real kdotl = 0.0;
    if (has_periodicity) {
      int aindx = d->HzbComplex.indxHO_to_col[k] / d->NBasisStates;
      kdotl = d->geo.phaseBloch(d->kxL, d->kyL, d->kzL, d->Atom_to_lm[aindx][0],
                                d->Atom_to_lm[aindx][1],
                                d->HzbComplex.indxHO_nbr[k]);
    }

    if (has_periodicity && kdotl != 0.0) {
      real cs = cos(kdotl);
      real sn = sin(kdotl);

      for (i = 0; i < Nb_2_20; i++) {
        for (j = 0; j < Nb_2_20; j++) {
          int indx = Nb_2_20 * (Nb_2_20 * k + i) + j;
          int I = d->HzbComplex.indxHO_to_row[k] + i;
          int J = d->HzbComplex.indxHO_to_col[k] + j;
          register double Ho_indx_r =
              cs * d->HzbComplex.Ho[indx].r - sn * d->HzbComplex.Ho[indx].i;
          register double Ho_indx_i =
              sn * d->HzbComplex.Ho[indx].r + cs * d->HzbComplex.Ho[indx].i;

          yc[I].r += Ho_indx_r * x[J].r - Ho_indx_i * x[J].i;
          yc[I].i += Ho_indx_r * x[J].i + Ho_indx_i * x[J].r;
        }
      }
    } else {
      for (i = 0; i < Nb_2_20; i++) {
        for (j = 0; j < Nb_2_20; j++) {
          int indx = Nb_2_20 * (Nb_2_20 * k + i) + j;
          int I = d->HzbComplex.indxHO_to_row[k] + i;
          int J = d->HzbComplex.indxHO_to_col[k] + j;

          yc[I].r += d->HzbComplex.Ho[indx].r * x[J].r -
                     d->HzbComplex.Ho[indx].i * x[J].i;
          yc[I].i += d->HzbComplex.Ho[indx].r * x[J].i +
                     d->HzbComplex.Ho[indx].i * x[J].r;
        }
      }
    }
  }
}

void Hv_stored_10_HzbComplex(cvectr y, cvectr yc, cvectr x, qd_struct d) {
  const int Nb_2_10 = 5;
  int atom, indx, i, j, k, thislength;
  bool has_periodicity = d->geo.hasPeriodicity();
  ham_mem_real_ptr H_IJ;
  ham_mem_complex_ptr Hc_IJ;

  /* y += Hdd * x */
  indx = 0;
  for (atom = 0; atom < d->HzbComplex.Natom; atom++) {
    for (i = 0; i < d->NBasisStates; i++) {
      y[indx].r += d->HzbComplex.Hdd[indx] * x[indx].r;
      y[indx].i += d->HzbComplex.Hdd[indx] * x[indx].i;
      indx++;
    }
  }

  /* y += HBxy * x */
  if (d->opt.ExecParam.Phys.MagneticFieldOn) {
    int nb = d->NBasisStates / 2;
    for (atom = 0; atom < d->HzbComplex.Natom; atom++) {
      int offset = atom * d->NBasisStates;
      for (int i = 0; i < nb; i++) {
        y[offset + i].r += d->HzbComplex.HBxy.r * x[offset + i + nb].r -
                           d->HzbComplex.HBxy.i * x[offset + i + nb].i;
        y[offset + i].i += d->HzbComplex.HBxy.r * x[offset + i + nb].i +
                           d->HzbComplex.HBxy.i * x[offset + i + nb].r;
        y[offset + i + nb].r += d->HzbComplex.HBxy.r * x[offset + i].r +
                                d->HzbComplex.HBxy.i * x[offset + i].i;
        y[offset + i + nb].i += d->HzbComplex.HBxy.r * x[offset + i].i -
                                d->HzbComplex.HBxy.i * x[offset + i].r;
      }
    }
  }

  /* y += Hdu * x */
  Hc_IJ = &d->HzbComplex.Hdu[0];
  for (atom = 0; atom < d->HzbComplex.Natom; atom++) {
    int offset = atom * d->NBasisStates;
    int I, J;

    I = offset + 2;
    J = offset + 3;
    y[I].r += Hc_IJ->r * x[J].r - Hc_IJ->i * x[J].i;
    y[I].i += Hc_IJ->r * x[J].i + Hc_IJ->i * x[J].r;
    y[J].r += Hc_IJ->r * x[I].r + Hc_IJ->i * x[I].i;
    y[J].i += Hc_IJ->r * x[I].i - Hc_IJ->i * x[I].r;
    Hc_IJ++;

    I = offset + 7;
    J = offset + 8;
    y[I].r += Hc_IJ->r * x[J].r - Hc_IJ->i * x[J].i;
    y[I].i += Hc_IJ->r * x[J].i + Hc_IJ->i * x[J].r;
    y[J].r += Hc_IJ->r * x[I].r + Hc_IJ->i * x[I].i;
    y[J].i += Hc_IJ->r * x[I].i - Hc_IJ->i * x[I].r;
    Hc_IJ++;

    I = offset + 4;
    J = offset + 7;
    y[I].r += Hc_IJ->r * x[J].r - Hc_IJ->i * x[J].i;
    y[I].i += Hc_IJ->r * x[J].i + Hc_IJ->i * x[J].r;
    y[J].r += Hc_IJ->r * x[I].r + Hc_IJ->i * x[I].i;
    y[J].i += Hc_IJ->r * x[I].i - Hc_IJ->i * x[I].r;
    Hc_IJ++;

    I = offset + 4;
    J = offset + 8;
    y[I].r += Hc_IJ->r * x[J].r - Hc_IJ->i * x[J].i;
    y[I].i += Hc_IJ->r * x[J].i + Hc_IJ->i * x[J].r;
    y[J].r += Hc_IJ->r * x[I].r + Hc_IJ->i * x[I].i;
    y[J].i += Hc_IJ->r * x[I].i - Hc_IJ->i * x[I].r;
    Hc_IJ++;

    I = offset + 2;
    J = offset + 9;
    y[I].r += Hc_IJ->r * x[J].r - Hc_IJ->i * x[J].i;
    y[I].i += Hc_IJ->r * x[J].i + Hc_IJ->i * x[J].r;
    y[J].r += Hc_IJ->r * x[I].r + Hc_IJ->i * x[I].i;
    y[J].i += Hc_IJ->r * x[I].i - Hc_IJ->i * x[I].r;
    Hc_IJ++;

    I = offset + 3;
    J = offset + 9;
    y[I].r += Hc_IJ->r * x[J].r - Hc_IJ->i * x[J].i;
    y[I].i += Hc_IJ->r * x[J].i + Hc_IJ->i * x[J].r;
    y[J].r += Hc_IJ->r * x[I].r + Hc_IJ->i * x[I].i;
    y[J].i += Hc_IJ->r * x[I].i - Hc_IJ->i * x[I].r;
    Hc_IJ++;
  }

  /* y += Hds * x */
  H_IJ = &d->HzbComplex.Hds[0];
  for (int Zsurf = 0; Zsurf < d->HzbComplex.Natom_surf; Zsurf++) {
    int I, J;
    int offset = d->surfaceAtoms[Zsurf] * d->NBasisStates;

    I = offset + 1;
    J = offset + 2;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;

    I = offset + 6;
    J = offset + 7;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;
    H_IJ++;

    I = offset + 1;
    J = offset + 3;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;

    I = offset + 6;
    J = offset + 8;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;
    H_IJ++;

    I = offset + 1;
    J = offset + 4;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;

    I = offset + 6;
    J = offset + 9;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;
    H_IJ++;

    I = offset + 4;
    J = offset + 2;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;

    I = offset + 9;
    J = offset + 7;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;
    H_IJ++;

    I = offset + 4;
    J = offset + 3;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;

    I = offset + 9;
    J = offset + 8;
    y[I].r += (*H_IJ) * x[J].r;
    y[I].i += (*H_IJ) * x[J].i;
    y[J].r += (*H_IJ) * x[I].r;
    y[J].i += (*H_IJ) * x[I].i;
    H_IJ++;
  }

  /* y += Hu * x */
  indx = 0;
  thislength = nml_iv_extent(d->HzbComplex.indxHU_to_row);
  for (k = 0; k < thislength; k++) {

    real kdotl = 0.0;
    if (has_periodicity) {
      int aindx = d->HzbComplex.indxHU_to_col[k] / d->NBasisStates;
      kdotl = d->geo.phaseBloch(d->kxL, d->kyL, d->kzL, d->Atom_to_lm[aindx][0],
                                d->Atom_to_lm[aindx][1],
                                d->HzbComplex.indxHU_nbr[k]);
    }

    if (has_periodicity && kdotl != 0.0) {
      real cs = cos(kdotl);
      real sn = sin(kdotl);

      for (i = 0; i < Nb_2_10; i++) {
        int row = d->HzbComplex.indxHU_to_row[k] + i;
        complex *x_I = x + row;
        complex *y_I = y + row;
        complex *x_I2 = x + row + Nb_2_10;
        complex *y_I2 = y + row + Nb_2_10;
        for (j = 0; j < Nb_2_10; j++) {
          int col = d->HzbComplex.indxHU_to_col[k] + j;
          complex *x_J = x + col;
          complex *y_J = y + col;
          complex *x_J2 = x + col + Nb_2_10;
          complex *y_J2 = y + col + Nb_2_10;
          register double Hu_indx_r =
              cs * d->HzbComplex.Hu[indx].r - sn * d->HzbComplex.Hu[indx].i;
          register double Hu_indx_i =
              sn * d->HzbComplex.Hu[indx].r + cs * d->HzbComplex.Hu[indx].i;

          /* do lower (upper?) triangular part */
          y_I->r += Hu_indx_r * x_J->r - Hu_indx_i * x_J->i;
          y_I->i += Hu_indx_r * x_J->i + Hu_indx_i * x_J->r;
          y_I2->r += Hu_indx_r * x_J2->r - Hu_indx_i * x_J2->i;
          y_I2->i += Hu_indx_r * x_J2->i + Hu_indx_i * x_J2->r;

          /* do upper (lower?) triangular part */
          y_J->r += Hu_indx_r * x_I->r + Hu_indx_i * x_I->i;
          y_J->i += Hu_indx_r * x_I->i - Hu_indx_i * x_I->r;
          y_J2->r += Hu_indx_r * x_I2->r + Hu_indx_i * x_I2->i;
          y_J2->i += Hu_indx_r * x_I2->i - Hu_indx_i * x_I2->r;

          indx++;
        }
      }
    } else {
      for (i = 0; i < Nb_2_10; i++) {
        int row = d->HzbComplex.indxHU_to_row[k] + i;
        complex *x_I = x + row;
        complex *y_I = y + row;
        complex *x_I2 = x + row + Nb_2_10;
        complex *y_I2 = y + row + Nb_2_10;
        for (j = 0; j < Nb_2_10; j++) {
          int col = d->HzbComplex.indxHU_to_col[k] + j;
          complex *x_J = x + col;
          complex *y_J = y + col;
          complex *x_J2 = x + col + Nb_2_10;
          complex *y_J2 = y + col + Nb_2_10;
          register double Hu_indx_r = d->HzbComplex.Hu[indx].r;
          register double Hu_indx_i = d->HzbComplex.Hu[indx].i;

          /* do lower (upper?) triangular part */
          y_I->r += Hu_indx_r * x_J->r - Hu_indx_i * x_J->i;
          y_I->i += Hu_indx_r * x_J->i + Hu_indx_i * x_J->r;
          y_I2->r += Hu_indx_r * x_J2->r - Hu_indx_i * x_J2->i;
          y_I2->i += Hu_indx_r * x_J2->i + Hu_indx_i * x_J2->r;

          /* do upper (lower?) triangular part */
          y_J->r += Hu_indx_r * x_I->r + Hu_indx_i * x_I->i;
          y_J->i += Hu_indx_r * x_I->i - Hu_indx_i * x_I->r;
          y_J2->r += Hu_indx_r * x_I2->r + Hu_indx_i * x_I2->i;
          y_J2->i += Hu_indx_r * x_I2->i - Hu_indx_i * x_I2->r;

          indx++;
        }
      }
    }
  }

  /* yc += Ho * x */
  thislength = nml_iv_extent(d->HzbComplex.indxHO_to_row);
  for (k = 0; k < thislength; k++) {

    real kdotl = 0.0;
    if (has_periodicity) {
      int aindx = d->HzbComplex.indxHO_to_col[k] / d->NBasisStates;
      kdotl = d->geo.phaseBloch(d->kxL, d->kyL, d->kzL, d->Atom_to_lm[aindx][0],
                                d->Atom_to_lm[aindx][1],
                                d->HzbComplex.indxHO_nbr[k]);
    }

    if (has_periodicity && kdotl != 0.0) {
      real cs = cos(kdotl);
      real sn = sin(kdotl);

      for (i = 0; i < Nb_2_10; i++) {
        int I = d->HzbComplex.indxHO_to_row[k] + i;
        for (j = 0; j < Nb_2_10; j++) {
          int indx = Nb_2_10 * (Nb_2_10 * k + i) + j;
          int J = d->HzbComplex.indxHO_to_col[k] + j;
          register double Ho_indx_r =
              cs * d->HzbComplex.Ho[indx].r - sn * d->HzbComplex.Ho[indx].i;
          register double Ho_indx_i =
              sn * d->HzbComplex.Ho[indx].r + cs * d->HzbComplex.Ho[indx].i;

          yc[I].r += Ho_indx_r * x[J].r - Ho_indx_i * x[J].i;
          yc[I].i += Ho_indx_r * x[J].i + Ho_indx_i * x[J].r;

          yc[I + Nb_2_10].r +=
              Ho_indx_r * x[J + Nb_2_10].r - Ho_indx_i * x[J + Nb_2_10].i;
          yc[I + Nb_2_10].i +=
              Ho_indx_r * x[J + Nb_2_10].i + Ho_indx_i * x[J + Nb_2_10].r;
        }
      }
    } else {
      for (i = 0; i < Nb_2_10; i++) {
        int I = d->HzbComplex.indxHO_to_row[k] + i;
        for (j = 0; j < Nb_2_10; j++) {
          int indx = Nb_2_10 * (Nb_2_10 * k + i) + j;
          int J = d->HzbComplex.indxHO_to_col[k] + j;

          yc[I].r += d->HzbComplex.Ho[indx].r * x[J].r -
                     d->HzbComplex.Ho[indx].i * x[J].i;
          yc[I].i += d->HzbComplex.Ho[indx].r * x[J].i +
                     d->HzbComplex.Ho[indx].i * x[J].r;

          yc[I + Nb_2_10].r += d->HzbComplex.Ho[indx].r * x[J + Nb_2_10].r -
                               d->HzbComplex.Ho[indx].i * x[J + Nb_2_10].i;
          yc[I + Nb_2_10].i += d->HzbComplex.Ho[indx].r * x[J + Nb_2_10].i +
                               d->HzbComplex.Ho[indx].i * x[J + Nb_2_10].r;
        }
      }
    }
  }
}
