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
$Header: /repo/nemo3d/src/base/Hmult_stored__Hv_stored_20_Hzb_V1.c,v 1.5
2005/02/04 20:13:06 gekco Exp $
*****************************************************************************/

/* Matrix-vector multiplication customized for the
   spin-up and spin-down sp3d5s* model Bands_20_sp3d5ss_spin */
void Hv_stored_20_Hzb_V1(cvectr y, cvectr yc, cvectr x, qd_struct d) {
  const int Nb_2 = 10;
  int atom, indx, i, j, k, thislength;
  bool has_periodicity = d->geo.hasPeriodicity();
  ham_mem_real_ptr H_IJ;
  ham_mem_complex_ptr Hc_IJ;
#ifdef MPI_TIMEmeasure
  MPI_TIME_INIT(Tstart);
  MPI_TIC(Tstart);
#endif;

  /* y += Hdd * x */
  indx = 0;
  for (atom = 0; atom < d->Hzb.Natom; atom++) {
    for (i = 0; i < d->NBasisStates; i++) {
      y[indx].r += d->Hzb.Hdd[indx] * x[indx].r;
      y[indx].i += d->Hzb.Hdd[indx] * x[indx].i;
      indx++;
    }
  }

#ifdef MPI_TIMEmeasure
  MPI_TOC(T_V1_diag, Tstart);
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
  MPI_TOC(T_V1_mag, Tstart);
  MPI_TIC(Tstart);
#endif

  /* y += Hdu * x */
  Hc_IJ = &d->Hzb.Hdu[0];
  for (atom = 0; atom < d->Hzb.Natom; atom++) {
    int offset = atom * d->NBasisStates;
    int I, J;

    register double hc_ij1_r = Hc_IJ->r;
    register double hc_ij1_i = Hc_IJ->i;

    I = offset + 2;
    J = offset + 3;
    register double x_i_r = x[I].r;
    register double x_i_i = x[I].i;
    register double x_j_r = x[J].r;
    register double x_j_i = x[J].i;
    y[I].r += hc_ij1_r * x_j_r - hc_ij1_i * x_j_i;
    y[I].i += hc_ij1_r * x_j_i + hc_ij1_i * x_j_r;
    y[J].r += hc_ij1_r * x_i_r + hc_ij1_i * x_i_i;
    y[J].i += hc_ij1_r * x_i_i - hc_ij1_i * x_i_r;
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
  MPI_TOC(T_V1_Hdu, Tstart);
  MPI_TIC(Tstart);
#endif

  /* y += Hds * x */
  H_IJ = &d->Hzb.Hds[0];
  for (int Zsurf = 0; Zsurf < d->Hzb.Natom_surf; Zsurf++) {
    int I, J;
    int offset = d->surfaceAtoms[Zsurf] * d->NBasisStates;
    register double h_ij1 = (*H_IJ);
    ;
    H_IJ++;
    register double h_ij2 = (*H_IJ);
    ;
    H_IJ++;
    register double h_ij3 = (*H_IJ);
    ;
    H_IJ++;
    register double h_ij4 = (*H_IJ);
    ;
    H_IJ++;
    register double h_ij5 = (*H_IJ);
    ;
    H_IJ++;

    I = offset + 1;
    register double x_I_r = x[I].r;
    register double x_I_i = x[I].i;
    J = offset + 2;
    y[I].r += h_ij1 * x[J].r;
    y[I].i += h_ij1 * x[J].i;
    y[J].r += h_ij1 * x_I_r;
    y[J].i += h_ij1 * x_I_i;
    J = offset + 3;
    y[I].r += h_ij2 * x[J].r;
    y[I].i += h_ij2 * x[J].i;
    y[J].r += h_ij2 * x_I_r;
    y[J].i += h_ij2 * x_I_i;
    J = offset + 4;
    y[I].r += h_ij3 * x[J].r;
    y[I].i += h_ij3 * x[J].i;
    y[J].r += h_ij3 * x_I_r;
    y[J].i += h_ij3 * x_I_i;

    I = offset + 11;
    x_I_r = x[I].r;
    x_I_i = x[I].i;
    J = offset + 12;
    y[I].r += h_ij1 * x[J].r;
    y[I].i += h_ij1 * x[J].i;
    y[J].r += h_ij1 * x_I_r;
    y[J].i += h_ij1 * x_I_i;
    J = offset + 13;
    y[I].r += h_ij2 * x[J].r;
    y[I].i += h_ij2 * x[J].i;
    y[J].r += h_ij2 * x_I_r;
    y[J].i += h_ij2 * x_I_i;
    J = offset + 14;
    y[I].r += h_ij3 * x[J].r;
    y[I].i += h_ij3 * x[J].i;
    y[J].r += h_ij3 * x_I_r;
    y[J].i += h_ij3 * x_I_i;

    I = offset + 4;
    J = offset + 2;
    y[I].r += h_ij4 * x[J].r;
    y[I].i += h_ij4 * x[J].i;
    y[J].r += h_ij4 * x[I].r;
    y[J].i += h_ij4 * x[I].i;

    I = offset + 4;
    J = offset + 3;
    y[I].r += h_ij5 * x[J].r;
    y[I].i += h_ij5 * x[J].i;
    y[J].r += h_ij5 * x[I].r;
    y[J].i += h_ij5 * x[I].i;

    I = offset + 14;
    J = offset + 12;
    y[I].r += h_ij4 * x[J].r;
    y[I].i += h_ij4 * x[J].i;
    y[J].r += h_ij4 * x[I].r;
    y[J].i += h_ij4 * x[I].i;

    I = offset + 14;
    J = offset + 13;
    y[I].r += h_ij5 * x[J].r;
    y[I].i += h_ij5 * x[J].i;
    y[J].r += h_ij5 * x[I].r;
    y[J].i += h_ij5 * x[I].i;
  }

#ifdef MPI_TIMEmeasure
  MPI_TOC(T_V1_Hds, Tstart);
  MPI_TIC(Tstart);
#endif

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

#ifdef MPI_TIMEmeasure
  MPI_TOC(T_V1_Hu, Tstart);
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
  MPI_TOC(T_V1_Ho, Tstart);
  /* printf("Vector-time:  T_V1_diag=%g T_V1_mag=%g T_V1_Hdu=%g T_V1_Hds=%g
   * T_V1_Hu=%g T_V1_Ho=%g\n",T_V1_diag, T_V1_mag, T_V1_Hdu, T_V1_Hds, T_V1_Hu,
   * T_VV1Ho), fflush(stdout); */
#endif
}
