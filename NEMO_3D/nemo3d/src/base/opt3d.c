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
$Header: /repo/nemo3d/src/base/opt3d.c,v 1.7 2008/06/21 00:12:21 gekco Exp $
*****************************************************************************/

#include "opt3d.h"

real calc_rate(qd_struct d, rvectr init, rvectr finl, real px, real py,
               real pz) {
  int i, j, k, l, m, n;
  real Tf, mxm, mxp, sumx, sumy, sumz, sumi, sumf;

  sumx = sumy = sumz = 0;

  Tf = -HBAR_MASSFACTOR / d->geo.lattice_x / d->geo.lattice_x;

  /* Normalize initial and final states */
  sumi = sumf = 0;
  for (l = 0; l < d->geo.N_Cell; l++) {
    sumi += init[l] * init[l];
    sumf += finl[l] * finl[l];
  }
  sumi = sqrt(sumi);
  sumf = sqrt(sumf);
  for (l = 0; l < d->geo.N_Cell; l++) {
    init[l] /= sumi;
    finl[l] /= sumf;
  }

  sumi = 0;
  sumf = 0;
  for (l = 0; l < d->geo.N_Cell; l++) {
    sumi += finl[l] * finl[l];
    sumf += finl[l] * init[l];
    i = d->geo.cell__ijk[l][0];
    j = d->geo.cell__ijk[l][1];
    k = d->geo.cell__ijk[l][2];

    if (px != 0) {
      m = d->geo.ijk__cell[i + 1][j][k];
      n = d->geo.ijk__cell[i - 1][j][k];

      if (m >= 0) {
        mxp = (d->mstar[d->mat[i + 1][j][k]] + d->mstar[d->mat[i][j][k]]) / 2.0;
        sumx += finl[l] * Tf / mxp * init[m];
      }

      if (n >= 0) {
        mxm = (d->mstar[d->mat[i - 1][j][k]] + d->mstar[d->mat[i][j][k]]) / 2.0;
        sumx -= finl[l] * Tf / mxm * init[n];
      }
    }
    if (py != 0) {
      m = d->geo.ijk__cell[i][j + 1][k];
      n = d->geo.ijk__cell[i][j - 1][k];

      if (m >= 0) {
        mxp = (d->mstar[d->mat[i][j + 1][k]] + d->mstar[d->mat[i][j][k]]) / 2.0;
        sumy += finl[l] * Tf / mxp * init[m];
      }

      if (n >= 0) {
        mxm = (d->mstar[d->mat[i][j - 1][k]] + d->mstar[d->mat[i][j][k]]) / 2.0;
        sumy -= finl[l] * Tf / mxm * init[n];
      }
    }
    if (pz != 0) {
      m = d->geo.ijk__cell[i][j][k + 1];
      n = d->geo.ijk__cell[i][j][k - 1];

      if (m >= 0) {
        mxp = (d->mstar[d->mat[i][j][k + 1]] + d->mstar[d->mat[i][j][k]]) / 2.0;
        sumz += finl[l] * Tf / mxp * init[m];
      }

      if (n >= 0) {
        mxm = (d->mstar[d->mat[i][j][k - 1]] + d->mstar[d->mat[i][j][k]]) / 2.0;
        sumz -= finl[l] * Tf / mxm * init[n];
      }
    }
  }

  return px * px * sumx * sumx + py * py * sumy * sumy + pz * pz * sumz * sumz;
}

/* Calculates dipole optical matrix elements.  Polarization is given
   by fractions px, py, pz.  Generates file that contains transition
   energies in ascending order with associated state indices and
   optical matrix elements.  Also prints file to plot bright and
   dark transitions vs. energy.  */
int calc_opt_mat_old(qd_struct d, real px, real py, real pz) {
  int i, j, n = d->neigv, nt, i0, i1;
  real a;
  rmatrix chi;
  FILE *f_tran, *f_tran1;

  if (d->opmat == NULL) {
    d->opmat = Rmatrix(n, n);
    d->te = Rvectr((n * n - n) / 2);
    d->tindx = Imatrix((n * n - n) / 2, 2);
  }
  chi = d->wf;

  nt = 0;
  for (i = 0; i < n - 1; i++)
    for (j = i + 1; j < n; j++) {
      d->opmat[i][j] = calc_rate(d, chi[i], chi[j], px, py, pz);
      d->te[nt] = d->eigv[j].r - d->eigv[i].r;
      d->tindx[nt][0] = i;
      d->tindx[nt][1] = j;
      nt++;
    }

  /* Sort transition energies in ascending order */
  for (j = 1; j < nt; j++) {
    a = d->te[j];
    i0 = d->tindx[j][0];
    i1 = d->tindx[j][1];
    i = j - 1;
    while (i >= 0 && d->te[i] > a) {
      d->te[i + 1] = d->te[i];
      d->tindx[i + 1][0] = d->tindx[i][0];
      d->tindx[i + 1][1] = d->tindx[i][1];
      i--;
    }
    d->te[i + 1] = a;
    d->tindx[i + 1][0] = i0;
    d->tindx[i + 1][1] = i1;
  }

  f_tran = fopen("transitions.out", "w");
  f_tran1 = fopen("tran_vs_e.out", "w");
  /* Print out transition energies with corresponding states.  Also
     print file containing delta functions at each transition */
  for (i = 0; i < nt; i++) {
    fprintf(f_tran, "%i %e %i %i %e\n", i, d->te[i], d->tindx[i][0],
            d->tindx[i][1], d->opmat[d->tindx[i][0]][d->tindx[i][1]]);
    fprintf(f_tran1, "%e %e\n", d->te[i] - 1e-12, 1e-10);
    if (d->opmat[d->tindx[i][0]][d->tindx[i][1]] > 1e-6)
      fprintf(f_tran1, "%e %e\n", d->te[i], 1.0);
    else
      fprintf(f_tran1, "%e %e\n", d->te[i], 0.1);
    fprintf(f_tran1, "%e %e\n", d->te[i] + 1e-12, 1e-10);
  }
  fclose(f_tran);
  fclose(f_tran1);

  return (0);
}

/* Calculates dipole optical matrix elements.  Polarization is input
   in spherical coordinates with angles in degrees.  Generates file
   that contains transition energies in ascending order with associated
   state indices and optical matrix elements.  Also prints file to plot
   bright and dark transitions vs. energy.  */
int calc_opt_mat(qd_struct d, real er, real ephi, real etheta) {
  int i, j, n = d->neigv, nt, i0, i1;
  real a, px, py, pz, theta, phi;
  rmatrix chi;
  FILE *f_tran, *f_tran1;

  if (d->opt.Dev.band_model != Dev_struct::Bands_1_s_nospin) {
    return 0;
  }
  if (mpi_n3d_numprocs != 1) {
    return 0;
  }

  if (d->wf == NULL)
    d->wf = Rmatrix(n, d->geo.N_Cell);
  if (d->eigv == NULL)
    d->eigv = Cvectr(n);

  {
    int i, j;
    cvectr cvec = NULL;
    rvectr rvec = NULL;
    for (i = 0; i < n; i++) {
      cvec = d->wfc[i];
      rvec = d->wf[i];
      d->eigv[i].r = d->WF[i].E * HAMILTONIAN_SCALE_VALUE;

      for (j = 0; j < d->geo.N_Cell; j++) {
        // rvec[j] = sqrt(cvec[j].r*cvec[j].r + cvec[j].i*cvec[j].i);
        rvec[j] = cvec[j].r;
      }
    }
  }

  if (d->opmat == NULL) {
    d->opmat = Rmatrix(n, n);
    d->te = Rvectr((n * n - n) / 2);
    d->tindx = Imatrix((n * n - n) / 2, 2);
  }
  chi = d->wf;

  theta = etheta / 180 * Pi;
  phi = ephi / 180 * Pi;

  /*
      px = er*sin(theta)*cos(phi);
      py = er*sin(theta)*sin(phi);
      pz = er*cos(theta);
  */

  px = er * cos(theta) * cos(phi);
  py = er * cos(theta) * sin(phi);
  pz = er * sin(theta);

  nt = 0;
  for (i = 0; i < n - 1; i++)
    for (j = i + 1; j < n; j++) {
      d->opmat[i][j] = calc_rate(d, chi[i], chi[j], px, py, pz);
      d->te[nt] = d->eigv[j].r - d->eigv[i].r;
      d->tindx[nt][0] = i;
      d->tindx[nt][1] = j;
      nt++;
    }

  /* Sort transition energies in ascending order */
  for (j = 1; j < nt; j++) {
    a = d->te[j];
    i0 = d->tindx[j][0];
    i1 = d->tindx[j][1];
    i = j - 1;
    while (i >= 0 && d->te[i] > a) {
      d->te[i + 1] = d->te[i];
      d->tindx[i + 1][0] = d->tindx[i][0];
      d->tindx[i + 1][1] = d->tindx[i][1];
      i--;
    }
    d->te[i + 1] = a;
    d->tindx[i + 1][0] = i0;
    d->tindx[i + 1][1] = i1;
  }

  f_tran = fopen("transitions.out", "w");
  f_tran1 = fopen("tran_vs_e.out", "w");
  /* Print out transition energies with corresponding states.  Also
     print file containing delta functions at each transition */
  for (i = 0; i < nt; i++) {
    fprintf(f_tran, "%i %e %i %i %e\n", i, d->te[i], d->tindx[i][0],
            d->tindx[i][1], d->opmat[d->tindx[i][0]][d->tindx[i][1]]);
    fprintf(f_tran1, "%e %e\n", d->te[i] - 1e-12, 1e-10);
    if (d->opmat[d->tindx[i][0]][d->tindx[i][1]] > 1e-6)
      fprintf(f_tran1, "%e %e\n", d->te[i], 1.0);
    else
      fprintf(f_tran1, "%e %e\n", d->te[i], 0.1);
    fprintf(f_tran1, "%e %e\n", d->te[i] + 1e-12, 1e-10);
  }
  fclose(f_tran);
  fclose(f_tran1);

  rm_rmatrix(&d->wf);

  return (0);
}

/* Calculates dipole optical absorption rates between the first n
   states.  Polarization is given by fractions px, py, pz. */
int calc_opt_rates(qd_struct d, real er, real ephi, real etheta) {
  int i, j, ie, nt, n = d->neigv, cnt;
  real Gamma, e, de, den, sum, temp, ef, f, mf;
  rmatrix f1mf;
  static int p = 0;

  char line[80];
  FILE *f_rve, *f_rvec, *f_rveci;

  /* Check output flag */
  if (d->opt.ExecParam.Absorp.AbsCoeff) {
    if (p == 0) {
      p = 1;

      f_rve = fopen("AbsCoeff_2d.m", "w");

      fprintf(f_rve, "hold on;\ntitle('Absorption Coefficient');\n");
      fprintf(f_rve, "clear x%d y%d;\n", p, p);
    } else {
      p++;
      f_rve = fopen("AbsCoeff_2d.m", "a");
      fprintf(f_rve, "clear x%d y%d;\n", p, p);
    }
  }

  /* If optical matrix elements haven't been calculated, do so */
  if (d->opmat == NULL)
    calc_opt_mat(d, er, ephi, etheta);
  if (d->opmat == NULL)
    return 0;

  f1mf = Rmatrix(n, n);

  Gamma = d->opt.ExecParam.Absorp.HbarOverTau;
  temp = d->opt.ExecParam.Abs_Ramp.Etemp;
  ef = d->opt.ExecParam.Abs_Ramp.Efermi;
  /* measure the fermi level from the bottom-most state */
  ef = d->eigv[0].r + ef;
  nt = 0;
  for (i = 0; i < n - 1; i++)
    for (j = i + 1; j < n; j++) {
      f = 1 / (1 + exp((d->eigv[i].r - ef) / (Kboltzmann * temp)));
      mf = 1.0 - 1 / (1 + exp((d->eigv[j].r - ef) / (Kboltzmann * temp)));
      f1mf[i][j] = f * mf;
      nt++;
    }

  /* Calculate transition rate vs. energy */
  cnt = 1;
  ie = 0;
  de = Gamma / 2.0;
  e = d->te[ie] - 10.0 * Gamma;
  do {
    if (ie < nt) {
      if (d->te[ie] < e + de) {
        e = d->te[ie];
        ie++;
      } else
        e += de;
    } else
      e += de;

    sum = 0;
    for (i = 0; i < n - 1; i++)
      for (j = i + 1; j < n; j++) {
        den = pow(d->eigv[j].r - d->eigv[i].r - e, 2) + Gamma * Gamma;
        sum += Gamma * d->opmat[i][j] * f1mf[i][j] / den;
      }
    if (d->opt.ExecParam.Absorp.AbsCoeff)
      fprintf(f_rve, "x%d(%i)=%e; y%d(%i)=%e;\n", p, cnt, e, p, cnt, sum);
    cnt++;
  } while (e < d->te[nt - 1] + 10 * Gamma);

  if (d->opt.ExecParam.Absorp.AbsCoeff) {
    fprintf(f_rve, "plot (x%d,y%d,'.');\n", p, p);
    fprintf(f_rve, "xlabel ('Energy (eV)');\n");
    fprintf(f_rve, "ylabel ('Log Absorption (arb units)');\n");
    fclose(f_rve);
  }

  rm_rmatrix(&f1mf);

  return (0);
}

#define swapReal(g, h)                                                         \
  {                                                                            \
    real y = (g);                                                              \
    (g) = (h);                                                                 \
    (h) = y;                                                                   \
  }

static void quick_rsort(real *x, int beg, int end) {
  int p, i;

  if (beg >= end)
    return;

  swapReal(x[beg], x[(beg + end) / 2]);

  p = beg;
  for (i = beg + 1; i <= end; i++)
    if (x[i] < x[beg]) {
      ++p;
      swapReal(x[i], x[p]);
    }

  swapReal(x[beg], x[p]);

  quick_rsort(x, beg, p - 1);
  quick_rsort(x, p + 1, end);
}

rvectr make_egrid_lineshape(rvectr trans, int n_trans, real Gamma, int *ne) {
#define GAMMARES_FACTOR 10.0 /* 10 */
#define GAMMA_RANGE 5.0      /* 50 */
  real de_small = Gamma / GAMMARES_FACTOR;
  real de_range = GAMMA_RANGE * Gamma;
  real emin = trans[0] - de_range;
  real emax = trans[n_trans - 1] + de_range;
  int n_lt_rt = (int)(de_range / de_small) + 1;
  int n_homog = 500;
  int n_total = n_homog + n_trans * (2 * n_lt_rt + 1);
  rvectr egrid = NULL;
  rvectr egrid_out = NULL;
  int i_trans = 0;
  int i_e = 0, i;
  real de_large = 0;
  real erange_ext = 0.1 * (emax - emin);

  /* gekco note
   * June 19 2008
   * Really the issue is here that the number of transitions goes up as the
   * square of the number of states that are computed and not all of these
   * transitions are bright, yet there is an energy grid that will be associated
   * with it to resolve all possible transitions.
   *
   * A better solution woulf be to improve the code to really just
   * resolve the bright transitions and not all possible transitions.
  printf("\n\n****** make_egrid_lineshape   optical absorptiongrid
  n_total=%d\n",n_total);
   */

  emin = emin - erange_ext;
  emax = emax + erange_ext;

  de_large = (emax - emin) / (n_homog + 1);

  egrid = Rvectr(n_total);
  i_trans = 0;
  i_e = 0;
  for (i = -n_lt_rt; i <= n_lt_rt; i++) {
    egrid[i_e] = trans[i_trans] + i * de_small;
    i_e++;
  }
  i_trans++;

  while (i_trans < n_trans) {
    if ((trans[i_trans] - trans[i_trans - 1]) > de_small / 2) {
      for (i = -n_lt_rt; i <= n_lt_rt; i++) {
        egrid[i_e] = trans[i_trans] + i * de_small;
        i_e++;
      }
    }
    i_trans++;
  }

  for (i = 1; i <= n_homog; i++) {
    egrid[i_e] = emin + (i * de_large);
    i_e++;
  }

  quick_rsort(&egrid[0], 0, i_e - 1);

  if (i_e < n_total) {
    egrid_out = Rvectr(i_e);
    for (i = 0; i < i_e; i++) {
      egrid_out[i] = egrid[i];
    }
    rm_rvectr(&egrid);
  } else {
    egrid_out = egrid;
  }

  *ne = i_e;

  return egrid_out;
}

/* Calculates dipole optical absorption rates between the first n
    states.  Polarization is given by fractions px, py, pz. */
real calc_opt_rates_int(qd_struct d, real er, real ephi, real etheta) {
  int i, j, ie, nt, n = d->neigv, cnt, ne;
  real Gamma, e, de, den, sum, sum_int, temp, ef, f, mf, em1, summ1;
  rmatrix f1mf;
  char line[80];
  /* char s1, s2, s3, chdum[20]; */
  FILE *f_rve = NULL, *f_rvec = NULL, *f_rveci = NULL, *f_xmgr = NULL;
  static int p = 0;
  rvectr egrid = NULL;
  char label[100], *filename = NULL;
  real px, py, pz;

  Gamma = d->opt.ExecParam.Absorp.HbarOverTau;
  temp = d->opt.ExecParam.Abs_Ramp.Etemp;
  ef = d->opt.ExecParam.Abs_Ramp.Efermi;
  /* measure the fermi level from the bottom-most state */
  ef = d->eigv[0].r + ef;
  /* Check output flag */
  if (d->opt.ExecParam.Absorp.AbsCoeff) {

    sprintf(label, "nd_absor_%d", p);
    filename = n3d_strdup_n(d->inputfile);
    n3d_FileTypeSet(&filename, label, TRUE);
    f_rve = fopen(filename, "w");
    str_free(filename);
    fprintf(f_rve, "# Ef_rel=%g Ef_abs=%g T=%g Phi=%g Theta=%g\n",
            ef - d->eigv[0].r, ef, temp, ephi, etheta);
    p++;
  }

  f1mf = Rmatrix(n, n);

  nt = 0;
  for (i = 0; i < n - 1; i++)
    for (j = i + 1; j < n; j++) {
      f = 1 / (1 + exp((d->eigv[i].r - ef) / (Kboltzmann * temp)));
      mf = 1.0 - 1 / (1 + exp((d->eigv[j].r - ef) / (Kboltzmann * temp)));
      f1mf[i][j] = f * mf;
      nt++;
    }

  /* Get matrix elements */
  if (d->opmat_pxyz == NULL)
    calc_opt_mat_cmplx(d, ephi, etheta);
  else {
    real theta, phi;
    theta = etheta / 180 * Pi;
    phi = ephi / 180 * Pi;

    /*
           px = sin(theta)*cos(phi);
           py = sin(theta)*sin(phi);
           pz = cos(theta);
    */

    px = cos(theta) * cos(phi);
    py = cos(theta) * sin(phi);
    pz = sin(theta);

    for (i = 0; i < n - 1; i++) {
      for (j = i + 1; j < n; j++) {
        d->opmat_pxyz[i][j][3] =
            px * px * d->opmat_pxyz[i][j][0] * d->opmat_pxyz[i][j][0] +
            py * py * d->opmat_pxyz[i][j][1] * d->opmat_pxyz[i][j][1] +
            pz * pz * d->opmat_pxyz[i][j][2] * d->opmat_pxyz[i][j][2];
      }
    }
  }

  /* Calculate transition rate vs. energy */
  egrid = make_egrid_lineshape(d->te, nt, Gamma, &ne);
  cnt = 1;
  sum_int = 0;

  e = egrid[0];
  em1 = e;
  summ1 = 0.0;

  for (ie = 0; ie < ne; ie++) {
    e = egrid[ie];
    sum = 0;
    for (i = 0; i < n - 1; i++)
      for (j = i + 1; j < n; j++) {
        den = pow(d->eigv[j].r - d->eigv[i].r - e, 2) + Gamma * Gamma;
        sum += Gamma * d->opmat_pxyz[i][j][3] * f1mf[i][j] / den;
      }

    if (d->opt.ExecParam.Absorp.AbsCoeff)
      fprintf(f_rve, "%e %e\n", e, sum);

    if (fabs(e) > 1e-15)
      sum_int += 0.5 * (e - em1) * (sum + summ1);

    em1 = e;
    summ1 = sum;
  }

  rm_rmatrix(&f1mf);

  if (d->opt.ExecParam.Absorp.AbsCoeff) {
    fclose(f_rve);
  }

  rm_rvectr(&egrid);

  return sum_int;
}

int calc_opt_rates_int_ramp(qd_struct d, real *target, real start, real end,
                            int n) {
  int i, p;
  real aint;
  /* real px, py, pz; */
  char line[80];
  FILE *f_rve = NULL, *f_rvec = NULL, *f_rveci = NULL, *f_abs = NULL;
  char label[100], *filename = NULL;

  if (d->opt.ExecParam.Absorp.IntegrAbs) {
    sprintf(label, "nd_int_absor");
    filename = n3d_strdup_n(d->inputfile);
    n3d_FileTypeSet(&filename, label, TRUE);
    f_rve = fopen(filename, "w");
    str_free(filename);
    switch (d->opt.ExecParam.Absorp.ramp) {
    case Absorp_struct::fermi:
      fprintf(f_rve, "# Fermi Level (eV)\n");
      break;
    case Absorp_struct::temp:
      fprintf(f_rve, "# Temperature (K)\n");
      break;
    case Absorp_struct::strength:
      fprintf(f_rve, "# Estrength (V/nm)\n");
      break;
    case Absorp_struct::phi:
      fprintf(f_rve, "# Phi (deg)\n");
      break;
    case Absorp_struct::theta:
      fprintf(f_rve, "# Theta (deg)\n");
      break;
    default:
      die("Error in absorption ramp\n");
      break;
    }
    fprintf(f_rve, "# Absorption (arb units)\n");
  }

  for (i = 0; i < n; i++) {
#undef LOG_RAMP
#ifdef LOG_RAMP
    if (i == 0)
      *target = start;
    else
      *target = start + pow(pow(end - start, 1.0 / (n - 1)), i);
#else
    *target = i * (end - start) / (n - 1) + start;
#endif
    aint = calc_opt_rates_int(d, d->opt.ExecParam.Abs_Ramp.Estrength,
                              d->opt.ExecParam.Abs_Ramp.Ephi,
                              d->opt.ExecParam.Abs_Ramp.Etheta);
    if (d->opt.ExecParam.Absorp.IntegrAbs)
      fprintf(f_rve, "%e %e\n", *target, aint);
  }

  if (d->opt.ExecParam.Absorp.IntegrAbs) {
    fclose(f_rve);
  }
  return (0);
}

/* Calculates dipole optical matrix elements.  Polarization is input
   in spherical coordinates with angles in degrees.  Generates file
   that contains transition energies in ascending order with associated
   state indices and optical matrix elements.  Also prints file to plot
   bright and dark transitions vs. energy.
   calculation is based on complex wavefunctions instead of the original
   code which was based on real wavefunctions
*/
int calc_opt_mat_cmplx(qd_struct d, real phi, real theta) {
  int i, j, n = d->neigv, nt, i0, i1, l;
  real a, px, py, pz;
  cmatrix chi;
  cvectr chivec;
  real sumi;
  FILE *f_tran, *f_tran_px, *f_tran_py, *f_tran_pz, *f_tran_pa;
  real erange_ext;
  char *filename = NULL;

  if (d->opt.Dev.band_model != Dev_struct::Bands_1_s_nospin) {
    return 0;
  }
  if (mpi_n3d_numprocs != 1) {
    return 0;
  }

  if (d->wf == NULL)
    d->wf = Rmatrix(n, d->geo.N_Cell);
  if (d->eigv == NULL)
    d->eigv = Cvectr(n);

  {
    int i, j;
    cvectr cvec = NULL;
    rvectr rvec = NULL;
    for (i = 0; i < n; i++) {
      cvec = d->wfc[i];
      rvec = d->wf[i];
      d->eigv[i].r = d->WF[i].E * HAMILTONIAN_SCALE_VALUE;

      for (j = 0; j < d->geo.N_Cell; j++) {
        rvec[j] = sqrt(cvec[j].r * cvec[j].r + cvec[j].i * cvec[j].i);
      }
    }
  }

  if (d->opmat_pxyz == NULL) {
    d->opmat_pxyz = R3tensor(0, n - 1, 0, n - 1, 0, 3);
    d->te = Rvectr((n * n - n) / 2);
    d->tindx = Imatrix((n * n - n) / 2, 2);
  }
  chi = d->wfc;

  for (i = 0; i < n; i++) {
    chivec = chi[i];
    sumi = 0;
    for (l = 0; l < d->geo.N_Cell; l++) {
      sumi += chivec[l].r * chivec[l].r + chivec[l].i * chivec[l].i;
    }
    sumi = sqrt(sumi);
    for (l = 0; l < d->geo.N_Cell; l++) {
      chivec[l].r /= sumi;
      chivec[l].i /= sumi;
    }
  }

  theta = theta / 180 * Pi;
  phi = phi / 180 * Pi;

  /*
      px = sin(theta)*cos(phi);
      py = sin(theta)*sin(phi);
      pz = cos(theta);
  */

  px = cos(theta) * cos(phi);
  py = cos(theta) * sin(phi);
  pz = sin(theta);

  nt = 0;
  for (i = 0; i < n - 1; i++)
    for (j = i + 1; j < n; j++) {
      calc_rate_cmplx(&d->opmat_pxyz[i][j][0], d, chi[i], chi[j], px, py, pz);
      d->te[nt] = d->eigv[j].r - d->eigv[i].r;
      d->tindx[nt][0] = i;
      d->tindx[nt][1] = j;
      nt++;
    }

  /* Sort transition energies in ascending order */
  for (j = 1; j < nt; j++) {
    a = d->te[j];
    i0 = d->tindx[j][0];
    i1 = d->tindx[j][1];
    i = j - 1;
    while (i >= 0 && d->te[i] > a) {
      d->te[i + 1] = d->te[i];
      d->tindx[i + 1][0] = d->tindx[i][0];
      d->tindx[i + 1][1] = d->tindx[i][1];
      i--;
    }
    d->te[i + 1] = a;
    d->tindx[i + 1][0] = i0;
    d->tindx[i + 1][1] = i1;
  }

  filename = n3d_strdup_n(d->inputfile);
  n3d_FileTypeSet(&filename, "nd_trans_list", TRUE);
  f_tran = fopen(filename, "w");
  str_free(filename);

  filename = n3d_strdup_n(d->inputfile);
  n3d_FileTypeSet(&filename, "nd_trans_vs_e_px", TRUE);
  f_tran_px = fopen(filename, "w");
  str_free(filename);

  filename = n3d_strdup_n(d->inputfile);
  n3d_FileTypeSet(&filename, "nd_trans_vs_e_py", TRUE);
  f_tran_py = fopen(filename, "w");
  str_free(filename);

  filename = n3d_strdup_n(d->inputfile);
  n3d_FileTypeSet(&filename, "nd_trans_vs_e_pz", TRUE);
  f_tran_pz = fopen(filename, "w");
  str_free(filename);

  filename = n3d_strdup_n(d->inputfile);
  n3d_FileTypeSet(&filename, "nd_trans_vs_e_pangle", TRUE);
  f_tran_pa = fopen(filename, "w");
  str_free(filename);

  fprintf(f_tran, "# index energy initial final strength_x strength_y "
                  "strength_z strength_angle\n");
  fprintf(f_tran_px, "# symolic plot of light and dark spectral lines - light "
                     "polarixed in x\n");
  fprintf(f_tran_py, "# symolic plot of light and dark spectral lines - light "
                     "polarixed in y\n");
  fprintf(f_tran_pz, "# symolic plot of light and dark spectral lines - light "
                     "polarixed in z\n");
  fprintf(f_tran_pa, "# symolic plot of light and dark spectral lines - light "
                     "polarixed at angle\n");
  /* Print out transition energies with corresponding states.  Also
     print file containing delta functions at each transition */
  erange_ext = 0.1 * (d->te[nt - 1] - d->te[0]);
  fprintf(f_tran_px, "%e %e\n", d->te[0] - erange_ext, 1e-10);
  fprintf(f_tran_py, "%e %e\n", d->te[0] - erange_ext, 1e-10);
  fprintf(f_tran_pz, "%e %e\n", d->te[0] - erange_ext, 1e-10);
  fprintf(f_tran_pa, "%e %e\n", d->te[0] - erange_ext, 1e-10);
  for (i = 0; i < nt; i++) {
#define LOWVAL 1e-30
    fprintf(f_tran, "%i %e %i %i %e %e %e %e\n", i, d->te[i], d->tindx[i][0],
            d->tindx[i][1], d->opmat_pxyz[d->tindx[i][0]][d->tindx[i][1]][0],
            d->opmat_pxyz[d->tindx[i][0]][d->tindx[i][1]][1],
            d->opmat_pxyz[d->tindx[i][0]][d->tindx[i][1]][2],
            d->opmat_pxyz[d->tindx[i][0]][d->tindx[i][1]][3]);
    fprintf(f_tran_px, "%e %e\n", d->te[i] - 1e-12, LOWVAL);
    fprintf(f_tran_py, "%e %e\n", d->te[i] - 1e-12, LOWVAL);
    fprintf(f_tran_pz, "%e %e\n", d->te[i] - 1e-12, LOWVAL);
    fprintf(f_tran_pa, "%e %e\n", d->te[i] - 1e-12, LOWVAL);

    fprintf(f_tran_px, "%e %e\n", d->te[i],
            d->opmat_pxyz[d->tindx[i][0]][d->tindx[i][1]][0]);
    fprintf(f_tran_py, "%e %e\n", d->te[i],
            d->opmat_pxyz[d->tindx[i][0]][d->tindx[i][1]][1]);
    fprintf(f_tran_pz, "%e %e\n", d->te[i],
            d->opmat_pxyz[d->tindx[i][0]][d->tindx[i][1]][2]);
    fprintf(f_tran_pa, "%e %e\n", d->te[i],
            d->opmat_pxyz[d->tindx[i][0]][d->tindx[i][1]][3]);
    /*
    if ( d->opmat_pxyz[d->tindx[i][0]][d->tindx[i][1]][0] > 1e-20 )
        fprintf( f_tran_px,"%e %e\n", d->te[i], 1.0);
    else
        fprintf( f_tran_px,"%e %e\n", d->te[i], 0.1);

    if ( d->opmat_pxyz[d->tindx[i][0]][d->tindx[i][1]][1] > 1e-20 )
        fprintf( f_tran_py,"%e %e\n", d->te[i], 1.0);
    else
        fprintf( f_tran_py,"%e %e\n", d->te[i], 0.1);

    if ( d->opmat_pxyz[d->tindx[i][0]][d->tindx[i][1]][2] > 1e-20 )
        fprintf( f_tran_pz,"%e %e\n", d->te[i], 1.0);
    else
        fprintf( f_tran_pz,"%e %e\n", d->te[i], 0.1);

    if ( d->opmat_pxyz[d->tindx[i][0]][d->tindx[i][1]][3] > 1e-20 )
        fprintf( f_tran_pa,"%e %e\n", d->te[i], 1.0);
    else
        fprintf( f_tran_pa,"%e %e\n", d->te[i], 0.1);

*/

    fprintf(f_tran_px, "%e %e\n", d->te[i] + 1e-12, LOWVAL);
    fprintf(f_tran_py, "%e %e\n", d->te[i] + 1e-12, LOWVAL);
    fprintf(f_tran_pz, "%e %e\n", d->te[i] + 1e-12, LOWVAL);
    fprintf(f_tran_pa, "%e %e\n", d->te[i] + 1e-12, LOWVAL);
  }
  fprintf(f_tran_px, "%e %e\n", d->te[nt - 1] + erange_ext, LOWVAL);
  fprintf(f_tran_py, "%e %e\n", d->te[nt - 1] + erange_ext, LOWVAL);
  fprintf(f_tran_pz, "%e %e\n", d->te[nt - 1] + erange_ext, LOWVAL);
  fprintf(f_tran_pa, "%e %e\n", d->te[nt - 1] + erange_ext, LOWVAL);

  fclose(f_tran);
  fclose(f_tran_px);
  fclose(f_tran_py);
  fclose(f_tran_pz);
  fclose(f_tran_pa);

  return 0;
}

int calc_rate_cmplx(rvectr optmatvec, qd_struct d, cvectr init, cvectr finl,
                    real px, real py, real pz) {
  int i, j, k, l, m, n;
  real Tf, mxm, mxp, sumx, sumy, sumz, sumi, sumf;

  sumx = sumy = sumz = 0;

  Tf = -HBAR_MASSFACTOR / d->geo.lattice_x / d->geo.lattice_x;

  sumi = 0;
  sumf = 0;
  for (l = 0; l < d->geo.N_Cell; l++) {
    i = d->geo.cell__ijk[l][0];
    j = d->geo.cell__ijk[l][1];
    k = d->geo.cell__ijk[l][2];

    m = d->geo.ijk__cell[i + 1][j][k];
    n = d->geo.ijk__cell[i - 1][j][k];
    if (m >= 0) {
      mxp = (d->mstar[d->mat[i + 1][j][k]] + d->mstar[d->mat[i][j][k]]) / 2.0;
      sumx += Tf / mxp * nml_dcreal(cmul(finl[l], cconj(init[m])));
    }
    if (n >= 0) {
      mxm = (d->mstar[d->mat[i - 1][j][k]] + d->mstar[d->mat[i][j][k]]) / 2.0;
      sumx -= Tf / mxm * nml_dcreal(cmul(finl[l], cconj(init[n])));
    }

    m = d->geo.ijk__cell[i][j + 1][k];
    n = d->geo.ijk__cell[i][j - 1][k];
    if (m >= 0) {
      mxp = (d->mstar[d->mat[i][j + 1][k]] + d->mstar[d->mat[i][j][k]]) / 2.0;
      sumy += Tf / mxp * nml_dcreal(cmul(finl[l], cconj(init[m])));
    }
    if (n >= 0) {
      mxm = (d->mstar[d->mat[i][j - 1][k]] + d->mstar[d->mat[i][j][k]]) / 2.0;
      sumy -= Tf / mxm * nml_dcreal(cmul(finl[l], cconj(init[n])));
    }

    m = d->geo.ijk__cell[i][j][k + 1];
    n = d->geo.ijk__cell[i][j][k - 1];
    if (m >= 0) {
      mxp = (d->mstar[d->mat[i][j][k + 1]] + d->mstar[d->mat[i][j][k]]) / 2.0;
      sumz += Tf / mxp * nml_dcreal(cmul(finl[l], cconj(init[m])));
    }
    if (n >= 0) {
      mxm = (d->mstar[d->mat[i][j][k - 1]] + d->mstar[d->mat[i][j][k]]) / 2.0;
      sumz -= Tf / mxm * nml_dcreal(cmul(finl[l], cconj(init[n])));
    }
  }
  optmatvec[0] = sumx * sumx;
  optmatvec[1] = sumy * sumy;
  optmatvec[2] = sumz * sumz;
  optmatvec[3] =
      px * px * sumx * sumx + py * py * sumy * sumy + pz * pz * sumz * sumz;

  return 0;
}
