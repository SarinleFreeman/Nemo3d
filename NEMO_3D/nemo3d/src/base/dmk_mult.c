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
  Olga L. Lazarenkova (LOL@jpl.nasa.gov)

Written by:  Chris Bowen
             Gerhard Klimeck
             Fabiano Oyafuso
             Seungwon Lee
             Olga Lazarenkova
             Hook Hua

This product includes software developed by the Apache Software Foundation
(http://www.apache.org/).

*****************************************************************************
$Header: /repo/nemo3d/src/base/dmk_mult.c,v 1.4 2007/08/10 15:03:32 baeh Exp $
*****************************************************************************/

#include "dmk_mult.h"

void av(qd_struct d, int row_init, int row_final, cvectr v, cvectr w)
// Computes w<==DM*v
{
  int j, n, n0, jmax;

  for (n = row_init, n0 = 0; n <= row_final; n++, n0++) {
    jmax = d->phon.Ind_DM[n][0];
    w[n0].r = 0.0;
    w[n0].i = 0.0;
    for (j = 1; j <= jmax; j++) {
      w[n0].r += (d->phon.DM[n][j - 1].r) * (v[d->phon.Ind_DM[n][j]].r) -
                 (d->phon.DM[n][j - 1].i) * (v[d->phon.Ind_DM[n][j]].i);
      w[n0].i += (d->phon.DM[n][j - 1].r) * (v[d->phon.Ind_DM[n][j]].i) +
                 (d->phon.DM[n][j - 1].i) * (v[d->phon.Ind_DM[n][j]].r);
    }
  }
}

void av_mpi(qd_struct d, cvectr v, cvectr w)
// Calls av to compute (1) residential w<==DM*v,
//                     (2) wR, wL to be sent to right/left neighbors
// Communicates with neighbors (1) send wR, wL
//                             (2) receive vR, vL
// Computes the final value of w
{
  int i; // dummy
  // Residential part
  av(d, d->phon.seg_col_s[mpi_n3d_id] - d->phon.min_row[mpi_n3d_id],
     (d->phon.seg_col_s[mpi_n3d_id] + d->phon.seg_col_ln[mpi_n3d_id]) -
         d->phon.min_row[mpi_n3d_id] - 1,
     v, w);
  //////////////////////////////

  if (d->phon.szCommL[mpi_n3d_id]) { // Block to be sent to the left processor
    if (mpi_n3d_numprocs != 1 && !mpi_n3d_id && d->geo.StrainIsPeriodicX()) {
      av(d, d->phon.N_rows[mpi_n3d_id] - d->phon.szCommL[mpi_n3d_id],
         d->phon.N_rows[mpi_n3d_id] - 1, v, d->phon.wL);
    } else {
      av(d,
         d->phon.seg_col_s[mpi_n3d_id] - d->phon.szCommL[mpi_n3d_id] -
             d->phon.min_row[mpi_n3d_id],
         d->phon.seg_col_s[mpi_n3d_id] - d->phon.min_row[mpi_n3d_id] - 1, v,
         d->phon.wL);
    }
  }
  if (d->phon.szCommR[mpi_n3d_id]) { // Block to be sent to the right processor
    if (mpi_n3d_numprocs != 1 && mpi_n3d_id == mpi_n3d_numprocs - 1 &&
        d->geo.StrainIsPeriodicX()) {
      av(d, d->phon.N_rows[mpi_n3d_id] - d->phon.szCommR[mpi_n3d_id],
         d->phon.N_rows[mpi_n3d_id] - 1, v, d->phon.wR);
    } else {
      av(d,
         d->phon.seg_col_s[mpi_n3d_id] + d->phon.seg_col_ln[mpi_n3d_id] -
             d->phon.min_row[mpi_n3d_id],
         d->phon.seg_col_s[mpi_n3d_id] + d->phon.seg_col_ln[mpi_n3d_id] +
             d->phon.szCommR[mpi_n3d_id] - d->phon.min_row[mpi_n3d_id] - 1,
         v, d->phon.wR);
    }
  }

  if (mpi_n3d_numprocs > 1) {
    //////////////////////////////////
    // Communication with neighbors //
    //////////////////////////////////
    if (mpi_n3d_id % 2) {
      if (d->phon.mpi_left[mpi_n3d_id] >= 0) {
        MPI_Sendrecv(d->phon.wL, 2 * d->phon.szCommL[mpi_n3d_id], MPI_DOUBLE,
                     d->phon.mpi_left[mpi_n3d_id], 1, d->phon.vL,
                     2 * d->phon.szCommR[d->phon.mpi_left[mpi_n3d_id]],
                     MPI_DOUBLE, d->phon.mpi_left[mpi_n3d_id], 0,
                     MPI_COMM_WORLD, &mpi_n3d_status);
      }
      if (d->phon.mpi_right[mpi_n3d_id] >= 0) {
        MPI_Sendrecv(d->phon.wR, 2 * d->phon.szCommR[mpi_n3d_id], MPI_DOUBLE,
                     d->phon.mpi_right[mpi_n3d_id], 0, d->phon.vR,
                     2 * d->phon.szCommL[d->phon.mpi_right[mpi_n3d_id]],
                     MPI_DOUBLE, d->phon.mpi_right[mpi_n3d_id], 1,
                     MPI_COMM_WORLD, &mpi_n3d_status);
      }
    } else {
      if (d->phon.mpi_right[mpi_n3d_id] >= 0) {
        MPI_Sendrecv(d->phon.wR, 2 * d->phon.szCommR[mpi_n3d_id], MPI_DOUBLE,
                     d->phon.mpi_right[mpi_n3d_id], 0, d->phon.vR,
                     2 * d->phon.szCommL[d->phon.mpi_right[mpi_n3d_id]],
                     MPI_DOUBLE, d->phon.mpi_right[mpi_n3d_id], 1,
                     MPI_COMM_WORLD, &mpi_n3d_status);
      }
      if (d->phon.mpi_left[mpi_n3d_id] >= 0) {
        MPI_Sendrecv(d->phon.wL, 2 * d->phon.szCommL[mpi_n3d_id], MPI_DOUBLE,
                     d->phon.mpi_left[mpi_n3d_id], 1, d->phon.vL,
                     2 * d->phon.szCommR[d->phon.mpi_left[mpi_n3d_id]],
                     MPI_DOUBLE, d->phon.mpi_left[mpi_n3d_id], 0,
                     MPI_COMM_WORLD, &mpi_n3d_status);
      }
    }
  } // mpi communication
  /////////////////////////////////////////////////////////
  // Addition of the blocks received from the neighbours //
  /////////////////////////////////////////////////////////
  if (d->phon.szCommR[d->phon.mpi_left[mpi_n3d_id]]) { // Addition from the left
                                                       // processor vL
    for (i = 0; i < d->phon.szCommR[d->phon.mpi_left[mpi_n3d_id]]; i++) {
      w[i].r += d->phon.vL[i].r;
      w[i].i += d->phon.vL[i].i;
    }
  }
  if (d->phon.szCommL[d->phon.mpi_right[mpi_n3d_id]]) { // Addition from the
                                                        // right processor vR
    for (i = 0; i < d->phon.szCommL[d->phon.mpi_right[mpi_n3d_id]]; i++) {
      w[d->phon.seg_col_ln[mpi_n3d_id] -
        d->phon.szCommL[d->phon.mpi_right[mpi_n3d_id]] + i]
          .r += d->phon.vR[i].r;
      w[d->phon.seg_col_ln[mpi_n3d_id] -
        d->phon.szCommL[d->phon.mpi_right[mpi_n3d_id]] + i]
          .i += d->phon.vR[i].i;
    }
  }
} // end of void av_mpi

#define swapReal(g, h)                                                         \
  {                                                                            \
    real y = (g);                                                              \
    (g) = (h);                                                                 \
    (h) = y;                                                                   \
  }

// static void qsortEph(cvectr Eph,//[nconv==iparam[4]] vector of converged
// eigenvalues
void qsortEph(cvectr Eph, //[nconv==iparam[4]] vector of converged eigenvalues
              int beg, int end // lower and upper bounds of the arrays
              )
/*Modification of quick_sort to
 sort obtained phonon eigenvalues
 in ascending order of its real part.
*/
{
  int p, i;

  if (beg >= end)
    return;

  swapReal(Eph[beg].r, Eph[(beg + end) / 2].r);
  swapReal(Eph[beg].i, Eph[(beg + end) / 2].i);

  p = beg;
  for (i = beg + 1; i <= end; i++)
    if (Eph[i].r < Eph[beg].r) {
      ++p;
      swapReal(Eph[i].r, Eph[p].r);
      swapReal(Eph[i].i, Eph[p].i);
    }

  swapReal(Eph[beg].r, Eph[p].r);
  swapReal(Eph[beg].i, Eph[p].i);

  qsortEph(Eph, beg, p - 1);
  qsortEph(Eph, p + 1, end);
}

int find_min_ivec_mpi(qd_struct d, ivectr iv) {
  int mini = iv[0];
  for (int i = 1; i < mpi_n3d_numprocs; i++) {
    if (iv[i] < mini)
      mini = iv[i];
  }
  return (mini);
}
