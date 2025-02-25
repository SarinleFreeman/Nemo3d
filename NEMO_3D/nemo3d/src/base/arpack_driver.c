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
$Header: /repo/nemo3d/src/base/arpack_driver.c,v 1.16 2006/10/17 14:10:14
hoonryu Exp $
*****************************************************************************/

#ifdef NO_PARPACK

#include "arpack_driver.h"

int arpack_driver(qd_struct d) {

  // verify_ham(d);
  if (d->rc == NULL)
    d->rc = Cvectr(d->nvmap[mpi_n3d_id] * d->NBasisStates);

  int nev = d->opt.ExecParam.ElCalc.ElAlg.Parpack.NumEigVal;
  int ncv = nev * d->opt.ExecParam.ElCalc.ElAlg.Parpack.NumArnoldiVector;
  int rvec = d->opt.ExecParam.ElCalc.ElAlg.Parpack.CalEigVec;
  real refE = d->opt.ExecParam.ElCalc.ElAlg.Parpack.RefEigVal;
  real tol = d->opt.ExecParam.ElCalc.ElAlg.Parpack.Tolerance;

  int ido = 0;          /* for reverse comm.; initialized to 0 */
  char bmat[] = "I";    /* standard eigval problem */
  int n = d->n_ham_tot; /* dimension of eigenproblem */
  int nloc = d->seg_ln[mpi_n3d_id];
  char which[] = "SM"; /* smallest eigenvalues */
  int ldv = nloc;      /* leading dimension of v */
  int iparam[11];
  int ipntr[14];
  int lworkl = (3 * ncv + 5) * ncv; /* size of work space */
  cvectr resid = Cvectr(nloc); /* random if info==0; contains final resid */
  cvectr vloc = Cvectr(ncv * nloc); /* ncv column vectors of size nloc */
  cvectr workd = Cvectr(3 * nloc);  /* work array for reverse comm */
  cvectr workl = Cvectr(lworkl);    /* work array for reverse comm */
  rvectr rwork = Rvectr(ncv);       /* Private (replicated) array on each PE
                                       or array allocated on the front end*/
  int info = 0;                     /* use random starting resid */
  char c = 'A';
  int *select = new int[ncv];
  complex sigma;
  sigma.r = sigma.i = 0.0;
  cvectr d1 = Cvectr(ncv * 2);
  cvectr workev = Cvectr(3 * ncv);

  iparam[0] = 1; /* shift method */
  iparam[2] = d->opt.ExecParam.ElCalc.ElAlg.Parpack
                  .MaxIter; /* max # arnoldi update iter */
  iparam[3] = 1;            /* block size (must be 1) */
  iparam[4] = 0;            /* on output: # converged ritz values */
  iparam[6] = 1;            /* type of eigenproblem to be solved */
  iparam[7] = 0;            /* used only if ido==3 */
  iparam[8] = 1;            /* total number of OP*x operations */
  iparam[9] = 1;            /* total number of B*x operations if BMAT='G' */
  iparam[10] = 1;           /* total number of steps of re-orthogonalization */

  double Tstart = 0.0;
  do {
    MPI_TIC(Tstart);
    znaupd(&ido, bmat, &nloc, which, &nev, &tol, resid, &ncv, vloc, &ldv,
           iparam, ipntr, workd, workl, &lworkl, rwork, &info);
    MPI_TOC(mpiTiming.parpack_function, Tstart);
    Hmult_spds_complete(&workd[ipntr[1] - 1], d, refE, &workd[ipntr[0] - 1]);
  } while (ido == -1 || ido == 1);

  if (info < 0) {
    printf("ERROR in znaupd of ARPACK_SOLVER %d  (# eig=%d) (n=%d)\n", info,
           iparam[4], n);
  } else {
    if (mpi_n3d_id == 0) {
      cout << endl;
      cout << " ====== Results of ZNAUPD ===== " << endl << endl;
      cout << " The number of Ritz values requested is " << nev << endl;
      cout << " The number of Arnoldi vectors generated  is " << ncv << endl;
      cout << " What portion of the spectrum: " << which << endl;
      cout << " The number of converged Ritz values is " << iparam[4] << endl;
      cout << " The number of Implicit Arnoldi update iterations taken is "
           << iparam[2] << endl;
      cout << " The number of OP*x is " << iparam[8] << endl;
      cout << " The number of reorthogonalization steps is " << iparam[10]
           << endl;
      cout << " The convergence criterion is " << tol << endl;
    }
    MPI_TIC(Tstart);
    zneupd(&rvec, &c, select, d1, vloc, &ldv, &sigma, workev, bmat, &nloc,
           which, &nev, &tol, resid, &ncv, vloc, &ldv, iparam, ipntr, workd,
           workl, &lworkl, rwork, &info);
    MPI_TOC(mpiTiming.parpack_function, Tstart);
  }

  if (info < 0) {
    printf("ERROR in zneupd of ARPACK SOLVER %d \n", info);
  }

  else if (iparam[4] > 0 && !mpi_n3d_id) {
    cout << endl << "Found eigval." << endl;
    for (int i = 0; i < iparam[4]; i++)
      cout << i << " " << d1[i].r + refE << " " << d1[i].i << endl;
    cout << endl;
  }

  rvectr orderedE;
  ivectr orderedI;

  if (info >= 0 && iparam[4] > 0) {
    d->neigv = iparam[4];
    d->WF = vector<waveFunction>(d->neigv);
    orderedE = Rvectr(d->neigv);
    orderedI = Ivectr(d->neigv);
    for (int ni = 0; ni < d->neigv; ni++) {
      orderedE[ni] = d1[ni].r + refE;
      orderedI[ni] = ni;
    }
    for (int a = 0; a < d->neigv - 1; a++)
      for (int b = a + 1; b < d->neigv; b++)
        if (orderedE[a] > orderedE[b]) {
          double tmp1 = orderedE[a];
          orderedE[a] = orderedE[b];
          orderedE[b] = tmp1;
          int tmp2 = orderedI[a];
          orderedI[a] = orderedI[b];
          orderedI[b] = tmp2;
        }
    for (int i = 0; i < d->neigv; i++) {
      waveFunction &wf = d->WF[i];
      wf.indx = i;
      wf.E = orderedE[i];
      wf.valid = 1;
    }
    writeEigenvalues(d, 'a');
  }

  if (info >= 0 && iparam[4] > 0 &&
      d->opt.ExecParam.ElCalc.ElAlg.Parpack.CalEigVec) {
    d->wfc = Cmatrix(ncv, d->seg_ln[mpi_n3d_id]);
    for (int ii = 0; ii < iparam[4]; ii++) {
      // int index_for_ii = orderedI[ii];
      int index_for_vloc = orderedI[ii] * d->seg_ln[mpi_n3d_id];
      for (int jj = 0; jj < d->seg_ln[mpi_n3d_id]; jj++) {
        d->wfc[ii][jj].r = vloc[index_for_vloc + jj].r;
        d->wfc[ii][jj].i = vloc[index_for_vloc + jj].i;
      }
    }

    if (d->opt.Dev.band_model != Dev_struct::Bands_1_s_nospin) {
      writePsiSqr(d->wfc, d->neigv, d, 0);
      for (int ie = 0; ie < d->neigv; ie++)
        writeEigenvector(d, d->wfc[ie], ie);
    } else {
      // in the single band case the data file is small enough
      //  enough that we write out the full wavefunction
      //  print_psisq_gnu_pdb(d->wfc,  d->neigv,d,0);
    }

    compute_ritz_residual_norm_deviation(d);
#if 0
     compute_ritz_value(d);   
     compute_residual(d);
     compute_eigenvector_norm(d);
#endif
  }

  // write_cubic_struct_wavefunctions(  d );

  delete[] select;
  rm_cvectr(&d1);
  rm_cvectr(&workev);
  rm_cvectr(&resid);
  rm_cvectr(&vloc);
  rm_cvectr(&workd);
  rm_cvectr(&workl);
  rm_rvectr(&rwork);
  rm_rvectr(&orderedE);
  rm_ivectr(&orderedI);

  return d->neigv;
}

#endif // NO_PARPACK
