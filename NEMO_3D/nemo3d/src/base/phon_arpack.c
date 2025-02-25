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
$Header: /repo/nemo3d/src/base/phon_arpack.c,v 1.1 2004/02/27 04:09:05 lol Exp $
*****************************************************************************/

#ifdef NO_PARPACK

#include "phon_arpack.h"
#include "dmk_mult.h"
#include "phon_Keating.h"

#define EnergyScale                                                            \
  0.0065821222 //=hbar[6.5821222*10^-16 eV*s^2]*sqrt{VFFC[1 N/m]/M_atom[10^-26
               //kg]}

// #define PR_EIGENVALUES_CONVERGED
// #define PR_ARPACK_PARAMETERS

int phon_arpack(qd_struct d,
                int iq_ph // number of the point in the reciprocal space
) {
  int num_states =
      3 * d->n_atom_tot; // d->phon.N_cols[mpi_n3d_id];//total number of states
  int minn = find_min_ivec_mpi(d, d->phon.N_cols);
  int nev = (d->opt.ExecParam.PhonCalc.N_branches >= num_states)
                ? num_states - 1
                : d->opt.ExecParam.PhonCalc.N_branches;
  int ncv = ((nev * 3 >= num_states) ? num_states : 3 * nev);
  if (ncv >
      minn) { // Ensure the same size of the Krylov's subspace on all processors
    ncv = minn;
    if (nev > ncv - 2) {
      nev = ncv - 2;
    }
  }
  int ido = 0;        /* for reverse comm.; initialized to 0 */
  char bmat[] = "I";  /* standard eigval problem */
  int n = num_states; /* dimension of eigenproblem */
  int nloc =
      d->phon.N_cols[mpi_n3d_id]; // d->phon.seg_col_ln[mpi_n3d_id];//The number
                                  // of elements of the vector V that result in
                                  // nonzero product on the processor
  char *which;
  if (d->opt.ExecParam.PhonCalc.PS == PhonCalc_struct::Higher)
    which = "LR";
  else
    which = "SR";
  real tol = d->opt.ExecParam.PhonCalc.PhonParpack.Tolerance; //!!!!!!!!!!!!!!
  // stopping criterion: BOUNDS(I) <= TOL*ABS(RITZ(I))
  int ldv = nloc; // num_states;//nloc; /* leading dimension of v */
  int iparam[11];
  int ipntr[14];
  int lworkl = (3 * ncv + 5) * ncv; // size of work space
  cvectr resid = Cvectr(nloc);      // random if info==0; contains final resid
  cvectr vloc = Cvectr(ncv * nloc); // ncv column vectors of size nloc
  cvectr workd = Cvectr(3 * nloc);  // work array for reverse comm
  cvectr workl = Cvectr(lworkl);    // work array for reverse comm
  rvectr rwork = Rvectr(ncv);       /* Private (replicated) array on each PE
                                       or array allocated on the front end*/
  int info = 0;                     /* use random starting resid */
  int rvec = d->opt.ExecParam.PhonCalc.PhonOut.PhonDV;
  char c = 'A';
  int *select = new int[ncv];
  complex sigma;
  sigma.r = sigma.i = 0.0;
  cvectr d1 = Cvectr(ncv * 2);
  cvectr workev = Cvectr(3 * ncv);

  iparam[0] = 1; /* shift method */
  iparam[2] = d->opt.ExecParam.PhonCalc.PhonParpack
                  .MaxIter; /* max # arnoldi update iter */
  iparam[3] = 1;            /* block size (must be 1) */
  iparam[4] = 0;            /* on output: # converged ritz values */
  iparam[6] = 1;            /* type of eigenproblem to be solved */
  iparam[7] = 0;            /* used only if ido==3 */
  iparam[8] = 1;            /* total number of OP*x operations */
  iparam[9] = 1;            /* total number of B*x operations if BMAT='G' */
  iparam[10] = 1;           /* total number of steps of re-orthogonalization */

  d->phon
      .comm_vectors_allocate(); // Allocates memory for parts of the eigenvector
                                // that participate in communication

  int dummy = 0,
      NiterLooked =
          20; // NiterLooked determines the # of Arnoldy itteration print-outs
  double Tstart = 0.0;

  div_t result;
  real precision = ((d->opt.ExecParam.PhonCalc.PhonParpack.Tolerance > 0.0 &&
                     d->opt.ExecParam.PhonCalc.PhonParpack.Tolerance < 1e-10)
                        ? 1000 * d->opt.ExecParam.PhonCalc.PhonParpack.Tolerance
                        : 1e-10); // precision for check of
                                  // reasonability of eigenvalues
  int i;     // dummy
  bool flag; //=true for troubles
#ifdef PR_PARPACK_PARAMETERS
  printf("\nmpi_n3d_id=%d : pznaupd & pzneupd parameters:"
         "\n                comm=%d ido=%d bmat=%s which=%s tol=%g rvec=%d"
         "\n                nev=%d ncv=%d nloc=%d ldv=%d lworkl=%d\n",
         mpi_n3d_id, comm, ido, bmat, which, tol, rvec, nev, ncv, nloc, ldv,
         lworkl);
#endif // PR_PARPACK_PARAMETERS
  while (TRUE) {
    /* Repeatedly call the routine PZNAUPD and take actions
       indicated by parameter IDO until either convergence
       is indicated or maxitr has been exceeded.  */
    MPI_TIC(Tstart);
    znaupd(&ido, bmat, &nloc, which, &nev, &tol, resid, &ncv, vloc, &ldv,
           iparam, ipntr, workd, workl, &lworkl, rwork, &info);
    MPI_TOC(mpiTiming.parpack_function, Tstart);
    if (ido == -1 || ido == 1) {
      dummy++;
      result = div(dummy, NiterLooked);
      if (!result.rem)
        masterPrint("phon_arpack : %d-th time going to multiply...", dummy);
      av_mpi(d, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
    } else
      break;
  }

  if (info < 0) {
    printf("ERROR in znaupd of ARPACK_SOLVER %d  (# eig=%d) (n=%d)\n", info,
           iparam[4], n);
  } else {
    if (mpi_n3d_id == 0) {
#ifdef PR_ARPACK_PARAMETERS
      cout << endl;
      cout << " ===== Results of ZNAUPD in phon_arpack===== " << endl << endl;
      cout << " The number of Ritz values requested is " << nev << endl;
      cout << " The number of Arnoldi vectors generated  is " << ncv << endl;
      cout << " The portion of the spectrum: " << which << endl;
      cout << " The number of converged Ritz values is " << iparam[4] << endl;
      cout << " The number of Implicit Arnoldi update iterations taken is "
           << iparam[2] << endl;
      cout << " The number of OP*x is " << iparam[8] << endl;
      cout << " The number of reorthogonalization steps is " << iparam[10]
           << endl;
      cout << " The convergence criterion is " << tol << endl;
#endif // PR_ARPACK_PARAMETERS
#ifndef PR_ARPACK_PARAMETERS
      cout << "Converged " << iparam[4] << " from requested " << nev << " "
           << which << " eigenvalues" << endl;
#endif // PR_ARPACK_PARAMETERS
    }
    MPI_TIC(Tstart);
    zneupd(&rvec, &c, select, d1, vloc, &ldv, &sigma, workev, bmat, &nloc,
           which, &nev, &tol, resid, &ncv, vloc, &ldv, iparam, ipntr, workd,
           workl, &lworkl, rwork, &info);
    MPI_TOC(mpiTiming.parpack_function, Tstart);
  }

  if (info < 0) {
    printf("ERROR in pzneupd of ARPACK SOLVER %d \n", info);
  } else if (info >= 0 && iparam[4] > 0) {
    qsortEph(d1, 0, iparam[4] - 1);
    if (!mpi_n3d_id) {
      flag = false;
      for (i = 0; i < iparam[4]; i++) { // loop over eigenvalues with
                                        // Check if eigenvalues are reasonable
        if (d1[i].i > precision) {
          masterPrint("!!!!!!!Large imaginary part of %d-th eigenvalue!!!!!!!",
                      i);
          flag = true;
        }
        if (d1[i].r < 0) {
          if (d1[i].r < -precision) {
            masterPrint(
                "!!!!!!!Large negative real part of %d-th eigenvalue!!!!!!!",
                i);
            flag = true;
          } else
            d->phon.E_ph[iq_ph][i + 3] =
                0.0; // Substitute small negative eigenvalue with zero
        } else
          d->phon.E_ph[iq_ph][i + 3] =
              EnergyScale * sqrt(d1[i].r); // write the scaled to eV
                                           // real part of the eigenvalues
                                           // to the structure
      } // loop over eigenvalues
#ifdef PR_EIGENVALUES_CONVERGED
      flag = true;
#endif            // PR_EIGENVALUES_CONVERGED
      if (flag) { // if something is wrong with eigenvalues
        masterPrint("Parpack thinks that there are %d eigenvalues converged:",
                    iparam[4]);
        for (i = 0; i < iparam[4]; i++) { // print loop
          masterPrint("%d : %g+i%g", i, d1[i].r, d1[i].i);
        }
      } // print if something is wrong with eigenvalues

    } // on the master processor
  } // if there are converged eigenvalues

  /* Output of the eigenvectors, i.e. atomic displacement vector */
  /*
  if(info>=0 && iparam[4]>0) {
    d->wfc = Cmatrix(ncv, d->seg_ln[mpi_n3d_id]);
    for ( int ii=0; ii < iparam[4]; ii++ ) {
      int index_for_ii = orderedI[ii];
      int index_for_vloc = orderedI[ii]*d->seg_ln[mpi_n3d_id];
      for ( int jj= 0; jj <d->seg_ln[mpi_n3d_id]; jj++ ){
           d->wfc[ii][jj].r = vloc[index_for_vloc+jj].r;
           d->wfc[ii][jj].i = vloc[index_for_vloc+jj].i;
      }
    }
    if(d->opt.ExecParam.Output.PsiSqrGnuPlot) writePsiSqr(d->wfc, d->neigv, d,
  0); if(d->opt.ExecParam.Output.FullEigvec) for(int ie=0; ie<d->neigv; ie++)
  writeEigenvector(d, d->wfc[ie], ie);

    compute_ritz_value(d);
    compute_residual(d);
    compute_eigenvector_norm(d);
  }
  */

  rm_cvectr(&d1);
  rm_cvectr(&workev);
  rm_cvectr(&resid);
  rm_cvectr(&vloc);
  rm_cvectr(&workd);
  rm_cvectr(&workl);
  rm_rvectr(&rwork);
  d->phon.comm_vectors_deallocate(); // Allocates memory for parts of the
                                     // eigenvector that participate in
                                     // communication

  return 0;
}

#endif // NO_PARPACK
