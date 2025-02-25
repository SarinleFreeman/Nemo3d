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
$Header: /repo/nemo3d/src/base/eigsys3d_par.c,v 1.78 2008/09/24 02:09:07 hoonryu
Exp $
*****************************************************************************/

#include "eigsys3d_par.h"
#include "CubicProj.h"

#undef MOMENTUM_MATRIX
#undef COMPUTE_Px
#undef COMPUTE_Py
#undef COMPUTE_Pz
#undef COMPUTE_TRANSP
#undef DUMP_MOMENTUM
#undef DUMP_TRANSP
#undef DUMP_OVLAP
#undef UNFLD_PROJ

// #define MOMENTUM_MATRIX
// #define COMPUTE_Px
// #define COMPUTE_Py
// #define COMPUTE_Pz
// #define COMPUTE_TRANSP
// #define DUMP_MOMENTUM
// #define DUMP_TRANSP
// #define DUMP_OVLAP
// #define UNFLD_PROJ

static void check_residual(qd_struct d) {
  int seglen = d->seg_ln[mpi_n3d_id];
  int **argList = (int **)nml_malloc(sizeof(int *) * 2);

  cout_master << "Computing residual of eigenvector.\n";

  cvectr Hpsi = Cvectr(seglen);
  cvectr psi = Cvectr(seglen);
  for (int i = 0; i < seglen; i++) {
    Hpsi[i].r = Hpsi[i].i = 0.0;
    psi[i].r = d->wfc[0][i].r;
    psi[i].i = d->wfc[0][i].i;
  }

  // double shift = d->eigv[0].r;
  double shift = 0.0;
  argList[1] = (int *)&shift;
  matvect_mult((const int **)argList, Hpsi, psi);

  // compute Hpsi^2
  double sum = 0.0;
  for (int i = 0; i < seglen; i++) {
    sum += Hpsi[i].r * Hpsi[i].r + Hpsi[i].i * Hpsi[i].i;
  }

  char filename[100];
  sprintf(filename, "jnk_%d", mpi_n3d_id);
  FILE *fp = fopen(filename, "w");
  for (int i = 0; i < seglen; i++) {
    fprintf(fp, "%e %e   %e %e\n", psi[i].r, psi[i].i, Hpsi[i].r, Hpsi[i].i);
  }
  fclose(fp);

  double tot;
  MPI_Allreduce(&sum, &tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  // compute Hpsi^2
  double sum2 = 0.0;
  for (int i = 0; i < seglen; i++) {
    sum2 += psi[i].r * psi[i].r + psi[i].i * psi[i].i;
  }
  double tot2;
  MPI_Allreduce(&sum2, &tot2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  cout_master << "eigval=" << d->eigv[0].r << "  Hpsi2=" << tot
              << "  psi2=" << tot2 << "\n";

  nml_free(argList);
}

int lanczos_initialize(qd_struct d) {
  MPI_TIME_INIT(Tstart);

  real sum_local = 0.0;
  real sum = 0.0;
  cvectr r0 = d->r0;
  cvectr r0_start = d->r0_start;
  char filename[40];
  int j;

  switch (d->opt.ExecParam.ElCalc.ElAlg.Lanczos.Start) {
  case Lanczos_struct::Random:
    for (j = 0; j < d->seg_ln[mpi_n3d_id]; j++) {
      r0[j].r = 1.0 * rand() / RAND_MAX;
      r0[j].i = 0.0;
      sum_local += r0[j].r * r0[j].r + r0[j].i * r0[j].i;
    }
    print_com3("CPU_id=%d call MPI_Allreduce wavefunction norm, sum_local=%g\n",
               mpi_n3d_id, sum_local);
    break;
  case Lanczos_struct::CyclicRandom: {
    int sz = d->opt.ExecParam.ElCalc.ElAlg.Lanczos.CyclicString.size();
    for (j = 0; j < d->seg_ln[mpi_n3d_id]; j++) {
      int indx = j + d->seg_s[mpi_n3d_id];
      r0[j].r =
          d->opt.ExecParam.ElCalc.ElAlg.Lanczos.CyclicString[(2 * indx) % sz];
      r0[j].i = d->opt.ExecParam.ElCalc.ElAlg.Lanczos
                    .CyclicString[(2 * indx + 1) % sz];
      sum_local += r0[j].r * r0[j].r + r0[j].i * r0[j].i;
    }
  } break;
  case Lanczos_struct::Const1:
    for (j = 0; j < d->seg_ln[mpi_n3d_id]; j++) {
      r0[j].r = 1.0;
      r0[j].i = 1.0;
      sum_local += r0[j].r * r0[j].r + r0[j].i * r0[j].i;
    }
    print_com3("CPU_id=%d call MPI_Allreduce wavefunction norm, sum_local=%g\n",
               mpi_n3d_id, sum_local);
    break;

  case Lanczos_struct::ReadFile:
    if (d->opt.ExecParam.ElCalc.ElAlg.Lanczos.StartFileRead == "NULL") {
      die("ERROR in lanc_r0_setup_par\n"
          "The input request that the Lanczos vector is read in from a file.\n"
          "However, no file name was specified!\n");
    }
    {
      int i, proc;
      cvectr r0_full = Cvectr(d->n_ham_tot);

      if (mpi_n3d_id == mpi_n3d_masterid) {
        sprintf(filename, "%s",
                d->opt.ExecParam.ElCalc.ElAlg.Lanczos.StartFileRead.c_str());
        initStartVect(filename, d, r0_full);
        for (proc = 1; proc < mpi_n3d_numprocs; proc++) {
          /* might use non-blocking send? */
          MPI_Send(&r0_full[d->seg_s[proc]], 2 * d->seg_ln[proc], MPI_DOUBLE,
                   proc, mpi_n3d_masterid, MPI_COMM_WORLD);
        }
        for (i = 0; i < d->seg_ln[mpi_n3d_id]; i++) {
          r0[i] = r0_full[i];
        }
      } else {
        MPI_Recv(r0, 2 * d->seg_ln[mpi_n3d_id], MPI_DOUBLE, mpi_n3d_masterid,
                 MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_n3d_status);
      }
      MPI_Barrier(MPI_COMM_WORLD);

      rm_cvectr(&r0_full);

      for (j = 0; j < d->seg_ln[mpi_n3d_id]; j++) {
        sum_local += r0[j].r * r0[j].r + r0[j].i * r0[j].i;
      }
    }
    break;
  default:
    die("ERROR in lanc_r0_setup_par\n"
        "Non-implemented lanczos startvector initialization\n");
    break;
  }

  Allreduce_MPI_sp(&sum_local, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  print_com3("CPU_id=%d done MPI_Allreduce wavefunction norm, sum = %g\n",
             mpi_n3d_id, sum);

  sum = sqrt(sum);

  if (r0_start) {
    for (j = 0; j < d->seg_ln[mpi_n3d_id]; j++) {
      r0[j].r /= sum;
      r0[j].i /= sum;
      r0_start[j] = r0[j];
    }
  } else {
    for (j = 0; j < d->seg_ln[mpi_n3d_id]; j++) {
      r0[j].r /= sum;
      r0[j].i /= sum;
    }
  }

  return 0;
}

int lanczos_driver(qd_struct d) {
  char *fil_trace = NULL, *fil_tridiag = NULL, *fil_eigval = NULL;

  int ConvCheckStartIter =
      d->opt.ExecParam.ElCalc.ElAlg.Lanczos.ConvCheckStartIter;
  int ConvCheckSkipRate =
      d->opt.ExecParam.ElCalc.ElAlg.Lanczos.ConvCheckSkipRate;
  int neigv1 = d->opt.ExecParam.ElCalc.ElAlg.Lanczos.NumEigReq_vb;
  int neigv2 = d->opt.ExecParam.ElCalc.ElAlg.Lanczos.NumEigReq_cb;
  int max_iter;    // = d->opt.ExecParam.ElCalc.ElAlg.Lanczos.MaxIter;
  double conv_tol; // = d->opt.ExecParam.ElCalc.ElAlg.Lanczos.ConvTolerance;

  double emin1 = 0.0; // d->opt.ExecParam.ElCalc.ElAlg.Lanczos.Erange_vb[0];
  double emax1 = 0.0; // d->opt.ExecParam.ElCalc.ElAlg.Lanczos.Erange_vb[1];
  double emin2 = 0.0; // d->opt.ExecParam.ElCalc.ElAlg.Lanczos.Erange_cb[0];
  double emax2 = 0.0; // d->opt.ExecParam.ElCalc.ElAlg.Lanczos.Erange_cb[1];

  int seglen = d->seg_ln[mpi_n3d_id];
  int **argList = (int **)nml_malloc(sizeof(int *) * 3);

  // hard-code this for now
  const int verbose = 1;
  double shift = 0.0;
  int Nproj = 0;

  // Block Lanczos parameters
  int p = d->opt.ExecParam.ElCalc.ElAlg.BLanczos.BL_size; // 6;
  double deflation_tolerance =
      d->opt.ExecParam.ElCalc.ElAlg.BLanczos.DefTolerance;          // 1e-7;
  int eig_dump = d->opt.ExecParam.ElCalc.ElAlg.BLanczos.Eig_dump;   // 0;
  int num_vb = d->opt.ExecParam.ElCalc.ElAlg.BLanczos.num_vb;       // 0;
  int num_cb = d->opt.ExecParam.ElCalc.ElAlg.BLanczos.num_cb;       // 15;
  int req_vb = d->opt.ExecParam.ElCalc.ElAlg.BLanczos.NumEigReq_vb; // 0;
  int req_cb = d->opt.ExecParam.ElCalc.ElAlg.BLanczos.NumEigReq_cb; // 15;
  int max_no_success =
      d->opt.ExecParam.ElCalc.ElAlg.BLanczos.Max_no_success; // 2000;
  int save_BT_band = d->opt.ExecParam.ElCalc.ElAlg.BLanczos
                         .Dump_blockTridiagonal_Sparse_Matrix_Format; // 0;
  int save_BT_SMF = d->opt.ExecParam.ElCalc.ElAlg.BLanczos
                        .Dump_blockTridiagonal_Band_Matrix_Format; // 0;
  int read_BT_band = d->opt.ExecParam.ElCalc.ElAlg.BLanczos
                         .Read_blockTridiagonal_Band_Matrix_Format;      // 0;
  int deflation_on = d->opt.ExecParam.ElCalc.ElAlg.BLanczos.DeflationOn; // 1;
  int conv_method = d->opt.ExecParam.ElCalc.ElAlg.BLanczos.conv_method;  // 1;
  int conv_vec_method =
      d->opt.ExecParam.ElCalc.ElAlg.BLanczos.conv_vec_method; // 2;

  int do_lancz_eigvect =
      (d->opt.ExecParam.ElCalc.ElOut.PsiSqr ||
       d->opt.ExecParam.ElCalc.ElOut.PsiSqr_cell ||
       d->opt.ExecParam.ElCalc.ElOut.PsiSqrPdb ||
       d->opt.ExecParam.ElCalc.ElOut.Psi4_IPR ||
       d->opt.ExecParam.ElCalc.ElOut.Eigvect) &&
      ((d->opt.ExecParam.ElCalc.ElAlg.ResFind == ElAlg_struct::Lancz) ||
       (d->opt.ExecParam.ElCalc.ElAlg.ResFind == ElAlg_struct::BlockLanczos));

  if (d->opt.ExecParam.ElCalc.ElAlg.ResFind == ElAlg_struct::Lancz) {
    max_iter = d->opt.ExecParam.ElCalc.ElAlg.Lanczos.MaxIter;
    conv_tol = d->opt.ExecParam.ElCalc.ElAlg.Lanczos.ConvTolerance;
    emin1 = d->opt.ExecParam.ElCalc.ElAlg.Lanczos.Erange_vb[0];
    emax1 = d->opt.ExecParam.ElCalc.ElAlg.Lanczos.Erange_vb[1];
    emin2 = d->opt.ExecParam.ElCalc.ElAlg.Lanczos.Erange_cb[0];
    emax2 = d->opt.ExecParam.ElCalc.ElAlg.Lanczos.Erange_cb[1];
    argList[0] = (int *)&d;
    argList[1] = (int *)&shift;
    argList[2] = (int *)&Nproj;
  } else if (d->opt.ExecParam.ElCalc.ElAlg.ResFind ==
             ElAlg_struct::BlockLanczos) {
    do_lancz_eigvect = 1;
    max_iter = d->opt.ExecParam.ElCalc.ElAlg.BLanczos.MaxIter;
    conv_tol = d->opt.ExecParam.ElCalc.ElAlg.BLanczos.ConvTolerance;
    emin1 = d->opt.ExecParam.ElCalc.ElAlg.BLanczos.Erange_vb[0];
    emax1 = d->opt.ExecParam.ElCalc.ElAlg.BLanczos.Erange_vb[1];
    emin2 = d->opt.ExecParam.ElCalc.ElAlg.BLanczos.Erange_cb[0];
    emax2 = d->opt.ExecParam.ElCalc.ElAlg.BLanczos.Erange_cb[1];

    argList[0] = (int *)&d;
    argList[1] = (int *)&shift;

#ifdef BLANCZ_MATMULT
    argList[2] = (int *)&p;
    argList[3] = (int *)&Nproj;
#else
    argList[2] = (int *)&Nproj;
#endif
  }

  nml_memory_report(
      "beginning of lanczos_driver \n************\n**************** from here "
      "on all memory should be distributed\n");

  rvectr ev = Rvectr(max_iter);
  nml_memory_report("allocated ev in lanczos_driver\n");

  if (d->opt.ExecParam.ElCalc.ElOut.TraceEigval) {
    fil_trace = n3d_strdup_n(d->inputfile);
    n3d_FileTypeSet(&fil_trace, "nd_eval_trace", TRUE);
  }

  if (d->opt.ExecParam.ElCalc.ElOut.SaveTridiagInfo) {
    fil_tridiag = n3d_strdup_n(d->inputfile);
    n3d_FileTypeSet(&fil_tridiag, "nd_tridiag_bin", TRUE);

    fil_eigval = n3d_strdup_n(d->inputfile);
    n3d_FileTypeSet(&fil_eigval, "nd_eval_bin", TRUE);
  }

  print_com3("CPU=%d d->seg_ln[mpi_n3d_id]=%d\n", mpi_n3d_id,
             d->seg_ln[mpi_n3d_id]);
  print_com3("CPU=%d d->nvmap[mpi_n3d_id]=%d\n", mpi_n3d_id,
             d->nvmap[mpi_n3d_id]);
  print_com3("CPU=%d maxcom=%d\n", mpi_n3d_id, maxcom);

  d->r0 = Cvectr(d->seg_ln[mpi_n3d_id]);
  d->rc =
      Cvectr(d->nvmap[mpi_n3d_id] * d->NBasisStates); /* what does this do? */
  d->q0 = Cvectr(d->seg_ln[mpi_n3d_id]);
  d->q_m10 = Cvectr(d->seg_ln[mpi_n3d_id]);
  d->r0_start = (do_lancz_eigvect ? Cvectr(d->seg_ln[mpi_n3d_id]) : NULL);
  nml_memory_report("allocated auxiliary vectors in lanczos_driver\n");

  lanczos_initialize(d);
  nml_memory_report("after lanczos_initialize in lanczos_driver\n");

  cout_master << "-----------------------------------------------\n";
  cout_master << "  beginning Lanczos iteration\n";
  cout_master << "-----------------------------------------------\n";

#undef FILEIO_BARRIER_TEST
#ifdef FILEIO_BARRIER_TEST
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  double Tstart;
  MPI_TIC(Tstart);
#ifdef MPIP
  MPI_Pcontrol(2); // discard previous profile data
  MPI_Pcontrol(1); // turn on profiling
#endif

#ifdef MPI_TIMING
  mpiTiming.section = mpiTiming.electr_mpi;
#endif

  if (d->opt.ExecParam.ElCalc.ElAlg.ResFind == ElAlg_struct::Lancz) {
    d->neigv = eigensolve_lanczos(
        &ev, ConvCheckStartIter, ConvCheckSkipRate, max_iter, neigv1, neigv2,
        seglen, (int **)argList, emin1, emax1, emin2, emax2, conv_tol, verbose,
        fil_trace, fil_tridiag, fil_eigval, (cvectr)(d->r0),
        (cvectr)(d->r0_start), (cvectr)(d->q0), (cvectr)(d->q_m10),
        (cvectr *)&(d->eigv), // beware mem leak
        (cmatrix *)&(d->wfc), // beware mem leak
        &matvect_mult, d->opt.ExecParam.ElCalc.ElAlg.Lanczos.ReadInTridiag);
  } else if (d->opt.ExecParam.ElCalc.ElAlg.ResFind ==
             ElAlg_struct::BlockLanczos) {
    d->neigv = eigensolve_blk_lanczos(
        &ev, ConvCheckStartIter, ConvCheckSkipRate, max_iter, req_vb, req_cb,
        seglen, (int **)argList, emin1, emax1, emin2, emax2, conv_tol, verbose,
        fil_trace, fil_tridiag, fil_eigval, (cvectr)(d->r0),
        (cvectr)(d->r0_start), (cvectr)(d->q0), (cvectr)(d->q_m10),
        (cvectr *)&(d->eigv), // beware mem leak
        (cmatrix *)&(d->wfc), // beware mem leak

#ifndef BLANCZ_MATMULT
        &matvect_mult,
#else
        &matmat_mult,
#endif

        d->opt.ExecParam.ElCalc.ElAlg.Lanczos.ReadInTridiag, p,
        deflation_tolerance, eig_dump, num_vb, num_cb, max_no_success,
        save_BT_band, save_BT_SMF, read_BT_band, deflation_on, conv_method,
        conv_vec_method);
  }

#ifdef MPI_TIMING
  mpiTiming.section = mpiTiming.other_mpi;
#endif
#ifdef MPIP
  MPI_Pcontrol(3); // generate verbose report
  MPI_Pcontrol(0); // turn off profiling
  MPI_Pcontrol(2); // discard everything
#endif

  MPI_TOC(mpiTiming.lanczos_function, Tstart);

  cout_master << "----------------------\nFound " << d->neigv
              << " eigvalues\n----------------------\n\n";

  if (d->opt.ExecParam.ElCalc.ElAlg.ResFind == ElAlg_struct::Lancz) {

    d->WF = vector<waveFunction>(d->neigv);

#ifdef MOMENTUM_MATRIX

    // define Constant
    // double m0   		 = 5.11*pow(10,8);
    // double c    		 = 3.0*pow(10,17);
    // double hbar 	 	 = 6.5822*pow(10,-16);
    // double Momentum_Constant = m0/hbar/c/c;
    double Momentum_Constant = 1;

    int isx = d->cell_s[mpi_n3d_id];
    int iex = isx + d->cell_ln[mpi_n3d_id];
    int aNumber = (iex - isx) * d->geo.AtomsPerCellMax();

    complex *psiDegenInit = new complex[seglen];
    complex *psiDegenFine = new complex[seglen];
    complex *Hpsi = new complex[seglen];
    complex *HpsiDegen = new complex[seglen];

    complex *xpsi = new complex[seglen];
    complex *xpsiDegen = new complex[seglen];
    complex *xHpsi = new complex[seglen];
    complex *xHpsiDegen = new complex[seglen];
    complex *Hxpsi = new complex[seglen];
    complex *HxpsiDegen = new complex[seglen];

    cmatrix pXmatrix = Cmatrix(d->neigv, 8 * d->neigv);

#ifdef COMPUTE_TRANSP

    double *TransProbX = new real[2 * (d->neigv) * (d->neigv - 1)];
    int PointIndx = 0;

#ifdef DUMP_OVLAP

    complex *OVERLAP = new complex[(d->neigv) * (d->neigv - 1)];
    int OVLAPIndx = 0;

#endif // DUMP_OVLAP

#endif // COMPUTE_TRANSP

    real atomPos[3];
    complex temp1, temp2;

#ifdef COMPUTE_Px

    cout << endl
         << "[OPTICAL MATRIX] " << "Computing Momentum Matrix Px" << endl;

    for (int fstate = 0; fstate < d->neigv; fstate++) {

      int localIndx1 = 0;

      // writing Kramer's Degenerate State for final State
      for (int aIndx = 0; aIndx < aNumber; aIndx++) {
        for (int oIndx = 0; oIndx < 10; oIndx++) {
          psiDegenFine[aIndx * d->NBasisStates + oIndx + 10].r =
              d->wfc[fstate][aIndx * d->NBasisStates + oIndx].r;
          psiDegenFine[aIndx * d->NBasisStates + oIndx + 10].i =
              (-1) * d->wfc[fstate][aIndx * d->NBasisStates + oIndx].i;
          psiDegenFine[aIndx * d->NBasisStates + oIndx].r =
              (-1) * d->wfc[fstate][aIndx * d->NBasisStates + oIndx + 10].r;
          psiDegenFine[aIndx * d->NBasisStates + oIndx].i =
              d->wfc[fstate][aIndx * d->NBasisStates + oIndx + 10].i;
        }
      }

      // Computing :
      //		H*psiFinal
      //		H*psidegenFinal

      Hmult_spds_complete(Hpsi, d, 0.0, &d->wfc[fstate][0]);
      Hmult_spds_complete(HpsiDegen, d, 0.0, psiDegenFine);

      // Computing :
      //		x*H*psiFinal
      //		x*H*PsidegenFinal
      //		x*psiFinal
      //		x*psidegenFinal

      for (int cellIndx = isx; cellIndx < iex; cellIndx++) {
        for (int atomIndx = 0; atomIndx < d->geo.AtomsPerCellMax();
             atomIndx++) {

          d->geo.getPosition_strained(atomPos, cellIndx, atomIndx);
          for (int basisIndx = 0; basisIndx < d->NBasisStates; basisIndx++) {

            xHpsi[localIndx1].r = atomPos[0] * Hpsi[localIndx1].r;
            xHpsi[localIndx1].i = atomPos[0] * Hpsi[localIndx1].i;
            xHpsiDegen[localIndx1].r = atomPos[0] * HpsiDegen[localIndx1].r;
            xHpsiDegen[localIndx1].i = atomPos[0] * HpsiDegen[localIndx1].i;

            xpsi[localIndx1].r = atomPos[0] * d->wfc[fstate][localIndx1].r;
            xpsi[localIndx1].i = atomPos[0] * d->wfc[fstate][localIndx1].i;
            xpsiDegen[localIndx1].r = atomPos[0] * psiDegenFine[localIndx1].r;
            xpsiDegen[localIndx1].i = atomPos[0] * psiDegenFine[localIndx1].i;

            localIndx1++;

          } // end of basisIndx
        } // end of atomIndx
      } // end of cellIndx

      // Computing :
      //		H*x*psiFinal
      //		H*x*psidegenFinal

      Hmult_spds_complete(Hxpsi, d, 0.0, xpsi);
      Hmult_spds_complete(HxpsiDegen, d, 0.0, xpsiDegen);

      for (int istate = 0; istate < d->neigv; istate++) {

        temp1.r = 0.0;
        temp1.i = 0.0;
        temp2.r = 0.0;
        temp2.i = 0.0;

        // writing Kramer's Degenerate State for Initial State.
        for (int aIndx = 0; aIndx < aNumber; aIndx++) {
          for (int oIndx = 0; oIndx < 10; oIndx++) {
            psiDegenInit[aIndx * d->NBasisStates + oIndx + 10].r =
                d->wfc[istate][aIndx * d->NBasisStates + oIndx].r;
            psiDegenInit[aIndx * d->NBasisStates + oIndx + 10].i =
                (-1) * d->wfc[istate][aIndx * d->NBasisStates + oIndx].i;
            psiDegenInit[aIndx * d->NBasisStates + oIndx].r =
                (-1) * d->wfc[istate][aIndx * d->NBasisStates + oIndx + 10].r;
            psiDegenInit[aIndx * d->NBasisStates + oIndx].i =
                d->wfc[istate][aIndx * d->NBasisStates + oIndx + 10].i;
          }
        }

        // original(initial) vs. original(final)
        temp1 = vect_dot_vect(&d->wfc[istate][0], xHpsi, d);
        temp2 = vect_dot_vect(&d->wfc[istate][0], Hxpsi, d);

        pXmatrix[istate][8 * fstate].r = (double)(istate + 1);
        pXmatrix[istate][8 * fstate].i = (double)(fstate + 1);
        pXmatrix[istate][8 * fstate + 1].i =
            Momentum_Constant * (temp2.r - temp1.r);
        pXmatrix[istate][8 * fstate + 1].r =
            Momentum_Constant * (temp1.i - temp2.i);

        // original(initial) vs. degenerate(final)
        temp1.r = 0.0;
        temp1.i = 0.0;
        temp2.r = 0.0;
        temp2.i = 0.0;

        temp1 = vect_dot_vect(&d->wfc[istate][0], xHpsiDegen, d);
        temp2 = vect_dot_vect(&d->wfc[istate][0], HxpsiDegen, d);

        pXmatrix[istate][8 * fstate + 2].r = (double)(istate + 1);
        pXmatrix[istate][8 * fstate + 2].i = (-1.0) * (double)(fstate + 1);
        pXmatrix[istate][8 * fstate + 3].i =
            Momentum_Constant * (temp2.r - temp1.r);
        pXmatrix[istate][8 * fstate + 3].r =
            Momentum_Constant * (temp1.i - temp2.i);

        // degenerate(initial) vs. original(final)
        temp1.r = 0.0;
        temp1.i = 0.0;
        temp2.r = 0.0;
        temp2.i = 0.0;

        temp1 = vect_dot_vect(psiDegenInit, xHpsi, d);
        temp2 = vect_dot_vect(psiDegenInit, Hxpsi, d);

        pXmatrix[istate][8 * fstate + 4].r = (-1.0) * (double)(istate + 1);
        pXmatrix[istate][8 * fstate + 4].i = (double)(fstate + 1);
        pXmatrix[istate][8 * fstate + 5].i =
            Momentum_Constant * (temp2.r - temp1.r);
        pXmatrix[istate][8 * fstate + 5].r =
            Momentum_Constant * (temp1.i - temp2.i);

        // degenerate(initial) vs. degenerate(final)
        temp1.r = 0.0;
        temp1.i = 0.0;
        temp2.r = 0.0;
        temp2.i = 0.0;

        temp1 = vect_dot_vect(psiDegenInit, xHpsiDegen, d);
        temp2 = vect_dot_vect(psiDegenInit, HxpsiDegen, d);

        pXmatrix[istate][8 * fstate + 6].r = (-1.0) * (double)(istate + 1);
        pXmatrix[istate][8 * fstate + 6].i = (-1.0) * (double)(fstate + 1);
        pXmatrix[istate][8 * fstate + 7].i =
            Momentum_Constant * (temp2.r - temp1.r);
        pXmatrix[istate][8 * fstate + 7].r =
            Momentum_Constant * (temp1.i - temp2.i);

#ifdef COMPUTE_TRANSP

        if (istate > fstate) {

#ifdef DUMP_OVLAP

          OVERLAP[OVLAPIndx].r = (double)istate + 1.0;
          OVERLAP[OVLAPIndx++].i = (double)fstate + 1.0;
          OVERLAP[OVLAPIndx++] = vect_dot_vect(&d->wfc[istate][0], xpsi, d);
          //	OVERLAP[OVLAPIndx++] = vect_dot_vect(&d->wfc[istate][0], xpsi,
          //d);

#endif // DUMP_OVLAP

          TransProbX[PointIndx++] = (double)istate + 1.0;
          TransProbX[PointIndx++] = (double)fstate + 1.0;
          TransProbX[PointIndx++] = ev[istate] - ev[fstate];
          TransProbX[PointIndx] = pXmatrix[istate][8 * fstate + 1].r *
                                      pXmatrix[istate][8 * fstate + 1].r +
                                  pXmatrix[istate][8 * fstate + 1].i *
                                      pXmatrix[istate][8 * fstate + 1].i;
          TransProbX[PointIndx] += pXmatrix[istate][8 * fstate + 3].r *
                                       pXmatrix[istate][8 * fstate + 3].r +
                                   pXmatrix[istate][8 * fstate + 3].i *
                                       pXmatrix[istate][8 * fstate + 3].i;
          TransProbX[PointIndx] += pXmatrix[istate][8 * fstate + 5].r *
                                       pXmatrix[istate][8 * fstate + 5].r +
                                   pXmatrix[istate][8 * fstate + 5].i *
                                       pXmatrix[istate][8 * fstate + 5].i;
          TransProbX[PointIndx++] += pXmatrix[istate][8 * fstate + 7].r *
                                         pXmatrix[istate][8 * fstate + 7].r +
                                     pXmatrix[istate][8 * fstate + 7].i *
                                         pXmatrix[istate][8 * fstate + 7].i;
        }

#endif // COMPUTE_TRANSP

      } // end of istate
    } // end of fstate

    cout << endl
         << "[OPTICAL MATRIX] " << "Finished : Computing Momentum Matrix Px"
         << endl;

#ifdef DUMP_MOMENTUM
    cout << endl << "[OPTICAL MATRIX] " << "Writing Momentum Matrix Px" << endl;
    writeMomentumMatrix(d, &pXmatrix[0][0], 0, "Px");
#endif // DUMP_MOMENTUM

#ifdef COMPUTE_TRANSP

#ifdef DUMP_TRANSP
    cout << endl
         << "[OPTICAL MATRIX] " << "Writing Transition Probability X" << endl;
    writeTransitionRate(d, &TransProbX[0], 0, "TransX");

#ifdef DUMP_OVLAP
    cout << endl
         << "[OPTICAL MATRIX] " << "Writing Overlap Matrix for X" << endl;
    writeOverLap(d, &OVERLAP[0], 0, "ovlapX");
#endif // DUMP_OVLAP

#endif // DUMP_TRANSP

    PointIndx = 0;

#endif // COMPUTE_TRANSP

#endif // COMPUTE_Px

#ifdef COMPUTE_Py

    cout << endl
         << "[OPTICAL MATRIX] " << "Computing Momentum Matrix Py" << endl;

    for (int fstate = 0; fstate < d->neigv; fstate++) {

      int localIndx1 = 0;

      // writing Kramer's Degenerate State for final State
      for (int aIndx = 0; aIndx < aNumber; aIndx++) {
        for (int oIndx = 0; oIndx < 10; oIndx++) {
          psiDegenFine[aIndx * d->NBasisStates + oIndx + 10].r =
              d->wfc[fstate][aIndx * d->NBasisStates + oIndx].r;
          psiDegenFine[aIndx * d->NBasisStates + oIndx + 10].i =
              (-1) * d->wfc[fstate][aIndx * d->NBasisStates + oIndx].i;
          psiDegenFine[aIndx * d->NBasisStates + oIndx].r =
              (-1) * d->wfc[fstate][aIndx * d->NBasisStates + oIndx + 10].r;
          psiDegenFine[aIndx * d->NBasisStates + oIndx].i =
              d->wfc[fstate][aIndx * d->NBasisStates + oIndx + 10].i;
        }
      }

      // Computing :
      //		H*psiFinal
      //		H*psidegenFinal

      Hmult_spds_complete(Hpsi, d, 0.0, &d->wfc[fstate][0]);
      Hmult_spds_complete(HpsiDegen, d, 0.0, psiDegenFine);

      // Computing :
      //		x*H*psiFinal
      //		x*H*PsidegenFinal
      //		x*psiFinal
      //		x*psidegenFinal

      for (int cellIndx = isx; cellIndx < iex; cellIndx++) {
        for (int atomIndx = 0; atomIndx < d->geo.AtomsPerCellMax();
             atomIndx++) {

          d->geo.getPosition_strained(atomPos, cellIndx, atomIndx);
          for (int basisIndx = 0; basisIndx < d->NBasisStates; basisIndx++) {

            xHpsi[localIndx1].r = atomPos[1] * Hpsi[localIndx1].r;
            xHpsi[localIndx1].i = atomPos[1] * Hpsi[localIndx1].i;
            xHpsiDegen[localIndx1].r = atomPos[1] * HpsiDegen[localIndx1].r;
            xHpsiDegen[localIndx1].i = atomPos[1] * HpsiDegen[localIndx1].i;

            xpsi[localIndx1].r = atomPos[1] * d->wfc[fstate][localIndx1].r;
            xpsi[localIndx1].i = atomPos[1] * d->wfc[fstate][localIndx1].i;
            xpsiDegen[localIndx1].r = atomPos[1] * psiDegenFine[localIndx1].r;
            xpsiDegen[localIndx1].i = atomPos[1] * psiDegenFine[localIndx1].i;

            localIndx1++;

          } // end of basisIndx
        } // end of atomIndx
      } // end of cellIndx

      // Computing :
      //		H*x*psiFinal
      //		H*x*psidegenFinal

      Hmult_spds_complete(Hxpsi, d, 0.0, xpsi);
      Hmult_spds_complete(HxpsiDegen, d, 0.0, xpsiDegen);

      for (int istate = 0; istate < d->neigv; istate++) {

        temp1.r = 0.0;
        temp1.i = 0.0;
        temp2.r = 0.0;
        temp2.i = 0.0;

        // writing Kramer's Degenerate State for Initial State.
        for (int aIndx = 0; aIndx < aNumber; aIndx++) {
          for (int oIndx = 0; oIndx < 10; oIndx++) {
            psiDegenInit[aIndx * d->NBasisStates + oIndx + 10].r =
                d->wfc[istate][aIndx * d->NBasisStates + oIndx].r;
            psiDegenInit[aIndx * d->NBasisStates + oIndx + 10].i =
                (-1) * d->wfc[istate][aIndx * d->NBasisStates + oIndx].i;
            psiDegenInit[aIndx * d->NBasisStates + oIndx].r =
                (-1) * d->wfc[istate][aIndx * d->NBasisStates + oIndx + 10].r;
            psiDegenInit[aIndx * d->NBasisStates + oIndx].i =
                d->wfc[istate][aIndx * d->NBasisStates + oIndx + 10].i;
          }
        }

        // original(initial) vs. original(final)
        temp1 = vect_dot_vect(&d->wfc[istate][0], xHpsi, d);
        temp2 = vect_dot_vect(&d->wfc[istate][0], Hxpsi, d);

        pXmatrix[istate][8 * fstate].r = (double)(istate + 1);
        pXmatrix[istate][8 * fstate].i = (double)(fstate + 1);
        pXmatrix[istate][8 * fstate + 1].i =
            Momentum_Constant * (temp2.r - temp1.r);
        pXmatrix[istate][8 * fstate + 1].r =
            Momentum_Constant * (temp1.i - temp2.i);

        // original(initial) vs. degenerate(final)
        temp1.r = 0.0;
        temp1.i = 0.0;
        temp2.r = 0.0;
        temp2.i = 0.0;

        temp1 = vect_dot_vect(&d->wfc[istate][0], xHpsiDegen, d);
        temp2 = vect_dot_vect(&d->wfc[istate][0], HxpsiDegen, d);

        pXmatrix[istate][8 * fstate + 2].r = (double)(istate + 1);
        pXmatrix[istate][8 * fstate + 2].i = (-1.0) * (double)(fstate + 1);
        pXmatrix[istate][8 * fstate + 3].i =
            Momentum_Constant * (temp2.r - temp1.r);
        pXmatrix[istate][8 * fstate + 3].r =
            Momentum_Constant * (temp1.i - temp2.i);

        // degenerate(initial) vs. original(final)
        temp1.r = 0.0;
        temp1.i = 0.0;
        temp2.r = 0.0;
        temp2.i = 0.0;

        temp1 = vect_dot_vect(psiDegenInit, xHpsi, d);
        temp2 = vect_dot_vect(psiDegenInit, Hxpsi, d);

        pXmatrix[istate][8 * fstate + 4].r = (-1.0) * (double)(istate + 1);
        pXmatrix[istate][8 * fstate + 4].i = (double)(fstate + 1);
        pXmatrix[istate][8 * fstate + 5].i =
            Momentum_Constant * (temp2.r - temp1.r);
        pXmatrix[istate][8 * fstate + 5].r =
            Momentum_Constant * (temp1.i - temp2.i);

        // degenerate(initial) vs. degenerate(final)
        temp1.r = 0.0;
        temp1.i = 0.0;
        temp2.r = 0.0;
        temp2.i = 0.0;

        temp1 = vect_dot_vect(psiDegenInit, xHpsiDegen, d);
        temp2 = vect_dot_vect(psiDegenInit, HxpsiDegen, d);

        pXmatrix[istate][8 * fstate + 6].r = (-1.0) * (double)(istate + 1);
        pXmatrix[istate][8 * fstate + 6].i = (-1.0) * (double)(fstate + 1);
        pXmatrix[istate][8 * fstate + 7].i =
            Momentum_Constant * (temp2.r - temp1.r);
        pXmatrix[istate][8 * fstate + 7].r =
            Momentum_Constant * (temp1.i - temp2.i);

#ifdef COMPUTE_TRANSP

        if (istate > fstate) {

          TransProbX[PointIndx++] = (double)istate + 1.0;
          TransProbX[PointIndx++] = (double)fstate + 1.0;
          TransProbX[PointIndx++] = ev[istate] - ev[fstate];
          TransProbX[PointIndx] = pXmatrix[istate][8 * fstate + 1].r *
                                      pXmatrix[istate][8 * fstate + 1].r +
                                  pXmatrix[istate][8 * fstate + 1].i *
                                      pXmatrix[istate][8 * fstate + 1].i;
          TransProbX[PointIndx] += pXmatrix[istate][8 * fstate + 3].r *
                                       pXmatrix[istate][8 * fstate + 3].r +
                                   pXmatrix[istate][8 * fstate + 3].i *
                                       pXmatrix[istate][8 * fstate + 3].i;
          TransProbX[PointIndx] += pXmatrix[istate][8 * fstate + 5].r *
                                       pXmatrix[istate][8 * fstate + 5].r +
                                   pXmatrix[istate][8 * fstate + 5].i *
                                       pXmatrix[istate][8 * fstate + 5].i;
          TransProbX[PointIndx++] += pXmatrix[istate][8 * fstate + 7].r *
                                         pXmatrix[istate][8 * fstate + 7].r +
                                     pXmatrix[istate][8 * fstate + 7].i *
                                         pXmatrix[istate][8 * fstate + 7].i;
        }

#endif // COMPUTE_TRANSP

      } // end of istate
    } // end of fstate

    cout << endl
         << "[OPTICAL MATRIX] " << "Finished : Computing Momentum Matrix Py"
         << endl;

#ifdef DUMP_MOMENTUM
    cout << endl << "[OPTICAL MATRIX] " << "Writing Momentum Matrix Py" << endl;
    writeMomentumMatrix(d, &pXmatrix[0][0], 0, "Py");
#endif // DUMP_MOMENTUM

#ifdef COMPUTE_TRANSP

#ifdef DUMP_TRANSP
    writeTransitionRate(d, &TransProbX[0], 0, "TransY");
    cout << endl
         << "[OPTICAL MATRIX] " << "Writing Transition Probability Y" << endl;
#endif // DUMP_TRANSP

    PointIndx = 0;

#endif // COMPUTE_TRANSP

#endif // COMPUTE_Py

#ifdef COMPUTE_Pz

    cout << endl
         << "[OPTICAL MATRIX] " << "Computing Momentum Matrix Pz" << endl;

    for (int fstate = 0; fstate < d->neigv; fstate++) {

      int localIndx1 = 0;

      // writing Kramer's Degenerate State for final State
      for (int aIndx = 0; aIndx < aNumber; aIndx++) {
        for (int oIndx = 0; oIndx < 10; oIndx++) {
          psiDegenFine[aIndx * d->NBasisStates + oIndx + 10].r =
              d->wfc[fstate][aIndx * d->NBasisStates + oIndx].r;
          psiDegenFine[aIndx * d->NBasisStates + oIndx + 10].i =
              (-1) * d->wfc[fstate][aIndx * d->NBasisStates + oIndx].i;
          psiDegenFine[aIndx * d->NBasisStates + oIndx].r =
              (-1) * d->wfc[fstate][aIndx * d->NBasisStates + oIndx + 10].r;
          psiDegenFine[aIndx * d->NBasisStates + oIndx].i =
              d->wfc[fstate][aIndx * d->NBasisStates + oIndx + 10].i;
        }
      }

      // Computing :
      //		H*psiFinal
      //		H*psidegenFinal

      Hmult_spds_complete(Hpsi, d, 0.0, &d->wfc[fstate][0]);
      Hmult_spds_complete(HpsiDegen, d, 0.0, psiDegenFine);

      // Computing :
      //		x*H*psiFinal
      //		x*H*PsidegenFinal
      //		x*psiFinal
      //		x*psidegenFinal

      for (int cellIndx = isx; cellIndx < iex; cellIndx++) {
        for (int atomIndx = 0; atomIndx < d->geo.AtomsPerCellMax();
             atomIndx++) {

          d->geo.getPosition_strained(atomPos, cellIndx, atomIndx);
          for (int basisIndx = 0; basisIndx < d->NBasisStates; basisIndx++) {

            xHpsi[localIndx1].r = atomPos[2] * Hpsi[localIndx1].r;
            xHpsi[localIndx1].i = atomPos[2] * Hpsi[localIndx1].i;
            xHpsiDegen[localIndx1].r = atomPos[2] * HpsiDegen[localIndx1].r;
            xHpsiDegen[localIndx1].i = atomPos[2] * HpsiDegen[localIndx1].i;

            xpsi[localIndx1].r = atomPos[2] * d->wfc[fstate][localIndx1].r;
            xpsi[localIndx1].i = atomPos[2] * d->wfc[fstate][localIndx1].i;
            xpsiDegen[localIndx1].r = atomPos[2] * psiDegenFine[localIndx1].r;
            xpsiDegen[localIndx1].i = atomPos[2] * psiDegenFine[localIndx1].i;

            localIndx1++;

          } // end of basisIndx
        } // end of atomIndx
      } // end of cellIndx

      // Computing :
      //		H*x*psiFinal
      //		H*x*psidegenFinal

      Hmult_spds_complete(Hxpsi, d, 0.0, xpsi);
      Hmult_spds_complete(HxpsiDegen, d, 0.0, xpsiDegen);

      for (int istate = 0; istate < d->neigv; istate++) {

        temp1.r = 0.0;
        temp1.i = 0.0;
        temp2.r = 0.0;
        temp2.i = 0.0;

        // writing Kramer's Degenerate State for Initial State.
        for (int aIndx = 0; aIndx < aNumber; aIndx++) {
          for (int oIndx = 0; oIndx < 10; oIndx++) {
            psiDegenInit[aIndx * d->NBasisStates + oIndx + 10].r =
                d->wfc[istate][aIndx * d->NBasisStates + oIndx].r;
            psiDegenInit[aIndx * d->NBasisStates + oIndx + 10].i =
                (-1) * d->wfc[istate][aIndx * d->NBasisStates + oIndx].i;
            psiDegenInit[aIndx * d->NBasisStates + oIndx].r =
                (-1) * d->wfc[istate][aIndx * d->NBasisStates + oIndx + 10].r;
            psiDegenInit[aIndx * d->NBasisStates + oIndx].i =
                d->wfc[istate][aIndx * d->NBasisStates + oIndx + 10].i;
          }
        }

        // original(initial) vs. original(final)
        temp1 = vect_dot_vect(&d->wfc[istate][0], xHpsi, d);
        temp2 = vect_dot_vect(&d->wfc[istate][0], Hxpsi, d);

        pXmatrix[istate][8 * fstate].r = (double)(istate + 1);
        pXmatrix[istate][8 * fstate].i = (double)(fstate + 1);
        pXmatrix[istate][8 * fstate + 1].i =
            Momentum_Constant * (temp2.r - temp1.r);
        pXmatrix[istate][8 * fstate + 1].r =
            Momentum_Constant * (temp1.i - temp2.i);

        // original(initial) vs. degenerate(final)
        temp1.r = 0.0;
        temp1.i = 0.0;
        temp2.r = 0.0;
        temp2.i = 0.0;

        temp1 = vect_dot_vect(&d->wfc[istate][0], xHpsiDegen, d);
        temp2 = vect_dot_vect(&d->wfc[istate][0], HxpsiDegen, d);

        pXmatrix[istate][8 * fstate + 2].r = (double)(istate + 1);
        pXmatrix[istate][8 * fstate + 2].i = (-1.0) * (double)(fstate + 1);
        pXmatrix[istate][8 * fstate + 3].i =
            Momentum_Constant * (temp2.r - temp1.r);
        pXmatrix[istate][8 * fstate + 3].r =
            Momentum_Constant * (temp1.i - temp2.i);

        // degenerate(initial) vs. original(final)
        temp1.r = 0.0;
        temp1.i = 0.0;
        temp2.r = 0.0;
        temp2.i = 0.0;

        temp1 = vect_dot_vect(psiDegenInit, xHpsi, d);
        temp2 = vect_dot_vect(psiDegenInit, Hxpsi, d);

        pXmatrix[istate][8 * fstate + 4].r = (-1.0) * (double)(istate + 1);
        pXmatrix[istate][8 * fstate + 4].i = (double)(fstate + 1);
        pXmatrix[istate][8 * fstate + 5].i =
            Momentum_Constant * (temp2.r - temp1.r);
        pXmatrix[istate][8 * fstate + 5].r =
            Momentum_Constant * (temp1.i - temp2.i);

        // degenerate(initial) vs. degenerate(final)
        temp1.r = 0.0;
        temp1.i = 0.0;
        temp2.r = 0.0;
        temp2.i = 0.0;

        temp1 = vect_dot_vect(psiDegenInit, xHpsiDegen, d);
        temp2 = vect_dot_vect(psiDegenInit, HxpsiDegen, d);

        pXmatrix[istate][8 * fstate + 6].r = (-1.0) * (double)(istate + 1);
        pXmatrix[istate][8 * fstate + 6].i = (-1.0) * (double)(fstate + 1);
        pXmatrix[istate][8 * fstate + 7].i =
            Momentum_Constant * (temp2.r - temp1.r);
        pXmatrix[istate][8 * fstate + 7].r =
            Momentum_Constant * (temp1.i - temp2.i);

#ifdef COMPUTE_TRANSP

        if (istate > fstate) {

          TransProbX[PointIndx++] = (double)istate + 1.0;
          TransProbX[PointIndx++] = (double)fstate + 1.0;
          TransProbX[PointIndx++] = ev[istate] - ev[fstate];
          TransProbX[PointIndx] = pXmatrix[istate][8 * fstate + 1].r *
                                      pXmatrix[istate][8 * fstate + 1].r +
                                  pXmatrix[istate][8 * fstate + 1].i *
                                      pXmatrix[istate][8 * fstate + 1].i;
          TransProbX[PointIndx] += pXmatrix[istate][8 * fstate + 3].r *
                                       pXmatrix[istate][8 * fstate + 3].r +
                                   pXmatrix[istate][8 * fstate + 3].i *
                                       pXmatrix[istate][8 * fstate + 3].i;
          TransProbX[PointIndx] += pXmatrix[istate][8 * fstate + 5].r *
                                       pXmatrix[istate][8 * fstate + 5].r +
                                   pXmatrix[istate][8 * fstate + 5].i *
                                       pXmatrix[istate][8 * fstate + 5].i;
          TransProbX[PointIndx++] += pXmatrix[istate][8 * fstate + 7].r *
                                         pXmatrix[istate][8 * fstate + 7].r +
                                     pXmatrix[istate][8 * fstate + 7].i *
                                         pXmatrix[istate][8 * fstate + 7].i;
        }

#endif // COMPUTE_TRANSP

      } // end of istate
    } // end of fstate

    cout << endl
         << "[OPTICAL MATRIX] " << "Finished : Computing Momentum Matrix Pz"
         << endl;

#ifdef DUMP_MOMENTUM
    cout << endl << "[OPTICAL MATRIX] " << "Writing Momentum Matrix Pz" << endl;
    writeMomentumMatrix(d, &pXmatrix[0][0], 0, "Pz");
#endif // DUMP_MOMENTUM

#ifdef COMPUTE_TRANSP

#ifdef DUMP_TRANSP
    writeTransitionRate(d, &TransProbX[0], 0, "TransZ");
    cout << endl
         << "[OPTICAL MATRIX] " << "Writing Transition Probability Z" << endl;
#endif // DUMP_TRANSP

    PointIndx = 0;

#endif // COMPUTE_TRANSP

#endif // COMPUTE_Pz

    MPI_Barrier(MPI_COMM_WORLD);

    delete[] psiDegenInit;
    delete[] psiDegenFine;
    delete[] Hpsi;
    delete[] HpsiDegen;
    delete[] TransProbX;

    delete[] xpsi;
    delete[] xpsiDegen;
    delete[] xHpsi;
    delete[] xHpsiDegen;
    delete[] Hxpsi;
    delete[] HxpsiDegen;
#ifdef DUMP_OVLAP
    delete[] OVERLAP;
#endif // DUMP_OVLAP
    cout << endl << "[OPTICAL MATRIX] " << "Dynamic Memory Freed" << endl;

#endif // MOMENTUM_MATRIX

    // copy to waveFunction data structure
    for (int i = 0; i < d->neigv; i++) {

      waveFunction &wf = d->WF[i];
      wf.indx = i;
      wf.E = ev[i];
      wf.psi = vector<complex>(seglen);

      if (do_lancz_eigvect) {
        wf.valid = 1;

        for (int j = 0; j < seglen; j++) {

          wf.psi[j].r = d->wfc[i][j].r;
          wf.psi[j].i = d->wfc[i][j].i;
        }
      }

      else
        wf.valid = 0;
    }

    writeEigenvalues(d, 'a');

    if (do_lancz_eigvect) {

      // write_cubic_struct_wavefunctions(  d );
      static int EigenvectorCount = 0;
      for (int i = 0; i < d->neigv; i++) {
        EigenvectorCount++;
        writeEigenvector(d, d->wfc[i], EigenvectorCount - 1 /*i*/);
      }

      static int WaveFunctionCount = 0;
      writePsiSqr(d->wfc, d->neigv, d, WaveFunctionCount /*0*/);
      WaveFunctionCount = WaveFunctionCount + d->neigv;

      if (d->opt.ExecParam.ElCalc.ElOut.Visualization_3D) {
        // if(!d->opt.ExecParam.Strain.StrainOut.TensionOutput){
        writePos(d, "shape_dx_small");
        //}
      }

      if (d->opt.ExecParam.Output.AtomPosEquil)
        writePos(d, "equil");

      writePos(d, "equil_cell");

      if (d->opt.Dev.band_model == Dev_struct::Bands_1_s_nospin) {

        // in the single band case the data file is small enough
        // enough that we write out the full wavefunction
        print_psisq_gnu_pdb(d->wfc, d->neigv, d, 0);
      }

#ifdef CHECK_RESIDUAL
      check_residual(d);
#endif

      compute_ritz_residual_norm_deviation(d);

#if 0
            	compute_ritz_value(d);
            	compute_residual(d);
            	compute_eigenvector_norm(d);
#endif
    }

  }

  else if (d->opt.ExecParam.ElCalc.ElAlg.ResFind ==
           ElAlg_struct::BlockLanczos) {

    // Filter  out valid eigenpairs

    int num_eig = d->neigv;
    cvectr Hx = Cvectr(seglen);
    // char mask[d->neigv];
    char *mask = (char *)nml_malloc(sizeof(char) * d->neigv);
    memset(mask, 0, sizeof(char) * d->neigv);

    // If deviation is larger than tolerance(TBD) discard
    for (int i = 0; i < d->neigv; i++) {

      Hmult_spds_complete(Hx, d, 0.0, &d->wfc[i][0]);
      complex xHx = vect_dot_vect(&d->wfc[i][0], Hx, d);
      complex xx = vect_dot_vect(&d->wfc[i][0], &d->wfc[i][0], d);
      real ritz_value = xHx.r / xx.r;

      double E = ev[i];
      real deviation = fabs(E - ritz_value);

      if (E != 0)
        deviation /= fabs(E);

      // Discard the pair if the deviation is larger than 0.01 percent
      // if(deviation > (conv_tol * 1))
      if (deviation > 0.0001) {

        printf(
            "[BLKLANCZ] E = %1.15e, ritz_value = %1.15e, deviation = %1.15e\n",
            E, ritz_value, deviation);
        mask[i] = 1;
        num_eig--;
      }
    }

    // Degenerate eigenvalues with same eigenvectors, discard the one
    for (int i = 0; i < d->neigv; i++) {

      if (mask[i] == 0) {

        double ref_ev = ev[i];
        int vector_err;

        for (int j = i + 1; j < d->neigv; j++) {

          if (mask[j] == 0) {
            // if it seems to be degenerate, compare the eigenvectors
            // It should be orthonormal to each other, otherwise just discard
            // it! TBD TBD TBD TBD
            if (fabs(ref_ev - ev[j]) < conv_tol) {

              complex inner_prod_ij =
                  vect_dot_vect(&d->wfc[i][0], &d->wfc[j][0], d);
              // printf("i = %d, j = %d, inner_prod_ij = %1.10e\n", i, j,
              // inner_prod_ij.r);
            }
          }
        }
      }
    }

    // complex wf_tmp[num_eig][seglen];
    cmatrix wf_tmp = Cmatrix(num_eig, seglen);
    //		double ev_tmp[num_eig];
    double *ev_tmp = (double *)nml_malloc(sizeof(double) * num_eig);
    int tmp_indx = 0;

    for (int i = 0; i < d->neigv; i++) {

      if (mask[i] == 0) {

        for (int j = 0; j < seglen; j++) {

          wf_tmp[tmp_indx][j].r = d->wfc[i][j].r;
          wf_tmp[tmp_indx][j].i = d->wfc[i][j].i;
        }

        ev_tmp[tmp_indx] = ev[i];
        tmp_indx++;
      }
    }

    d->neigv = num_eig;

    for (int i = 0; i < d->neigv; i++) {
      for (int j = 0; j < seglen; j++) {

        d->wfc[i][j].r = wf_tmp[i][j].r;
        d->wfc[i][j].i = wf_tmp[i][j].i;
      }

      ev[i] = ev_tmp[i];
    }

    rm_cmatrix(&wf_tmp);
    nml_free(ev_tmp);

    printf(
        "[BLKLANCZ] Valid number of Eigenvalues after residual filtering: %d\n",
        num_eig);

    if (d->neigv > 0) {

      d->WF = vector<waveFunction>(num_eig);

      for (int i = 0; i < d->neigv; i++) {

        waveFunction &wf = d->WF[i];
        wf.indx = i;
        wf.E = ev[i];
        wf.psi = vector<complex>(seglen);

        if (do_lancz_eigvect) {

          wf.valid = 1;

          for (int j = 0; j < seglen; j++) {

            wf.psi[j].r = d->wfc[i][j].r;
            wf.psi[j].i = d->wfc[i][j].i;
          }
        }

        else
          wf.valid = 0;
      }

      writeEigenvalues(d, 'a');
    }

    if (do_lancz_eigvect) {

      if (d->neigv > 0) {

        // write_cubic_struct_wavefunctions(  d );

        for (int i = 0; i < d->neigv; i++)
          writeEigenvector(d, d->wfc[i], i);

        writePsiSqr(d->wfc, d->neigv, d, 0);

        if (d->opt.ExecParam.ElCalc.ElOut.Visualization_3D) {
          // if(!d->opt.ExecParam.Strain.StrainOut.TensionOutput){
          writePos(d, "shape_dx_small");
          //}
        }

        if (d->opt.ExecParam.Output.AtomPosEquil)
          writePos(d, "equil");

        writePos(d, "equil_cell");

        if (d->opt.Dev.band_model == Dev_struct::Bands_1_s_nospin) {
          // in the single band case the data file is small enough
          //  enough that we write out the full wavefunction
          print_psisq_gnu_pdb(d->wfc, d->neigv, d, 0);
        }

#ifdef CHECK_RESIDUAL
        check_residual(d);
#endif

        compute_ritz_residual_norm_deviation(d);
      }
    }

    nml_free(mask);
  }

  if (d->opt.ExecParam.ElCalc.ElOut.TraceEigval)
    str_free(fil_trace);

  if (d->opt.ExecParam.ElCalc.ElOut.SaveTridiagInfo) {
    str_free(fil_tridiag);
    str_free(fil_eigval);
  }

  nml_free(argList);
  rm_rvectr(&ev);

  return d->neigv;
}

void generate_wave_vectors(qd_struct d, rmatrix *kr, int *Nk) {
  real kxL_min = d->opt.ExecParam.ElCalc.ElAlg.k0[0];
  real kyL_min = d->opt.ExecParam.ElCalc.ElAlg.k0[1];
  real kzL_min = d->opt.ExecParam.ElCalc.ElAlg.k0[2];
  real kxL_max = d->opt.ExecParam.ElCalc.ElAlg.kf[0];
  real kyL_max = d->opt.ExecParam.ElCalc.ElAlg.kf[1];
  real kzL_max = d->opt.ExecParam.ElCalc.ElAlg.kf[2];

  *Nk = d->opt.ExecParam.ElCalc.ElAlg.Nk;

  /* CONSTRAINT:  if we're doing Raleigh-Ritz, allow only k=0 */
  if (d->opt.ExecParam.ElCalc.ElAlg.ResFind == ElAlg_struct::Rayleigh) {
    *Nk = 1;
    kxL_min = kyL_min = kzL_min = 0.0;
  }

  /* generate kr ... */
  *kr = Rmatrix(*Nk, 3);
  (*kr)[0][0] = kxL_min;
  (*kr)[0][1] = kyL_min;
  (*kr)[0][2] = kzL_min;

  for (int i = 1; i < *Nk; i++) {
    (*kr)[i][0] = kxL_min + (real)i * (kxL_max - kxL_min) / (real)((*Nk) - 1);
    (*kr)[i][1] = kyL_min + (real)i * (kyL_max - kyL_min) / (real)((*Nk) - 1);
    (*kr)[i][2] = kzL_min + (real)i * (kzL_max - kzL_min) / (real)((*Nk) - 1);
  }

  /* ... then check that definition is consistent with periodicity */
  if (!d->geo.isPeriodicX() && (kxL_min != 0.0 || kxL_max != 0.0))
    die("   No periodicity in x => k_x must vanish!\n");
  if (!d->geo.isPeriodicY() && (kyL_min != 0.0 || kyL_max != 0.0))
    die("   No periodicity in y => k_y must vanish!\n");
  if (!d->geo.isPeriodicZ() && (kzL_min != 0.0 || kzL_max != 0.0))
    die("   No periodicity in z => k_z must vanish!\n");
}

/* Calculates neigv eigenvalues and eigenvectors of the hamiltonian in the
   energy range [emin:emax].  max_iter is the maximum number of iterations. */
void compute_electronic_structure(qd_struct d) {
  real emin = d->opt.ExecParam.ElCalc.ElAlg.Lanczos.Erange_vb[0];
  real emax = d->opt.ExecParam.ElCalc.ElAlg.Lanczos.Erange_cb[1];

  int i, n_converged;
  rmatrix kr = NULL;
  int Nk = 0;

  if (d->opt.ExecParam.ElCalc.ElAlg.HamSize == ElAlg_struct::Full3d) {
    if (d->opt.ExecParam.ElCalc.ElOut.LocalBands) {
      d->_primitiveCell_granularity = 1;
      d->_bandstruct_on_full_domain = 0;
      d->_bandstruct_on_local_domain = 1;
    } else {
      d->_primitiveCell_granularity = 0;
      d->_bandstruct_on_full_domain = 1;
      d->_bandstruct_on_local_domain = 0;
    }
  } else if (d->opt.ExecParam.ElCalc.ElAlg.HamSize ==
             ElAlg_struct::VCA_1_cell) {
    if (d->opt.ExecParam.ElCalc.ElOut.LocalBands) {
      d->_primitiveCell_granularity = 1;
      d->_bandstruct_on_full_domain = 1;
      d->_bandstruct_on_local_domain = 1;
    } else {
      d->_primitiveCell_granularity = 1;
      d->_bandstruct_on_full_domain = 1;
      d->_bandstruct_on_local_domain = 0;
    }
  } else {
    die("compute_electronic_structure: Unimplemented Option\n");
  }

  /* override some options for consistency */
  d->_bandstruct_on_local_domain =
      d->_bandstruct_on_local_domain & d->_primitiveCell_granularity;

  /* perform conistency checks */
  if (d->opt.ExecParam.ElCalc.ElAlg.ResFind == ElAlg_struct::Direct &&
      mpi_n3d_numprocs != 1)
    die("   ERROR: Can't use 'direct' method with multiple processors!\n");

  generate_wave_vectors(d, &kr, &Nk);

/* Added for Zone-unfolding */
#ifdef UNFLD_PROJ
  int is_unitary;
  double k[3], tmp_vec[3], **kpGvecs, **Gvecs, **Rv, **Ur, **Ui;
  int NBASIS = 0, n_nanostruct[3], Nc_nano;
  /* For simple cubic sm cell size is conv unit cell cube edge */
  const double a_sm = A_SM_SC;
  const double bz_max = M_PI / a_sm;
  double *vr, *vi, *wvr, *wvi, **Ekmat;
  FILE *fpk, *fpbin;
  int eix;
  int Norb = d->NBasisStates;                 // NOTE
  int Natom = d->n_ham_tot / d->NBasisStates; // NOTE
  int max_ks = Nk;
  int n_nemo[3];
  n_nemo[0] = d->geo.cell_xmax - d->geo.cell_xmin + 1;
  n_nemo[1] = d->geo.cell_ymax - d->geo.cell_ymin + 1;
  n_nemo[2] = d->geo.cell_zmax - d->geo.cell_zmin + 1;

  int do_lancz_eigvect =
      (d->opt.ExecParam.ElCalc.ElOut.PsiSqr ||
       d->opt.ExecParam.ElCalc.ElOut.PsiSqr_cell ||
       d->opt.ExecParam.ElCalc.ElOut.PsiSqrPdb ||
       d->opt.ExecParam.ElCalc.ElOut.Psi4_IPR ||
       d->opt.ExecParam.ElCalc.ElOut.Eigvect) &&
      ((d->opt.ExecParam.ElCalc.ElAlg.ResFind == ElAlg_struct::Lancz) ||
       (d->opt.ExecParam.ElCalc.ElAlg.ResFind == ElAlg_struct::BlockLanczos));

  if (n_nemo[0] == 0 || !do_lancz_eigvect) {
    printf("Unfolding not needed because supercell has same size as the small "
           "cell. \n");
    printf("Unfolding can not be done if eigenvectros are not calculated. \n");
  } else
  //  if(n_nemo[0] > 0 && do_lancz_eigvect)
  { /* unfold */
    printf("\nNEERAVS TEST 2\n");
    int periodic_y = d->geo.isPeriodicY() ? 1 : 0;
    int periodic_z = d->geo.isPeriodicZ() ? 1 : 0;
    int primitive_p_cell; // NOTE

    if (d->geo.isCellZB_110_small())
      primitive_p_cell = 2;
    else if (d->geo.isCellZB_110())
      primitive_p_cell = 8;
    else if (d->geo.isCellZB_111())
      primitive_p_cell = 6;
    else
      primitive_p_cell = 4;

    printf("\nprimitive_p_cell = %d\n", primitive_p_cell);

    int natoms, nrec;

    if (mpi_n3d_id == mpi_n3d_masterid) {
      /* Adjust the unit cell size and basis size according to the periodicity
       * of supercell */
      /* n_nanostruct[i]'s are the number of unit cells in direction i, these
         are defined taking into account the periodicity of the lattice */
      define_small_cell_of_nanostructure(n_nemo, n_nanostruct, &NBASIS,
                                         &periodic_y, &periodic_z,
                                         primitive_p_cell);
      Nc_nano = n_nanostruct[0] * n_nanostruct[1] * n_nanostruct[2];

      /* Compute length of eigenvector */
      natoms = Nc_nano * NBASIS;
      nrec = natoms * Norb; /* d->n_ham_tot*/ /* length of eigenvector */

      /* Test print */
      printf("natoms = %3d  nrec = %3d NBASIS = %3d Nc_nano = %3d\n", natoms,
             nrec, NBASIS, Nc_nano);
      printf("n_nanostruct = (%d %d %d)\n", n_nanostruct[0], n_nanostruct[1],
             n_nanostruct[2]);

      /* Allocate storage */
      kpGvecs = R2tensor(0, Nc_nano, 0, 2); // allocate_real_matrix(0,Nc-1,0,2);
      Gvecs = R2tensor(0, Nc_nano, 0, 2);
      Rv = R2tensor(0, Nc_nano, 0, 2); // allocate_real_matrix(0,Nc-1,0,2);
      Ur = R2tensor(0, Nc_nano, 0,
                    Nc_nano); // allocate_real_matrix(0,Nc-1,0,Nc-1);
      Ui = R2tensor(0, Nc_nano, 0,
                    Nc_nano); // allocate_real_matrix(0,Nc-1,0,Nc-1);

      /* Generate R vectors */
      gen_Rvecs_cubic(Rv, n_nanostruct);

      /* Generate small cell Gvecs  */
      gen_SmCell_Gvecs(Gvecs, n_nanostruct, bz_max);

      /* Generate projection matrix */
      projection_matrix(Gvecs, Rv, Ur, Ui, Nc_nano);

      /* Check unitarity and report */
      is_unitary = check_unitarity(Ur, Ui, Nc_nano, EPS_UNITARY);

      if (is_unitary)
        printf("Unitarity check passed.\n");
      else
        printf("Unitarity check failed.\n");

      /* Allocate storage */
      vr = Rvectr(nrec);  // allocate_real_vector(0,nrec-1);
      vi = Rvectr(nrec);  // allocate_real_vector(0,nrec-1);
      wvr = Rvectr(nrec); // allocate_real_vector(0,Nc-1);
      wvi = Rvectr(nrec); // allocate_real_vector(0,Nc-1);

      /* Open file for printing available small cell k to file */
      fpk = fopen("SmCell_k.txt", "w");
      if (fpk == NULL)
        system_error("Cannot open k file for write.");

      /* Save Ekmatrix, k-vectors, and energies to binary file */
      char *filename = n3d_strdup_n(d->inputfile);
      n3d_FileTypeSet(&filename, "nd_proj_unfold", TRUE);

      fpbin = fopen(filename, "wb");
      //     fpbin = fopen("Projected_states.bin","wb");
      if (fpbin == NULL)
        system_error("Cannot open output binary file for write");

      fwrite(&max_ks, sizeof(int), 1, fpbin);
      fwrite(&Nc_nano, sizeof(int), 1, fpbin);

      /* compute # of bands, actual number of bands found from projection
       * algorithm will always be less than maxbands */
      int maxbands = NBASIS * Norb;

      fwrite(&maxbands, sizeof(int), 1, fpbin);
    }
  } /* if(n_nemo[0] == 1) */
#endif
  /* End: Added for Zone-unfolding */

  // initialize file that will contain values of discovered resonances
  writeEigenvalues(d, 'i');

  if (d->_primitiveCell_granularity) {
    for (i = 0; i < Nk; i++) {
      d->kxL = kr[i][0];
      d->kyL = kr[i][1];
      d->kzL = kr[i][2];
      cout << "calling local_bandstruct_zincBlende_cubic" << endl;
#define TETRAHEDRAL_LOCAL // do 5 atom calculation
#ifdef TETRAHEDRAL_LOCAL
      local_bandstruct_zincBlende_cubic(d);
#else
      local_bandstruct_cubic(d);
#endif
    }
  } else if (d->opt.ExecParam.ElCalc.ElAlg.ResFind == ElAlg_struct::Rayleigh) {
    min_rayleigh(d);
  } else if (d->opt.ExecParam.ElCalc.ElAlg.ResFind ==
             ElAlg_struct::Readstates) {
    readstates_driver_c(d);
  } else if (d->opt.ExecParam.ElCalc.ElAlg.ResFind ==
             ElAlg_struct::PCTraceminSolver) {
    ///////////////////////////////////////////////
    // Tracemin Driver for Chebyshev Mapping
#if (defined USE_TRACEMIN_PAR || defined USE_TRACEMIN_SER)
    ctracemin_driver(d); // First Try
#else
    die("This option is not compiled into this particular executable!\n");
#endif
  } else if (d->opt.ExecParam.ElCalc.ElAlg.ResFind ==
             ElAlg_struct::PQTraceminSolver) {
    // Tracemin Driver for Quadratic Mapping (using conjugate gradient)
#if (defined USE_TRACEMIN_PAR || defined USE_TRACEMIN_SER)
    qtracemin_driver(d); // Second Try
#else
    die("This option is not compiled into this particular executable!\n");
#endif
  }
  ///////////////////////////////////////////////
  else if (d->opt.ExecParam.ElCalc.ElAlg.ResFind ==
           ElAlg_struct::ParpackSolver) {
#ifdef NO_PARPACK
    cout << "arpack_driver is called" << endl;
    arpack_driver(d);
#else
    cout << "parpack_driver is called" << endl;
    parpack_driver(d);
#endif
  } else if (d->opt.ExecParam.ElCalc.ElAlg.ResFind ==
             ElAlg_struct::LanczAndRayleigh) {
    d->kxL = kr[0][0];
    d->kyL = kr[0][1];
    d->kzL = kr[0][2];

    int N = d->opt.ExecParam.ElCalc.ElAlg.RayleighRitz.EigValGuess.size();

    if (N == 0) {
      n_converged =
          (d->opt.ExecParam.ElCalc.ElAlg.ResFind == ElAlg_struct::Direct
               ? h_eigval_spds(d, emin, emax)
               : lanczos_driver(d));
    } else {
      d->WF = vector<waveFunction>(N);
      for (int i = 0; i < N; i++) {
        waveFunction &wf = d->WF[i];
        wf.indx = i;
        wf.E = d->opt.ExecParam.ElCalc.ElAlg.RayleighRitz.EigValGuess[i];
        wf.psi = vector<complex>(d->seg_ln[mpi_n3d_id]);
        wf.valid = 1;
      }
    }
    cout << "performing Rayleigh Quotient Iteration" << endl;
    double Tstart;
    MPI_TIC(Tstart);
    rayleighQuotientIter(d);
    MPI_TOC(mpiTiming.rayleigh_function, Tstart);
  } else { // either Direct or Lanczos Solver
    for (i = 0; i < Nk; i++) {
      d->kxL = kr[i][0];
      d->kyL = kr[i][1];
      d->kzL = kr[i][2];

      n_converged =
          (d->opt.ExecParam.ElCalc.ElAlg.ResFind == ElAlg_struct::Direct
               ? h_eigval_spds(d, emin, emax)
               : lanczos_driver(d));

/* Added for Zone-unfolding */
#ifdef UNFLD_PROJ
      if (!(n_nemo[0] == 0) && do_lancz_eigvect)
      //         if(n_nemo[0] > 1 && do_lancz_eigvect)
      {
        /* Allocate storage */
        Ekmat = R2tensor(0, n_converged, 0, Nc_nano);

        /* scale k according to Nc in x,y,z */
        for (int ii = 0; ii < 3; ii++)
          k[ii] = kr[i][ii] * M_PI / n_nanostruct[ii];

        /* Prepare data structures for MPI_Gatherv */
        rvectr sbuf, rbuf;
        ivectr displs, rcounts;
        int sum = 0;
        rbuf = Rvectr(2 * d->n_ham_tot);
        sbuf = Rvectr(2 * d->seg_ln[mpi_n3d_id]);
        displs = Ivectr(mpi_n3d_numprocs);
        rcounts = Ivectr(mpi_n3d_numprocs);

        rcounts[0] = 2 * d->seg_ln[mpi_n3d_masterid];
        displs[0] = 0;
        for (int proc = 1; proc < mpi_n3d_numprocs; proc++) {
          rcounts[proc] = 2 * d->seg_ln[proc];
          sum += rcounts[proc - 1];
          displs[proc] = sum;
        }

        eix = 0;
        for (int j = 0; j < n_converged; j++, eix++) {
          /* Copy j th eigenvector to sbuf array */
          int jj;
          for (int ii = 0; ii < d->seg_ln[mpi_n3d_id]; ii++) {
            jj = 2 * ii;
            sbuf[jj] = d->wfc[j][ii].r;
            jj++;
            sbuf[jj] = d->wfc[j][ii].i;
          }

          /* Gather it into rbuf array */
          MPI_Gatherv(sbuf, 2 * d->seg_ln[mpi_n3d_id], MPI_DOUBLE, rbuf,
                      rcounts, displs, MPI_DOUBLE, mpi_n3d_masterid,
                      MPI_COMM_WORLD);

          if (mpi_n3d_id == mpi_n3d_masterid) {
            jj = 0;
            for (int ii = 0; ii < d->n_ham_tot; ii++) {
              vr[ii] = rbuf[jj];
              jj++;
              vi[ii] = rbuf[jj];
              jj++;
              // printf("%3d  %12.10le  %12.10le\n",ii,vr[ii],vi[ii]);
            } /* for ii */

            /* Project out this energy eigenvector */
            project_probs(Ekmat, Ur, Ui, vr, vi, wvr, wvi, NBASIS, Norb,
                          Nc_nano, n_converged, eix, k, n_nanostruct);

            /* Display progress on the command line */
            printf("Supercell k %d of %d, Eigenvector %d of %d\n", i + 1,
                   max_ks, eix + 1, n_converged);
          } /* if master_id */
        } /* for j: n_converged */

        /* Free allocated memory */
        rm_rvectr(&sbuf);
        rm_rvectr(&rbuf);
        rm_ivectr(&displs);
        rm_ivectr(&rcounts);

        if (mpi_n3d_id == mpi_n3d_masterid) {
          /* Generate small cell kpGvecs */
          gen_SmCell_kpGvecs(Gvecs, k, kpGvecs, bz_max, Nc_nano);

          /* Print available small cell k to file */
          fprintf(fpk, "\n%12.10le %12.10le %12.10le  \n", k[0], k[1], k[2]);
          for (int ii = 0; ii < Nc_nano; ii++) {
            for (int jj = 0; jj < 3; jj++)
              tmp_vec[jj] = kpGvecs[ii][jj] / M_PI;
            fprintf(fpk, "%6d   %12.5le   %12.5le   %12.5le\n", ii, tmp_vec[0],
                    tmp_vec[1], tmp_vec[2]);
          } /* for ii */

          /* Write data in a binary file */
          fwrite(&(n_converged), sizeof(int), 1, fpbin);

          for (int ii = 0; ii < n_converged; ii++)
            fwrite(&(d->WF[ii].E), sizeof(double), 1, fpbin);

          /* Write out matrix of k-vectors */
          for (int ii = 0; ii < Nc_nano; ii++)
            for (int jj = 0; jj < 3; jj++)
              fwrite(&(kpGvecs[ii][jj]), sizeof(double), 1, fpbin);

          /* Now write out Ekmatrix */
          for (int ii = 0; ii < n_converged; ii++)
            for (int jj = 0; jj < Nc_nano; jj++)
              fwrite(&(Ekmat[ii][jj]), sizeof(double), 1, fpbin);

        } /* if master_id */
        rm_r2tensor(&Ekmat, 0);
      } /* if(!n_nemo[0] == 1) */
#endif
      /* End: Added for Zone-unfolding */

      // verify_ham ( d );
      // MPI_Finalize();
      // exit(1);
    }

/* Added for Zone-unfolding */
#ifdef UNFLD_PROJ
    if (mpi_n3d_id == mpi_n3d_masterid) {
      if (!(n_nemo[0] == 0) && do_lancz_eigvect)
      //         if(n_nemo[0] > 1 && do_lancz_eigvect)
      {
        /* Free up storage */
        rm_r2tensor(&kpGvecs, 0); // free_real_matrix(Gv,0,Nc-1,0);
        rm_r2tensor(&Gvecs, 0);
        rm_r2tensor(&Rv, 0); // free_real_matrix(Rv,0,Nc-1,0);
        rm_r2tensor(&Ur, 0); // free_real_matrix(Ur,0,Nc-1,0);
        rm_r2tensor(&Ui, 0); // free_real_matrix(Ui,0,Nc-1,0);
        rm_rvectr(&vr);      // free_real_vector(vr,0);
        rm_rvectr(&vi);      // free_real_vector(vi,0);
        rm_rvectr(&wvr);     // free_real_vector(wvr,0);
        rm_rvectr(&wvi);     // free_real_vector(wvi,0);
      } /* if(!n_nemo[0] == 1) */
    }
#endif
    /* End: Added for Zone-unfolding */
  }

  rm_rmatrix(&kr);
}

int h_eigval_spds(qd_struct d, real emin, real emax) {
  cmatrix wfc;

  cvectr x = Cvectr(d->n_ham_tot);
  cvectr y = Cvectr(d->n_ham_tot);
  cmatrix h = Coperator(d->n_ham_tot);

  if (mpi_n3d_numprocs != 1 && mpi_n3d_id == mpi_n3d_masterid)
    die("   ERROR: Can't use 'direct' method with multiple processors!\n");

  for (int i = 0; i < d->n_ham_tot; i++) {
    int kk;

    x[i].r = 1.0;
    if ((i % 100) == 0)
      printf("%d ", i), fflush(stdout);
    if (i)
      x[i - 1].r = 0;

    for (kk = 0; kk < d->n_ham_tot; kk++) {
      x[kk].r = x[kk].i = y[kk].r = y[kk].i = 0.0;
    }
    x[i].r = 1.0;

    Hmult_spds_complete(y, d, 0, x);

    for (int j = 0; j < d->n_ham_tot; j++) {
      h[j][i] = y[j];
    }
  }

#undef DEBUG_PRINT_AND_EXIT
#ifdef DEBUG_PRINT_AND_EXIT
  {
    FILE *fp = fopen("Hind", "w");
    int r, c;

    for (r = 0; r < d->n_ham_tot; r++) {
      for (c = 0; c < d->n_ham_tot; c++) {
        if (h[r][c].r != 0.0 || h[r][c].i != 0.0) {
          double norm = sqrt(h[r][c].r * h[r][c].r + h[r][c].i * h[r][c].i);
          fprintf(fp, "%d %d %e %e %e\n", r, c, h[r][c].r, h[r][c].i, norm);
        }
      }
    }
    fclose(fp);
  }
  cout_master << "Printing H and exiting!\n";
  exit(1);
#endif

  d->WF = vector<waveFunction>(d->n_ham_tot);

  if (!d->opt.ExecParam.ElCalc.ElOut.PsiSqr &&
      !d->opt.ExecParam.ElCalc.ElOut.PsiSqr_cell &&
      !d->opt.ExecParam.ElCalc.ElOut.Psi4_IPR &&
      !d->opt.ExecParam.ElCalc.ElOut.PsiSqrPdb &&
      !d->opt.ExecParam.ElCalc.ElOut.Eigvect) {

    ceigval_ns_full(x, h);

    // copy to waveFunction data structure
    for (int i = 0; i < d->n_ham_tot; i++) {
      waveFunction &wf = d->WF[i];
      wf.indx = i;
      wf.E = x[i].r;
      wf.valid = 0;
    }

  } else {
    wfc = Coperator(d->n_ham_tot);

    ceigvec_ns_full(x, wfc, h);

    wfc = ctranspose_nocopy(wfc);
    eigvec_srt(x, wfc, d->n_ham_tot);

    print_psisq_gnu_pdb(wfc, d->n_ham_tot, d, 0);

    // copy to waveFunction data structure
    for (int i = 0; i < d->n_ham_tot; i++) {
      waveFunction &wf = d->WF[i];
      wf.indx = i;
      wf.E = x[i].r;
      wf.psi = vector<complex>(d->n_ham_tot);
      if (d->opt.ExecParam.ElCalc.ElAlg.ResFind == ElAlg_struct::Lancz) {
        wf.valid = 1;
        for (int j = 0; j < d->n_ham_tot; j++) {
          wf.psi[j].r = d->wfc[i][j].r;
          wf.psi[j].i = d->wfc[i][j].i;
        }
      } else {
        wf.valid = 0;
      }
    }

    for (int i = 0; i < d->n_ham_tot; i++)
      writeEigenvector(d, wfc[i], i);

    rm_cmatrix(&wfc);
  }

  writeEigenvalues(d, 'a');

  rm_cvectr(&x);
  rm_cvectr(&y);
  rm_cmatrix(&h);

  return d->n_ham_tot;
}

void verify_ham(qd_struct d) {
  int ii, jj, kk;
  char dum[100];
  FILE *fp_loc;
  int s0 = d->seg_s[mpi_n3d_id];
  int sf = d->seg_s[mpi_n3d_id] + d->seg_ln[mpi_n3d_id];
  cvectr yy = Cvectr(sf - s0);
  cvectr xx = Cvectr(sf - s0);

  cout_master << "DEBUG:  printing hamiltonian\n\n";

  sprintf(dum, "H%3.2f_%d\n", d->kxL, mpi_n3d_id);
  dum[7] = '\0';
  fp_loc = fopen(dum, "w");

  for (ii = 0; ii < d->n_ham_tot; ii++) {
    cout_master << ii << "...";

    for (kk = 0; kk < sf - s0; kk++) {
      xx[kk].r = xx[kk].i = yy[kk].r = yy[kk].i = 0.0;
    }
    if (ii >= s0 && ii < sf)
      xx[ii - s0].r = 1.0;

    Hmult_spds_complete(yy, d, 0.0, xx);

    for (jj = 0; jj < sf - s0; jj++) {
      if (yy[jj].r != 0.0 || yy[jj].i != 0.0)
        fprintf(fp_loc, "%8d %8d %25.12e %25.12e\n", jj + s0, ii, yy[jj].r,
                yy[jj].i);
    }
    fprintf(fp_loc, "\n");
  }

  cout_master << "\ndone\n\n";

  fclose(fp_loc);
  rm_cvectr(&yy);
  rm_cvectr(&xx);
}

void verify_ham_col(qd_struct d, int col0, int colf) {
  int ii, jj, kk;
  char dum[100];
  FILE *fp_loc;
  int s0 = d->seg_s[mpi_n3d_id];
  int sf = d->seg_s[mpi_n3d_id] + d->seg_ln[mpi_n3d_id];
  cvectr yy = Cvectr(sf - s0);
  cvectr xx = Cvectr(sf - s0);

  cout_master << "DEBUG:  printing a hamiltonian column\n\n";

  sprintf(dum, "Hcol(%d,%d)_%3.2f_%d\n", col0, colf - 1, d->kxL, mpi_n3d_id);
  fp_loc = fopen(dum, "w");

  for (ii = col0; ii < colf; ii++) {
    cout_master << ii << "...";

    for (kk = 0; kk < sf - s0; kk++) {
      xx[kk].r = xx[kk].i = yy[kk].r = yy[kk].i = 0.0;
    }
    if (ii >= s0 && ii < sf)
      xx[ii - s0].r = 1.0;

    Hmult_spds_complete(yy, d, 0.0, xx);

    for (jj = 0; jj < sf - s0; jj++) {
      if (yy[jj].r != 0.0 || yy[jj].i != 0.0)
        fprintf(fp_loc, "%8d %8d %25.12e %25.12e\n", jj + s0, ii, yy[jj].r,
                yy[jj].i);
    }
    fprintf(fp_loc, "\n");
  }

  cout_master << "\ndone\n\n";

  fclose(fp_loc);
  rm_cvectr(&yy);
  rm_cvectr(&xx);
}

/* Check unitarity of matrix Ur, Ui[0..n-1][0..n-1].  eps is effectively 0 */
/* Return 1 for unitary, 0 for not unitary (to within eps, of course). */
int check_unitarity(double **Ur, double **Ui, int n, double eps) {
  int i, j, m, pass;
  double sumr, sumi;

  /* Set pass test to 1 and compute matrix (U+)*U:  should be 1 */
  i = 0;
  pass = 1;
  while ((i < n) && (pass)) {
    j = 0;
    while ((j < n) && (pass)) {
      sumr = 0.0;
      sumi = 0.0;

      /* Compute (U+)*U[i][j] */
      for (m = 0; m < n; m++) {
        sumr += (Ur[m][i] * Ur[m][j] + Ui[m][i] * Ui[m][j]);
        sumi += (Ur[m][i] * Ui[m][j] - Ui[m][i] * Ur[m][j]);
      } /* for m */

      /* Different checks for diagonal and off-diag elements */
      if (j == i)
        pass = ((fabs(sumi) < eps) && (fabs(1.0 - sumr) < eps));
      else
        pass = ((fabs(sumi) < eps) && (fabs(sumr) < eps));
      j++;
    } /* while j */
    // printf("\ni = %d",i);
    i++;
  } /* while i */

  return (pass);
} /* check_unitarity */

void define_small_cell_of_nanostructure(int *n_nemo, int *n_nanostruct,
                                        int *NBASIS, int *periodic_y,
                                        int *periodic_z, int primitive_p_cell) {
  if (*periodic_y == 1 && *periodic_z == 1) /* Bulk */
  {
    n_nanostruct[0] = n_nemo[0];
    n_nanostruct[1] = n_nemo[1];
    n_nanostruct[2] = n_nemo[2];
    *NBASIS = 2 * primitive_p_cell;
  }

  if (*periodic_y == 1 &&
      *periodic_z == 0) /* Quantum well confined in z direction */
  {
    n_nanostruct[0] = n_nemo[0];
    n_nanostruct[1] = n_nemo[1];
    n_nanostruct[2] = 1;
    *NBASIS = 2 * primitive_p_cell * n_nemo[2];
  }

  if (*periodic_y == 0 &&
      *periodic_z == 1) /* Quantum well confined in y direction */
  {
    n_nanostruct[0] = n_nemo[0];
    n_nanostruct[1] = 1;
    n_nanostruct[2] = n_nemo[2];
    *NBASIS = 2 * primitive_p_cell * n_nemo[1];
  }

  if (*periodic_y == 0 &&
      *periodic_z == 0) /* Nanowire confined in y and z directions */
  {
    n_nanostruct[0] = n_nemo[0];
    n_nanostruct[1] = 1;
    n_nanostruct[2] = 1;
    *NBASIS = 2 * primitive_p_cell * n_nemo[1] * n_nemo[2];
  }
} /* define_small_cell_of_nanostructure */

void gen_Rvecs_cubic(double **R_a, int Nrect[]) /* renamed */
/* Note that first index of R_a starts from 0 */
{
  int i, nemo[3];
  int nemo_ix, aix;

  /* Now loop through the NEMO cubes */
  for (nemo[0] = 0; nemo[0] < Nrect[0]; nemo[0]++)
    for (nemo[1] = 0; nemo[1] < Nrect[1]; nemo[1]++)
      for (nemo[2] = 0; nemo[2] < Nrect[2]; nemo[2]++) {
        /* Compute NEMO index */
        nemo_ix = ((nemo[0] * Nrect[1] + nemo[1]) * Nrect[2] +
                   nemo[2]); //*FCC_P_CUBE; // NO NEED TO MULTIPLY BY FCC_P_CUBE
        aix = nemo_ix;

        for (i = 0; i < 3; i++)
          R_a[aix][i] = nemo[i];

      } /* for nemo[2] */

} /* gen_Rvecs */

void gen_SmCell_Gvecs(double **Gvecs, int *Nrect, double bz_max) {
  int i, j, Nmin[3], Nmax[3], ix, iy, iz, iix, iiy, iiz, vindx;
  double Gv[3], kpG[3];
  double **Bvec;

  Bvec = R2tensor(0, 2, 0, 2);

  /* Assign Supercell Recip Lattice Vectors; for SC they  */
  /* are one-component each and orthogonal.                               */
  /* Each vector is stored as a COLUMN of Bvec.                   */
  double twoPI = 2.0 * M_PI;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++)
      Bvec[j][i] = 0.0;

    /* Units are inverse small cell cube edge */
    Bvec[i][i] = twoPI / ((double)Nrect[i]);
  } /* for i */

  /* Compute min and max integers for Supercell Recip vects */
  for (i = 0; i < 3; i++) {
    if (Nrect[i] % 2) {
      Nmax[i] = (Nrect[i] - 1) / 2;
      Nmin[i] = -Nmax[i];
    } else {
      Nmin[i] = -(Nrect[i] - 2) / 2;
      Nmax[i] = Nrect[i] / 2;
    }
  } /* for i */

  /* Now compute bulk vectors for each recip latt vect */
  for (ix = Nmin[0]; ix <= Nmax[0]; ix++)
    for (iy = Nmin[1]; iy <= Nmax[1]; iy++)
      for (iz = Nmin[2]; iz <= Nmax[2]; iz++) {
        iix = ix - Nmin[0];
        iiy = iy - Nmin[1];
        iiz = iz - Nmin[2];

        /* For SC dot prod is trivial so do it manually */
        Gv[0] = ((double)ix) * Bvec[0][0]; /* + iy*Bvec[1][0] + iz*Bvec[2][0] */
        Gv[1] = ((double)iy) * Bvec[1][1]; /* + ix*Bvec[0][1] + iz*Bvec[2][1] */
        Gv[2] = ((double)iz) * Bvec[2][2]; /* + ix*Bvec[0][2] + iy*Bvec[1][2] */

        vindx = ((iix * Nrect[1] + iiy) * Nrect[2] + iiz);

        for (i = 0; i < 3; i++)
          Gvecs[vindx][i] = Gv[i]; /// M_PI;
      } /* for iz*/
  rm_r2tensor(&Bvec, 0);
} /* gen_SmCell_Gvecs */

void gen_SmCell_kpGvecs(double **Gvecs, double *k, double **kpGvecs,
                        double bz_max, int Nc) {
  double kpG[3], ktemp[3];

  for (int ii = 0; ii < Nc; ii++) {
    for (int i = 0; i < 3; i++) {
      kpG[i] = Gvecs[ii][i] + k[i];
      /* Test if out of1st BZ; shift back if it is */
      if (kpG[i] > bz_max)
        kpG[i] -= 2.0 * bz_max;
      else if (kpG[i] < -bz_max)
        kpG[i] += 2.0 * bz_max;
    } /* for i */

    for (int i = 0; i < 3; i++)
      kpGvecs[ii][i] = kpG[i];

  } /* for ii */
} /* gen_SmCell_kpGvecs */

/* Construct unitary projection matrix Ur, Ui[0..Nc-1][0..Nc-1] using G vectors
 */
/* Gv[0..Nc-1][0..2] and FCC prim cell vectors Rv[0..Nc-1][0..2] where row is */
/* vector and cols are components, 0->x etc. */
void projection_matrix(double **Gv, double **Rv, double **Ur, double **Ui,
                       int Nc) {
  int r, c, j;
  double norm, RGdot;

  /* Overall normalization */
  norm = sqrt((double)Nc);

  /* Loop over rows, then cols of Ur, Ui */
  for (r = 0; r < Nc; r++)
    for (c = 0; c < Nc; c++) {
      /* Compute Rv[r].Gv[c]  */
      RGdot = 0.0;
      for (j = 0; j < 3; j++)
        RGdot += Rv[r][j] * Gv[c][j];

      /* Compute exp(I*Rv[r].Gv[c])/sqrt(Nc), store in Ur, Ui */
      Ur[r][c] = cos(RGdot) / norm;
      Ui[r][c] = sin(RGdot) / norm;
    } /* for c */
} /* projection_matrix */

/* Project out probabilities for wavevectors (as ordered by gen_kpGvec) from
 * vector     */
/* vr, vi[0..nrec], where nrec = 2*Nc*FCC_P_CUBE*Norb. Work-space */
/* is the vector wvr, wvi[0..nrec-1].  Store result (after projecting out and
 * summing   */
/* over all orbitals) as row eix of matrix Ekmat[0..Negy-1][0..Nc-1] where Negy
 * =       */
/* #energy eigenvectors we read in from the nemo3d output.  Cols of Ekmat are */
/* prim cell k's (ie supercell kpGvecs's) with order as determined by
 * gen_kpGvecs.      */
/* Use unitary projection matrix Ur, Ui[0..Nc-1][0..Nc-1] to do projection. */
void project_probs(double **Ekmat, double **Ur, double **Ui, double vr[],
                   double vi[], double wvr[], double wvi[], int NBASIS,
                   int Norb, int Nc_nano, int Negy, int eix, double k[],
                   int *n_nanostruct) {
  int i, j, orb, ip;
  double norm, kRdot, exr, exi, tr, ti;
  int n, ii, jj, Gix, indx;
  int ix, iy, iz, Rv[3], row, row_rv, atom;
  int Norbnano, CNorbnano;

  Norbnano = NBASIS * Norb;

  /* Zero row eix of Ekmat since we deal with eix-th energy only */
  for (j = 0; j < Nc_nano; j++)
    Ekmat[eix][j] = 0.0;

  /* First dephase the matrix W for this state all orbs */
  for (ix = 0; ix < n_nanostruct[0]; ix++) {
    Rv[0] = ix;
    for (iy = 0; iy < n_nanostruct[1]; iy++) {
      Rv[1] = iy;
      for (iz = 0; iz < n_nanostruct[2]; iz++) {
        Rv[2] = iz;
        //                                row_rv = Rv[2] + n_nanostruct[2]*Rv[1]
        //                                +
        //                                n_nanostruct[2]*n_nanostruct[1]*Rv[0];
        indx =
            ix * n_nanostruct[2] * n_nanostruct[1] + iy * n_nanostruct[2] + iz;

        /* Compute k.R */
        kRdot = 0.0;
        for (i = 0; i <= 2; i++)
          kRdot += ((double)Rv[i]) * k[i];

        /* Dephase; replace entry in matrix by dephased entry */
        for (int atom = 0; atom < NBASIS; atom++) {

          /* Compute exp(-i*k.R); remember cosine is EVEN! */
          exr = cos(kRdot);
          exi = -sin(kRdot);

          /* Dephase; replace entry in matrix by dephased entry */
          for (int orb = 0; orb < Norb; orb++) {
            //                                                        row =
            //                                                        CNorbnano[ix]
            //                                                        + atomorb;
            row = indx * NBASIS * Norb + atom * Norb + orb;
            tr = exr * vr[row] - exi * vi[row];
            ti = exr * vi[row] + exi * vr[row];
            vr[row] = tr;
            vi[row] = ti;
          } /* for atomorb */
        } /* for atom */
      } /* for iz */
    } /* for iy */
  } /* for ix */

  // printf("\n\nTEST: Dephasing done.\n\n");

  /* Do the projection and normalization one state at a time */
  /* Proceed by orbital type on projections */
  /* Project out with herm conj of U */
  for (i = 0; i < Nc_nano; i++) {
    for (int atomorb = 0; atomorb < Norbnano; atomorb++) {
      /* ii is index in v, wv */
      //                        ii = CNorbnano[i] + atomorb;
      ii = NBASIS * Norb * i + atomorb;
      wvr[ii] = 0.0;
      wvi[ii] = 0.0;

      for (j = 0; j < Nc_nano; j++) {
        /* jj is index in v */
        //                                jj = CNorbnano[j] + atomorb;
        jj = NBASIS * Norb * j + atomorb;
        wvr[ii] += (Ur[j][i] * vr[jj] + Ui[j][i] * vi[jj]);
        wvi[ii] += (Ur[j][i] * vi[jj] - Ui[j][i] * vr[jj]);
      } /* for j */
    } /* for atomorb */
  } /* for i */

  // printf("\n\nTEST: Projection done.\n\n");

  /* Now normalize bulk coeffs for each G_n for this state */
  for (i = 0; i < Nc_nano; i++) {
    norm = 0.0;
    for (int atomorb = 0; atomorb < Norbnano; atomorb++) {
      //                        ii = CNorbnano[i] + atomorb;
      ii = NBASIS * Norb * i + atomorb;
      norm += wvr[ii] * wvr[ii] + wvi[ii] * wvi[ii];
    } /* for atomorb */
    Ekmat[eix][i] = norm;
  } /* for i */

  // printf("\n\nTEST: Normalization done.\n\n");

} /* project_probs */
