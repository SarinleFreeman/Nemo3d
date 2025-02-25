#include "ctracemin_driver.h"
#include "eigsys3d_par.h"

// #define MATRIX_MULT

#ifdef MATRIX_MULT
#include "h_cmat_mult.h"
#endif

int ctracemin_driver(qd_struct d) {

#ifdef USE_TRACEMIN_PAR /* USE_TRACEMIN_PAR */
  if (d->rc == NULL)
    d->rc = Cvectr(d->nvmap[mpi_n3d_id] *
                   d->NBasisStates); // kept this line for parpack_driver.c, not
                                     // sure what it's for???
  int neig = d->opt.ExecParam.ElCalc.ElAlg.PCTracemin
                 .NumEigVal; // number of eigenvalues wanted
  real refE = d->opt.ExecParam.ElCalc.ElAlg.PCTracemin
                  .RefEigVal; // shift the matrix A for it to be hermitian
                              // positive definite (A - refE*I)
  real tol =
      d->opt.ExecParam.ElCalc.ElAlg.PCTracemin.Tolerance; // accuracy wanted
  real seig =
      d->opt.ExecParam.ElCalc.ElAlg.PCTracemin
          .WantEigRange[0]; // min. of interval where eigenvalues are wanted
  real leig =
      d->opt.ExecParam.ElCalc.ElAlg.PCTracemin
          .WantEigRange[1]; // max. of interval where eigenvalues are wanted
  real saeig = d->opt.ExecParam.ElCalc.ElAlg.PCTracemin
                   .AllEigRange[0]; // min. all eigenvalues
  real laeig = d->opt.ExecParam.ElCalc.ElAlg.PCTracemin
                   .AllEigRange[1]; // max all eigenvalues
  int maxi = d->opt.ExecParam.ElCalc.ElAlg.PCTracemin.maxi;
  int chebdegree = d->opt.ExecParam.ElCalc.ElAlg.PCTracemin
                       .ChebyDeg; // degree of chebyshev polynomial
  int n = d->n_ham_tot;           // dimension of eigenproblem
  int nloc =
      d->seg_ln[mpi_n3d_id]; // number of vector elements on current processor
  int info;                  // flag for reverse communication; initialized to 0
  // declare important MPI parameters
  MPI_Fint comm = MPI_Comm_c2f(
      MPI_COMM_WORLD); // transform C into Fortran communicator (currently
                       // tracemin uses MPI_COMM_WORLD and comm doesn't matter)
  int rank; // rank of the current processor
  int np;   // number of processors
  // allocate arrays containing result of tracemin
  rvectr E = Rvectr(neig); // real eigenvalues
  // allocate necessary temporary arrays
  cvectr H = Cvectr(neig * neig);
  cvectr Hloc = Cvectr(neig * neig);
  cvectr G = Cvectr(neig * neig);
  cvectr ZWork = Cvectr(max(1, 2 * neig - 1));
  rvectr RWork = Rvectr(max(1, 3 * neig - 2));
  ivectr IWork = Ivectr(neig);

#ifdef MATRIX_MULT
  cmatrix R = Cmatrix(neig, nloc); // estimate of residual
  cmatrix Y = Cmatrix(neig, nloc); // complex eigenvectors
  cmatrix Q = Cmatrix(neig, nloc);
  cmatrix Qt = Cmatrix(neig, nloc);
  cmatrix S = Cmatrix(neig, nloc);
  cmatrix Vn = Cmatrix(neig, nloc);
  cmatrix Vp = Cmatrix(neig, nloc);
  cmatrix T = Cmatrix(neig, nloc);
#else
  cvectr R = Cvectr(nloc * neig); // estimate of residual
  cvectr Y = Cvectr(nloc * neig); // complex eigenvectors
  cvectr Q = Cvectr(nloc * neig);
  cvectr Qt = Cvectr(nloc * neig);
  cvectr S = Cvectr(nloc * neig);
  cvectr Vn = Cvectr(nloc * neig);
  cvectr Vp = Cvectr(nloc * neig);
  cvectr T = Cvectr(nloc * neig);
#endif

  int i, j, k, mcsinfo, einfo;
  real gamma; // parameter used by tracemin
  // real seig; //smallest eigenvalue you are looking for real leig; //largest
  // eigenvalue you are looking for real saeig; //smallest of alleigenvalues
  // real laeig; //largest of all eigenvalues
  int nloct; // total number of elements on this processor: nloc*neig
  // int chebdegree; //degree of chebyshev polynomial int maxi; //maximum number
  // of iterations of tracemin (there are neig*2*6*2 matvecs per iteration)
  complex ytloc; // temporary variables used in computation of eigenvalues after
                 // tracemin (next 4 lines)
  complex yyloc;
  real yt;
  real yy;
  rvectr Rnrm = Rvectr(neig); // norm of the residual of of individual
                              // eigenpairs
  // variables used in sorting of the eigenvalues after tracemin
  rvectr orderedE;
  ivectr orderedI;
  double dtmp;
  int itmp;
  // declare some constants needed in the code
  int IONE = 1;
  real DTWO = 2.0;
  int iseed[4];

  printf("ENTERING PCTRACEMIN DRIVER\n");
  printf("setup some variables\n");
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); // find out my rank
  MPI_Comm_size(MPI_COMM_WORLD, &np);   // find out number of processors
  ////if(SMALLER_THEN_ALL_EIG <= 0){
  ////  refE = (real)(SMALLER_THEN_ALL_EIG - 1);
  ////}
  ////else{
  refE = 0.0;
  ////}
  // seig = emin1;                ////emin1 - refE; ///-30-refE; //emin1 -
  // refE;//???                   //define an interval where to look for the
  // eigenvalues leig = emax2;                ////emax2 - refE; ///30-refE;
  // //emax2 - refE;//???                   //THIS PARAMETER AND ABOVE ONE CAN
  // BE AN OPTION IN THE INPUT FILE!!! saeig= SMALLER_THEN_ALL_EIG;
  // ////SMALLER_THEN_ALL_EIG - refE; //define an interval in which all
  // eigenvalues are contained laeig= LARGER_THEN_ALL_EIG;
  // ////LARGER_THEN_ALL_EIG  - refE; //THIS PARAMETER AND ABOVE ONE CAN BE AN
  // OPTION IN THE INPUT FILE!!!
  gamma = min(DTWO / ((saeig - seig) * (saeig - leig)),
              DTWO / ((laeig - seig) * (laeig - leig)));
  nloct = neig *
          nloc; // number of elements of all vectors combined on this processor
  // maxi       = 40;  //50;   (comment by marta) //THIS PARAMETER CAN BE AN
  // OPTION IN THE INPUT FILE!!! chebdegree = 240; //6; //THIS PARAMETER CAN BE
  // AN OPTION IN THE INPUT FILE!!!
  iseed[0] = 103 + rank; // set up seed for the random number generator (making
                         // sure Y is random across processor)
  iseed[1] = 101 * rank;
  iseed[2] = 7;
  iseed[3] = 57;

#ifdef MATRIX_MULT
  zlarnv(&IONE, iseed, &nloct, &Y[0][0]); // generate random Y
#else
  zlarnv(&IONE, iseed, &nloct, Y);
#endif

  printf("----- start zgtracemin -----\n");
  printf("parameters: n=%d, nloc=%d, neig=%d, [a,b]=[%f,%f], [c,d]=[%f,%f], "
         "chebdegree=%d, maxi=%d, tol=%E\n",
         n, nloc, neig, seig, leig, saeig, laeig, chebdegree, maxi, tol);
  printf("some other stuff: seig=%f,leig=%f; gamma=%f; refE=%f\n", seig, leig,
         gamma, refE);

  printf("perform initial map_chebyshev_shift multiply\n");

#ifdef MATRIX_MULT
  zcopy(&nloct, &Y[0][0], &IONE, &Q[0][0], &IONE);
  zcopy(&nloct, &Q[0][0], &IONE, &S[0][0], &IONE);
  Hmatmult_spds_complete(Vn, d, refE, Q, neig);
#else
  zcopy(&nloct, Y, &IONE, Q, &IONE);
  zcopy(&nloct, Q, &IONE, S, &IONE);
  for (k = 0; k < neig; k++) {
    Hmult_spds_complete(&Vn[k * nloc], d, refE, &Q[k * nloc]);
  }
#endif

  for (j = 1; j <= 4 * chebdegree; j++) {

#ifdef MATRIX_MULT
    map_chebyshev_shift_wrapper(&n, &nloc, &neig, &Vn[0][0], &Q[0][0], &R[0][0],
                                &Vp[0][0], &S[0][0], &chebdegree, &gamma, &seig,
                                &leig, &j, &mcsinfo);
#else
    map_chebyshev_shift_wrapper(&n, &nloc, &neig, Vn, Q, R, Vp, S, &chebdegree,
                                &gamma, &seig, &leig, &j, &mcsinfo);
#endif

    if (mcsinfo == 0) {
      // printf("--- map chebyshev shift done ---\n");
      break;
    }
    if ((mcsinfo == 2) || (mcsinfo == 1)) {
#ifdef MATRIX_MULT
      Hmatmult_spds_complete(T, d, refE, Vn, neig);
      zcopy(&nloct, &T[0][0], &IONE, &Vn[0][0], &IONE);
#else
      for (k = 0; k < neig; k++) {
        Hmult_spds_complete(&T[k * nloc], d, refE, &Vn[k * nloc]);
      }
      zcopy(&nloct, T, &IONE, Vn, &IONE);
#endif
    }
    if (mcsinfo < 0) {
      // printf("--- Error in the map chebyshev shift ---");
      info = mcsinfo;
      break;
    }
  }
  printf("start the main for loop\n");
  for (i = 1; i <= 4 * maxi;
       i++) { // every iteration requires two matrix-vector multiplies (the
              // algorithm will stop when needed by itself)

#ifdef MATRIX_MULT
    zgtracemin_wrapper(&n, &neig, &Y[0][0], E, &Q[0][0], &Qt[0][0], H, Hloc, G,
                       ZWork, RWork, IWork, &R[0][0], &tol, &maxi, &i, &nloc,
                       &rank, &np, &comm, &info);
#else
    zgtracemin_wrapper(&n, &neig, Y, E, Q, Qt, H, Hloc, G, ZWork, RWork, IWork,
                       R, &tol, &maxi, &i, &nloc, &rank, &np, &comm, &info);
#endif

    if (info == 0) {
      printf("--- Algorithm Converged --- (iteration=%d)\n", i / 2);
      break;
    }
    if ((info == 1) || (info == 2)) {
      printf("perform map_chebyshev_shift multiply, iteration %d\n", i);

#ifdef MATRIX_MULT
      zcopy(&nloct, &Q[0][0], &IONE, &S[0][0], &IONE);
      Hmatmult_spds_complete(Vn, d, refE, Q, neig);
#else
      zcopy(&nloct, Q, &IONE, S, &IONE);
      for (k = 0; k < neig; k++) {
        Hmult_spds_complete(&Vn[k * nloc], d, refE, &Q[k * nloc]);
      }
#endif

      for (j = 1; j <= 4 * chebdegree; j++) {

#ifdef MATRIX_MULT
        map_chebyshev_shift_wrapper(&n, &nloc, &neig, &Vn[0][0], &Q[0][0],
                                    &R[0][0], &Vp[0][0], &S[0][0], &chebdegree,
                                    &gamma, &seig, &leig, &j, &mcsinfo);
#else
        map_chebyshev_shift_wrapper(&n, &nloc, &neig, Vn, Q, R, Vp, S,
                                    &chebdegree, &gamma, &seig, &leig, &j,
                                    &mcsinfo);
#endif

        if (mcsinfo == 0) {
          // printf("--- map chebyshev shift done ---\n");
          break;
        }
        if ((mcsinfo == 1) || (mcsinfo == 2)) {

#ifdef MATRIX_MULT
          Hmatmult_spds_complete(T, d, refE, Vn, neig);
          zcopy(&nloct, &T[0][0], &IONE, &Vn[0][0], &IONE);
#else
          for (k = 0; k < neig; k++) {
            Hmult_spds_complete(&T[k * nloc], d, refE, &Vn[k * nloc]);
          }
          zcopy(&nloct, T, &IONE, Vn, &IONE);
#endif
        }
        if (mcsinfo < 0) {
          printf("--- Error in the map chebyshev shift ---\n");
          info = mcsinfo;
          break;
        }
      }
    }
    if (info == 3) {
      printf(
          "--- Algorithm reached maximum # of iterations --- (iteration=%d)\n",
          i / 2);
      break;
    }
    if (info < 0) {
      printf(
          "--- Error in the zgtracemin routine --- (iteration=%d; info=%d)\n",
          i / 2, info);
      break;
    }
  }
  maxi = i / 2;

  printf("A orthonormalize Y\n");
  einfo = 1;

#ifdef MATRIX_MULT
  pcextract_eigenpairs_wrapper(&Y[0][0], E, &neig, &nloc, H, Hloc, ZWork, RWork,
                               &T[0][0], &rank, &comm, &einfo);
  Hmatmult_spds_complete(T, d, refE, Y, neig);
  pcextract_eigenpairs_wrapper(&Y[0][0], E, &neig, &nloc, H, Hloc, ZWork, RWork,
                               &T[0][0], &rank, &comm, &einfo);
#else
  pcextract_eigenpairs_wrapper(Y, E, &neig, &nloc, H, Hloc, ZWork, RWork, T,
                               &rank, &comm, &einfo);
  // pcfmv(M,Y,T,neig,nloc,ne,rank,np)
  for (k = 0; k < neig; k++) {
    Hmult_spds_complete(&T[k * nloc], d, refE, &Y[k * nloc]);
  }
  pcextract_eigenpairs_wrapper(Y, E, &neig, &nloc, H, Hloc, ZWork, RWork, T,
                               &rank, &comm, &einfo);
#endif

// recompute eigenvalues through lambda=x'*A*x/x'*x
#ifdef MATRIX_MULT
  Hmatmult_spds_complete(T, d, refE, Y, neig);
  for (i = 0; i < neig; i++) {
    zdotc(&ytloc, &nloc, Y[i], &IONE, T[i], &IONE);
    zdotc(&yyloc, &nloc, Y[i], &IONE, Y[i], &IONE);
    MPI_Allreduce(&ytloc.r, &yt, 1, MPI_DOUBLE_PRECISION, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&yyloc.r, &yy, 1, MPI_DOUBLE_PRECISION, MPI_SUM,
                  MPI_COMM_WORLD);
    E[i] = yt / yy;
  }
#else
  for (k = 0; k < neig; k++) {
    Hmult_spds_complete(&T[k * nloc], d, refE, &Y[k * nloc]);
  }
  for (i = 0; i < neig; i++) {
    zdotc(&ytloc, &nloc, &Y[i * nloc], &IONE, &T[i * nloc], &IONE);
    zdotc(&yyloc, &nloc, &Y[i * nloc], &IONE, &Y[i * nloc], &IONE);
    MPI_Allreduce(&ytloc.r, &yt, 1, MPI_DOUBLE_PRECISION, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&yyloc.r, &yy, 1, MPI_DOUBLE_PRECISION, MPI_SUM,
                  MPI_COMM_WORLD);
    E[i] = yt / yy;
  }
#endif

  printf(
      "print the eigenvalues and the norm of the estimate of the residual\n");
  for (i = 0; i < neig; i++) {

#ifdef MATRIX_MULT
    Rnrm[i] = dznrm2(&nloc, R[i], &IONE);
#else
    Rnrm[i] = dznrm2(&nloc, &R[i * nloc], &IONE);
#endif

    printf("eigenvalue= %E, norm of its residual=%E\n", E[i], Rnrm[i]);
  }

  printf("print the eigenvalues and the norm of the real residual\n");
  for (i = 0; i < neig; i++) {
    yyloc.r = -E[i];
    yyloc.i = 0.0;
    yy = 0.0;

#ifdef MATRIX_MULT
    zaxpy(&nloc, &yyloc, Y[i], &IONE, T[i], &IONE);
    zdotc(&yyloc, &nloc, T[i], &IONE, T[i], &IONE);
#else
    zaxpy(&nloc, &yyloc, &Y[i * nloc], &IONE, &T[i * nloc], &IONE);
    zdotc(&yyloc, &nloc, &T[i * nloc], &IONE, &T[i * nloc], &IONE);
#endif

    MPI_Allreduce(&yyloc.r, &yy, 1, MPI_DOUBLE_PRECISION, MPI_SUM,
                  MPI_COMM_WORLD);
    printf("eigenvalue= %E, norm of its residual=%E\n", E[i], sqrt(yy));
  }

  printf("deallocate temporary variables\n");

#ifdef MATRIX_MULT
  rm_cmatrix(&Q);
  rm_cmatrix(&Qt);
  rm_cmatrix(&S);
  rm_cmatrix(&Vn);
  rm_cmatrix(&Vp);
  rm_cmatrix(&T);
#else
  rm_cvectr(&Q);
  rm_cvectr(&Qt);
  rm_cvectr(&S);
  rm_cvectr(&Vn);
  rm_cvectr(&Vp);
  rm_cvectr(&T);
#endif

  rm_cvectr(&H);
  rm_cvectr(&Hloc);
  rm_cvectr(&G);
  rm_cvectr(&ZWork);
  rm_rvectr(&RWork);
  rm_ivectr(&IWork);

  printf("put results back into qd_struct d\n");
  /*
  //assuming eigenvalues are already sorted
  if(info >=0) {
    //put back the eigenvalues
    d->neigv = neig;
    d->WF = vector<waveFunction>(d->neigv);
    for (i=0; i<d->neigv; i++) {                                        //wf
  points to d->WF[i] (C++ syntax) waveFunction& wf = d->WF[i]; wf.indx = i; wf.E
  = E[i] + refE; wf.valid         = 1;
    }
    writeEigenvalues(d, 'a');                                           //write
  eigenvalues into a file

    //put back the eigenvectors
    if(d->opt.ExecParam.ElCalc.ElAlg.PCTracemin.CalEigVec) {
      d->wfc = Cmatrix(d->neigv, nloc); //d->neigv was ncv, but I think it was a
  mistake and should've been nev in the original version??? for (i=0;
  i<d->neigv; i++){ for (j=0; j<nloc; j++){ #ifdef MATRIX_MULT d->wfc[i][j].r =
  Y[i][j].r; d->wfc[i][j].i = Y[i][j].i; #else d->wfc[i][j].r = Y[i*nloc+j].r;
              d->wfc[i][j].i = Y[i*nloc+j].i;
           #endif
        }
      }

      writePsiSqr(d->wfc,d->neigv,d,0);                                 //not
  sure what this does??? if
  (d->opt.Dev.band_model==Dev_struct::Bands_1_s_nospin){         //in the single
  band case the data file is small enough so that we write out the full
  wavefunction print_psisq_gnu_pdb(d->wfc,d->neigv,d,0);
      }

      for(i=0; i<d->neigv; i++){
        writeEigenvector(d, d->wfc[i], i);                              //write
  eigenvector to a file
      }
    }
  }
  */

  // assuming needs to sort the eigenvalues
  if (info >= 0) {
    // put back sorted eigenvalues
    d->neigv = neig;
    d->WF = vector<waveFunction>(d->neigv);
    orderedE = Rvectr(d->neigv);
    orderedI = Ivectr(d->neigv);
    for (i = 0; i < d->neigv; i++) {
      orderedE[i] =
          E[i] + refE; // assign shifted back eigenvalues to another array
      orderedI[i] = i; // keep track of the original index of the eigenvalue
    }
    for (i = 0; i < d->neigv - 1; i++) { // do bubble sort
      for (j = i + 1; j < d->neigv; j++) {
        if (orderedE[i] > orderedE[j]) {
          dtmp = orderedE[i];
          orderedE[i] = orderedE[j];
          orderedE[j] = dtmp;
          itmp = orderedI[i];
          orderedI[i] = orderedI[j];
          orderedI[j] = itmp;
        }
      }
    }
    for (i = 0; i < d->neigv; i++) { // wf points to d->WF[i] (C++ syntax)
      waveFunction &wf = d->WF[i];
      wf.indx = i;
      wf.E = orderedE[i];
      wf.valid = 1;
      printf("shifted back eig:%E\n", orderedE[i]);
    }
    writeEigenvalues(d, 'a'); // write eigenvalues into a file

    // put back the eigenvectors accounting for previous sort of the eigenvalues
    if (d->opt.ExecParam.ElCalc.ElAlg.PCTracemin.CalEigVec) {
      d->wfc = Cmatrix(
          d->neigv, nloc); // d->neigv was ncv, but I think it was a mistake and
                           // should've been nev in the original version???
      for (i = 0; i < d->neigv; i++) {
        for (j = 0; j < nloc; j++) {
#ifdef MATRIX_MULT
          d->wfc[i][j].r = Y[orderedI[i]][j].r;
          d->wfc[i][j].i = Y[orderedI[i]][j].i;
#else
          d->wfc[i][j].r = Y[orderedI[i] * nloc + j].r;
          d->wfc[i][j].i = Y[orderedI[i] * nloc + j].i;
#endif
        }
      }

      writePsiSqr(d->wfc, d->neigv, d, 0); // not sure what this does???
      if (d->opt.Dev.band_model ==
          Dev_struct::Bands_1_s_nospin) { // in the single band case the data
                                          // file is small enough so that we
                                          // write out the full wavefunction
        print_psisq_gnu_pdb(d->wfc, d->neigv, d, 0);
      }
      for (i = 0; i < d->neigv; i++) {
        writeEigenvector(d, d->wfc[i], i); // write eigenvector to a file
      }
    }
  }

  compute_ritz_residual_norm_deviation(d);
  // compute_ritz_value(d);
  // compute_residual(d);
  // compute_eigenvector_norm(d);

  printf("EXITING PCTRACEMIN DRIVER\n");
  return d->neigv;

#else
#ifdef USE_TRACEMIN_SER
  // if(d->rc==NULL) d->rc = Cvectr(d->nvmap[mpi_n3d_id]*d->NBasisStates);
  // //kept this line for parpack_driver.c, not sure what it's for???
  int neig = d->opt.ExecParam.ElCalc.ElAlg.PCTracemin
                 .NumEigVal; // number of eigenvalues wanted
  real refE = d->opt.ExecParam.ElCalc.ElAlg.PCTracemin
                  .RefEigVal; // shift the matrix A for it to be hermitian
                              // positive definite (A - refE*I)
  real tol =
      d->opt.ExecParam.ElCalc.ElAlg.PCTracemin.Tolerance; // accuracy wanted
  real seig =
      d->opt.ExecParam.ElCalc.ElAlg.PCTracemin
          .WantEigRange[0]; // min. of interval where eigenvalues are wanted
  real leig =
      d->opt.ExecParam.ElCalc.ElAlg.PCTracemin
          .WantEigRange[1]; // max. of interval where eigenvalues are wanted
  real saeig = d->opt.ExecParam.ElCalc.ElAlg.PCTracemin
                   .AllEigRange[0]; // min. all eigenvalues
  real laeig = d->opt.ExecParam.ElCalc.ElAlg.PCTracemin
                   .AllEigRange[1]; // max all eigenvalues
  int maxi = d->opt.ExecParam.ElCalc.ElAlg.PCTracemin.maxi;
  int chebdegree = d->opt.ExecParam.ElCalc.ElAlg.PCTracemin
                       .ChebyDeg;          // degree of chebyshev polynomial
  int n = d->n_ham_tot;                    // dimension of eigenproblem
  int nloc = 0; /*d->seg_ln[mpi_n3d_id];*/ // number of vector elements on
                                           // current processor
  int info; // flag for reverse communication; initialized to 0
  // declare important MPI parameters
  // MPI_Fint comm = MPI_Comm_c2f(MPI_COMM_WORLD); //transform C into Fortran
  // communicator (currently tracemin uses MPI_COMM_WORLD and comm doesn't
  // matter)
  int comm = 0;
  int rank = 0; // rank of the current processor
  int np = 0;   // number of processors
  // allocate arrays containing result of tracemin
  rvectr E = Rvectr(neig); // real eigenvalues
  // allocate necessary temporary arrays
  cvectr H = Cvectr(neig * neig);
  cvectr Hloc = Cvectr(neig * neig);
  cvectr G = Cvectr(neig * neig);
  cvectr ZWork = Cvectr(max(1, 2 * neig - 1));
  rvectr RWork = Rvectr(max(1, 3 * neig - 2));
  ivectr IWork = Ivectr(neig);

#ifdef MATRIX_MULT
  cmatrix R = Cmatrix(neig, n); // estimate of residual
  cmatrix Y = Cmatrix(neig, n); // complex eigenvectors
  cmatrix Q = Cmatrix(neig, n);
  cmatrix Qt = Cmatrix(neig, n);
  cmatrix S = Cmatrix(neig, n);
  cmatrix Vn = Cmatrix(neig, n);
  cmatrix Vp = Cmatrix(neig, n);
  cmatrix T = Cmatrix(neig, n);
#else
  cvectr R = Cvectr(n * neig); // estimate of residual
  cvectr Y = Cvectr(n * neig); // complex eigenvectors
  cvectr Q = Cvectr(n * neig);
  cvectr Qt = Cvectr(n * neig);
  cvectr S = Cvectr(n * neig);
  cvectr Vn = Cvectr(n * neig);
  cvectr Vp = Cvectr(n * neig);
  cvectr T = Cvectr(n * neig);
#endif

  int i, j, k, mcsinfo, einfo;
  real gamma; // parameter used by tracemin
  // real seig; //smallest eigenvalue you are looking for real leig; //largest
  // eigenvalue you are looking for real saeig; //smallest of alleigenvalues
  // real laeig; //largest of all eigenvalues
  int nloct; // total number of elements on this processor: nloc*neig
  // int chebdegree; //degree of chebyshev polynomial int maxi; //maximum number
  // of iterations of tracemin (there are neig*2*6*2 matvecs per iteration)
  complex ytloc; // temporary variables used in computation of eigenvalues after
                 // tracemin (next 4 lines)
  complex yyloc;
  real yt;
  real yy;
  rvectr Rnrm = Rvectr(neig); // norm of the residual of of individual
                              // eigenpairs
  // variables used in sorting of the eigenvalues after tracemin
  rvectr orderedE;
  ivectr orderedI;
  double dtmp;
  int itmp;
  // declare some constants needed in the code
  int IONE = 1;
  real DTWO = 2.0;
  int iseed[4];

  printf("ENTERING SCTRACEMIN DRIVER\n");
  printf("setup some variables\n");
  // MPI_Comm_rank(MPI_COMM_WORLD, &rank); //find out my rank
  // MPI_Comm_size(MPI_COMM_WORLD, &np); //find out number of processors
  ////if(SMALLER_THEN_ALL_EIG <= 0){
  ////  refE = (real)(SMALLER_THEN_ALL_EIG - 1);
  ////}
  ////else{
  refE = 0.0;
  ////}
  // seig = emin1;                ////emin1 - refE; ///-30-refE; //emin1 -
  // refE;//???                   //define an interval where to look for the
  // eigenvalues leig = emax2;                ////emax2 - refE; ///30-refE;
  // //emax2 - refE;//???                   //THIS PARAMETER AND ABOVE ONE CAN
  // BE AN OPTION IN THE INPUT FILE!!! saeig= SMALLER_THEN_ALL_EIG;
  // ////SMALLER_THEN_ALL_EIG - refE; //define an interval in which all
  // eigenvalues are contained laeig= LARGER_THEN_ALL_EIG;
  // ////LARGER_THEN_ALL_EIG  - refE; //THIS PARAMETER AND ABOVE ONE CAN BE AN
  // OPTION IN THE INPUT FILE!!!
  gamma = min(DTWO / ((saeig - seig) * (saeig - leig)),
              DTWO / ((laeig - seig) * (laeig - leig)));
  nloct =
      neig * n; // number of elements of all vectors combined on this processor
  // maxi       = 40;  //50;   (comment by marta) //THIS PARAMETER CAN BE AN
  // OPTION IN THE INPUT FILE!!! chebdegree = 240; //6; //THIS PARAMETER CAN BE
  // AN OPTION IN THE INPUT FILE!!!
  iseed[0] = 103 + rank; // set up seed for the random number generator (making
                         // sure Y is random across processor)
  iseed[1] = 101 * rank;
  iseed[2] = 7;
  iseed[3] = 57;
#ifdef MATRIX_MULT
  zlarnv(&IONE, iseed, &nloct, &Y[0][0]); // generate random Y
#else
  zlarnv(&IONE, iseed, &nloct, Y); // generate random Y
#endif

  printf("----- start zgtracemin -----\n");
  printf("parameters: n=%d, nloc=%d, neig=%d, [a,b]=[%f,%f], [c,d]=[%f,%f], "
         "chebdegree=%d, maxi=%d, tol=%E\n",
         n, nloc, neig, seig, leig, saeig, laeig, chebdegree, maxi, tol);
  printf("some other stuff: seig=%f,leig=%f; gamma=%f; refE=%f\n", seig, leig,
         gamma, refE);

  printf("perform initial map_chebyshev_shift multiply\n");

#ifdef MATRIX_MULT
  zcopy(&nloct, &Y[0][0], &IONE, &Q[0][0], &IONE);
  zcopy(&nloct, &Q[0][0], &IONE, &S[0][0], &IONE);
  Hmatmult_spds_complete(Vn, d, refE, Q, neig);
#else
  zcopy(&nloct, Y, &IONE, Q, &IONE);
  zcopy(&nloct, Q, &IONE, S, &IONE);
  for (k = 0; k < neig; k++) {
    Hmult_spds_complete(&Vn[k * n], d, refE, &Q[k * n]);
  }
#endif

  for (j = 1; j <= 4 * chebdegree; j++) {
#ifdef MATRIX_MULT
    map_chebyshev_shift_wrapper(&n, &nloc, &neig, &Vn[0][0], &Q[0][0], &R[0][0],
                                &Vp[0][0], &S[0][0], &chebdegree, &gamma, &seig,
                                &leig, &j, &mcsinfo);
#else
    map_chebyshev_shift_wrapper(&n, &nloc, &neig, Vn, Q, R, Vp, S, &chebdegree,
                                &gamma, &seig, &leig, &j, &mcsinfo);
#endif

    if (mcsinfo == 0) {
      // printf("--- map chebyshev shift done ---\n");
      break;
    }
    if ((mcsinfo == 2) || (mcsinfo == 1)) {
#ifdef MATRIX_MULT
      Hmatmult_spds_complete(T, d, refE, Vn, neig);
      zcopy(&nloct, &T[0][0], &IONE, &Vn[0][0], &IONE);
#else
      for (k = 0; k < neig; k++) {
        Hmult_spds_complete(&T[k * n], d, refE, &Vn[k * n]);
      }
      zcopy(&nloct, T, &IONE, Vn, &IONE);
#endif
    }
    if (mcsinfo < 0) {
      // printf("--- Error in the map chebyshev shift ---");
      info = mcsinfo;
      break;
    }
  }
  printf("start the main for loop\n");
  for (i = 1; i <= 4 * maxi;
       i++) { // every iteration requires two matrix-vector multiplies (the
              // algorithm will stop when needed by itself)

#ifdef MATRIX_MULT
    zgtracemin_wrapper(&n, &neig, &Y[0][0], E, &Q[0][0], &Qt[0][0], H, Hloc, G,
                       ZWork, RWork, IWork, &R[0][0], &tol, &maxi, &i, &nloc,
                       &rank, &np, &comm, &info);
#else
    zgtracemin_wrapper(&n, &neig, Y, E, Q, Qt, H, Hloc, G, ZWork, RWork, IWork,
                       R, &tol, &maxi, &i, &nloc, &rank, &np, &comm, &info);
#endif

    if (info == 0) {
      printf("--- Algorithm Converged --- (iteration=%d)\n", i / 2);
      break;
    }
    if ((info == 1) || (info == 2)) {
      printf("perform map_chebyshev_shift multiply, iteration %d\n", i);

#ifdef MATRIX_MULT
      zcopy(&nloct, &Q[0][0], &IONE, &S[0][0], &IONE);
      Hmatmult_spds_complete(Vn, d, refE, Q, neig);
#else
      zcopy(&nloct, Q, &IONE, S, &IONE);
      for (k = 0; k < neig; k++) {
        Hmult_spds_complete(&Vn[k * n], d, refE, &Q[k * n]);
      }
#endif

      for (j = 1; j <= 4 * chebdegree; j++) {

#ifdef MATRIX_MULT
        map_chebyshev_shift_wrapper(&n, &nloc, &neig, &Vn[0][0], &Q[0][0],
                                    &R[0][0], &Vp[0][0], &S[0][0], &chebdegree,
                                    &gamma, &seig, &leig, &j, &mcsinfo);
#else
        map_chebyshev_shift_wrapper(&n, &nloc, &neig, Vn, Q, R, Vp, S,
                                    &chebdegree, &gamma, &seig, &leig, &j,
                                    &mcsinfo);
#endif

        if (mcsinfo == 0) {
          // printf("--- map chebyshev shift done ---\n");
          break;
        }
        if ((mcsinfo == 1) || (mcsinfo == 2)) {
#ifdef MATRIX_MULT
          Hmatmult_spds_complete(T, d, refE, Vn, neig);
          zcopy(&nloct, &T[0][0], &IONE, &Vn[0][0], &IONE);
#else
          for (k = 0; k < neig; k++) {
            Hmult_spds_complete(&T[k * n], d, refE, &Vn[k * n]);
          }
          zcopy(&nloct, T, &IONE, Vn, &IONE);
#endif
        }
        if (mcsinfo < 0) {
          printf("--- Error in the map chebyshev shift ---\n");
          info = mcsinfo;
          break;
        }
      }
    }
    if (info == 3) {
      printf(
          "--- Algorithm reached maximum # of iterations --- (iteration=%d)\n",
          i / 2);
      break;
    }
    if (info < 0) {
      printf(
          "--- Error in the zgtracemin routine --- (iteration=%d; info=%d)\n",
          i / 2, info);
      break;
    }
  }
  maxi = i / 2;

  printf("A orthonormalize Y\n");
  einfo = 1;

#ifdef MATRIX_MULT
  pcextract_eigenpairs_wrapper(&Y[0][0], E, &neig, &n, H, Hloc, ZWork, RWork,
                               &T[0][0], &rank, &comm, &einfo);
  Hmatmult_spds_complete(T, d, refE, Y, neig);
  pcextract_eigenpairs_wrapper(&Y[0][0], E, &neig, &n, H, Hloc, ZWork, RWork,
                               &T[0][0], &rank, &comm, &einfo);
#else
  pcextract_eigenpairs_wrapper(Y, E, &neig, &n, H, Hloc, ZWork, RWork, T, &rank,
                               &comm, &einfo);
  // pcfmv(M,Y,T,neig,nloc,ne,rank,np)
  for (k = 0; k < neig; k++) {
    Hmult_spds_complete(&T[k * n], d, refE, &Y[k * n]);
  }
  pcextract_eigenpairs_wrapper(Y, E, &neig, &n, H, Hloc, ZWork, RWork, T, &rank,
                               &comm, &einfo);
#endif

// recompute eigenvalues through lambda=x'*A*x/x'*x
#ifdef MATRIX_MULT
  Hmatmult_spds_complete(T, d, refE, Y, neig);
  for (i = 0; i < neig; i++) {
    zdotc(&ytloc, &n, Y[i], &IONE, T[i], &IONE);
    zdotc(&yyloc, &n, Y[i], &IONE, Y[i], &IONE);
    yt = ytloc.r;
    yy = yyloc.r;
    E[i] = yt / yy;
  }
#else
  for (k = 0; k < neig; k++) {
    Hmult_spds_complete(&T[k * n], d, refE, &Y[k * n]);
  }
  for (i = 0; i < neig; i++) {
    zdotc(&ytloc, &n, &Y[i * n], &IONE, &T[i * n], &IONE);
    zdotc(&yyloc, &n, &Y[i * n], &IONE, &Y[i * n], &IONE);
    // MPI_Allreduce(&ytloc.r, &yt, 1, MPI_DOUBLE_PRECISION, MPI_SUM,
    // MPI_COMM_WORLD); MPI_Allreduce(&yyloc.r, &yy, 1, MPI_DOUBLE_PRECISION,
    // MPI_SUM, MPI_COMM_WORLD);
    yt = ytloc.r;
    yy = yyloc.r;
    E[i] = yt / yy;
  }
#endif

  printf(
      "print the eigenvalues and the norm of the estimate of the residual\n");
  for (i = 0; i < neig; i++) {

#ifdef MATRIX_MULT
    Rnrm[i] = dznrm2(&n, R[i], &IONE);
#else
    Rnrm[i] = dznrm2(&n, &R[i * n], &IONE);
#endif

    printf("eigenvalue= %E, norm of its residual=%E\n", E[i], Rnrm[i]);
  }

  printf("print the eigenvalues and the norm of the real residual\n");
  for (i = 0; i < neig; i++) {
    yyloc.r = -E[i];
    yyloc.i = 0.0;
    yy = 0.0;

#ifdef MATRIX_MULT
    zaxpy(&n, &yyloc, Y[i], &IONE, T[i], &IONE);
    zdotc(&yyloc, &n, T[i], &IONE, T[i], &IONE);
#else
    zaxpy(&n, &yyloc, &Y[i * n], &IONE, &T[i * n], &IONE);
    zdotc(&yyloc, &n, &T[i * n], &IONE, &T[i * n], &IONE);
#endif

    // MPI_Allreduce(&yyloc.r, &yy, 1, MPI_DOUBLE_PRECISION, MPI_SUM,
    // MPI_COMM_WORLD);
    yy = yyloc.r;
    printf("eigenvalue= %E, norm of its residual=%E\n", E[i], sqrt(yy));
  }

  printf("deallocate temporary variables\n");

#ifdef MATRIX_MULT
  rm_cmatrix(&Q);
  rm_cmatrix(&Qt);
  rm_cmatrix(&S);
  rm_cmatrix(&Vn);
  rm_cmatrix(&Vp);
  rm_cmatrix(&T);
#else
  rm_cvectr(&Q);
  rm_cvectr(&Qt);
  rm_cvectr(&S);
  rm_cvectr(&Vn);
  rm_cvectr(&Vp);
  rm_cvectr(&T);
#endif

  rm_cvectr(&H);
  rm_cvectr(&Hloc);
  rm_cvectr(&G);
  rm_cvectr(&ZWork);
  rm_rvectr(&RWork);
  rm_ivectr(&IWork);

  printf("put results back into qd_struct d\n");
  /*
  //assuming eigenvalues are already sorted
  if(info >=0) {
    //put back the eigenvalues
    d->neigv = neig;
    d->WF = vector<waveFunction>(d->neigv);
    for (i=0; i<d->neigv; i++) {                                        //wf
  points to d->WF[i] (C++ syntax) waveFunction& wf = d->WF[i]; wf.indx = i; wf.E
  = E[i] + refE; wf.valid         = 1;
    }
    writeEigenvalues(d, 'a');                                           //write
  eigenvalues into a file

    //put back the eigenvectors
    if(d->opt.ExecParam.ElCalc.ElAlg.PCTracemin.CalEigVec) {
      d->wfc = Cmatrix(d->neigv, nloc); //d->neigv was ncv, but I think it was a
  mistake and should've been nev in the original version??? for (i=0;
  i<d->neigv; i++){ for (j=0; j<n; j++){ #ifdef MATRIX_MULT d->wfc[i][j].r =
  Y[i][j].r; d->wfc[i][j].i = Y[i][j].i; #else d->wfc[i][j].r = Y[i*n+j].r;
              d->wfc[i][j].i = Y[i*n+j].i;
           #endif
        }
      }

      writePsiSqr(d->wfc,d->neigv,d,0);                                 //not
  sure what this does??? if
  (d->opt.Dev.band_model==Dev_struct::Bands_1_s_nospin){         //in the single
  band case the data file is small enough so that we write out the full
  wavefunction print_psisq_gnu_pdb(d->wfc,d->neigv,d,0);
      }

      for(i=0; i<d->neigv; i++){
        writeEigenvector(d, d->wfc[i], i);                              //write
  eigenvector to a file
      }
    }
  }
  */

  // assuming needs to sort the eigenvalues
  if (info >= 0) {
    // put back sorted eigenvalues
    d->neigv = neig;
    d->WF = vector<waveFunction>(d->neigv);
    orderedE = Rvectr(d->neigv);
    orderedI = Ivectr(d->neigv);
    for (i = 0; i < d->neigv; i++) {
      orderedE[i] =
          E[i] + refE; // assign shifted back eigenvalues to another array
      orderedI[i] = i; // keep track of the original index of the eigenvalue
    }
    for (i = 0; i < d->neigv - 1; i++) { // do bubble sort
      for (j = i + 1; j < d->neigv; j++) {
        if (orderedE[i] > orderedE[j]) {
          dtmp = orderedE[i];
          orderedE[i] = orderedE[j];
          orderedE[j] = dtmp;
          itmp = orderedI[i];
          orderedI[i] = orderedI[j];
          orderedI[j] = itmp;
        }
      }
    }
    for (i = 0; i < d->neigv; i++) { // wf points to d->WF[i] (C++ syntax)
      waveFunction &wf = d->WF[i];
      wf.indx = i;
      wf.E = orderedE[i];
      wf.valid = 1;
      printf("shifted back eig:%E\n", orderedE[i]);
    }
    writeEigenvalues(d, 'a'); // write eigenvalues into a file

    // put back the eigenvectors accounting for previous sort of the eigenvalues
    if (d->opt.ExecParam.ElCalc.ElAlg.PCTracemin.CalEigVec) {
      d->wfc = Cmatrix(d->neigv,
                       n); // d->neigv was ncv, but I think it was a mistake and
                           // should've been nev in the original version???
      for (i = 0; i < d->neigv; i++) {
        for (j = 0; j < n; j++) {
#ifdef MATRIX_MULT
          d->wfc[i][j].r = Y[orderedI[i]][j].r;
          d->wfc[i][j].i = Y[orderedI[i]][j].i;
#else
          d->wfc[i][j].r = Y[orderedI[i] * n + j].r;
          d->wfc[i][j].i = Y[orderedI[i] * n + j].i;
#endif
        }
      }

      writePsiSqr(d->wfc, d->neigv, d, 0); // not sure what this does???
      if (d->opt.Dev.band_model ==
          Dev_struct::Bands_1_s_nospin) { // in the single band case the data
                                          // file is small enough so that we
                                          // write out the full wavefunction
        print_psisq_gnu_pdb(d->wfc, d->neigv, d, 0);
      }
      for (i = 0; i < d->neigv; i++) {
        writeEigenvector(d, d->wfc[i], i); // write eigenvector to a file
      }
    }
  }

  compute_ritz_residual_norm_deviation(d);
  // compute_ritz_value(d);
  // compute_residual(d);
  // compute_eigenvector_norm(d);

  printf("EXITING SCTRACEMIN DRIVER\n");
  return d->neigv;
#endif
  return 0;
#endif // USE_TRACEMIN_SER
}
