#include "eigsys3d_par.h"
#include "qtracemin_driver.h"

//#define MATRIX_MULT
#ifdef MATRIX_MULT
   #include "h_cmat_mult.h"
#endif

int qtracemin_driver(qd_struct d){

#ifdef USE_TRACEMIN_PAR
  if(d->rc==NULL) d->rc = Cvectr(d->nvmap[mpi_n3d_id]*d->NBasisStates); //kept this line from parpack_driver.c, not sure what it's for???
  int  neig  = d->opt.ExecParam.ElCalc.ElAlg.PQTracemin.NumEigVal;         //number of eigenvalues wanted                                                          
  real refE  = d->opt.ExecParam.ElCalc.ElAlg.PQTracemin.RefEigVal;         //shift the matrix A so that you are finding eigenvalues closest to the shift
  real tol   = d->opt.ExecParam.ElCalc.ElAlg.PQTracemin.Tolerance;         //accuracy wanted 
  real cgtol  = d->opt.ExecParam.ElCalc.ElAlg.PQTracemin.cgtol;         //initial tol of GC
  int maxi = d->opt.ExecParam.ElCalc.ElAlg.PQTracemin.maxi;  		//this was set to 40, in line 48 
  int cgmaxi  = d->opt.ExecParam.ElCalc.ElAlg.PQTracemin.cgmaxim;         // Number of iterations of CG 
  int n      = d->n_ham_tot;                                            //dimension of eigenproblem                          
  int nloc   = d->seg_ln[mpi_n3d_id];                                   //number of vector elements on current processor                                         
  int info;                                                             //flag for reverse communication; initialized to 0                                            
  //declare important MPI parameters
  MPI_Fint comm = MPI_Comm_c2f(MPI_COMM_WORLD);                         //transform C into Fortran communicator (currently tracemin uses MPI_COMM_WORLD and comm doesn't matter)
  int rank;                                                             //rank of the current processor               
  int np;                                                               //number of processors 
  //allocate arrays containing result of tracemin
  rvectr E = Rvectr(neig);                                              //real eigenvalues                                                                   

  //allocate necessary temporary arrays  
  #ifdef MATRIX_MULT
     cmatrix Y    = Cmatrix(neig,nloc);                                         //complex eigenvectors 
     cmatrix Q    = Cmatrix(neig,nloc);                                             
     cmatrix T    = Cmatrix(neig,nloc);
     cmatrix Qcg  = Cmatrix(neig,nloc);
     cmatrix P    = Cmatrix(neig,nloc);
     cmatrix X    = Cmatrix(neig,nloc);
     cmatrix Rcg  = Cmatrix(neig,nloc);
     cmatrix R    = Cmatrix(neig,nloc);                                         //estimate of residual                                                           
  #else
     cvectr Y = Cvectr(nloc*neig);                                         //complex eigenvectors 
     cvectr Q    = Cvectr(nloc*neig);                                             
     cvectr T    = Cvectr(nloc*neig);
     cvectr Qcg  = Cvectr(nloc*neig);
     cvectr P    = Cvectr(nloc*neig);
     cvectr X    = Cvectr(nloc*neig);
     cvectr Rcg  = Cvectr(nloc*neig);
     cvectr R = Cvectr(nloc*neig);                                     //estimate of residual                                  
  #endif

  cvectr H    = Cvectr(neig*neig);  
  cvectr Hloc = Cvectr(neig*neig);
  cvectr ZWork= Cvectr(max(1,2*neig-1)); 
  rvectr RWork= Rvectr(max(1,3*neig-2));
  ivectr IWork= Ivectr(neig);
  cvectr G    = Cvectr(neig*neig);
  cvectr Gloc = Cvectr(neig*neig);
  cvectr gamma_old = Cvectr(neig);
  rvectr rcgtol    = Rvectr(neig);
//commented by marta: defs of cgmaxi and cgtol to be defined inpt deck; 
  int  i,j,k,cginfo;//,cgmaxi;
//real cgtol;                                                           //parameter used by tracemin
  int nloct;                                                            //total number of elements on this processor: nloc*neig
  //int maxi;                                                             //maximum number of iterations of tracemin (there are neig*2*6*2 matvecs per iteration)
  complex ytloc;                                                        //temporary variables used in computation of eigenvalues after tracemin (next 4 lines)
  complex yyloc;                
  real yt;                      
  real yy;                      
  rvectr Rnrm = Rvectr(neig);                                           //norm of the residual of of individual eigenpairs
  //variables used in sorting of the eigenvalues after tracemin
  rvectr orderedE;
  ivectr orderedI;
  double dtmp;
  int    itmp;
  //declare some constants needed in the code
  int  IONE = 1;
  real DONE = 1.0;
  real DTWO = 2.0;  
  int iseed[4];

  printf("ENTERING PQTRACEMIN DRIVER\n");
  printf("Setup some variables\n");
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);                                 //find out my rank
  MPI_Comm_size(MPI_COMM_WORLD, &np);                                   //find out number of processors  
  nloct   =neig*nloc;                                                   //number of elements of all vectors combined on this processor 
  //maxi    =40;  (commented by marta)                                           //maximum number of tracemin iterations (SHOULD BE AN OPTION IN THE INPUT FILE) 
  //cgtol   = DONE;    (commented by marta)                                                   //initial cg tolerance is set to 1.0    (SHOULD BE AN OPTION IN THE INPUT FILE)
  iseed[0]=103+rank;                                                    //set up seed for the random number generator (making sure Y is random across processor)
  iseed[1]=101*rank;                                                      
  iseed[2]=7;
  iseed[3]=57;

  #ifdef MATRIX_MULT
     zlarnv(&IONE,iseed,&nloct,&Y[0][0]);                                         //generate random Y 
  #else
     zlarnv(&IONE,iseed,&nloct,Y);                                         //generate random Y 
  #endif

  printf("%d,----- start ztracemin -----\n",rank);
  printf("%d,parameters: n=%d, nloc=%d, neig=%d, maxi=%d, tol=%E\n",rank,n,nloc,neig,maxi,tol);
  printf("%d, refE=%f\n",rank,refE);

  printf("%d,perform initial two matrix-vector multiplies\n",rank);
  //call pcfmv(M,Y,T,neig,nloc,ne,rank,np)
  //call pcfmv(M,T,Q,neig,nloc,ne,rank,np)

  #ifdef MATRIX_MULT
     Hmatmult_spds_complete(T, d, refE, Y, neig);
	 Hmatmult_spds_complete(Q, d, refE, T, neig);
  #else
     for (k=0; k<neig; k++){
        Hmult_spds_complete(&T[k*nloc], d, refE, &Y[k*nloc]);
     }
     for (k=0; k<neig; k++){
        Hmult_spds_complete(&Q[k*nloc], d, refE, &T[k*nloc]);
     }
  #endif

  printf("%d,start the main for loop\n",rank);
  for(i=1; i<=2*maxi+1; i++){ //every iteration requires one matrix vector multiply and one system solve with matrix M=A^2 
    #ifdef MATRIX_MULT
       ztracemin_wrapper(&n,&neig,&Y[0][0],E,&Q[0][0],H,Hloc,ZWork,RWork,IWork,&R[0][0],&tol,&maxi,&i,&nloc,&rank,&np,&comm,&info);
	#else
       ztracemin_wrapper(&n,&neig,Y,E,Q,H,Hloc,ZWork,RWork,IWork,R,&tol,&maxi,&i,&nloc,&rank,&np,&comm,&info);
	#endif
    if(info == 1){
      printf("%d,--- ztracemin has converged --- (iteration=%d; info=%d)\n",rank,i/2,info);
      break;
    }
    if (info == 2){
      printf("%d,Perform a solve\n",rank);
      //setup CG parameters
      //call pcfmv(M,Y,T,  neig,nloc,ne,rank,np)
      //call pcfmv(M,T,Qcg,neig,nloc,ne,rank,np)

	  #ifdef MATRIX_MULT
	     Hmatmult_spds_complete(T, d, refE, Y, neig);
	     Hmatmult_spds_complete(Qcg, d, refE, T, neig);
	  #else
         for (k=0; k<neig; k++){
	        Hmult_spds_complete(&T[k*nloc], d, refE, &Y[k*nloc]);
         }
         for (k=0; k<neig; k++){
	        Hmult_spds_complete(&Qcg[k*nloc], d, refE, &T[k*nloc]);
         }
	  #endif

      //cgmaxi=min(n,1000);
      if (i%10 == 0){
	cgtol = cgtol/10;
	if (cgtol < tol){
	  cgtol=tol;
	}
      }
      printf("CG: tolerance=%f,maxi=%d\n",cgtol,cgmaxi);

      //start CG
      for (j=1; j<=cgmaxi+2; j++){ //notice that one iteration is for setup
	  #ifdef MATRIX_MULT
	     modified_conjugate_gradient_wrapper(&n,&neig,&X[0][0],&Y[0][0],&Qcg[0][0],G,Gloc,&P[0][0],&Rcg[0][0],gamma_old,rcgtol,&j,&cgmaxi,&cgtol,&tol,&nloc,&rank,&np,&comm,&cginfo);
	  #else
	     modified_conjugate_gradient_wrapper(&n,&neig,X,Y,Qcg,G,Gloc,P,Rcg,gamma_old,rcgtol,&j,&cgmaxi,&cgtol,&tol,&nloc,&rank,&np,&comm,&cginfo);
	  #endif
	if (cginfo == 1){
	  printf("%d,modified CG converged at iteration %d\n",rank,j);
	  break;
	}
	if (cginfo == 2){
	  //printf("%d,Perform two matrix-vector multiplies (CG)\n",rank);
	  //call pcfmv(M,P,T,  neig,nloc,ne,rank,np)
	  //call pcfmv(M,T,Qcg,neig,nloc,ne,rank,np)
	  #ifdef MATRIX_MULT
	     Hmatmult_spds_complete(T, d, refE, P,neig);
	     Hmatmult_spds_complete(Qcg, d, refE, T,neig);
	  #else
	     for (k=0; k<neig; k++){
	        Hmult_spds_complete(&T[k*nloc], d, refE, &P[k*nloc]);
	     }
	     for (k=0; k<neig; k++){
	        Hmult_spds_complete(&Qcg[k*nloc], d, refE, &T[k*nloc]);
	     }
	  #endif
	}
	if (cginfo == 3){
	  printf("%d,modified CG has reached maximum number of iterations at iteration %d\n",rank,j);
	  break;
	}
	if (cginfo < 0){
	  printf("%d,ERROR in modified CG, at iteration %d\n",rank,j);
	  break;
	}
      }
      //copy the solution X=Y-D back into Y and compute Q=A*Y
	  #ifdef MATRIX_MULT
         zcopy(&nloct,&X[0][0],&IONE,&Y[0][0],&IONE);
	  #else
         zcopy(&nloct,X,&IONE,Y,&IONE);
	  #endif

      //call pcfmv(M,Y,T,neig,nloc,ne,rank,np)
      //call pcfmv(M,T,Q,neig,nloc,ne,rank,np)

	  #ifdef MATRIX_MULT
	     Hmatmult_spds_complete(T, d, refE, Y,neig);
	     Hmatmult_spds_complete(Q, d, refE, T,neig);
	  #else
         for (k=0; k<neig; k++){
	        Hmult_spds_complete(&T[k*nloc], d, refE, &Y[k*nloc]);
         }
         for (k=0; k<neig; k++){
	        Hmult_spds_complete(&Q[k*nloc], d, refE, &T[k*nloc]);
         }
	  #endif

    }
    if (info == 3){
      printf("%d,Perform two matrix-vector multiplies\n",rank);
      //call pcfmv(M,Y,T,neig,nloc,ne,rank,np)
      //call pcfmv(M,T,Q,neig,nloc,ne,rank,np)
	  #ifdef MATRIX_MULT
	     Hmatmult_spds_complete(T, d, refE, Y,neig);
	     Hmatmult_spds_complete(Q, d, refE, T,neig);
	  #else
         for (k=0; k<neig; k++){
	        Hmult_spds_complete(&T[k*nloc], d, refE, &Y[k*nloc]);
         }
         for (k=0; k<neig; k++){
	        Hmult_spds_complete(&Q[k*nloc], d, refE, &T[k*nloc]);
         }
	  #endif
    }
    if (info == 4){
      printf("%d,--- Maximum number of iterations reached ---(iteration=%d; info=%d)\n",rank,i/2,info);
      break;
    }
    if (info < 0){
      printf("%d,--- Error in the ztracemin routine --- (iteration=%d; info=%d)\n",rank,i/2,info);
      break; 
    }
    printf("%d,done with do loop iteration -->%d\n",rank,i);
  }
  maxi=i/2; 

  printf("%d,find the eigenvalues lambda using lambda=xt*A*x/xt*x\n",rank);
  //call pcfmv(M,Y,T,neig,nloc,ne,rank,np)
  #ifdef MATRIX_MULT
     Hmatmult_spds_complete(T, d, refE, Y,neig);
  #else
     for (k=0; k<neig; k++){
        Hmult_spds_complete(&T[k*nloc], d, refE, &Y[k*nloc]);
     }
  #endif

  for (i=0; i<neig; i++){
    #ifdef MATRIX_MULT
       zdotc(&ytloc,&nloc,Y[i],&IONE,T[i],&IONE);   
       zdotc(&yyloc,&nloc,Y[i],&IONE,Y[i],&IONE);   
	#else
       zdotc(&ytloc,&nloc,&Y[i*nloc],&IONE,&T[i*nloc],&IONE);   
       zdotc(&yyloc,&nloc,&Y[i*nloc],&IONE,&Y[i*nloc],&IONE);   
	#endif

    MPI_Allreduce(&ytloc.r, &yt, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&yyloc.r, &yy, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
    E[i]=yt/yy;
  }

  printf("%d,print the eigenvalues and the norm of the estimate of the residual\n",rank);
  for (i=0; i<neig; i++){
    #ifdef MATRIX_MULT
       Rnrm[i]=dznrm2(&nloc, R[i], &IONE);
	#else
       Rnrm[i]=dznrm2(&nloc, &R[i*nloc], &IONE);
	#endif
    printf("eigenvalue= %E, norm of its residual=%E\n",E[i],Rnrm[i]);
  }  
  printf("%d,print the eigenvalues and the norm of the real residual\n",rank);
  for (i=0; i<neig; i++){
    yyloc.r =-E[i];
    yyloc.i = 0.0;
    yy      = 0.0;

	#ifdef MATRIX_MULT
       zaxpy(&nloc,&yyloc,Y[i],&IONE,T[i],&IONE);
       zdotc(&yyloc,&nloc,T[i],&IONE,T[i],&IONE);
	#else
       zaxpy(&nloc,&yyloc,&Y[i*nloc],&IONE,&T[i*nloc],&IONE);
       zdotc(&yyloc,&nloc,&T[i*nloc],&IONE,&T[i*nloc],&IONE);
	#endif

    MPI_Allreduce(&yyloc.r, &yy, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
    printf("eigenvalue= %E, norm of its residual=%E\n",E[i],sqrt(yy));
  }
    
  printf("%d,deallocate temporary variables\n",rank);
  rm_cvectr(&H);
  rm_cvectr(&Hloc);
  rm_cvectr(&ZWork);
  rm_rvectr(&RWork);
  rm_ivectr(&IWork);
  rm_cvectr(&G); 
  rm_cvectr(&Gloc);
  rm_cvectr(&gamma_old);
  rm_rvectr(&rcgtol);  

  #ifdef MATRIX_MULT
     rm_cmatrix(&Rcg);
     rm_cmatrix(&X);
     rm_cmatrix(&P); 
     rm_cmatrix(&Q);
     rm_cmatrix(&Qcg);
     rm_cmatrix(&T); 
  #else
     rm_cvectr(&Rcg);
     rm_cvectr(&X);
     rm_cvectr(&P); 
     rm_cvectr(&Q);
     rm_cvectr(&Qcg);
     rm_cvectr(&T); 
  #endif

  printf("put results back into qd_struct d\n");
  /*
  //assuming eigenvalues are already sorted  
  if(info >=0) {                                                   
    //put back the eigenvalues
    d->neigv = neig;
    d->WF = vector<waveFunction>(d->neigv);
    for (i=0; i<d->neigv; i++) {                                        //wf points to d->WF[i] (C++ syntax)
      waveFunction& wf = d->WF[i];
      wf.indx          = i;
      wf.E             = E[i] + refE; 
      wf.valid         = 1;
    }
    writeEigenvalues(d, 'a');                                           //write eigenvalues into a file
 
    //put back the eigenvectors
    if(d->opt.ExecParam.ElCalc.ElAlg.PQTracemin.CalEigVec) { 
      d->wfc = Cmatrix(d->neigv, nloc);                                 //d->neigv was ncv, but I think it was a mistake and should've been nev in the original version???
      for (i=0; i<d->neigv; i++){
	for (j=0; j<nloc; j++){
	  #ifdef MATRIX_MULT
	     d->wfc[i][j].r = Y[i][j].r;
	     d->wfc[i][j].i = Y[i][j].i;
	  #else
	     d->wfc[i][j].r = Y[i*nloc+j].r;
	     d->wfc[i][j].i = Y[i*nloc+j].i;
	  #endif
	}
      }  

      writePsiSqr(d->wfc,d->neigv,d,0);                                 //not sure what this does???      
      if (d->opt.Dev.band_model==Dev_struct::Bands_1_s_nospin){         //in the single band case the data file is small enough so that we write out the full wavefunction 
	print_psisq_gnu_pdb(d->wfc,d->neigv,d,0);                  
      }

      for(i=0; i<d->neigv; i++){
	writeEigenvector(d, d->wfc[i], i);                              //write eigenvector to a file  
      }
    }
  }
  */
  
  //assuming needs to sort the eigenvalues
  if(info>=0) {
    //put back sorted eigenvalues
    d->neigv = neig;
    d->WF    = vector<waveFunction>(d->neigv);
    orderedE = Rvectr(d->neigv);
    orderedI = Ivectr(d->neigv);
    for(i=0; i<d->neigv; i++) {                     
      orderedE[i]=E[i] + refE;                                          //assign shifted back eigenvalues to another array 
      orderedI[i]=i;                                                    //keep track of the original index of the eigenvalue
    }    
    for(i=0; i<d->neigv-1;i++){                                         //do bubble sort
      for(j=i+1; j<d->neigv; j++){
	if(orderedE[i]>orderedE[j]) {
	  dtmp       =orderedE[i];
	  orderedE[i]=orderedE[j];
	  orderedE[j]=dtmp;
	  itmp       = orderedI[i];
	  orderedI[i]=orderedI[j];
	  orderedI[j]=itmp;
	}
      }
    }
    for (i=0; i<d->neigv; i++) {                                        //wf points to d->WF[i] (C++ syntax)
      waveFunction& wf = d->WF[i];
      wf.indx          = i;
      wf.E             = orderedE[i];
      wf.valid         = 1;
      printf("shifted back eig:%E\n",orderedE[i]);
    }
    writeEigenvalues(d, 'a');                                           //write eigenvalues into a file 
    
    //put back the eigenvectors accounting for previous sort of the eigenvalues  
    if(d->opt.ExecParam.ElCalc.ElAlg.PQTracemin.CalEigVec) {
      d->wfc = Cmatrix(d->neigv, nloc);                                 //d->neigv was ncv, but I think it was a mistake and should've been nev in the original version???       
      for (i=0; i<d->neigv; i++){
	for (j=0; j<nloc; j++){
	   #ifdef MATRIX_MULT
	      d->wfc[i][j].r = Y[orderedI[i]][j].r;
	      d->wfc[i][j].i = Y[orderedI[i]][j].i;
	   #else
	      d->wfc[i][j].r = Y[orderedI[i]*nloc+j].r;
	      d->wfc[i][j].i = Y[orderedI[i]*nloc+j].i;
	   #endif
	}
      }
  
      writePsiSqr(d->wfc, d->neigv, d, 0);                              //not sure what this does???       
      if (d->opt.Dev.band_model==Dev_struct::Bands_1_s_nospin){         //in the single band case the data file is small enough so that we write out the full wavefunction
	print_psisq_gnu_pdb(d->wfc,d->neigv,d,0);
      }
      for(i=0; i<d->neigv; i++){
	writeEigenvector(d, d->wfc[i], i);                              //write eigenvector to a file
      }
    }
  }
  

  compute_ritz_residual_norm_deviation(d);
  //compute_ritz_value(d);   
  //compute_residual(d);
  //compute_eigenvector_norm(d);

  printf("EXITING PQTRACEMIN DRIVER\n");
  return d->neigv;  

#else
#ifdef USE_TRACEMIN_SER
//  if(d->rc==NULL) d->rc = Cvectr(d->nvmap[mpi_n3d_id]*d->NBasisStates); //kept this line from parpack_driver.c, not sure what it's for???
  int  neig  = d->opt.ExecParam.ElCalc.ElAlg.PQTracemin.NumEigVal;         //number of eigenvalues wanted                                                          
  real refE  = d->opt.ExecParam.ElCalc.ElAlg.PQTracemin.RefEigVal;         //shift the matrix A so that you are finding eigenvalues closest to the shift
  real tol   = d->opt.ExecParam.ElCalc.ElAlg.PQTracemin.Tolerance;         //accuracy wanted 
  real cgtol  = d->opt.ExecParam.ElCalc.ElAlg.PQTracemin.cgtol;         //initial tol of GC
  int maxi = d->opt.ExecParam.ElCalc.ElAlg.PQTracemin.maxi;  		//this was set to 40, in line 48 
  int cgmaxi  = d->opt.ExecParam.ElCalc.ElAlg.PQTracemin.cgmaxim;         // Number of iterations of CG 
  int n      = d->n_ham_tot;                                            //dimension of eigenproblem                          
  int nloc   = d->seg_ln[mpi_n3d_id];                                   //number of vector elements on current processor                                         
  int info;                                                             //flag for reverse communication; initialized to 0                                            
  //declare important MPI parameters
  //MPI_Fint comm = MPI_Comm_c2f(MPI_COMM_WORLD);                         //transform C into Fortran communicator (currently tracemin uses MPI_COMM_WORLD and comm doesn't matter)
  int comm = 0;
  int rank = 0;                                                             //rank of the current processor               
  int np = 1;                                                               //number of processors 
  //allocate arrays containing result of tracemin
  rvectr E = Rvectr(neig);                                              //real eigenvalues                                                                   
  //allocate necessary temporary arrays  
  #ifdef MATRIX_MULT
     cmatrix Y    = Cmatrix(neig,n);                                         //complex eigenvectors 
     cmatrix Q    = Cmatrix(neig,n);                                             
     cmatrix T    = Cmatrix(neig,n);
     cmatrix Qcg  = Cmatrix(neig,n);
     cmatrix P    = Cmatrix(neig,n);
     cmatrix X    = Cmatrix(neig,n);
     cmatrix Rcg  = Cmatrix(neig,n);
     cmatrix R    = Cmatrix(neig,n);                                         //estimate of residual                                                           
  #else
     cvectr Y = Cvectr(n*neig);                                         //complex eigenvectors 
     cvectr Q    = Cvectr(n*neig);                                             
     cvectr T    = Cvectr(n*neig);
     cvectr Qcg  = Cvectr(n*neig);
     cvectr P    = Cvectr(n*neig);
     cvectr X    = Cvectr(n*neig);
     cvectr Rcg  = Cvectr(n*neig);
     cvectr R = Cvectr(n*neig);                                     //estimate of residual                                  
  #endif

  cvectr H    = Cvectr(neig*neig);  
  cvectr Hloc = Cvectr(neig*neig);
  cvectr ZWork= Cvectr(max(1,2*neig-1)); 
  rvectr RWork= Rvectr(max(1,3*neig-2));
  ivectr IWork= Ivectr(neig);
  cvectr G    = Cvectr(neig*neig);
  cvectr Gloc = Cvectr(neig*neig);
  cvectr gamma_old = Cvectr(neig);
  rvectr rcgtol    = Rvectr(neig);
//commented by marta: defs of cgmaxi and cgtol to be defined inpt deck; 
  int  i,j,k,cginfo;//,cgmaxi;
//real cgtol;                                                           //parameter used by tracemin
  int nloct;                                                            //total number of elements on this processor: nloc*neig
  //ina maxi;                                                             //maximum number of iterations of tracemin (there are neig*2*6*2 matvecs per iteration)
  complex ytloc;                                                        //temporary variables used in computation of eigenvalues after tracemin (next 4 lines)
  complex yyloc;                
  real yt;                      
  real yy;                      
  rvectr Rnrm = Rvectr(neig);                                           //norm of the residual of of individual eigenpairs
  //variables used in sorting of the eigenvalues after tracemin
  rvectr orderedE;
  ivectr orderedI;
  double dtmp;
  int    itmp;
  //declare some constants needed in the code
  int  IONE = 1;
  real DONE = 1.0;
  real DTWO = 2.0;  
  int iseed[4];

  printf("ENTERING SQTRACEMIN DRIVER\n");
  printf("Setup some variables\n");
  //MPI_Comm_rank(MPI_COMM_WORLD, &rank);                                 //find out my rank
  //MPI_Comm_size(MPI_COMM_WORLD, &np);                                   //find out number of processors  
  nloct   =neig*n;                                                   //number of elements of all vectors combined on this processor 
  //maxi    =40;  (commented by marta)                                           //maximum number of tracemin iterations (SHOULD BE AN OPTION IN THE INPUT FILE) 
  //cgtol   = DONE;    (commented by marta)                                                   //initial cg tolerance is set to 1.0    (SHOULD BE AN OPTION IN THE INPUT FILE)
  iseed[0]=103+rank;                                                    //set up seed for the random number generator (making sure Y is random across processor)
  iseed[1]=101*rank;                                                      
  iseed[2]=7;
  iseed[3]=57;

  #ifdef MATRIX_MULT
     zlarnv(&IONE,iseed,&nloct,&Y[0][0]);                                         //generate random Y 
  #else
     zlarnv(&IONE,iseed,&nloct,Y);                                         //generate random Y 
  #endif

  printf("%d,----- start ztracemin -----\n",rank);
  printf("%d,parameters: n=%d, nloc=%d, neig=%d, maxi=%d, tol=%E\n",rank,n,nloc,neig,maxi,tol);
  printf("%d, refE=%f\n",rank,refE);

  printf("%d,perform initial two matrix-vector multiplies\n",rank);
  //call pcfmv(M,Y,T,neig,nloc,ne,rank,np)
  //call pcfmv(M,T,Q,neig,nloc,ne,rank,np)

  #ifdef MATRIX_MULT
     Hmatmult_spds_complete(T, d, refE, Y, neig);
	 Hmatmult_spds_complete(Q, d, refE, T, neig);
  #else
     for (k=0; k<neig; k++){
        Hmult_spds_complete(&T[k*n], d, refE, &Y[k*n]);
     }
     for (k=0; k<neig; k++){
        Hmult_spds_complete(&Q[k*n], d, refE, &T[k*n]);
     }
  #endif

  printf("%d,start the main for loop\n",rank);
  for(i=1; i<=2*maxi+1; i++){ //every iteration requires one matrix vector multiply and one system solve with matrix M=A^2 
    #ifdef MATRIX_MULT
       ztracemin_wrapper(&n,&neig,&Y[0][0],E,&Q[0][0],H,Hloc,ZWork,RWork,IWork,&R[0][0],&tol,&maxi,&i,&nloc,&rank,&np,&comm,&info);
	#else
       ztracemin_wrapper(&n,&neig,Y,E,Q,H,Hloc,ZWork,RWork,IWork,R,&tol,&maxi,&i,&nloc,&rank,&np,&comm,&info);
	#endif
    if(info == 1){
      printf("%d,--- ztracemin has converged --- (iteration=%d; info=%d)\n",rank,i/2,info);
      break;
    }
    if (info == 2){
      printf("%d,Perform a solve\n",rank);
      //setup CG parameters
      //call pcfmv(M,Y,T,  neig,nloc,ne,rank,np)
      //call pcfmv(M,T,Qcg,neig,nloc,ne,rank,np)

	  #ifdef MATRIX_MULT
	     Hmatmult_spds_complete(T, d, refE, Y, neig);
	     Hmatmult_spds_complete(Qcg, d, refE, T, neig);
	  #else
         for (k=0; k<neig; k++){
	        Hmult_spds_complete(&T[k*n], d, refE, &Y[k*n]);
         }
         for (k=0; k<neig; k++){
	        Hmult_spds_complete(&Qcg[k*n], d, refE, &T[k*n]);
         }
	  #endif
      //cgmaxi=min(n,1000);
      if (i%10 == 0){
	cgtol = cgtol/10;
	if (cgtol < tol){
	  cgtol=tol;
	}
      }
      printf("CG: tolerance=%f,maxi=%d\n",cgtol,cgmaxi);

      //start CG
      for (j=1; j<=cgmaxi+2; j++){ //notice that one iteration is for setup
	  #ifdef MATRIX_MULT
	     modified_conjugate_gradient_wrapper(&n,&neig,&X[0][0],&Y[0][0],&Qcg[0][0],G,Gloc,&P[0][0],&Rcg[0][0],gamma_old,rcgtol,&j,&cgmaxi,&cgtol,&tol,&nloc,&rank,&np,&comm,&cginfo);
	  #else
	     modified_conjugate_gradient_wrapper(&n,&neig,X,Y,Qcg,G,Gloc,P,Rcg,gamma_old,rcgtol,&j,&cgmaxi,&cgtol,&tol,&nloc,&rank,&np,&comm,&cginfo);
	  #endif
	if (cginfo == 1){
	  printf("%d,modified CG converged at iteration %d\n",rank,j);
	  break;
	}
	if (cginfo == 2){
	  //printf("%d,Perform two matrix-vector multiplies (CG)\n",rank);
	  //call pcfmv(M,P,T,  neig,nloc,ne,rank,np)
	  //call pcfmv(M,T,Qcg,neig,nloc,ne,rank,np)

	  #ifdef MATRIX_MULT
	     Hmatmult_spds_complete(T, d, refE, P,neig);
	     Hmatmult_spds_complete(Qcg, d, refE, T,neig);
	  #else
	     for (k=0; k<neig; k++){
	        Hmult_spds_complete(&T[k*n], d, refE, &P[k*n]);
	     }
	     for (k=0; k<neig; k++){
	        Hmult_spds_complete(&Qcg[k*n], d, refE, &T[k*n]);
	     }
	  #endif

	}
	if (cginfo == 3){
	  printf("%d,modified CG has reached maximum number of iterations at iteration %d\n",rank,j);
	  break;
	}
	if (cginfo < 0){
	  printf("%d,ERROR in modified CG, at iteration %d\n",rank,j);
	  break;
	}
      }
      //copy the solution X=Y-D back into Y and compute Q=A*Y
	  #ifdef MATRIX_MULT
         zcopy(&nloct,&X[0][0],&IONE,&Y[0][0],&IONE);
	  #else
         zcopy(&nloct,X,&IONE,Y,&IONE);
	  #endif

      //call pcfmv(M,Y,T,neig,nloc,ne,rank,np)
      //call pcfmv(M,T,Q,neig,nloc,ne,rank,np)

	  #ifdef MATRIX_MULT
	     Hmatmult_spds_complete(T, d, refE, Y,neig);
	     Hmatmult_spds_complete(Q, d, refE, T,neig);
	  #else
         for (k=0; k<neig; k++){
	        Hmult_spds_complete(&T[k*n], d, refE, &Y[k*n]);
         }
         for (k=0; k<neig; k++){
	        Hmult_spds_complete(&Q[k*n], d, refE, &T[k*n]);
         }
	  #endif

    }
    if (info == 3){
      printf("%d,Perform two matrix-vector multiplies\n",rank);
      //call pcfmv(M,Y,T,neig,nloc,ne,rank,np)
      //call pcfmv(M,T,Q,neig,nloc,ne,rank,np)

	  #ifdef MATRIX_MULT
	     Hmatmult_spds_complete(T, d, refE, Y,neig);
	     Hmatmult_spds_complete(Q, d, refE, T,neig);
	  #else
         for (k=0; k<neig; k++){
	        Hmult_spds_complete(&T[k*n], d, refE, &Y[k*n]);
         }
         for (k=0; k<neig; k++){
	        Hmult_spds_complete(&Q[k*n], d, refE, &T[k*n]);
         }
	  #endif

    }
    if (info == 4){
      printf("%d,--- Maximum number of iterations reached ---(iteration=%d; info=%d)\n",rank,i/2,info);
      break;
    }
    if (info < 0){
      printf("%d,--- Error in the ztracemin routine --- (iteration=%d; info=%d)\n",rank,i/2,info);
      break; 
    }
    printf("%d,done with do loop iteration -->%d\n",rank,i);
  }
  maxi=i/2; 

  printf("%d,find the eigenvalues lambda using lambda=xt*A*x/xt*x\n",rank);
  //call pcfmv(M,Y,T,neig,nloc,ne,rank,np)

  #ifdef MATRIX_MULT
     Hmatmult_spds_complete(T, d, refE, Y,neig);
  #else
     for (k=0; k<neig; k++){
        Hmult_spds_complete(&T[k*n], d, refE, &Y[k*n]);
     }
  #endif

  for (i=0; i<neig; i++){
    #ifdef MATRIX_MULT
       zdotc(&ytloc,&n,Y[i],&IONE,T[i],&IONE);   
       zdotc(&yyloc,&n,Y[i],&IONE,Y[i],&IONE);   
	#else
       zdotc(&ytloc,&n,&Y[i*n],&IONE,&T[i*n],&IONE);   
       zdotc(&yyloc,&n,&Y[i*n],&IONE,&Y[i*n],&IONE);   
	#endif

    //MPI_Allreduce(&ytloc.r, &yt, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
    //MPI_Allreduce(&yyloc.r, &yy, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
	yt = ytloc.r;
	yy = yyloc.r;
    E[i]=yt/yy;
  }

  printf("%d,print the eigenvalues and the norm of the estimate of the residual\n",rank);
  for (i=0; i<neig; i++){
    #ifdef MATRIX_MULT
       Rnrm[i]=dznrm2(&n, R[i], &IONE);
	#else
       Rnrm[i]=dznrm2(&n, &R[i*n], &IONE);
	#endif

    printf("eigenvalue= %E, norm of its residual=%E\n",E[i],Rnrm[i]);
  }  
  printf("%d,print the eigenvalues and the norm of the real residual\n",rank);
  for (i=0; i<neig; i++){
    yyloc.r =-E[i];
    yyloc.i = 0.0;
    yy      = 0.0;

	#ifdef MATRIX_MULT
       zaxpy(&n,&yyloc,Y[i],&IONE,T[i],&IONE);
       zdotc(&yyloc,&n,T[i],&IONE,T[i],&IONE);
	#else
       zaxpy(&n,&yyloc,&Y[i*n],&IONE,&T[i*n],&IONE);
       zdotc(&yyloc,&n,&T[i*n],&IONE,&T[i*n],&IONE);
	#endif

    //MPI_Allreduce(&yyloc.r, &yy, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
	yy = yyloc.r;
    printf("eigenvalue= %E, norm of its residual=%E\n",E[i],sqrt(yy));
  }
    
  printf("%d,deallocate temporary variables\n",rank);
  rm_cvectr(&H);
  rm_cvectr(&Hloc);
  rm_cvectr(&ZWork);
  rm_rvectr(&RWork);
  rm_ivectr(&IWork);
  rm_cvectr(&G); 
  rm_cvectr(&Gloc);
  rm_cvectr(&gamma_old);
  rm_rvectr(&rcgtol);  

  #ifdef MATRIX_MULT
     rm_cmatrix(&Rcg);
     rm_cmatrix(&X);
     rm_cmatrix(&P); 
     rm_cmatrix(&Q);
     rm_cmatrix(&Qcg);
     rm_cmatrix(&T); 
  #else
     rm_cvectr(&Rcg);
     rm_cvectr(&X);
     rm_cvectr(&P); 
     rm_cvectr(&Q);
     rm_cvectr(&Qcg);
     rm_cvectr(&T); 
  #endif
  
  printf("put results back into qd_struct d\n");
  /*
  //assuming eigenvalues are already sorted  
  if(info >=0) {                                                   
    //put back the eigenvalues
    d->neigv = neig;
    d->WF = vector<waveFunction>(d->neigv);
    for (i=0; i<d->neigv; i++) {                                        //wf points to d->WF[i] (C++ syntax)
      waveFunction& wf = d->WF[i];
      wf.indx          = i;
      wf.E             = E[i] + refE; 
      wf.valid         = 1;
    }
    writeEigenvalues(d, 'a');                                           //write eigenvalues into a file
 
    //put back the eigenvectors
    if(d->opt.ExecParam.ElCalc.ElAlg.PQTracemin.CalEigVec) { 
      d->wfc = Cmatrix(d->neigv, n);                                 //d->neigv was ncv, but I think it was a mistake and should've been nev in the original version???
      for (i=0; i<d->neigv; i++){
	for (j=0; j<n; j++){
	  #ifdef MATRIX_MULT
	     d->wfc[i][j].r = Y[i][j].r;
	     d->wfc[i][j].i = Y[i][j].i;
	  #else
	     d->wfc[i][j].r = Y[i*n+j].r;
	     d->wfc[i][j].i = Y[i*n+j].i;
	  #endif
	}
      }  

      writePsiSqr(d->wfc,d->neigv,d,0);                                 //not sure what this does???      
      if (d->opt.Dev.band_model==Dev_struct::Bands_1_s_nospin){         //in the single band case the data file is small enough so that we write out the full wavefunction 
	print_psisq_gnu_pdb(d->wfc,d->neigv,d,0);                  
      }

      for(i=0; i<d->neigv; i++){
	writeEigenvector(d, d->wfc[i], i);                              //write eigenvector to a file  
      }
    }
  }
  */
  
  //assuming needs to sort the eigenvalues
  if(info>=0) {
    //put back sorted eigenvalues
    d->neigv = neig;
    d->WF    = vector<waveFunction>(d->neigv);
    orderedE = Rvectr(d->neigv);
    orderedI = Ivectr(d->neigv);
    for(i=0; i<d->neigv; i++) {                     
      orderedE[i]=E[i] + refE;                                          //assign shifted back eigenvalues to another array 
      orderedI[i]=i;                                                    //keep track of the original index of the eigenvalue
    }    
    for(i=0; i<d->neigv-1;i++){                                         //do bubble sort
      for(j=i+1; j<d->neigv; j++){
	if(orderedE[i]>orderedE[j]) {
	  dtmp       =orderedE[i];
	  orderedE[i]=orderedE[j];
	  orderedE[j]=dtmp;
	  itmp       = orderedI[i];
	  orderedI[i]=orderedI[j];
	  orderedI[j]=itmp;
	}
      }
    }
    for (i=0; i<d->neigv; i++) {                                        //wf points to d->WF[i] (C++ syntax)
      waveFunction& wf = d->WF[i];
      wf.indx          = i;
      wf.E             = orderedE[i];
      wf.valid         = 1;
      printf("shifted back eig:%E\n",orderedE[i]);
    }
    writeEigenvalues(d, 'a');                                           //write eigenvalues into a file 
    
    //put back the eigenvectors accounting for previous sort of the eigenvalues  
    if(d->opt.ExecParam.ElCalc.ElAlg.PQTracemin.CalEigVec) {
      d->wfc = Cmatrix(d->neigv, n);                                 //d->neigv was ncv, but I think it was a mistake and should've been nev in the original version???       
      for (i=0; i<d->neigv; i++){
	for (j=0; j<n; j++){
	   #ifdef MATRIX_MULT
	      d->wfc[i][j].r = Y[orderedI[i]][j].r;
	      d->wfc[i][j].i = Y[orderedI[i]][j].i;
	   #else
	      d->wfc[i][j].r = Y[orderedI[i]*n+j].r;
	      d->wfc[i][j].i = Y[orderedI[i]*n+j].i;
	   #endif
	}
      }
  
      writePsiSqr(d->wfc, d->neigv, d, 0);                              //not sure what this does???       
      if (d->opt.Dev.band_model==Dev_struct::Bands_1_s_nospin){         //in the single band case the data file is small enough so that we write out the full wavefunction
	print_psisq_gnu_pdb(d->wfc,d->neigv,d,0);
      }
      for(i=0; i<d->neigv; i++){
	writeEigenvector(d, d->wfc[i], i);                              //write eigenvector to a file
      }
    }
  }
  

  compute_ritz_residual_norm_deviation(d);
  //compute_ritz_value(d);   
  //compute_residual(d);
  //compute_eigenvector_norm(d);

  printf("EXITING SQTRACEMIN DRIVER\n");
  return d->neigv;  
#endif
  return 0;
#endif
}
