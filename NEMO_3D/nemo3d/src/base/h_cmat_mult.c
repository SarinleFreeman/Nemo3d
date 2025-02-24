#define print_com2(A,B) ;
#define print_com3(A,B,C) ;

//#define USE_BLAS
#undef USE_MKL
#ifdef USE_BLAS
   //#define TRANSPOSED
#endif

#ifdef USE_MKL
   #include "mkl.h"
#endif

// Uncomment this part if you want to send the matrix as one chunk, not the vector by number of column times
#define SENDRECV_ONE_CHUNK

#include "h_cvectr_mult.h"
#include "h_cmat_mult.h"


void matmat_mult(const int** argList, nml_dcmatrix* y, const nml_dcmatrix* x)
{
   qd_struct d = *((qd_struct*) argList[0]);
   double shift = *((double*) argList[1]);
   int num_col = *((int *) argList[2]);
   int Nproj = *((int*) argList[3]);

   Hmatmult_spds_complete((cmatrix) y, d, shift, (cmatrix) x, num_col);

   #if 0

   if (Nproj==0) return;

   int* eval_scoperti = *((int**) argList[3]);

   // compute ip_j = E_j*<psi_j|x>
   complex* ip_loc = new complex[Nproj];
   for (int j=0; j<Nproj; j++) {
      const waveFunction& wf = d->WF[ eval_scoperti[j] ];
      complex* ip_loc_j = ip_loc + j;
      ip_loc_j->r = ip_loc_j->i = 0.0;
      for (int i=0; i<d->seg_ln[mpi_n3d_id]; i++) {
         const complex& psi = wf.psi[i];
         const complex& xi = x[i];
         ip_loc_j->r += (psi.r*xi.r + psi.i*xi.i);
         ip_loc_j->i += (psi.r*xi.i - psi.i*xi.r);
      }
      ip_loc_j->r *= wf.E;
      ip_loc_j->i *= wf.E;
   }
   complex* ip = new complex[Nproj];
   MPI_Allreduce(ip_loc, ip, 2*Nproj, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   
   for (int j=0; j<Nproj; j++) {
      const waveFunction& wf = d->WF[ eval_scoperti[j] ];
      complex* ip_j = ip + j;
      printf("mpi_n3d_id=%d  j=%d  ip_loc=(%e %e) ip=(%e %e)\n",
             mpi_n3d_id,
             j,
             ip_loc[j].r, ip_loc[j].i,
             ip[j].r, ip[j].i);

      masterPrint("eval_scoperti[%d]=%d   (%lf %lf)  %lf\n", 
                  j, eval_scoperti[j],
                  ip_j->r, ip_j->i, wf.E);
   }

   // y -= sum_j { ip_j |psi1> }
   for (int j=0; j<Nproj; j++) {
      const waveFunction& wf = d->WF[ eval_scoperti[j] ];
      complex* ip_j = ip + j;
      for (int i=0; i<d->seg_ln[mpi_n3d_id]; i++) {
         const complex& psi = wf.psi[i];
         y[i].r -= (ip_j->r * psi.r - ip_j->i * psi.i);
         y[i].i -= (ip_j->r * psi.i + ip_j->i * psi.r);
      }
   }

   delete [] ip_loc;
   delete [] ip;
   #endif
}

/* Compute Y = (H - shift)*X st that each processor gets its appropriate
   portion of Y and no further communication is required. */
void Hmatmult_spds_complete( cmatrix Y, qd_struct d, real shift, cmatrix X, int num_col )
{
   real Tstart=0.0;
   static int matrix_multiply =0;
   MPI_TIC(Tstart);
   if (0 == matrix_multiply%100 && mpi_n3d_id== 0)
     fprintf(stdout,"# of matrix-matrix multiplies: %d num_col: %d total matrix-vector multiplies: %d\n",
	 matrix_multiply, num_col,matrix_multiply*num_col);
	 fflush(stdout);
   matrix_multiply ++;
   
   #define NEED_TO_USE_STORED_MATRIX_ALWAYS

   #ifndef MPI3d  // serial version doesn't need to do any fancy communication
      cmatrix YC=NULL; // a dummy placeholder that won't be used in matvect mult
   
      if (d->opt.ExecParam.ElCalc.ElAlg.HamCompute==ElAlg_struct::SparseStore && 
         (d->Hzb.isInitialized() || d->HzbComplex.isInitialized()))
         Hmatmult_spds_col_stored( Y, YC, d, shift, X, num_col, 0 ); 
      else if (d->opt.ExecParam.ElCalc.ElAlg.HamCompute==ElAlg_struct::FullStore && d->Ham)
         Hmatmult_spds_col_fullstored_par( Y, YC, d, shift, X, num_col, 0 ); 
      else{
         Hmatmult_spds_col( Y, YC, d, shift, X, num_col,  0 );

         #ifdef NEED_TO_USE_STORED_MATRIX_ALWAYS
         /* Here we have just computed the matrix vector product from an explicit Hamiltonian
         * construction.  If Storage of the Hamiltonian was requested, then the Hamiltonian
         * that was computed on the fly was indeed stored.
         *
         * Through some numerical experiments we have found that the accuracy that was 
         * obtained from the recomputed Hamiltonian and stored Hamiltonian is different
         * enough to cause problems in the Lanczos algorithm.
         * We therefore repeat the matrix-vector multiply here one more time 
         * to make sure we always use the stored matrix rather than the on-the fly
         * matrix. 
         */
            if (d->opt.ExecParam.ElCalc.ElAlg.HamCompute!=ElAlg_struct::Recompute){
               if (d->opt.ExecParam.ElCalc.ElAlg.HamCompute==ElAlg_struct::SparseStore && 
                  (d->Hzb.isInitialized() || d->HzbComplex.isInitialized()))
                  Hmatmult_spds_col_stored( Y, YC, d, shift, X, num_col, 0 ); 
               else if (d->opt.ExecParam.ElCalc.ElAlg.HamCompute==ElAlg_struct::FullStore && d->Ham)
                  Hmatmult_spds_col_fullstored_par( Y, YC, d, shift, X, num_col, 0 ); 
               else
                  die("ERROR in Hmatmult_spds_complete:\n"
	                  "Unidentified Hamiltonian storage system\n");
            }
         #endif

      }
   
   #else // The parallel version needs to set up the communication arrays
   
      static int send_max=0;
      static ivectr sflag=NULL;

	  #ifndef SENDRECV_ONE_CHUNK 
      static cmatrix receive_data=NULL;
      static cmatrix send_data=NULL;
	  #else
      static cvectr receive_data=NULL;
      static cvectr send_data=NULL;
	  #endif // SENDRECV_ONE_CHUNK

      static imatrix cmap_row, cmap_col;

      #ifdef USE_BLAS
         #ifdef TRANSPOSED
            if(d->mm==NULL) d->mm = Cmatrix(num_col, d->nvmap[mpi_n3d_id]*d->NBasisStates); 
         #else
            if(d->mm==NULL) d->mm = Cmatrix(d->nvmap[mpi_n3d_id]*d->NBasisStates,num_col); 
         #endif
      #else
         if(d->mm==NULL) d->mm = Cmatrix(num_col, d->nvmap[mpi_n3d_id]*d->NBasisStates); 
	  #endif // USE_BLAS

      cmatrix YC=d->mm;
      cmatrix source=NULL, target=NULL;
      int i, j, k, cnt, src_i, tgt_i, lo_indx, commstep;
	  static int maxcom;
      int nyc = 0;

      print_com2("CPU=%d start Hmatmult_spds_complete\n",mpi_n3d_id);
      
	  /* Allocate and initialize static stuff if required */

      if ( sflag == NULL ) {
         print_com2("CPU=%d initialize static vectors in Hmatmult_spds_complete\n", mpi_n3d_id);
         maxcom = 0;
         for ( i=0; i < d->nproc; i++ ) {
            for ( j=0; j < d->nproc; j++ ) {
               if ( d->commsize[i][j] > maxcom )
                  maxcom = d->commsize[i][j];
         }}
      
         cmap_row = Imatrix(d->nproc,maxcom);
         cmap_col = Imatrix(d->nproc,maxcom);
      
         maxcom *= d->NBasisStates;

		 printf("[MATRIX MULTIPLY] maximum communication size = %d\n", maxcom);
      
         for (tgt_i=0; tgt_i < d->nproc; tgt_i++){
	        if (d->commsize[mpi_n3d_id][tgt_i]>0){
               cnt = 0;
	           for ( i=0; i < d->nvmap[mpi_n3d_id]; i++ ) {
	              if ( d->vecmap[mpi_n3d_id][i] >= d->seg_s[tgt_i] && 
		               d->vecmap[mpi_n3d_id][i] < d->seg_s[tgt_i]+d->seg_ln[tgt_i] ){
		             cmap_row[tgt_i][cnt] = i;
		             cnt++;
	     }}}}
      
         for (tgt_i=0; tgt_i < d->nproc; tgt_i++){
	        if (d->commsize[mpi_n3d_id][tgt_i]>0){
               cnt = 0;
	           for ( i=0; i < d->nvmap[tgt_i]; i++ ){
	              if ( d->vecmap[tgt_i][i] >= d->seg_s[mpi_n3d_id] && 
		               d->vecmap[tgt_i][i] < d->seg_s[mpi_n3d_id]+d->seg_ln[mpi_n3d_id] ){
		             cmap_col[tgt_i][cnt] = i;
		             cnt++;
	     }}}}
      
	  #ifndef SENDRECV_ONE_CHUNK 
         receive_data = Cmatrix(num_col, maxcom);
         print_com2("CPU=%d did receive_data\n",mpi_n3d_id);
         send_data = Cmatrix(num_col, maxcom);
         print_com2("CPU=%d did send_data\n",mpi_n3d_id);
	  #else
         receive_data = Cvectr(num_col* maxcom);
         print_com2("CPU=%d did receive_data\n",mpi_n3d_id);
         send_data = Cvectr(num_col* maxcom);
         print_com2("CPU=%d did send_data\n",mpi_n3d_id);
	  #endif // SENDRECV_ONE_CHUNK
      }

      if (d->opt.ExecParam.ElCalc.ElAlg.HamCompute==ElAlg_struct::SparseStore && 
         (d->Hzb.isInitialized() || d->HzbComplex.isInitialized())){
         nyc = Hmatmult_spds_col_stored( Y, YC, d, shift, X, num_col, mpi_n3d_id ); 
      }
      else if (d->opt.ExecParam.ElCalc.ElAlg.HamCompute==ElAlg_struct::FullStore && d->Ham)
      {
         nyc = Hmatmult_spds_col_fullstored_par( Y, YC, d, shift, X, num_col, mpi_n3d_id ); 
	  }
      else{
         nyc = Hmatmult_spds_col( Y, YC, d, shift, X, num_col, mpi_n3d_id );

         #ifdef NEED_TO_USE_STORED_MATRIX_ALWAYS
         /* Here we have just computed the matrix vector product from an explicit Hamiltonian
          * construction.  If Storage of the Hamiltonian was requested, then the Hamiltonian
          * that was computed on the fly was indeed stored.
          *
          * Through some numerical experiments we have found that the accuracy that was 
          * obtained from the recomputed Hamiltonian and stored Hamiltonian is different
          * enough to cause problems in the Lanczos algorithm.
          * We therefore repeat the matrix-vector multiply here one more time 
          * to make sure we always use the stored matrix rather than the on-the fly
          * matrix. 
          */
             if (d->opt.ExecParam.ElCalc.ElAlg.HamCompute!=ElAlg_struct::Recompute){
                if (d->opt.ExecParam.ElCalc.ElAlg.HamCompute==ElAlg_struct::SparseStore && 
                   (d->Hzb.isInitialized() || d->HzbComplex.isInitialized()))
                   Hmatmult_spds_col_stored( Y, YC, d, shift, X, num_col, mpi_n3d_id ); 
                else if (d->opt.ExecParam.ElCalc.ElAlg.HamCompute==ElAlg_struct::FullStore && d->Ham)
                   Hmatmult_spds_col_fullstored_par( Y, YC, d, shift, X, num_col, mpi_n3d_id ); 
                else
                   die("ERROR in Hmatmult_spds_complete:\n"
	                   "Unidentified Hamiltonian storage system\n");
            }
         #endif
      }

      src_i = mpi_n3d_id;
   
      if (sflag==NULL) {
         sflag = Ivectr(d->nproc);
         for (tgt_i=0; tgt_i<d->nproc; tgt_i++){
	        if (d->commsize[src_i][tgt_i]>0){
	           for ( i=0; i < nyc; i++ ){
	              if ( d->vecmap[mpi_n3d_id][i] >= d->seg_s[tgt_i] && 
		               d->vecmap[mpi_n3d_id][i] < d->seg_s[tgt_i]+d->seg_ln[tgt_i] ){
		             sflag[tgt_i] = 1;
		             send_max++;
		             break;
      }}}}}
   
      // NB: see old version of code for old sendreceive
      for ( commstep=0; commstep < 2; commstep++ ) {
         tgt_i = d->commorder[src_i][commstep];
         if ( tgt_i != -1 ) {

	        //source = receive_data;

            #if 0 // I don't know why this part causes error so I commented out
            for(j = 0; j < num_col; j++) {
			   for(i = 0; i < maxcom; i++) {
			      send_data[j][i].r = 0.0;
				  send_data[j][i].i = 0.0;
			}}
            #endif

	        //target = send_data;
			
			//for ( k=0; k < num_col; k++)
			//{

		    #ifdef SENDRECV_ONE_CHUNK
			int comm_data_size = 0;
			//printf("src_i = %d, tgt_i = %d, commsize = %d %d\n", src_i, tgt_i, d->commsize[src_i][tgt_i], d->commsize[tgt_i][src_i]);
			if(d->commsize[src_i][tgt_i] >= d->commsize[tgt_i][src_i]) comm_data_size = d->commsize[src_i][tgt_i];
			else comm_data_size = d->commsize[tgt_i][src_i];
			#endif

			for ( k=0; k < num_col; k++)
	           for ( i=0; i < d->commsize[src_i][tgt_i]; i++ )
	              for ( j=0; j < d->NBasisStates; j++ ) {
	                 //target[k][i*d->NBasisStates+j] = YC[k][cmap_row[tgt_i][i]*d->NBasisStates+j];   
					 #ifdef USE_BLAS
					    #ifdef TRANSPOSED
	                       send_data[k][i*d->NBasisStates+j] = YC[k][cmap_row[tgt_i][i]*d->NBasisStates+j];   
					    #else
	                       send_data[k][i*d->NBasisStates+j] = YC[cmap_row[tgt_i][i]*d->NBasisStates+j][k];   
					    #endif
					 #else

					    #ifdef SENDRECV_ONE_CHUNK
	                    //send_data[k*d->NBasisStates*d->commsize[src_i][tgt_i]+i*d->NBasisStates+j] = YC[k][cmap_row[tgt_i][i]*d->NBasisStates+j];   
	                    send_data[k*d->NBasisStates*comm_data_size+i*d->NBasisStates+j] = YC[k][cmap_row[tgt_i][i]*d->NBasisStates+j];   
						//comm_data_size += d->NBasisStates*d->commsize[src_i][tgt_i];
						#else
	                    send_data[k][i*d->NBasisStates+j] = YC[k][cmap_row[tgt_i][i]*d->NBasisStates+j];   
						#endif // SENDRECV_ONE_CHUNK
					
					 #endif // USE_BLAS
				  }

            real matmult_Tstart = 0.0;
	        MPI_TIC(matmult_Tstart);

            #ifdef SENDRECV_ONE_CHUNK
	        MPI_Sendrecv(//target[k], 
			      &send_data[0],
				  d->NBasisStates*num_col*comm_data_size*2,
                  //2*d->commsize[src_i][tgt_i]*d->NBasisStates*num_col, 
                  //2*maxcom*num_col, 
		          MPI_DOUBLE,tgt_i, MPI_rc_distr, 
		          //source[k], 
				  &receive_data[0],
				  d->NBasisStates*num_col*comm_data_size*2,
                  //2*d->commsize[tgt_i][src_i]*d->NBasisStates*num_col, 
                  //2*maxcom*num_col, 
		          MPI_DOUBLE,tgt_i, MPI_rc_distr, 
				  MPI_COMM_WORLD, 
                  &mpi_n3d_status);
            #else
			for ( k=0; k < num_col; k++)
			{
	           MPI_Sendrecv(//target[k], 
			      send_data[k],
                  2*d->commsize[src_i][tgt_i]*d->NBasisStates, 
		          MPI_DOUBLE,tgt_i, MPI_rc_distr, 
		          //source[k], 
				  receive_data[k],
                  2*d->commsize[tgt_i][src_i]*d->NBasisStates, 
		          MPI_DOUBLE,tgt_i, MPI_rc_distr, 
				  MPI_COMM_WORLD, 
                  &mpi_n3d_status);
			}
            #endif // SENDRECV_ONE_CHUNK

	        MPI_TOC(mpiTiming.matmult_sendrecv, matmult_Tstart);

	        /* add the scattered packages to the main array */
	        //target = Y;	
	        lo_indx=d->seg_s[src_i];

            for ( k=0; k < num_col; k++) {
	           for ( i=0; i < d->commsize[tgt_i][src_i]; i++ ) {
	              for ( j=0; j < d->NBasisStates; j++ ) {
                     //target[k][d->vecmap[tgt_i][cmap_col[tgt_i][i]] - lo_indx+j].r += source[k][i*d->NBasisStates + j].r;
                     //target[k][d->vecmap[tgt_i][cmap_col[tgt_i][i]] - lo_indx+j].i += source[k][i*d->NBasisStates + j].i;
					 #ifdef USE_BLAS
					    #ifdef TRANSPOSED
                           Y[k][d->vecmap[tgt_i][cmap_col[tgt_i][i]] - lo_indx+j].r += receive_data[k][i*d->NBasisStates + j].r;
                           Y[k][d->vecmap[tgt_i][cmap_col[tgt_i][i]] - lo_indx+j].i += receive_data[k][i*d->NBasisStates + j].i;
			            #else
                           Y[d->vecmap[tgt_i][cmap_col[tgt_i][i]] - lo_indx+j][k].r += receive_data[k][i*d->NBasisStates + j].r;
                           Y[d->vecmap[tgt_i][cmap_col[tgt_i][i]] - lo_indx+j][k].i += receive_data[k][i*d->NBasisStates + j].i;
			            #endif
			         #else
					    #ifdef SENDRECV_ONE_CHUNK
                        //Y[k][d->vecmap[tgt_i][cmap_col[tgt_i][i]] - lo_indx+j].r += receive_data[k*d->NBasisStates*d->commsize[tgt_i][src_i]+i*d->NBasisStates + j].r;
                        //Y[k][d->vecmap[tgt_i][cmap_col[tgt_i][i]] - lo_indx+j].i += receive_data[k*d->NBasisStates*d->commsize[tgt_i][src_i]+i*d->NBasisStates + j].i;
                        Y[k][d->vecmap[tgt_i][cmap_col[tgt_i][i]] - lo_indx+j].r += receive_data[k*d->NBasisStates*comm_data_size+i*d->NBasisStates + j].r;
                        Y[k][d->vecmap[tgt_i][cmap_col[tgt_i][i]] - lo_indx+j].i += receive_data[k*d->NBasisStates*comm_data_size+i*d->NBasisStates + j].i;
						#else
                        Y[k][d->vecmap[tgt_i][cmap_col[tgt_i][i]] - lo_indx+j].r += receive_data[k][i*d->NBasisStates + j].r;
                        Y[k][d->vecmap[tgt_i][cmap_col[tgt_i][i]] - lo_indx+j].i += receive_data[k][i*d->NBasisStates + j].i;
						#endif
					 #endif // USE_BLAS
            }}}
			//}
         }
      }
   #endif /* MPI3d */

   MPI_TOC(mpiTiming.Hv_tot, Tstart);
}


/* Computes Y=(H-Is)X where H is the Hamiltonian, s is input shift,
   and X is input matrix.  In this function X is a submatrix starting at isx 
   and ending at iex.  isx and iex are atomic (not actual) indices.  
   On output Y contains the [isx:isy] submatrix.  The rest of Y is contained 
   in YC. It is arranged as follows:
   YC[i*d-nb] = Y[ycmap[i][0]*d->NBasisStates]. 
*/
int Hmatmult_spds_col_fullstored_par ( cmatrix Y, cmatrix YC, qd_struct d, real shift, 
				                       cmatrix X, int num_col, int proc )
{
   MPI_TIME_INIT(Tstart);

   int atom, nbr, i, j;

   #ifdef USE_BLAS
      #ifdef USE_MKL
         MKL_Complex16 ALPHA;
	     MKL_Complex16 BETA;
	  #else
	     complex ALPHA;
		 complex BETA;
	  #endif

	  int M = d->NBasisStates;
      int K = d->NBasisStates;
      int N = num_col;
      int offsetx;
      int offsety;
	  int vec_len = M*N;
	  int zz,yy;

      char TRANSA = 'T';

	  #ifdef TRANSPOSED
         char TRANSB = 'N';
	  #else
         char TRANSB = 'T';
	  #endif

      #ifndef TRANSPOSED
	     static cvectr TMP = NULL;

         if(TMP == NULL)
	        TMP = Cvectr(vec_len);
	  #else
	     static cvectr TMPx = NULL;
		 static cvectr TMPy = NULL;

		 if(TMPx == NULL)
		    TMPx = Cvectr(vec_len);
		 if(TMPy == NULL)
		    TMPy = Cvectr(vec_len);
      #endif

      static cvectr HAMT = NULL;
	  if(HAMT == NULL) HAMT = Cvectr(M*K);
	  
	  #ifdef USE_MKL
         ALPHA.real = 1;
         ALPHA.imag = 0;
         BETA.real = 0;
         BETA.imag = 0;
	  #else
         ALPHA.r = 1;
         ALPHA.i = 0;
         BETA.r = 0;
         BETA.i = 0;
	  #endif // USE_MKL
   #endif // USE_BLAS

   if ( d->geo.hasPeriodicity() )
      die("Can't do periodic BC with (waste) full storage.\n");

   #ifdef USE_BLAS
      #ifdef TRANSPOSED
         for ( j = 0; j < num_col; j++)
            for ( i=0; i < d->nvmap[proc]*d->NBasisStates; i++ ) 
               YC[j][i].r = YC[j][i].i = 0.0;

         for ( j = 0; j < num_col; j++)
            for ( i=0; i < d->seg_ln[proc]; i++ )
               Y[j][i].r = Y[j][i].i = 0.0;
      #else
         for ( j = 0; j < num_col; j++)
            for ( i=0; i < d->nvmap[proc]*d->NBasisStates; i++ ) 
               YC[i][j].r = YC[i][j].i = 0.0;

         for ( j = 0; j < num_col; j++)
            for ( i=0; i < d->seg_ln[proc]; i++ )
               Y[i][j].r = Y[i][j].i = 0.0;
      #endif
   #else
      for ( j = 0; j < num_col; j++)
         for ( i=0; i < d->nvmap[proc]*d->NBasisStates; i++ ) 
            YC[j][i].r = YC[j][i].i = 0.0;

      for ( j = 0; j < num_col; j++)
         for ( i=0; i < d->seg_ln[proc]; i++ )
            Y[j][i].r = Y[j][i].i = 0.0;
   #endif // USE_BLAS

   MPI_TIC(Tstart);
 
   #undef _PRINT_HAM
   #ifdef  _PRINT_HAM
      printHam(d, "pos");
      printHam(d, "ham"); die("");
   #endif	

   if (d->NBasisStates==20) {
      for (atom=0;atom<vdim1(d->Ham);atom++) {
         for (nbr=0;nbr<=d->Ham_neighbor[atom];nbr++) {

		    #ifdef USE_BLAS
			   offsetx = d->Ham_offset_x[atom][nbr];
			   offsety = d->Ham_offset_y[atom][nbr];
			#endif // USE_BLAS

            if (d->Ham_yc[atom][nbr]){
               /* the outside elements (outside this CPU are unique and do not
                * have a Hermitian conjugate associated with them.... */
               #ifdef USE_BLAS
			      #ifdef TRANSPOSED
			         for(zz = 0; zz < num_col; zz++)
				        for(yy = 0; yy < M; yy++)
					       TMPx[zz*M+yy] = X[zz][offsetx+yy];
			      
				     #ifdef USE_MKL
                        zgemm(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, (MKL_Complex16 *)(d->Ham[atom][nbr]), 
					          &K, (MKL_Complex16 *)(TMPx), &K, &BETA, (MKL_Complex16 *)(TMPy), &M);
					 #else
                        zgemm(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, (d->Ham[atom][nbr]), 
					          &K, (TMPx), &K, &BETA, (TMPy), &M);
					 #endif

				     for(zz = 0; zz < num_col; zz++)
				        for(yy = 0; yy < M; yy++)
					       YC[zz][offsety+yy] = TMPy[zz*M+yy];
			      #else
				     #ifdef USE_MKL
			            // Transpose first
                        zgemm(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, (MKL_Complex16 *)(d->Ham[atom][nbr]), 
					       &K, (MKL_Complex16 *)(X[offsetx]), &N, &BETA, /*&(YC[offsety][0])*/(MKL_Complex16 *)(TMP), &M);
					 #else
                        zgemm(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, (d->Ham[atom][nbr]), 
					       &K, (X[offsetx]), &N, &BETA, (TMP), &M);
					 #endif // USE_MKL

				     // Transpose the result
				     for(zz = 0; zz < num_col; zz++)
				        for(yy = 0; yy < M; yy++)
					       YC[offsety+yy][zz] = TMP[zz*M+yy];
			      #endif
			   #else
			      for( j=0; j < num_col; j++)   
                     cmatmul_spds_20(YC[j], d->Ham[atom][nbr], X[j], 
                                  d->Ham_offset_y[atom][nbr], d->Ham_offset_x[atom][nbr], 
								  d->NBasisStates);
			   #endif // USE_BLAS
            }
	        else{
	           if (d->Ham_transp[atom][nbr]){
		          /* HC = cadjoint(HC,d->Ham[atom][nbr]); */
				  #ifdef USE_BLAS
				     TRANSA = 'N';
					 for(zz = 0; zz < M; zz++)
					    for(yy = 0; yy < K; yy++)
						{
						   HAMT[zz*K+yy].r = (d->Ham[atom][nbr][zz][yy].r); 
						   HAMT[zz*K+yy].i = -(d->Ham[atom][nbr][zz][yy].i); 
						}

				     #ifdef TRANSPOSED
			            for(zz = 0; zz < num_col; zz++)
				           for(yy = 0; yy < M; yy++)
					          TMPx[zz*M+yy] = X[zz][offsetx+yy];
			      
				        #ifdef USE_MKL
						   // transpose conjugate
                           zgemm(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, (MKL_Complex16 *)(HAMT), 
					          &K, (MKL_Complex16 *)(TMPx), &K, &BETA, (MKL_Complex16 *)(TMPy), &M);
						#else
                           zgemm(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, (HAMT), 
					          &K, (TMPx), &K, &BETA, (TMPy), &M);
						#endif // USE_MKL

				        for(zz = 0; zz < num_col; zz++)
				           for(yy = 0; yy < M; yy++)
					          Y[zz][offsety+yy] = TMPy[zz*M+yy];
				     #else
					    #ifdef USE_MKL
				           // Transpose conjugate
                           zgemm(&TRANSA,&TRANSB, &M, &N, &K, &ALPHA, (MKL_Complex16 *)(HAMT), 
					          &K, (MKL_Complex16 *)(X[offsetx]), &N, &BETA, /*&(Y[offsety][0])*/(MKL_Complex16 *)(TMP), &M);
						#else
                           zgemm(&TRANSA,&TRANSB, &M, &N, &K, &ALPHA, HAMT, 
					          &K, (X[offsetx]), &N, &BETA, (TMP), &M);
						#endif

				        // Transpose the result
				        for(zz = 0; zz < num_col; zz++)
				           for(yy = 0; yy < M; yy++)
					          Y[offsety+yy][zz] = TMP[zz*M+yy];
				     #endif
				  #else
				     for (j = 0; j < num_col; j++)
		                cmatmul_spds_hc_20(Y[j], d->Ham[atom][nbr], X[j], 
				                        d->Ham_offset_y[atom][nbr], 
				                        d->Ham_offset_x[atom][nbr], 
                                        d->NBasisStates);
				  #endif
	           } 
			   else{
			      #ifdef USE_BLAS
				     #ifdef TRANSPOSED
			            for(zz = 0; zz < num_col; zz++)
				           for(yy = 0; yy < M; yy++)
					          TMPx[zz*M+yy] = X[zz][offsetx+yy];
			      
				        #ifdef USE_MKL
                           zgemm(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, (MKL_Complex16 *)(d->Ham[atom][nbr]), 
					       &K, (MKL_Complex16 *)(TMPx), &K, &BETA, (MKL_Complex16 *)(TMPy), &M);
						#else
                           zgemm(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, (d->Ham[atom][nbr]), 
					       &K, (TMPx), &K, &BETA, (TMPy), &M);
						#endif

				        for(zz = 0; zz < num_col; zz++)
				           for(yy = 0; yy < M; yy++)
					          Y[zz][offsety+yy] = TMPy[zz*M+yy];
				     #else
				        // Transpose first
					    #ifdef USE_MKL
                           zgemm(&TRANSA,&TRANSB, &M, &N, &K, &ALPHA, (MKL_Complex16 *)(d->Ham[atom][nbr]), 
					       &K, (MKL_Complex16 *)(X[offsetx]), &N, &BETA, /*&(Y[offsety][0])*/(MKL_Complex16 *)(TMP), &M);
					    #else
                           zgemm(&TRANSA,&TRANSB, &M, &N, &K, &ALPHA, (d->Ham[atom][nbr]), 
					       &K, (X[offsetx]), &N, &BETA, (TMP), &M);
					    #endif

				        // Transpose the result
				        for(zz = 0; zz < num_col; zz++)
				           for(yy = 0; yy < M; yy++)
					          Y[offsety+yy][zz] = TMP[zz*M+yy];
				     #endif
				  #else
			         for (j = 0; j < num_col; j++)
		                cmatmul_spds_20(Y[j], d->Ham[atom][nbr], X[j], 
				                     d->Ham_offset_y[atom][nbr], 
				                     d->Ham_offset_x[atom][nbr], d->NBasisStates);
				  #endif // USE_BLAS
   }}}}}	           
   else if (d->NBasisStates==10) {
      for (atom=0;atom<vdim1(d->Ham);atom++) {
	     for (nbr=0;nbr<=d->Ham_neighbor[atom];nbr++) {
	        if (d->Ham_yc[atom][nbr]) {
	           /* the outside elements (outside this CPU are unique and do not 
		        * have a Hermitian conjugate associated with them.... */
			   #ifdef USE_BLAS

			   #else
			      for ( j = 0; j < num_col; j++)
	                 cmatmul_spds_10(YC[j], d->Ham[atom][nbr], X[j], 
                                  d->Ham_offset_y[atom][nbr], 
			                      d->Ham_offset_x[atom][nbr], 10);
			   #endif // USE_BLAS
	        }
	        else{
	           if (d->Ham_transp[atom][nbr]){
		          /* HC = cadjoint(HC,d->Ham[atom][nbr]); */
				  #ifdef USE_BLAS

				  #else
				     for ( j = 0; j < num_col; j++)
		                cmatmul_spds_hc_10(Y[j], d->Ham[atom][nbr], X[j], 
				                        d->Ham_offset_y[atom][nbr], 
				                        d->Ham_offset_x[atom][nbr], 10);
				  #endif // USE_BLAS
	           } else{
			      #ifdef USE_BLAS

				  #else
			         for ( j = 0; j < num_col; j++)
		                cmatmul_spds_10(Y[j], d->Ham[atom][nbr], X[j], 
				                     d->Ham_offset_y[atom][nbr], 
				                     d->Ham_offset_x[atom][nbr], 10);
				  #endif // USE_BLAS
   }}}}}   
   else if (d->NBasisStates==1) {
       	for (atom=0;atom<vdim1(d->Ham);atom++) {
	   		for (nbr=0;nbr<=d->Ham_neighbor[atom];nbr++) {
	       	   if (d->Ham_yc[atom][nbr]) {
		   			/* the outside elements (outside this CPU are unique and do not 
		      		 * have a Hermitian conjugate associated with them.... */
			      for ( j = 0; j < num_col; j++)
		   		     cmatmul_spds_1(YC[j], d->Ham[atom][nbr], X[j], 
                                    d->Ham_offset_y[atom][nbr],
                                    d->Ham_offset_x[atom][nbr], 1);
	       	   }
	       	   else{
		          if (d->Ham_transp[atom][nbr]){
		             /* HC = cadjoint(HC,d->Ham[atom][nbr]); */
					 for ( j = 0; j < num_col; j++)
		                cmatmul_spds_hc_1(Y[j], d->Ham[atom][nbr], X[j],
                                          d->Ham_offset_y[atom][nbr],
                                          d->Ham_offset_x[atom][nbr], 1);
		   	      } else{
				     for ( j = 0; j < num_col; j++)
		       			cmatmul_spds_1(Y[j], d->Ham[atom][nbr], X[j],
                                       d->Ham_offset_y[atom][nbr],
                                       d->Ham_offset_x[atom][nbr], 1);
   }}}}} 

   MPI_TOC(mpiTiming.matmul_spds, Tstart);

   if (shift!=0.0) {
      for ( j = 0; j < num_col; j++ ) {
         for ( i=0; i < d->seg_ln[proc]; i++ ) {
	        Y[j][i].r -= shift*X[j][i].r;
            Y[j][i].i -= shift*X[j][i].i;
   }}}      

   return d->nvmap[proc];
}


/* Computes Y=(H-Is)X where H is the Hamiltonian, shift is input shift,
   and X is input matrix.  In this function X is a submatrix starting at isx 
   and ending at iex.  isx and iex are atomic (not actual) indices.  
   On output Y contains the [isx:isy] submatrix.  The rest of Y is contained 
   in YC. It is arranged as follows:
   YC[i*d-nb] = Y[ycmap[i][0]*d->NBasisStates]. */
int Hmatmult_spds_col ( cmatrix Y, cmatrix YC, qd_struct d, real shift, 
                        cmatrix X, int num_col, int proc )
{
   static int counter=0;
   
   MPI_TIME_INIT(Tstart);

   int n, is0, isn, is1;

   #ifdef MPI3d
      int p;
   #endif	

   int Zo=0, Zu=0;
   int imin, imax, isyc, yccnt;
   complex eik;
   int Zatom=0,Zatom_m1, Zsurf=-1;
   int Znbr;
   rvectr nnv,nnv0;
   int done, rel_cindx, up_elements=0, down_elements=0;
   int AtomType_l_m;
   real a_lattice[3];
   
   bool has_periodicity = d->geo.hasPeriodicity();
   int isx = d->cell_s[proc];
   int iex = isx + d->cell_ln[proc];

   #ifdef    MPI_TIMING
      static bool firstrun = true;
      double Tsmall=0.0;
      double start_time_lanczos_iter_once=0;

      if (firstrun)   MPI_TIC(start_time_lanczos_iter_once);
   
      mpiTiming.h3ddiag_small=0.0;
      mpiTiming.h3doffd_small=0.0;
   #endif

   #ifdef USE_BLAS
      #ifdef USE_MKL
         MKL_Complex16 ALPHA;
	     MKL_Complex16 BETA;
	  #else
	     complex ALPHA;
		 complex BETA;
	  #endif

	  int M = d->NBasisStates;
      int K = d->NBasisStates;
      int N = num_col;
      int offsetx;
      int offsety;
	  int vec_len = M*N;
	  int zz,yy;

      char TRANSA = 'T';

	  #ifdef TRANSPOSED
         char TRANSB = 'N';
	  #else
         char TRANSB = 'T';
	  #endif

      #ifndef TRANSPOSED
	     static cvectr TMP = NULL;

         if(TMP == NULL)
	        TMP = Cvectr(vec_len);
	  #else
	     static cvectr TMPx = NULL;
		 static cvectr TMPy = NULL;

		 if(TMPx == NULL)
		    TMPx = Cvectr(vec_len);
		 if(TMPy == NULL)
		    TMPy = Cvectr(vec_len);
      #endif

	  #ifdef USE_MKL
         ALPHA.real = 1;
         ALPHA.imag = 0;
         BETA.real = 0;
         BETA.imag = 0;
	  #else
         ALPHA.r = 1;
         ALPHA.i = 0;
         BETA.r = 0;
         BETA.i = 0;
	  #endif // USE_MKL
   #endif // USE_BLAS

   nnv = Rvectr(3); 
   nnv0= Rvectr(3); 
   
   a_lattice[0] = d->geo.lattice_x;
   a_lattice[1] = d->geo.lattice_y;
   a_lattice[2] = d->geo.lattice_z;
   
   for (int i=0; i < d->seg_ln[proc]; i++)
      for (int j=0; j < num_col; j++)
         Y[j][i].r = Y[j][i].i = 0.0;

   /* Initialize segment boundary indices for processor proc.  
      isx and iex contain first and last atom assigned to processor 
      proc.  imin and imax contain the first and last row assigned 
      to processor proc.  There are d->NBasisStates rows for each 
      atom. */
   isyc = 0;
   yccnt = 0;
   imin = d->seg_s[proc];
   imax = d->seg_s[proc] + d->seg_ln[proc];
   
   #ifdef MPI3d
      for (int j=0; j < num_col; j++)
         for (int i=0; i < d->nvmap[proc]*d->NBasisStates; i++)
            YC[j][i].r = YC[j][i].i = 0.0;
   #endif

   if (d->opt.ExecParam.ElCalc.ElAlg.HamCompute==ElAlg_struct::SparseStore || 
       d->opt.ExecParam.ElCalc.ElAlg.HamCompute==ElAlg_struct::FullStore)
      Hmatmult_spds_col_prep_store (  Y,  YC, d,  shift,  X, num_col,  proc );

   /* Loop through all unit cells in the structure */
   for (int l=isx; l < iex; l++) {

      /* These are real space indexes for the unit cell */
      int i = d->geo.cell__ijk[l][0];
      int j = d->geo.cell__ijk[l][1];
      int k = d->geo.cell__ijk[l][2];
      
      /* Loop through all the atoms in the unit cell */
      for (int m=0; m < d->geo.AtomsPerCellMax(); m++) {
	     AtomType_l_m = (int) d->geo.AtomType[l][m];

	     /* Check if atom is not in simulation domain, skip */ 
	     if ( !AtomType_l_m ) continue;
         
         Zatom_m1=Zatom;
         Zatom++;

         /* Get starting index for the part of vector corresponding
            to this atom.  */
         is0 = d->geo.l_m_2_offset_ham_tot[l][m];
         is1 = is0 - imin;
         
         /* Hamiltonian diagonal (before strain correction).  
            If anion send anion basis vector, otherwise use cation vector . */
         H_Diag_Wrapper(d, l, m, i, j, k, AtomType_l_m, 
                (d->geo.isCation(l,m) ? d->BasisCation : d->BasisAnion) );

         if ( d->opt.ExecParam.Phys.PotentialFile != "NULL" ||
              d->opt.ExecParam.Phys.ImpNumber>0 || 
              d->opt.ExecParam.Phys.EfieldON) {
            for (int ii=0; ii<d->NBasisStates; ii++)
               d->hd[ii][ii].r -= d->phi[l][m];
         }

         /* Adjust the diagonal block of Hamiltonian for nonzero 
            magnetic field */
         if (d->opt.ExecParam.Phys.MagneticFieldOn) 
            apply_magnetic_field(d, l, m);

         /* Hamiltonian off-diagonal:  Here we loop through all this atom's
            nearest neighbors, get the coupling hamiltonians, and do the matrix
            vector multiplications. */
         for ( n = 0; n < d->geo.Neighbors(m); n++ ) {

            // Get relative position of unit cell where this neighbor lives.
            int i_nbr = i + d->geo.NbrCell(m,n,0);
            int j_nbr = j + d->geo.NbrCell(m,n,1);
            int k_nbr = k + d->geo.NbrCell(m,n,2);

            // Get the unit cell and atomic index for this neighbor 
            int aindx = d->geo.NbrCell(m,n,3);
            int cindx = d->geo.ijk__cell[i_nbr][j_nbr][k_nbr];
			
            /* Calculate absolute nearest neighbor vector to be used 
			 * in calculation of Peierl's phase due to Magnetic field. nnv_B is adjusted for strain later */
            double posAtom[3], posNbr[3], nnv_B[3];
            d->geo.getPosition_strained(posAtom,l,m);
            double* rNbr = d->geo.PositionInCell(aindx); /* position of this neighbor in the unit cell */
            posNbr[0] = d->geo.lattice_x * ( i_nbr + rNbr[0] );
            posNbr[1] = d->geo.lattice_y * ( j_nbr + rNbr[1] );
            posNbr[2] = d->geo.lattice_z * ( k_nbr + rNbr[2] );

            for(int ii=0;ii<3;ii++)
                nnv_B[ii] = posNbr[ii] - posAtom[ii];
            /* End: Calculate absolute nearest neighbor vector */

            #ifdef ELIMINATE_SSMAP
	           if (cindx < 0 && has_periodicity && ( cindx + d->geo.N_Cell + 1) >= 0) {
	              real phase = d->geo.phasePeriodic(d->kxL, d->kyL, d->kzL, i_nbr, j_nbr, k_nbr);
	      
	              eik.r = cos(phase);
	              eik.i = sin(phase);            
	              cindx = cindx + d->geo.N_Cell + 1 ;
	           }
                   if(d->geo.is_tilted) {
                      real phase = d->geo.phasePeriodic(d->kxL, d->kyL, d->kzL, i_nbr, j_nbr, k_nbr);
                      eik.r = cos(phase);
                      eik.i = sin(phase);
                      cindx = d->geo.get_cindx_tilted(l,m,n);
                   }

            #else /* ELIMINATE_SSMAP */
	           if (/*cindx < 0 &&*/ has_periodicity && d->geo.ssmap[l][m][n]!=-1) { 
			   /* When crystal is tilted ssmap should be used for connecting atoms in cell with index = 0 as well */ 
	      
	              real phase = d->geo.phasePeriodic(d->kxL, d->kyL, d->kzL, i_nbr, j_nbr, k_nbr);
	      
	              eik.r = cos(phase);
	              eik.i = sin(phase);            
	              cindx = d->geo.ssmap[l][m][n];
	           }
            #endif /* ELIMINATE_SSMAP */
            else {
               eik.r = 1.0;
               eik.i = 0.0;
            }
            
            // if nbr's cell is not in simulation domain, skip it
            if ( cindx<0 || !d->geo.AtomType[cindx][aindx])  continue;
            
            /* Modify nnv_B according to strain */
            for (int ii=0;ii<3;ii++)
               nnv_B[ii] += d->geo.l_m_2_atomdsp[cindx][aindx][ii];
            /* End: Modify nnv_B according to strain */

            // Get starting point on vector for this neigbhor.
            isn = d->geo.l_m_2_offset_ham_tot[cindx][aindx];
            
            bool out = false;

            #ifdef MPI3d
               // Check if neighbor is outside the x vector range 
               if ( isn >= imax || isn < imin ) {
                  // This is just a flag indicating that we are outside 
                  out = true;
               
                  // Set the offset isyc to the appropriate place 
                  for ( p=0; p < d->nvmap[proc]; p++ )
                     if ( isn == d->vecmap[proc][p] ) {
                        isyc = p*d->NBasisStates;
                        break;
                     } 
               }
            #endif
            
            H_DiagStrainCorr_OffDiag_Wrapper(d, l, m, AtomType_l_m, n, cindx, aindx, nnv, nnv0, a_lattice);
            /* Peierl's phase should be stored in memory so that it can be used later in matrix-vector multiplication. 
               Following lines were commented because phase was not stored. Additional code to store the phase is added 
               below at appropriate places */
            /* Adjust the offdiagonal block of Hamiltonian for nonzero 
               magnetic field */
            /*if(d->opt.ExecParam.Phys.MagneticFieldOn && 
               (d->HzbComplex.isInitialized() || 
                d->opt.ExecParam.ElCalc.ElAlg.HamCompute==ElAlg_struct::FullStore
				))  
               apply_magnetic_field(d, l, m, nnv);
            */
            /* End Peierl's phase should be stored */
            
            /* Matrix-vector multiply */
            if (out) {
               if (d->opt.ExecParam.ElCalc.ElAlg.HamCompute==ElAlg_struct::SparseStore){ 
                  int ii, jj;
                  int Nb_2;

                  if (d->BandModel==BM_10_sp3d5ss_nospin) Nb_2 = d->NBasisStates;
                  else                                    Nb_2 = d->NBasisStates / 2;

                  /* This matrix will be unique due to vector splitting 
                   * across different CPUs.  Therefore we do not have to 
                   * check whether we have encountered this matrix already.
                   */
                  d->Ham_neighbor[Zatom_m1]++;
                  Znbr = d->Ham_neighbor[Zatom_m1];

                  if (d->Ham_use_hermiticity){
                     d->Ham_atom_nbr_2_l[Zatom_m1][Znbr] = cindx-isx;
                     d->Ham_atom_nbr_2_m[Zatom_m1][Znbr] = aindx;
                  }

                  d->Ham_offset_y[Zatom_m1][Znbr] = isyc;
                  d->Ham_offset_x[Zatom_m1][Znbr] = is1;
                  d->Ham_yc[Zatom_m1][Znbr] = 1;
                  
                  /* set up off-diagonal superblocks corresponding to 
                   * (this proc atom, neighoring proc atom) pair */
                  if(d->Hzb.isInitialized()) {
                     for (ii=0; ii<Nb_2; ii++) {
                     for (jj=0; jj<Nb_2; jj++) {
                        int indx = Nb_2*(Nb_2*Zo + ii) + jj;
                        d->Hzb.Ho[indx] = d->ho[ii][jj].r;
                     }}
                  } else if(d->HzbComplex.isInitialized()){ 
                     /* Calculate and store Peierl's phase in HphaseOut/In */
                     complex dummy;

                     if(d->opt.ExecParam.Phys.MagneticFieldOn) {
                        dummy = vector_potential_phase(d, l, m, nnv_B);
                        d->HzbComplex.HphaseOut[Zo].r = (float) dummy.r;
                        d->HzbComplex.HphaseOut[Zo].i = (float) dummy.i;
                     }

                     /* End: Calculate and store Peierl's phase in HphaseOut/In */
                     for (ii=0; ii<Nb_2; ii++) {
                        for (jj=0; jj<Nb_2; jj++) {
                           int indx = Nb_2*(Nb_2*Zo + ii) + jj;

                           #ifdef  HamMemorySave
                              d->HzbComplex.Ho[indx].r = d->ho[ii][jj].r;
                              d->HzbComplex.Ho[indx].i = d->ho[ii][jj].i;
                           #else
                              d->HzbComplex.Ho[indx] = d->ho[ii][jj];
                           #endif
                     }}
                  }
                  
                  if(d->Hzb.isInitialized()) {
                     if(d->opt.ExecParam.Phys.MagneticFieldOn){
                        #ifdef  HamMemorySave
                           complex dummy = vector_potential_phase(d, l, m, nnv);
                           d->Hzb.HphaseOut[Zo].r = (float) dummy.r;
                           d->Hzb.HphaseOut[Zo].i = (float) dummy.i;
                        #else
                           d->Hzb.HphaseOut[Zo] = vector_potential_phase(d, l, m, nnv);
                        #endif
                     }

                     d->Hzb.indxHO_to_row[Zo] = isyc;
                     d->Hzb.indxHO_to_col[Zo] = d->NBasisStates*Zatom_m1;
                     d->Hzb.indxHO_nbr[Zo] = n;
                  }
                  else if(d->HzbComplex.isInitialized()) {
                     d->HzbComplex.indxHO_to_row[Zo] = isyc;
                     d->HzbComplex.indxHO_to_col[Zo] = d->NBasisStates*Zatom_m1;
                     d->HzbComplex.indxHO_nbr[Zo] = n;
                  } 	
				  
                  #ifdef _PRINT_MAPPING_INFO
                     if(d->Hzb.isInitialized())
                        printf("proc %d -- indx %d -- row %d -- col %d -- nbr %d -- %d -- (lm) %d %d\n", 
                            mpi_n3d_id, Zo,
                            d->Hzb.indxHO_to_row[Zo], d->Hzb.indxHO_to_col[Zo], 
                            d->Hzb.indxHO_nbr[Zo],
                            d->geo.NbrCell
                            (d->Atom_to_lm[ d->Hzb.indxHO_to_col[Zo] / d->NBasisStates ][1],
                             d->Hzb.indxHO_nbr[Zo],
                             3),
                            cindx, aindx
                            );
                  #endif

                  Zo++;
                  up_elements++;
               } 
               else if (d->opt.ExecParam.ElCalc.ElAlg.HamCompute==ElAlg_struct::FullStore){ 
                  /* This matrix will be unique due to vector splitting 
                   * across different CPUs.  Therefore we do not have to 
                   * check whether we have encountered this matrix already.
                   */
                  d->Ham_neighbor[Zatom_m1]++;
                  Znbr = d->Ham_neighbor[Zatom_m1];

                  if (d->Ham_use_hermiticity){
                     d->Ham_atom_nbr_2_l[Zatom_m1][Znbr] = cindx-isx;
                     d->Ham_atom_nbr_2_m[Zatom_m1][Znbr] = aindx;
                  }

                  d->Ham[Zatom_m1][Znbr] = cp_cmatrix(NULL,d->ho);
                  
                  d->Ham_offset_y[Zatom_m1][Znbr] = isyc;
                  d->Ham_offset_x[Zatom_m1][Znbr] = is1;
                  d->Ham_yc[Zatom_m1][Znbr] = 1;
                  up_elements++;
               }
               
               MPI_TIC(Tstart);
               
                  if (eik.i!=0.0 || eik.r!=1.0)
				  {
				     #ifdef USE_BLAS
				        #ifdef TRANSPOSED
			               for(zz = 0; zz < num_col; zz++)
				              for(yy = 0; yy < M; yy++)
					             TMPx[zz*M+yy] = X[zz][is1+yy];
			      
				           #ifdef USE_MKL
						      ALPHA.real = eik.r;
							  ALPHA.imag = eik.i;

                              zgemm(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, (MKL_Complex16 *)(d->ho), 
					          &K, (MKL_Complex16 *)(TMPx), &K, &BETA, (MKL_Complex16 *)(TMPy), &M);
						   #else
                              zgemm(&TRANSA, &TRANSB, &M, &N, &K, &eik, (d->ho), 
					          &K, (TMPx), &K, &BETA, (TMPy), &M);
						   #endif

				           for(zz = 0; zz < num_col; zz++)
				              for(yy = 0; yy < M; yy++)
					             YC[zz][isyc+yy] = TMPy[zz*M+yy];
				        #else
				           // Transpose first
					       #ifdef USE_MKL
						      ALPHA.real = eik.r;
							  ALPHA.imag = eik.i;

                              zgemm(&TRANSA,&TRANSB, &M, &N, &K, &ALPHA, (MKL_Complex16 *)(d->ho), 
					          &K, (MKL_Complex16 *)(X[is1]), &N, &BETA, (MKL_Complex16 *)(TMP), &M);
					       #else
                              zgemm(&TRANSA,&TRANSB, &M, &N, &K, &eik, (d->ho), 
					          &K, (X[is1]), &N, &BETA, (TMP), &M);
					       #endif

				           // Transpose the result
				           for(zz = 0; zz < num_col; zz++)
				              for(yy = 0; yy < M; yy++)
					             YC[isyc+yy][zz] = TMP[zz*M+yy];
				        #endif
					 #else
				        for ( int j = 0; j < num_col; j++)
                           cmatmul_spds_sc(YC[j], d->ho, X[j], isyc, is1, d->NBasisStates, eik);
				     #endif
				  }
                  else
				  {
				     #ifdef USE_BLAS
				        #ifdef TRANSPOSED
			               for(zz = 0; zz < num_col; zz++)
				              for(yy = 0; yy < M; yy++)
					             TMPx[zz*M+yy] = X[zz][is1+yy];
			      
				           #ifdef USE_MKL
                              zgemm(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, (MKL_Complex16 *)(d->ho), 
					          &K, (MKL_Complex16 *)(TMPx), &K, &BETA, (MKL_Complex16 *)(TMPy), &M);
						   #else
                              zgemm(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, (d->ho), 
					          &K, (TMPx), &K, &BETA, (TMPy), &M);
						   #endif

				           for(zz = 0; zz < num_col; zz++)
				              for(yy = 0; yy < M; yy++)
					             YC[zz][isyc+yy] = TMPy[zz*M+yy];
				        #else
				           // Transpose first
					       #ifdef USE_MKL
                              zgemm(&TRANSA,&TRANSB, &M, &N, &K, &ALPHA, (MKL_Complex16 *)(d->ho), 
					          &K, (MKL_Complex16 *)(X[is1]), &N, &BETA, (MKL_Complex16 *)(TMP), &M);
					       #else
                              zgemm(&TRANSA,&TRANSB, &M, &N, &K, &ALPHA, (d->ho), 
					          &K, (X[is1]), &N, &BETA, (TMP), &M);
					       #endif

				           // Transpose the result
				           for(zz = 0; zz < num_col; zz++)
				              for(yy = 0; yy < M; yy++)
					             YC[isyc+yy][zz] = TMP[zz*M+yy];
				        #endif
					 #else
				        for ( int j = 0; j < num_col; j++)
                           cmatmul_spds(YC[j], d->ho, X[j], isyc, is1, d->NBasisStates); 
					 #endif
				  }
               
               MPI_TOC(Tsmall, Tstart);
            }  // we are not out
            else {
               if (d->opt.ExecParam.ElCalc.ElAlg.HamCompute==ElAlg_struct::SparseStore) {
                  int ii, jj;
                  int Nb_2;

                  if (d->BandModel==BM_10_sp3d5ss_nospin) Nb_2 = d->NBasisStates;
                  else                                    Nb_2 = d->NBasisStates / 2;

                  done = 0;

                  d->Ham_neighbor[Zatom_m1]++;
                  Znbr = d->Ham_neighbor[Zatom_m1];
                  rel_cindx = cindx-isx;

                  /* See whether we have encountered this neigbor pair already */
                  /* What is the (l,m) index of this neighbor? (cindx,aindex)  */
                  if (d->Ham_use_hermiticity &&
                     (rel_cindx>=0) && (rel_cindx<iex-isx) && (aindx>=0) &&
                     (Zatom > d->l_m_2_atom[rel_cindx][aindx])   ){
                     /* We should have computed this off-diagonal matrix 
                      * already. Find it in the list of previous atoms.  */
                     {
                        int i_nbr,old_atom=d->l_m_2_atom[rel_cindx][aindx];

                        for (i_nbr=1;i_nbr<=d->Ham_neighbor[old_atom];i_nbr++){
                           if (d->Ham_atom_nbr_2_l[old_atom][i_nbr]==(l-isx) &&
                               d->Ham_atom_nbr_2_m[old_atom][i_nbr]==m){
                              /* found the right neigbor */
                              done = 1;
                              d->Ham_transp[Zatom_m1][Znbr] = 1;
                              d->Ham_offset_y[Zatom_m1][Znbr] = isn-imin;
                              d->Ham_offset_x[Zatom_m1][Znbr] = is1;
                              
                              down_elements++;
                              break;
                           }
                        }
                     }
                  } 

                  if (done==0){
                     if (d->Ham_use_hermiticity){
                        d->Ham_atom_nbr_2_l[Zatom_m1][Znbr] = rel_cindx;
                        d->Ham_atom_nbr_2_m[Zatom_m1][Znbr] = aindx;
                     }

                     d->Ham_offset_y[Zatom_m1][Znbr] = isn-imin;
                     d->Ham_offset_x[Zatom_m1][Znbr] = is1;
                     
                     if(d->Hzb.isInitialized()) {
                        for (ii=0; ii<Nb_2; ii++) {
                           for (jj=0; jj<Nb_2; jj++) {
                              int indx = Nb_2*(Nb_2*Zu + ii) + jj;
                              d->Hzb.Hu[indx] = d->ho[ii][jj].r;
                        }}
                     } else if(d->HzbComplex.isInitialized()) {
                        /* Calculate and store Peierl's phase in HphaseOut/In */
                        complex dummy;

                        if(d->opt.ExecParam.Phys.MagneticFieldOn) {
                           dummy = vector_potential_phase(d, l, m, nnv_B);
                           d->HzbComplex.HphaseIn[Zu].r = (float) dummy.r;
                           d->HzbComplex.HphaseIn[Zu].i = (float) dummy.i;
                        }
                        /* End: Calculate and store Peierl's phase in HphaseOut/In */
						
                        for (ii=0; ii<Nb_2; ii++) {
                           for (jj=0; jj<Nb_2; jj++) {
                              int indx = Nb_2*(Nb_2*Zu + ii) + jj;

                              #ifdef  HamMemorySave
                                 d->HzbComplex.Hu[indx].r = d->ho[ii][jj].r;
                                 d->HzbComplex.Hu[indx].i = d->ho[ii][jj].i;
                              #else
                                 d->HzbComplex.Hu[indx] = d->ho[ii][jj];
                              #endif
                        }}
                     }
                     
                     if(d->Hzb.isInitialized()) {
                        if(d->opt.ExecParam.Phys.MagneticFieldOn){
                           #ifdef  HamMemorySave
                              complex dummy = vector_potential_phase(d, l, m, nnv);
                              d->Hzb.HphaseIn[Zu].r = (float) dummy.r;
                              d->Hzb.HphaseIn[Zu].i = (float) dummy.i;
                           #else
                              d->Hzb.HphaseIn[Zu] = vector_potential_phase(d, l, m, nnv);
                           #endif
                        }

                        d->Hzb.indxHU_to_row[Zu] = isn-imin;
                        d->Hzb.indxHU_to_col[Zu] = d->NBasisStates*Zatom_m1;
                        d->Hzb.indxHU_nbr[Zu] = n;
                     }
                     else if(d->HzbComplex.isInitialized()) {
                        d->HzbComplex.indxHU_to_row[Zu] = isn-imin;
                        d->HzbComplex.indxHU_to_col[Zu] = d->NBasisStates*Zatom_m1;
                        d->HzbComplex.indxHU_nbr[Zu] = n;
                     }

                     Zu++;			      
                     up_elements++;
                  }
               }
               else if (d->opt.ExecParam.ElCalc.ElAlg.HamCompute==ElAlg_struct::FullStore){
                  done = 0;
                  d->Ham_neighbor[Zatom_m1]++;
                  Znbr = d->Ham_neighbor[Zatom_m1];
                  rel_cindx = cindx-isx;
                  /* See whether we have encountered this neigbor pair already*/
                  /* What is the (l,m) index of this neighbor? (cindx,aindex)  */
                  if (d->Ham_use_hermiticity &&
                     (rel_cindx>=0) && (rel_cindx<iex-isx) && (aindx>=0) &&
                     (Zatom > d->l_m_2_atom[rel_cindx][aindx])   ){
                     /* We should have computed this off-diagonal matrix 
                      * already. Find it in the list of previous atoms.  */
                     {
                        int i_nbr,old_atom=d->l_m_2_atom[rel_cindx][aindx];

                        for (i_nbr=1;i_nbr<=d->Ham_neighbor[old_atom];i_nbr++){
                           if (d->Ham_atom_nbr_2_l[old_atom][i_nbr]==(l-isx) &&
                               d->Ham_atom_nbr_2_m[old_atom][i_nbr]==m){
                              /* found the right neigbor */
                              /* printf("Found another matrix element normal\n");  */
                              done = 1;
                              d->Ham[Zatom_m1][Znbr] = d->Ham[old_atom][i_nbr];
                              d->Ham_transp[Zatom_m1][Znbr] = 1;
                              d->Ham_offset_y[Zatom_m1][Znbr] = isn-imin;
                              d->Ham_offset_x[Zatom_m1][Znbr] = is1;
                              down_elements++;
                              break;
                           }
                        }
                     }
                  } 
                  if (done==0){
                     if (d->Ham_use_hermiticity){
                        d->Ham_atom_nbr_2_l[Zatom_m1][Znbr] = rel_cindx;
                        d->Ham_atom_nbr_2_m[Zatom_m1][Znbr] = aindx;
                     }

                     d->Ham[Zatom_m1][Znbr] = cp_cmatrix(NULL,d->ho);
                     d->Ham_offset_y[Zatom_m1][Znbr] = isn-imin;
                     d->Ham_offset_x[Zatom_m1][Znbr] = is1;
                     up_elements++;
                  }
               }
               
               MPI_TIC(Tstart);
               
                  if (eik.i!=0.0 || eik.r!=1.0)
				  {
				     #ifdef USE_BLAS
				        #ifdef TRANSPOSED
			               for(zz = 0; zz < num_col; zz++)
				              for(yy = 0; yy < M; yy++)
					             TMPx[zz*M+yy] = X[zz][is1+yy];
			      
				           #ifdef USE_MKL
						      ALPHA.real = eik.r;
							  ALPHA.imag = eik.i;

                              zgemm(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, (MKL_Complex16 *)(d->ho), 
					          &K, (MKL_Complex16 *)(TMPx), &K, &BETA, (MKL_Complex16 *)(TMPy), &M);
						   #else
                              zgemm(&TRANSA, &TRANSB, &M, &N, &K, &eik, (d->ho), 
					          &K, (TMPx), &K, &BETA, (TMPy), &M);
						   #endif

				           for(zz = 0; zz < num_col; zz++)
				              for(yy = 0; yy < M; yy++)
					             Y[zz][isn-imin+yy] = TMPy[zz*M+yy];
				        #else
				           // Transpose first
					       #ifdef USE_MKL
						      ALPHA.real = eik.r;
							  ALPHA.imag = eik.i;

                              zgemm(&TRANSA,&TRANSB, &M, &N, &K, &ALPHA, (MKL_Complex16 *)(d->ho), 
					          &K, (MKL_Complex16 *)(X[is1]), &N, &BETA, (MKL_Complex16 *)(TMP), &M);
					       #else
                              zgemm(&TRANSA,&TRANSB, &M, &N, &K, &eik, (d->ho), 
					          &K, (X[is1]), &N, &BETA, (TMP), &M);
					       #endif

				           // Transpose the result
				           for(zz = 0; zz < num_col; zz++)
				              for(yy = 0; yy < M; yy++)
					             Y[isn-imin+yy][zz] = TMP[zz*M+yy];
				        #endif
					 #else
				        for (int j = 0; j < num_col; j++)
                           cmatmul_spds_sc( Y[j], d->ho, X[j], isn-imin, is1, d->NBasisStates, eik );
					 #endif
				  }
                  else
			      {
				     #ifdef USE_BLAS
				        #ifdef TRANSPOSED
			               for(zz = 0; zz < num_col; zz++)
				              for(yy = 0; yy < M; yy++)
					             TMPx[zz*M+yy] = X[zz][is1+yy];
			      
				           #ifdef USE_MKL
                              zgemm(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, (MKL_Complex16 *)(d->ho), 
					          &K, (MKL_Complex16 *)(TMPx), &K, &BETA, (MKL_Complex16 *)(TMPy), &M);
						   #else
                              zgemm(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, (d->ho), 
					          &K, (TMPx), &K, &BETA, (TMPy), &M);
						   #endif

				           for(zz = 0; zz < num_col; zz++)
				              for(yy = 0; yy < M; yy++)
					             Y[zz][isn-imin+yy] = TMPy[zz*M+yy];
				        #else
				           // Transpose first
					       #ifdef USE_MKL
                              zgemm(&TRANSA,&TRANSB, &M, &N, &K, &ALPHA, (MKL_Complex16 *)(d->ho), 
					          &K, (MKL_Complex16 *)(X[is1]), &N, &BETA, (MKL_Complex16 *)(TMP), &M);
					       #else
                              zgemm(&TRANSA,&TRANSB, &M, &N, &K, &ALPHA, (d->ho), 
					          &K, (X[is1]), &N, &BETA, (TMP), &M);
					       #endif

				           // Transpose the result
				           for(zz = 0; zz < num_col; zz++)
				              for(yy = 0; yy < M; yy++)
					             Y[isn-imin+yy][zz] = TMP[zz*M+yy];
				        #endif
					 #else
			            #ifndef REDUCE_RECOMPUTE
					       for (int j = 0; j < num_col; j++)
                              cmatmul_spds( Y[j], d->ho, X[j], isn-imin, is1, d->NBasisStates ); 
				        #else
					       for (int j = 0; j < num_col; j++)
                              cmatmul_spds_offdiag( Y[j], d->ho, X[j], isn-imin, is1, d->NBasisStates ); 
				        #endif // REDUCE_RECOMPUTE
					 #endif
			      }

         
            #ifdef _DEBUGGING_HAM
               char dum[100];

               sprintf(dum,"Hdir_%d",counter);

               FILE* fp = fopen(dum,"a");
               for (int ii=0; ii<20; ii++)
                  for (int jj=0; jj<20; jj++)
                     fprintf(fp,"%4d %4d %e %e\n", 
                       1+20*d->geo.NbrCell(m,n,3)+ii, 
                       1+20*m+jj, 
                       d->ho[ii][jj].r, d->ho[ii][jj].i);
               fclose(fp);
            #endif
         
               MPI_TOC(Tsmall, Tstart);
            }
	     }
         
         if (d->opt.ExecParam.ElCalc.ElAlg.HamCompute==ElAlg_struct::SparseStore) {
            d->Ham_offset_x[Zatom_m1][0] = is1;
            d->Ham_offset_y[Zatom_m1][0] = is1;
            
            int I = d->Ham_offset_y[Zatom_m1][0];
            int J = d->NBasisStates*Zatom_m1;
            
            if( d->opt.ExecParam.Phys.SurfacePassivation == 
                Phys_struct::PassivateDanglingBond &&
                d->geo.isOnSurface(l,m) ) {
               Zsurf++;
               d->surfaceAtoms[Zsurf] = Zatom_m1;

               if(d->Hzb.isInitialized()) 
                  d->Hzb.setDiag(d->hd, Zatom_m1, I, J, Zsurf);
               else if(d->HzbComplex.isInitialized()) 
                  d->HzbComplex.setDiag(d->hd, Zatom_m1, I, J, Zsurf);
            }
            else {
               if(d->Hzb.isInitialized()) 
                  d->Hzb.setDiag(d->hd, Zatom_m1, I, J, -1);
               else if(d->HzbComplex.isInitialized()) 
                  d->HzbComplex.setDiag(d->hd, Zatom_m1, I, J, -1);
            }
         } else if (d->opt.ExecParam.ElCalc.ElAlg.HamCompute==ElAlg_struct::FullStore){
            d->Ham[Zatom_m1][0] = cp_cmatrix(NULL,d->hd);
            d->Ham_offset_x[Zatom_m1][0] = is1;
            d->Ham_offset_y[Zatom_m1][0] = is1;
         }
         
         #ifdef _DEBUGGING_HAM
            char dum[100];

         sprintf(dum,"Hdir_%d",counter);

         FILE* fp = fopen(dum,"a");

         for (int ii=0; ii<20; ii++)
            for (int jj=0; jj<20; jj++)
               fprintf(fp,"%4d %4d %e %e\n", 
                    1+20*m+ii, 
                    1+20*m+jj, 
                    d->hd[ii][jj].r, d->hd[ii][jj].i);
         fclose(fp);
         #endif

         MPI_TIC(Tstart);
         
		 #ifdef USE_BLAS
		    #ifdef TRANSPOSED
			   for(zz = 0; zz < num_col; zz++)
				  for(yy = 0; yy < M; yy++)
					 TMPx[zz*M+yy] = X[zz][is1+yy];
			      
			   #ifdef USE_MKL
                  zgemm(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, (MKL_Complex16 *)(d->hd), 
					          &K, (MKL_Complex16 *)(TMPx), &K, &BETA, (MKL_Complex16 *)(TMPy), &M);
	           #else
                  zgemm(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, (d->hd), 
					          &K, (TMPx), &K, &BETA, (TMPy), &M);
			   #endif

			   for(zz = 0; zz < num_col; zz++)
				  for(yy = 0; yy < M; yy++)
					 Y[zz][is1+yy] = TMPy[zz*M+yy];
			#else
			   // Transpose first
			   #ifdef USE_MKL
                  zgemm(&TRANSA,&TRANSB, &M, &N, &K, &ALPHA, (MKL_Complex16 *)(d->hd), 
					          &K, (MKL_Complex16 *)(X[is1]), &N, &BETA, (MKL_Complex16 *)(TMP), &M);
			   #else
                  zgemm(&TRANSA,&TRANSB, &M, &N, &K, &ALPHA, (d->hd), 
					          &K, (X[is1]), &N, &BETA, (TMP), &M);
			   #endif

			   // Transpose the result
			   for(zz = 0; zz < num_col; zz++)
				  for(yy = 0; yy < M; yy++)
			         Y[is1+yy][zz] = TMP[zz*M+yy];
			#endif
		 #else
            #ifndef REDUCE_RECOMPUTE
			   for(int j=0; j < num_col; j++)
                  cmatmul_spds( Y[j], d->hd, X[j], is1, is1, d->NBasisStates ); 
	        #else
			   for(int j=0; j < num_col; j++)
                  cmatmul_spds_diag( Y[j], d->hd, X[j], is1, is1, d->NBasisStates ); 
		    #endif // REDUCE_RECOMPUTE
         #endif

         MPI_TOC(Tsmall, Tstart);
      }
   }

   rm_rvectr(&nnv);
   rm_rvectr(&nnv0);
   
   #ifdef    MPI_TIMING
      mpiTiming.matmul_spds += Tsmall;
      mpiTiming.h3ddiag     += mpiTiming.h3ddiag_small;
      mpiTiming.h3doffd     += mpiTiming.h3doffd_small;

      if (firstrun) {
         MPI_TOC(mpiTiming.lanczos_setup_once, start_time_lanczos_iter_once);
         firstrun=false;
      }
   #endif
   
   if (d->opt.ExecParam.ElCalc.ElAlg.HamCompute==ElAlg_struct::SparseStore || 
       d->opt.ExecParam.ElCalc.ElAlg.HamCompute==ElAlg_struct::FullStore) {
      if (mpi_n3d_id==0 || mpi_n3d_id==mpi_n3d_numprocs-1)
	     printf("CPU %d  up_diag %d  low_diag %d\n",proc,up_elements,down_elements),fflush(stdout);
   }
   
   if(shift != 0.0) {
      for (int i=0; i < d->seg_ln[proc]; i++) {
         for (int j=0; j < num_col; j++) {
            Y[j][i].r -= shift*X[j][i].r;
            Y[j][i].i -= shift*X[j][i].i;
   }}}
   
   counter++;

   #ifdef MPI3d
      if (d->opt.ExecParam.Phys.Disorder1!=0.0){
         Randomize_H (d, proc);
         d->nvmap[proc] = Hmatmult_spds_col_fullstored_par(Y, YC, d, shift, X, num_col,  proc);
      }

      return d->nvmap[proc];
   #else
      return 0;
   #endif       
}


int Hmatmult_spds_col_prep_store ( cmatrix Y, cmatrix YC, qd_struct d, real shift, cmatrix X, int num_col, int proc )
{
   int l, m;
   int isx, iex;
   int isyc, yccnt;
   
   int Zatom=0, Zatom_surf=0;
   int INT_LARGE = (int) pow(2.0,30.0);
   
   const int maxNbr = d->geo.NeighborsMax();

   if (d->Ham || d->Hzb.isInitialized() || d->HzbComplex.isInitialized()) 
      return 0;
   
   nml_memory_report("before Allocation of Hamiltonian Storage - Hmatmult_spds_col_prep_store ");

   /* Initialize segment boundary indices for processor proc.  
      isx and iex contain first and last atom assigned to processor proc.  
      imin and imax contain the first and last row assigned to processor proc.
      There are d->NBasisStates rows for each atom. */
   isx = d->cell_s[proc];
   iex = isx + d->cell_ln[proc];
   isyc = 0;
   yccnt = 0;
   
   if (d->Ham_use_hermiticity){

      d->l_m_2_atom = Imatrix(iex-isx,d->geo.AtomsPerCellMax());

      for (l=0;l<iex-isx;l++){
	     for (m=0;m<d->geo.AtomsPerCellMax();m++){
	        d->l_m_2_atom[l][m] = INT_LARGE;
	  }}
   }
   
   //  Count the number of real atoms and number on surface
   for ( l=isx; l < iex; l++ ) {
      for ( m=0; m < d->geo.AtomsPerCellMax(); m++ ) {
         if ( !d->geo.AtomType[l][m] )  continue;
      
         if (d->Ham_use_hermiticity) 
            d->l_m_2_atom[l-isx][m] = Zatom;
         Zatom++;

         if (d->geo.isOnSurface(l,m))
            Zatom_surf++;
   }}

   d->Ham_offset_x = Imatrix(Zatom,maxNbr+1);

   /*
   double memory_estimate=0;
   memory_estimate = (Zatom * (maxNbr+1) *3 * sizeof(int) ) +   // Imatrix stuff) 
                     (Zatom * 3 * sizeof(int) ) ;               // pointers in matrix 
   if (d->Ham_use_hermiticity){
      memory_estimate *= 2.0;
   }
   printf("mpi_id=%d memory estimate Hamiltonian storage prep: %g MB Natom=%d \n", 
           mpi_n3d_id, memory_estimate/1.0e6,Zatom); fflush(stdout);
   */

   d->Ham_offset_y = Imatrix(Zatom,maxNbr+1);
   d->Ham_yc       = Imatrix(Zatom,maxNbr+1);
   d->Ham_neighbor = Ivectr(Zatom);
   if (d->Ham_use_hermiticity){
      d->Ham_transp       = Imatrix(Zatom,maxNbr+1);
      d->Ham_atom_nbr_2_l = Imatrix(Zatom,maxNbr+1);
      d->Ham_atom_nbr_2_m = Imatrix(Zatom,maxNbr+1);
   }

   if (d->opt.ExecParam.ElCalc.ElAlg.HamCompute==ElAlg_struct::FullStore) {
      d->Ham = Cmatrixmatrix(Zatom,maxNbr+1);
   }
   else if (d->opt.ExecParam.ElCalc.ElAlg.HamCompute==ElAlg_struct::SparseStore) {
      int Nout, Nin;

      count_nbrs_for_sparse_zincblende_ham(&Nout, &Nin, d);
      
      complex HBxy;
      double bohr_magneton=0.00005788381749;  /* in unit of eV/T */ 

      HBxy.r = d->opt.ExecParam.Phys.Bx*bohr_magneton;
      HBxy.i = -d->opt.ExecParam.Phys.By*bohr_magneton;
	  
      /* Following lines were commented for fixing bugs in Magnetic field */
      /*
      if(!d->opt.ExecParam.Phys.MagneticFieldOn)
         d->Hzb.Initialize(d->BandModel,d->NBasisStates, Zatom, Zatom_surf, Nout, Nin, false, HBxy);
      else if(d->opt.ExecParam.ElCalc.ElAlg.Efficiency==ElAlg_struct::Memory)
         d->Hzb.Initialize(d->BandModel,d->NBasisStates, Zatom, Zatom_surf, Nout, Nin, true,  HBxy);
      else if(d->opt.ExecParam.ElCalc.ElAlg.Efficiency==ElAlg_struct::Speed)
         d->HzbComplex.Initialize(d->BandModel,d->NBasisStates, Zatom, Zatom_surf, Nout, Nin, HBxy);
      */
      /* End: Following lines were commented for fixing bugs in Magnetic field */

      /* Following new conditions were addeded for fixing bugs in Magnetic field */
      if(!d->opt.ExecParam.Phys.MagneticFieldOn && d->opt.ExecParam.ElCalc.ElAlg.Efficiency==ElAlg_struct::Memory)
         d->Hzb.Initialize(d->BandModel,d->NBasisStates, Zatom, Zatom_surf, Nout, Nin, false, HBxy);
      else if(d->opt.ExecParam.Phys.MagneticFieldOn && d->opt.ExecParam.ElCalc.ElAlg.Efficiency==ElAlg_struct::Memory)
         d->Hzb.Initialize(d->BandModel,d->NBasisStates, Zatom, Zatom_surf, Nout, Nin, true,  HBxy);
      if(!d->opt.ExecParam.Phys.MagneticFieldOn && d->opt.ExecParam.ElCalc.ElAlg.Efficiency==ElAlg_struct::Speed)
         d->HzbComplex.Initialize(d->BandModel,d->NBasisStates, Zatom, Zatom_surf, Nout, Nin,
                                  false, // Added for Magnetic field (Value of MagneticFieldOn) 
								  HBxy);
      if(d->opt.ExecParam.Phys.MagneticFieldOn && d->opt.ExecParam.ElCalc.ElAlg.Efficiency==ElAlg_struct::Speed)
         d->HzbComplex.Initialize(d->BandModel,d->NBasisStates, Zatom, Zatom_surf, Nout, Nin,
                                  true, // Added for Magnetic field (Value of MagneticFieldOn)
                                  HBxy);
      /* End: Following new conditions were addeded for fixing bugs in Magnetic field */

      d->surfaceAtoms = Ivectr(Zatom_surf);
   }
   
   /* print storage info to stdout */
   if (mpi_n3d_id==mpi_n3d_masterid) {
      real storage_main=0.0, storage_help=0.0;

      if (d->opt.ExecParam.ElCalc.ElAlg.HamCompute==ElAlg_struct::FullStore) {
	     if (d->Ham_use_hermiticity){ 
	        storage_main = 16.0*Zatom*(maxNbr/2+1)*d->NBasisStates*d->NBasisStates/1000000.0;
	        storage_help = (6.0*4.0*Zatom*(maxNbr+1) + 4.0*(iex-isx)*maxNbr)/1000000.0;
	     }
	     else {
	        storage_main = 16.0*Zatom*(maxNbr+1)*d->NBasisStates*d->NBasisStates/1000000.0;
	        storage_help = 3.0*4.0*Zatom*(maxNbr+1)/1000000.0;
	     }

         masterPrint("CPU %d, atom_count=%d\n",proc,Zatom);
	     /*
         masterPrint("CPU %d:  Estimated Hamiltonian Storage %.1fMB (main)  %.1fMB(help)\n",
                mpi_n3d_id, storage_main, storage_help);
	     */
      }
      else {
         /* outdated estimators.  Have details in the beginning of the Initializers\n"
         if(d->Hzb.isInitialized()) d->Hzb.printStorageInfo();
         else if(d->HzbComplex.isInitialized()) d->HzbComplex.printStorageInfo(); 
	     */
      }
   }
   
   fflush(stdout);

   return 0;
}

