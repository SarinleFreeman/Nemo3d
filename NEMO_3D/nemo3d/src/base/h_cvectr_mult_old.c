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
$Header: /repo/nemo3d/src/base/h_cvectr_mult_old.c,v 1.2 2007/08/10 15:03:32 baeh Exp $
*****************************************************************************/

#undef PRINT_DEBUG_COM
//#define PRINT_DEBUG_COM
#ifdef PRINT_DEBUG_COM
#define print_com2(A,B) "printf(A,B), fflush(stdout);"
#define print_com3(A,B,C) "printf(A,B,C), fflush(stdout);"
#else
#define print_com2(A,B) ;
#define print_com3(A,B,C) ;
#endif

#include "h_cvectr_mult.h"

/* Include matmul directly into this file for inlining optimization. */
#ifndef MACOSX
#include "matmul.c"
#endif


void compute_ritz_residual_norm_deviation(qd_struct d)
{
  cvectr Hx = Cvectr(d->seg_ln[mpi_n3d_id]);
  for(int iv=0; iv<d->neigv; iv++) {
    Hmult_spds_complete(Hx, d, 0.0, &d->wfc[iv][0]);
    complex xHx = vect_dot_vect(&d->wfc[iv][0], Hx, d);
    complex  xx = vect_dot_vect(&d->wfc[iv][0], &d->wfc[iv][0], d);
    real ritz_value = xHx.r/xx.r;

    double E = d->WF[iv].E;
    Hmult_spds_complete(Hx, d, E, &d->wfc[iv][0]);
    complex HxHx = vect_dot_vect(Hx, Hx, d);
    real residual = sqrt(HxHx.r*HxHx.r+HxHx.i*HxHx.i);

    xx = vect_dot_vect(&d->wfc[iv][0], &d->wfc[iv][0], d);
    real norm = sqrt(xx.r*xx.r + xx.i*xx.i);

    real deviation = fabs(E-ritz_value);
    if (E!=0.0)
	    deviation/= E;

    masterPrint("state %d Eigval %-e Ritzval %-e residual %-e (Eigval-Ritz)/Eigval %-e norm %-e" , iv, E,ritz_value,residual,deviation,norm);
  }
  masterPrint(" ");
}

void compute_ritz_value(qd_struct d)
{
  cvectr Hx = Cvectr(d->seg_ln[mpi_n3d_id]);
  for(int iv=0; iv<d->neigv; iv++) {
    Hmult_spds_complete(Hx, d, 0.0, &d->wfc[iv][0]);
    complex xHx = vect_dot_vect(&d->wfc[iv][0], Hx, d);
    complex  xx = vect_dot_vect(&d->wfc[iv][0], &d->wfc[iv][0], d);
    real ritz_value = xHx.r/xx.r;
    masterPrint("Ritz values of state %d is %e", iv, ritz_value);
  }
  masterPrint(" ");
}

void compute_residual(qd_struct d)
{
  cvectr Hx = Cvectr(d->seg_ln[mpi_n3d_id]);
  for(int iv=0; iv<d->neigv; iv++) {
    double E = d->WF[iv].E;
    Hmult_spds_complete(Hx, d, E, &d->wfc[iv][0]);
    complex HxHx = vect_dot_vect(Hx, Hx, d);
    real residual = sqrt(HxHx.r*HxHx.r+HxHx.i+HxHx.i);
    masterPrint("Residual of state %d is %e", iv, residual);
  }  
  masterPrint(" ");
}


void compute_eigenvector_norm (qd_struct d)
{
  for(int iv=0; iv<d->neigv; iv++) {
    complex xx = vect_dot_vect(&d->wfc[iv][0], &d->wfc[iv][0], d);
    real norm = sqrt(xx.r*xx.r + xx.i*xx.i);
    masterPrint("Norm of state %d is %e", iv, norm);
  } 
  masterPrint(" ");
}

complex vect_dot_vect (complex* v1, complex* v2, qd_struct d)
{
  complex loc_product;
  loc_product.r=0.0; loc_product.i=0.0;
  for(int i=0; i<d->seg_ln[mpi_n3d_id]; i++) {
    loc_product.r+=(v1[i].r*v2[i].r+v1[i].i*v2[i].i);
    loc_product.i+=(v1[i].r*v2[i].i-v1[i].i*v2[i].r);
  }
  complex tot_product;
#if (defined MPI3d && !defined FAKE_MPI)
  MPI_Allreduce(&loc_product, &tot_product, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  tot_product = loc_product;
#endif
  return tot_product;
}

void matvect_mult2(const int** argList, nml_dcscalar* y, const nml_dcscalar* x)
{
   qd_struct d = *((qd_struct*) argList[0]);
   double shift = *((double*) argList[1]);
   int Nproj = *((int*) argList[2]);
   int i;

#if 0
   double shift1 = 1.48-0.1;
   double shift2 = 1.48+0.1;
   double factor = -1.0;
   double offset = 1400.0;
   if (!d->cvec_aux){
     d->cvec_aux = Cvectr(d->seg_ln[mpi_n3d_id]);
   }

   Hmult_spds_complete((complex*) d->cvec_aux, d, shift1, (complex*) x);
   Hmult_spds_complete((complex*) y, d, shift2, (complex*) d->cvec_aux);
   for(int i=0; i<d->seg_ln[mpi_n3d_id]; i++) {
     y[i].r = (y[i].r*factor + x[i].r*offset);     
     y[i].i = (y[i].i*factor + x[i].i*offset);
   }
#endif

#if 1
   double shift1 = 1.48-0.1;
   double shift2 = 1.48+0.1;
   double factor = -10.0;
   double offset = 0.0;
   if (!d->cvec_aux){
     d->cvec_aux = Cvectr(d->seg_ln[mpi_n3d_id]);
   }

   Hmult_spds_complete((complex*) d->cvec_aux, d, shift1, (complex*) x);
   Hmult_spds_complete((complex*) y, d, shift2, (complex*) d->cvec_aux);
   for(int i=0; i<d->seg_ln[mpi_n3d_id]; i++) {
     y[i].r = (y[i].r*factor + x[i].r*offset);     
     y[i].i = (y[i].i*factor + x[i].i*offset);
   }
#endif

   return;
}

 

/* NOTE: argList is an array of pointers structured as follows:
   argList[0] = ptr to qd_struct d;
   argList[1] = ptr to double shift;
*/
/* FIX!!!:  we assume that all pointers are all the same size */
void matvect_mult(const int** argList, nml_dcscalar* y, const nml_dcscalar* x)
{
   qd_struct d = *((qd_struct*) argList[0]);
   double shift = *((double*) argList[1]);
   int Nproj = *((int*) argList[2]);

   Hmult_spds_complete((complex*) y, d, shift, (complex*) x);

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
}



/* Compute y = (H - shift)*x st that each processor gets its appropriate
   portion of y and no further communication is required. */
void Hmult_spds_complete( cvectr y, qd_struct d, real shift, cvectr x )
{
   real Tstart=0.0;
   static int n_multiply =0;
   MPI_TIC(Tstart);
   if (0 == n_multiply%100 && mpi_n3d_id== 0)
     fprintf(stdout,"# of matrix-vector multiplies: %d\n",n_multiply),fflush(stdout);
   n_multiply++;
   
   
#define NEED_TO_USE_STORED_MATRIX_ALWAYS
#ifndef MPI3d  // serial version doesn't need to do any fancy communication
   
   cvectr yc=NULL; // a dummy placeholder that won't be used in matvect mult
   
   if (d->opt.ExecParam.ElCalc.ElAlg.HamCompute==ElAlg_struct::SparseStore && 
       (d->Hzb.isInitialized() || d->HzbComplex.isInitialized()))
      Hmult_spds_col_stored( y, yc, d, shift, x, 0 ); 
   else if (d->opt.ExecParam.ElCalc.ElAlg.HamCompute==ElAlg_struct::FullStore && d->Ham)
      Hmult_spds_col_fullstored_par( y, yc, d, shift, x, 0 ); 
   else{
      Hmult_spds_col( y, yc, d, shift, x, 0 );
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
             Hmult_spds_col_stored( y, yc, d, shift, x, 0 ); 
          else if (d->opt.ExecParam.ElCalc.ElAlg.HamCompute==ElAlg_struct::FullStore && d->Ham)
             Hmult_spds_col_fullstored_par( y, yc, d, shift, x, 0 ); 
          else{
             die("ERROR in Hmult_spds_complete:\n"
	         "Unidentified Hamiltonian storage system\n");
	  }
       }
#endif
   }

   
#else // The parallel version needs to set up the communication arrays
   
   static int send_max=0;
   static ivectr sflag=NULL;
   static cvectr receive_data=NULL;
   static cvectr send_data=NULL;
   static imatrix cmap_row, cmap_col;

   if(d->rc==NULL) d->rc = Cvectr(d->nvmap[mpi_n3d_id]*d->NBasisStates); 
   cvectr yc=d->rc;
   cvectr source=NULL, target=NULL;
   int i, j, cnt, src_i, tgt_i, maxcom, lo_indx, commstep;
   int nyc = 0;
   print_com2("CPU=%d start Hmult_spds_complete\n",mpi_n3d_id);
   /* Allocate and initialize static stuff if required */

   if ( sflag == NULL ) {
      print_com2("CPU=%d initialize static vectors in Hmult_spds_complete\n",
                 mpi_n3d_id);
      maxcom = 0;
      for ( i=0; i < d->nproc; i++ ) {
         for ( j=0; j < d->nproc; j++ ) {
            if ( d->commsize[i][j] > maxcom )
               maxcom = d->commsize[i][j];
         }}
      
      cmap_row = Imatrix(d->nproc,maxcom);
      cmap_col = Imatrix(d->nproc,maxcom);
      
      maxcom *= d->NBasisStates;
      
      for (tgt_i=0; tgt_i < d->nproc; tgt_i++){
	 if (d->commsize[mpi_n3d_id][tgt_i]>0){
            cnt = 0;
	    for ( i=0; i < d->nvmap[mpi_n3d_id]; i++ ) {
	       if ( d->vecmap[mpi_n3d_id][i] >= d->seg_s[tgt_i] && 
		    d->vecmap[mpi_n3d_id][i] < d->seg_s[tgt_i]+d->seg_ln[tgt_i] ){
		  cmap_row[tgt_i][cnt] = i;
		  cnt++;
	       }
	    }
	 }
      }
      
      for (tgt_i=0; tgt_i < d->nproc; tgt_i++){
	 if (d->commsize[mpi_n3d_id][tgt_i]>0){
            cnt = 0;
	    for ( i=0; i < d->nvmap[tgt_i]; i++ ){
	       if ( d->vecmap[tgt_i][i] >= d->seg_s[mpi_n3d_id] && 
		    d->vecmap[tgt_i][i] < d->seg_s[mpi_n3d_id]+d->seg_ln[mpi_n3d_id] ){
		  cmap_col[tgt_i][cnt] = i;
		  cnt++;
	       }
	    }
	 }
      }
      
      receive_data = Cvectr(maxcom);
      print_com2("CPU=%d did receive_data\n",mpi_n3d_id);
      send_data = Cvectr(maxcom);
      print_com2("CPU=%d did send_data\n",mpi_n3d_id);
   }

  //      static int eventcount=0;
   
   if (d->opt.ExecParam.ElCalc.ElAlg.HamCompute==ElAlg_struct::SparseStore && 
       (d->Hzb.isInitialized() || d->HzbComplex.isInitialized())){
       //eventcount++;
       //nyc = Hmult_spds_col( y, yc, d, shift, x, mpi_n3d_id );
       //{
             //printf("*****************original*********\n");
           //for (i=0;i<nml_fcv_extent(y);i++){
             //printf("%d %g %g\n",i,y[i].r,y[i].i);
           //}
           //fflush(stdout);
        //}
      nyc = Hmult_spds_col_stored( y, yc, d, shift, x, mpi_n3d_id ); 
        //{
             //printf("*****************modified*********\n");
           //for (i=0;i<nml_fcv_extent(y);i++){
             //printf("%d %g %g\n",i,y[i].r,y[i].i);
           //}
           //fflush(stdout);
           //if (eventcount>=2)
        //   exit(-1);
        //}
   }
   else if (d->opt.ExecParam.ElCalc.ElAlg.HamCompute==ElAlg_struct::FullStore && d->Ham)
      nyc = Hmult_spds_col_fullstored_par( y, yc, d, shift, x, mpi_n3d_id ); 
   else{
      nyc = Hmult_spds_col( y, yc, d, shift, x, mpi_n3d_id );
#ifdef NEED_TO_USE_STORED_MATRIX_ALWAYS
       //eventcount++;
        //{
             //printf("*****************original*********\n");
           //for (i=0;i<nml_fcv_extent(y);i++){
             //printf("%d %g %g\n",i,y[i].r,y[i].i);
           //}
           //fflush(stdout);
        //}
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
             Hmult_spds_col_stored( y, yc, d, shift, x, mpi_n3d_id ); 
          else if (d->opt.ExecParam.ElCalc.ElAlg.HamCompute==ElAlg_struct::FullStore && d->Ham)
             Hmult_spds_col_fullstored_par( y, yc, d, shift, x, mpi_n3d_id ); 
          else{
             die("ERROR in Hmult_spds_complete:\n"
	         "Unidentified Hamiltonian storage system\n");
	  }
       }
        //{
             //printf("*****************modified*********\n");
           //for (i=0;i<nml_fcv_extent(y);i++){
             //printf("%d %g %g\n",i,y[i].r,y[i].i);
           //}
           //fflush(stdout);
           //if (eventcount>=1)
           //exit(-1);
        //}
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
	       }
	    }
	 }
      }
   }
   
   // NB: see old version of code for old sendreceive
   for ( commstep=0; commstep < 2; commstep++ ) {
      tgt_i = d->commorder[src_i][commstep];
      if ( tgt_i != -1 ) {
	 source = receive_data;
	 send_data=cvectrClear(send_data);
	 target = send_data;     
	 for ( i=0; i < d->commsize[src_i][tgt_i]; i++ )
	    for ( j=0; j < d->NBasisStates; j++ )
	       target[i*d->NBasisStates+j] = 
                  yc[cmap_row[tgt_i][i]*d->NBasisStates+j];   
	 MPI_Sendrecv(target, 
                      2*d->commsize[src_i][tgt_i]*d->NBasisStates, 
		      MPI_DOUBLE,tgt_i, MPI_rc_distr, 
		      source, 
                      2*d->commsize[tgt_i][src_i]*d->NBasisStates, 
		      MPI_DOUBLE,tgt_i, MPI_rc_distr, MPI_COMM_WORLD, 
                      &mpi_n3d_status);
	 /* add the scattered packages to the main array */
	 target = y;	
	 lo_indx=d->seg_s[src_i];
	 for ( i=0; i < d->commsize[tgt_i][src_i]; i++ )
	 for ( j=0; j < d->NBasisStates; j++ ) {
            target[d->vecmap[tgt_i][cmap_col[tgt_i][i]] - lo_indx+j].r += 
               source[i*d->NBasisStates + j].r;
            target[d->vecmap[tgt_i][cmap_col[tgt_i][i]] - lo_indx+j].i += 
               source[i*d->NBasisStates + j].i;
         }}
      /* MPI_Barrier(MPI_COMM_WORLD); */
   }
   
   /* Receive packages and add them in as they show up */
   tgt_i = mpi_n3d_id;
   for ( src_i=0; src_i < d->nproc; src_i++ ){
      /* assume symmetry in the original sflag matrix */
      if (src_i!=tgt_i && sflag[src_i]){
      }
   }
   
#endif /* MPI3d */
   MPI_TOC(mpiTiming.Hv_tot, Tstart);
}


/* Computes y=(H-Is)x where H is the Hamiltonian, s is input shift,
   and x is input vector.  In this function x is a subvector starting at isx 
   and ending at iex.  isx and iex are atomic (not actual) indices.  
   On output y contains the [isx:isy] subvector.  The rest of y is contained 
   in yc. It is arranged as follows:
   yc[i*d-nb] = y[ycmap[i][0]*d->NBasisStates]. 
*/
int Hmult_spds_col_fullstored_par ( cvectr y, cvectr yc, qd_struct d, real shift, 
				   cvectr x, int proc )
{
   MPI_TIME_INIT(Tstart);

   int atom, nbr, i;

   /* cmatrix HC=NULL; */
   /* printf("Entered Hmult_spds_col_fullstored_par\n"); */
   /* if (d->Ham_use_hermiticity)
      HC=Coperator(d->NBasisStates); */

   if ( d->geo.hasPeriodicity() )
      die("Can't do periodic BC with (waste) full storage.\n");

   for ( i=0; i < d->nvmap[proc]*d->NBasisStates; i++ )
      yc[i].r = yc[i].i = 0.0;

   for ( i=0; i < d->seg_ln[proc]; i++ ) {
      y[i].r = y[i].i = 0.0;
   }

   MPI_TIC(Tstart);
 
#undef _PRINT_HAM
#ifdef  _PRINT_HAM
   printHam(d, "pos");
   printHam(d, "ham"); die("");
#endif	

   if (d->NBasisStates==20){
      for (atom=0;atom<vdim1(d->Ham);atom++) {
	 for (nbr=0;nbr<=d->Ham_neighbor[atom];nbr++) {
	    if (d->Ham_yc[atom][nbr]){
	       /* the outside elements (outside this CPU are unique and do not
                  have a Hermitian conjugate associated with them.... */
	       cmatmul_spds_20(yc, d->Ham[atom][nbr], x, 
                               d->Ham_offset_y[atom][nbr], 
			       d->Ham_offset_x[atom][nbr], d->NBasisStates);
	    }
	    else{
	       if (d->Ham_transp[atom][nbr]){
		  /* HC = cadjoint(HC,d->Ham[atom][nbr]); */
		  cmatmul_spds_hc_20(y, d->Ham[atom][nbr], x, 
				     d->Ham_offset_y[atom][nbr], 
				     d->Ham_offset_x[atom][nbr], 
                                     d->NBasisStates);
	       } else{
		  cmatmul_spds_20(y, d->Ham[atom][nbr], x, 
				  d->Ham_offset_y[atom][nbr], 
				  d->Ham_offset_x[atom][nbr], d->NBasisStates);
	       }
	    }
	 }
      }
   } else if (d->NBasisStates==10) {
      for (atom=0;atom<vdim1(d->Ham);atom++) {
	 for (nbr=0;nbr<=d->Ham_neighbor[atom];nbr++) {
	    if (d->Ham_yc[atom][nbr]) {
	       /* the outside elements (outside this CPU are unique and do not 
		  have a Hermitian conjugate associated with them.... */
	       cmatmul_spds_10(yc, d->Ham[atom][nbr], x, 
                               d->Ham_offset_y[atom][nbr], 
			       d->Ham_offset_x[atom][nbr], 10);
	    }
	    else{
	       if (d->Ham_transp[atom][nbr]){
		  /* HC = cadjoint(HC,d->Ham[atom][nbr]); */
		  cmatmul_spds_hc_10(y, d->Ham[atom][nbr], x, 
				     d->Ham_offset_y[atom][nbr], 
				     d->Ham_offset_x[atom][nbr], 10);
	       } else{
		  cmatmul_spds_10(y, d->Ham[atom][nbr], x, 
				  d->Ham_offset_y[atom][nbr], 
				  d->Ham_offset_x[atom][nbr], 10);
	       }
	    }
	 }
      }
   } else if (d->NBasisStates==1) {

       for (atom=0;atom<vdim1(d->Ham);atom++) {
	   for (nbr=0;nbr<=d->Ham_neighbor[atom];nbr++) {
	       if (d->Ham_yc[atom][nbr]) {
		   /* the outside elements (outside this CPU are unique and do not 
		      h	ave a Hermitian conjugate associated with them.... */
		   cmatmul_spds_1(yc, d->Ham[atom][nbr], x, 
                                  d->Ham_offset_y[atom][nbr],
                                  d->Ham_offset_x[atom][nbr], 1);
	       }
	       else{
		   if (d->Ham_transp[atom][nbr]){
		       /* HC = cadjoint(HC,d->Ham[atom][nbr]); */
		       cmatmul_spds_hc_1(y, d->Ham[atom][nbr], x,
                                         d->Ham_offset_y[atom][nbr],
                                         d->Ham_offset_x[atom][nbr], 1);
		   } else{
		       cmatmul_spds_1(y, d->Ham[atom][nbr], x,
                                      d->Ham_offset_y[atom][nbr],
                                      d->Ham_offset_x[atom][nbr], 1);
		   }
	       }
	   }
       }
       
   } else {
       for (atom=0;atom<vdim1(d->Ham);atom++){
	   for (nbr=0;nbr<=d->Ham_neighbor[atom];nbr++){
	       if (d->Ham_yc[atom][nbr]){
		   /* the outside elements (outside this CPU are unique and 
                      do not have a Hermitian conjugate associated with them */
		   cmatmul_spds(yc, d->Ham[atom][nbr], x, 
                                d->Ham_offset_y[atom][nbr], 
				d->Ham_offset_x[atom][nbr], 
                                d->NBasisStates);
	       }
	       else{
		   if (d->Ham_transp[atom][nbr]){
		       /* HC	 = cadjoint(HC,d->Ham[atom][nbr]); */
		       cmatmul_spds_hc(y, d->Ham[atom][nbr], x, 
				       d->Ham_offset_y[atom][nbr], 
				       d->Ham_offset_x[atom][nbr], 
                                       d->NBasisStates);
		   } else{
		       cmatmul_spds(y, d->Ham[atom][nbr], x, 
				    d->Ham_offset_y[atom][nbr], 
				    d->Ham_offset_x[atom][nbr], 
                                    d->NBasisStates);
		   }
	       }
	   }
       }
   }

   MPI_TOC(mpiTiming.matmul_spds, Tstart);

   if (shift!=0.0){
       for ( i=0; i < d->seg_ln[proc]; i++ ) {
	   /* y[i].r = y[i].i = 0.0 */
	   y[i].r -= shift*x[i].r;
	   y[i].i -= shift*x[i].i;
       }
   }
   return d->nvmap[proc];
}

// replacement for H_diag_cation H_diag_anion
void H_Diag_Wrapper(qd_struct d, int l, int m, int i, int j, int k, int AtomType_l_m,
            ivectr Basis) 
{
   MPI_TIME_INIT(Tstart);

   bool _symS3 = d->geo.hasPeriodicity();/* d->geo.isPeriodicXYZ();  Looks like this was a bug.
                                            This boolean should be true for 1D and 2D periodic 
                                            crystals as well */

   // I don't think this averaging is necessary for the onsite terms  -- F.O.
   /* we actually do since the As in GaAs and InAs might have different onsite matrix elements */
   /* We need a weighted average of the orbital energies */
   for (int p=0; p < sp3d5s_nparam; p++ ) d->param[p]=0;

   real strainShift=0.0;
   int nbrcnt=0;
   for (int n=0; n < d->geo.Neighbors(m); n++ ) {
      int id = d->geo.NbrCell(m,n,0);
      int jd = d->geo.NbrCell(m,n,1);
      int kd = d->geo.NbrCell(m,n,2);
      int aindx = d->geo.NbrCell(m,n,3);
      int cindx = d->geo.ijk__cell[i+id][j+jd][k+kd];

#ifdef ELIMINATE_SSMAP
      if (cindx < 0 && _symS3 && ( cindx + d->geo.N_Cell + 1) >= 0) {
         cindx = cindx + d->geo.N_Cell + 1 ;
      }
#else /* ELIMINATE_SSMAP */
      if (cindx < 0 && _symS3 && d->geo.ssmap[l][m][n]!=-1) {
         cindx = d->geo.ssmap[l][m][n];
      }
#endif /* ELIMINATE_SSMAP */

      if ( cindx >= 0 ) {
         /* The new data storage can crash on the dereferencing of the casted data */
         /* It appears to crash on the new gcc compilers                   */
         /* const int& AtomType_nbr = (int) d->geo.AtomType[cindx][aindx]; */
         int AtomType_nbr = (int) d->geo.AtomType[cindx][aindx];
	 if ( !AtomType_nbr )  continue;

         int cation, anion;
         if (d->geo.isCation(l,m)) {
            cation = AtomType_l_m;
            anion  = AtomType_nbr;
#ifdef STRAIN_SHIFT_ENABLE
            if (d->opt.ExecParam.ElCalc.include_strain_in_Hamiltonian)
               strainShift += d->strnmat[cation][anion][eta_cat_shift];
#endif
         }
         else {
            cation = AtomType_nbr;
            anion  = AtomType_l_m;
#ifdef STRAIN_SHIFT_ENABLE
            if (d->opt.ExecParam.ElCalc.include_strain_in_Hamiltonian)
               strainShift += d->strnmat[cation][anion][eta_an_shift];
#endif
         }

         rvectr param_local = d->parmat[cation][anion];
         for (int p=0; p < sp3d5s_nparam; p++)
            d->param[p] += param_local[p];
         nbrcnt++;
      }
   }

   strainShift /= nbrcnt;

   for (int p=0; p < sp3d5s_nparam; p++ ) 
      d->param[p] /= nbrcnt;

   for (int p=0; p < sp3d5s_nparam; p++ ) {
      //printf("%d %d(%s) param[%d]=%f:  ", m, AtomType_l_m, (m<4?"C":"A"), p, d->param[p]);
      for (int n=0; n < d->geo.Neighbors(m); n++ ) {
         int id = d->geo.NbrCell(m,n,0);
         int jd = d->geo.NbrCell(m,n,1);
         int kd = d->geo.NbrCell(m,n,2);
         int aindx = d->geo.NbrCell(m,n,3);
         int cindx = d->geo.ijk__cell[i+id][j+jd][k+kd];

#ifdef ELIMINATE_SSMAP
	 if (cindx < 0 && _symS3 && ( cindx + d->geo.N_Cell + 1) >= 0) {
	   cindx = cindx + d->geo.N_Cell + 1 ;
	 }
#else /* ELIMINATE_SSMAP */
	 if (cindx < 0 && _symS3 && d->geo.ssmap[l][m][n]!=-1) {
	   cindx = d->geo.ssmap[l][m][n];
	 }
#endif /* ELIMINATE_SSMAP */



         if ( cindx >= 0 ) {
            /* The new data storage can crash on the dereferencing of the casted data */
            /* It appears to crash on the new gcc compilers                   */
            /* const int& AtomType_nbr = (int) d->geo.AtomType[cindx][aindx]; */
            int AtomType_nbr = (int) d->geo.AtomType[cindx][aindx];
            if ( !AtomType_nbr )  continue;

            int cation, anion;
            if (d->geo.isCation(l,m)) {
               cation = AtomType_l_m;
               anion  = AtomType_nbr;
            }
            else {
               cation = AtomType_nbr;
               anion  = AtomType_l_m;
            }
            //printf("%f  (%d,%d)   ", d->parmat[cation][anion][p], AtomType_l_m,AtomType_nbr);
         }
      }
      //printf("\n");
   }
   
   MPI_TIC(Tstart);
   h3d_diag(d->hd, Basis, d->Basis_Spin, d->param, d->NBasisStates );
   MPI_TOC(mpiTiming.h3ddiag_small, Tstart);

   if (d->geo.isCellWurtzite()){
      double pz_addition_energy = -0.04;
      if(pz_addition_energy !=0.0)
       h3d_diag_pz_addition_for_wurtzite_crystal_splitting(d, pz_addition_energy);
   }

   if (d->geo.isOnSurface(l,m)) {
      if( d->opt.ExecParam.Phys.SurfacePassivation == 
          Phys_struct::PassivateOrbital )
         h3d_diag_surface_atom_shift(d);
      if( d->opt.ExecParam.Phys.SurfacePassivation == 
          Phys_struct::PassivateDanglingBond )
         h3d_diag_dangling_bond_shift(d, l, m);
   }

#ifdef STRAIN_SHIFT_ENABLE			       
   for (int ii=0;ii<d->NBasisStates;ii++){
      d->hd[ii][ii].r += strainShift;
   }
#endif

#ifdef HAMILTONIAN_SCALE
   if (HAMILTONIAN_SCALE != 1.0) {
      int i1,i2;
      for (i1=0;i1<d->NBasisStates;i1++){
      for (i2=0;i2<d->NBasisStates;i2++){
         d->hd[i1][i2].r /= HAMILTONIAN_SCALE_VALUE;
      }}
   }
#endif

}


void H_DiagStrainCorr_OffDiag_Wrapper(qd_struct d, int l, int m, int AtomType_l_m, 
                            int n, int cindx, int aindx, rvectr nnv, 
                            rvectr nnv0, real a_lattice[3]) 
{
   MPI_TIME_INIT(Tstart);

   /* The new data storage can crash on the dereferencing of the casted data */
   /* It appears to crash on the new gcc compilers                   */
   /* const int& AtomType_nbr = (int) d->geo.AtomType[cindx][aindx]; */
   int AtomType_nbr ;
   real a_lattice_unstrnd;

   AtomType_nbr      = (int) d->geo.AtomType[cindx][aindx];
   a_lattice_unstrnd = d->parmat[AtomType_l_m][AtomType_nbr][pV_unstr_latt];

   d->geo.getRelativePosVec(nnv0, nnv, a_lattice_unstrnd, cindx, aindx, l, m, n);

   ivectr basis_this, basis_nbr;
   rvectr mb_strain=NULL, param;
   if ( d->geo.isCation(l,m) ) {
      basis_this = d->BasisCation;
      basis_nbr  = d->BasisAnion;
      param = d->parmat[AtomType_l_m][AtomType_nbr];
      if (d->opt.ExecParam.ElCalc.include_strain_in_Hamiltonian)
         mb_strain = d->strnmat[AtomType_l_m][AtomType_nbr];

      for (int p=0; p < 3; p++) {
         nnv[p] *= -1;
         nnv0[p] *= -1;
      }
   }
   else {
      basis_this = d->BasisAnion;
      basis_nbr  = d->BasisCation;
      param = d->parmat[AtomType_nbr][AtomType_l_m];
      if (d->opt.ExecParam.ElCalc.include_strain_in_Hamiltonian)
         mb_strain = d->strnmat[AtomType_nbr][AtomType_l_m];
   }

   /* Generate off-diagonal Hamiltonian */
   /* returns correct d->ho (original contents destroyed), 
      d->hd (new) = d->hd (old) + update */
   MPI_TIC(Tstart);

   h3d_offdiag_diag_corr(d->ho, d->hd, nnv, nnv0, 
                         basis_this, basis_nbr, 
                         d->Basis_Spin, 
                         param, 
                         mb_strain,
                         a_lattice_unstrnd,
                         d->NBasisStates );
   MPI_TOC(mpiTiming.h3doffd_small, Tstart);

#ifdef HAMILTONIAN_SCALE
   if (HAMILTONIAN_SCALE != 1.0) {
      int i1,i2;
      for (i1=0;i1<d->NBasisStates;i1++){
      for (i2=0;i2<d->NBasisStates;i2++){
         d->h[i1][i2].r /= HAMILTONIAN_SCALE_VALUE;
      }}
   }
#endif
}


/* Computes y=(H-Is)x where H is the Hamiltonian, shift is input shift,
   and x is input vector.  In this function x is a subvector starting at isx 
   and ending at iex.  isx and iex are atomic (not actual) indices.  
   On output y contains the [isx:isy] subvector.  The rest of y is contained 
   in yc. It is arranged as follows:
   yc[i*d-nb] = y[ycmap[i][0]*d->NBasisStates]. */
int Hmult_spds_col ( cvectr y, cvectr yc, qd_struct d, real shift, cvectr x, 
			  int proc )
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

   nnv = Rvectr(3); 
   nnv0= Rvectr(3); 
   
   a_lattice[0] = d->geo.lattice_x;
   a_lattice[1] = d->geo.lattice_y;
   a_lattice[2] = d->geo.lattice_z;
   
   for (int i=0; i < d->seg_ln[proc]; i++)
      y[i].r = y[i].i = 0.0;

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
   for (int i=0; i < d->nvmap[proc]*d->NBasisStates; i++)
      yc[i].r = yc[i].i = 0.0;
#endif
   
   if (d->opt.ExecParam.ElCalc.ElAlg.HamCompute==ElAlg_struct::SparseStore || 
       d->opt.ExecParam.ElCalc.ElAlg.HamCompute==ElAlg_struct::FullStore)
      Hmult_spds_col_prep_store (  y,  yc, d,  shift,  x,  proc );

   /* Loop through all unit cells in the structure */
   for (int l=isx; l < iex; l++) {

      /* These are real space indexes for the unit cell */
      int i = d->geo.cell__ijk[l][0];
      int j = d->geo.cell__ijk[l][1];
      int k = d->geo.cell__ijk[l][2];
      
      /* Loop through all the atoms in the unit cell */
      for (int m=0; m < d->geo.AtomsPerCellMax(); m++) {

	 AtomType_l_m = (int) d->geo.AtomType[l][m];
	 /*  Check if atom is not in simulation domain, skip */ 
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
              d->opt.ExecParam.Phys.includeImpurity ) {
            for (int ii=0; ii<d->NBasisStates; ii++)
               d->hd[ii][ii].r -= d->phi[l][m];
         }

         /* Adjust the diagonal block of Hamiltonian for nonzero 
            magnetic field */
         if(d->opt.ExecParam.Phys.MagneticFieldOn) 
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



#ifdef ELIMINATE_SSMAP
	    if (cindx < 0 && has_periodicity && ( cindx + d->geo.N_Cell + 1) >= 0) {
	      
	      real phase = d->geo.phasePeriodic(d->kxL, d->kyL, d->kzL,
						i_nbr, j_nbr, k_nbr);
	      
	      eik.r = cos(phase);
	      eik.i = sin(phase);            
	      cindx = cindx + d->geo.N_Cell + 1 ;
	    }
#else /* ELIMINATE_SSMAP */
	    if (cindx < 0 && has_periodicity && d->geo.ssmap[l][m][n]!=-1) {
	      
	      real phase = d->geo.phasePeriodic(d->kxL, d->kyL, d->kzL,
						i_nbr, j_nbr, k_nbr);
	      
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
            
            H_DiagStrainCorr_OffDiag_Wrapper(d, l, m, AtomType_l_m, n, 
                                   cindx, aindx, nnv, nnv0, a_lattice);
            /* Adjust the offdiagonal block of Hamiltonian for nonzero 
               magnetic field */

            if(d->opt.ExecParam.Phys.MagneticFieldOn && 
               (d->HzbComplex.isInitialized() || 
                d->opt.ExecParam.ElCalc.ElAlg.HamCompute==ElAlg_struct::FullStore))  
               apply_magnetic_field(d, l, m, nnv);
            
            /* Matrix-vector multiply */
            if (out) {
               if (d->opt.ExecParam.ElCalc.ElAlg.HamCompute==ElAlg_struct::SparseStore){ 
                  int ii, jj;

                  int Nb_2;
                  if (d->BandModel==BM_10_sp3d5ss_nospin) Nb_2 = d->NBasisStates;
                  else                                    Nb_2 = d->NBasisStates / 2;

                  /* This matrix will be unique due to vector splitting 
                     across different CPUs.  Therefore we do not have to 
                     check whether we have encountered this matrix already.
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
                     (this proc atom, neighoring proc atom) pair */
                  if(d->Hzb.isInitialized()) {
                     for (ii=0; ii<Nb_2; ii++) {
                     for (jj=0; jj<Nb_2; jj++) {
                        int indx = Nb_2*(Nb_2*Zo + ii) + jj;
                        d->Hzb.Ho[indx] = d->ho[ii][jj].r;
                     }}
                  }else if(d->HzbComplex.isInitialized()){ 
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
                     across different CPUs.  Therefore we do not have to 
                     check whether we have encountered this matrix already.
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
               
#ifdef FORTRAN_MUL
               zgemvt_f( &(y[isyc]), &(d->ho[0][0]), &(cone), &(cone), 
                         &(x[is1]), &(d->NBasisStates) ); 
#else		        
               if (eik.i!=0.0 || eik.r!=1.0)
                  cmatmul_spds_sc(yc, d->ho, x, isyc, is1, d->NBasisStates, 
                                  eik);
               else
                  cmatmul_spds(yc, d->ho, x, isyc, is1, d->NBasisStates); 
#endif
               
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
                  /* See whether we have encountered this neigbor pair 
                     already*/
                  /* What is the (l,m) index of this neighbor?   
                     (cindx,aindex)  */
                  if (d->Ham_use_hermiticity &&
                      (rel_cindx>=0) && (rel_cindx<iex-isx) && (aindx>=0) &&
                      (Zatom > d->l_m_2_atom[rel_cindx][aindx])   ){
                     /* We should have computed this off-diagonal matrix 
                        already. Find it in the list of previous atoms.  */
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
                     }else if(d->HzbComplex.isInitialized()) {
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
                  /* See whether we have encountered this neigbor pair 
                     already*/
                  /* What is the (l,m) index of this neighbor?   
                     (cindx,aindex)  */
                  if (d->Ham_use_hermiticity &&
                      (rel_cindx>=0) && (rel_cindx<iex-isx) && (aindx>=0) &&
                      (Zatom > d->l_m_2_atom[rel_cindx][aindx])   ){
                     /* We should have computed this off-diagonal matrix 
                        already. Find it in the list of previous atoms.  */
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
               
#ifdef FORTRAN_MUL
               zgemvt_f( &(y[isn-imin]), &(d->ho[0][0]), &(cone), &(cone), 
                         &(x[is1]), &(d->NBasisStates) ); 
#else
               if (eik.i!=0.0 || eik.r!=1.0)
                  cmatmul_spds_sc( y, d->ho, x, isn-imin, is1, d->NBasisStates, eik );
               else
                  cmatmul_spds( y, d->ho, x, isn-imin, is1, d->NBasisStates ); 
#endif

         
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
         
         
#undef FORTRAN_MUL
#ifdef FORTRAN_MUL
         zgemvt_f( &(y[is1]), &(d->hd[0][0]), &(cone), &(cone), &(x[is1]), 
                   &(d->NBasisStates) ); 
#else
         cmatmul_spds( y, d->hd, x, is1, is1, d->NBasisStates ); 
#endif
         
         MPI_TOC(Tsmall, Tstart);
         
      }
   }
   
   /* printf("Hmult_spds_col - 2\n"); */

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
	 printf("CPU %d  up_diag %d  low_diag %d\n",proc,up_elements,
		down_elements),fflush(stdout);
   }
   
   for (int i=0; i < d->seg_ln[proc]; i++) {
      y[i].r -= shift*x[i].r;
      y[i].i -= shift*x[i].i;
   }
   
   counter++;


#ifdef MPI3d
   if (d->opt.ExecParam.Phys.Disorder1!=0.0){
       Randomize_H (d, proc);
       d->nvmap[proc] = Hmult_spds_col_fullstored_par(y, yc, d, shift, x, proc);
   }
   return d->nvmap[proc];
#else
   return 0;
#endif       

}


int Hmult_spds_col_prep_store ( cvectr y, cvectr yc, qd_struct d, real shift, 
				     cvectr x, int proc )
{
   int l, m;
   int isx, iex;
   int isyc, yccnt;
   
   int Zatom=0, Zatom_surf=0;
   int INT_LARGE = (int) pow(2.0,30.0);
   
   const int maxNbr = d->geo.NeighborsMax();

   if (d->Ham || d->Hzb.isInitialized() || d->HzbComplex.isInitialized() ) 
      return 0;
   
   nml_memory_report("before Allocation of Hamiltonian Storage - Hmult_spds_col_prep_store ");

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
	 }
      }
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
   printf("mpi_id=%d memory estimate Hamiltonian storage prep: %g MB Natom=%d \n", mpi_n3d_id, memory_estimate/1.0e6,Zatom); fflush(stdout);
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
      if(!d->opt.ExecParam.Phys.MagneticFieldOn)
         d->Hzb.Initialize(d->BandModel,d->NBasisStates, Zatom, Zatom_surf, Nout, Nin, false, HBxy);
      else if(d->opt.ExecParam.ElCalc.ElAlg.Efficiency==ElAlg_struct::Memory)
         d->Hzb.Initialize(d->BandModel,d->NBasisStates, Zatom, Zatom_surf, Nout, Nin, true,  HBxy);
      else if(d->opt.ExecParam.ElCalc.ElAlg.Efficiency==ElAlg_struct::Speed)
         d->HzbComplex.Initialize(d->BandModel,d->NBasisStates, Zatom, Zatom_surf, Nout, Nin, HBxy);
      d->surfaceAtoms = Ivectr(Zatom_surf);
   }
   
   /* print storage info to stdout */
   if (mpi_n3d_id==mpi_n3d_masterid) {
      real storage_main=0.0, storage_help=0.0;

      if (d->opt.ExecParam.ElCalc.ElAlg.HamCompute==ElAlg_struct::FullStore) {
	 if (d->Ham_use_hermiticity){ 
	    storage_main = 16.0*Zatom*(maxNbr/2+1)*d->NBasisStates*
	       d->NBasisStates/1000000.0;
	    storage_help = (6.0*4.0*Zatom*(maxNbr+1) +
			    4.0*(iex-isx)*maxNbr)/1000000.0;
	 }
	 else {
	    storage_main = 16.0*Zatom*(maxNbr+1)*d->NBasisStates*
	       d->NBasisStates/1000000.0;
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


#ifdef DOES_NOT_SEEM_TO_BE_USED
void Ham_to_Hzb_20(qd_struct d) 
{
   int i, j ,atom, nbr;

   int Zo=0, Zu=0;
   int Nb_2 = d->NBasisStates / 2;
   
   if ( d->geo.hasPeriodicity() )
      die("Can't use Ham_to_Hzb_20 with periodic BC yet.");

   /* count the number of off-diagonal and upper-diagonal blocks */
   for (atom=0; atom<d->Hzb.Natom; atom++) {
   for (nbr=0; nbr<=d->Ham_neighbor[atom]; nbr++) {
      int I = d->Ham_offset_y[atom][nbr];
      int J = d->NBasisStates*atom;

      if ( d->Ham_yc[atom][nbr] )
	 Zo++;
      else if ( I < J )
	 Zu++;
   }}

   d->Hzb.indxHU_to_row = Ivectr(Zu);
   d->Hzb.indxHU_to_col = Ivectr(Zu);
   d->Hzb.indxHO_to_row = Ivectr(Zo);
   d->Hzb.indxHO_to_col = Ivectr(Zo);

   Zo=Zu=0;
   for (atom=0; atom<d->Hzb.Natom; atom++) {
   for (nbr=0; nbr<=d->Ham_neighbor[atom]; nbr++) {

      cmatrix thisHam = d->Ham[atom][nbr];
      int I = d->Ham_offset_y[atom][nbr];
      int J = d->NBasisStates*atom;

      if ( I == J ) {
	 /* set up diagonal */
	 for (i=0; i<d->NBasisStates; i++) {
	    d->Hzb.Hdd[J+i] = d->Ham[atom][nbr][i][i].r;
	 }
	 /* set up upper triangular off-diagonal elements of diagonal block */
	 d->Hzb.Hdu[6*atom+0].r = thisHam[2][3].r;
	 d->Hzb.Hdu[6*atom+0].i = thisHam[2][3].i;
	 d->Hzb.Hdu[6*atom+1].r = thisHam[12][13].r;
	 d->Hzb.Hdu[6*atom+1].i = thisHam[12][13].i;
	 d->Hzb.Hdu[6*atom+2].r = thisHam[4][12].r;
	 d->Hzb.Hdu[6*atom+2].i = thisHam[4][12].i;
	 d->Hzb.Hdu[6*atom+3].r = thisHam[4][13].r;
	 d->Hzb.Hdu[6*atom+3].i = thisHam[4][13].i;
	 d->Hzb.Hdu[6*atom+4].r = thisHam[2][14].r;
	 d->Hzb.Hdu[6*atom+4].i = thisHam[2][14].i;
	 d->Hzb.Hdu[6*atom+5].r = thisHam[3][14].r;
	 d->Hzb.Hdu[6*atom+5].i = thisHam[3][14].i;
      }
      else if ( d->Ham_yc[atom][nbr] ) {
	 /* set up off-diagonal superblocks corresponding to 
	    (this proc atom, neighoring proc atom) pair */
	 for (i=0; i<Nb_2; i++) {
	 for (j=0; j<Nb_2; j++) {
	    int indx = Nb_2*(Nb_2*Zo + i) + j;
	    d->Hzb.Ho[indx] = d->Ham[atom][nbr][i][j].r;
	 }}
	 
	 d->Hzb.indxHO_to_row[Zo] = I;
	 d->Hzb.indxHO_to_col[Zo] = J;
	 Zo++;
      }
      else if ( !d->Ham_transp[atom][nbr] ) {
	 /* set up upper diagonal blocks within this diagonal superblock */
	 for (i=0; i<Nb_2; i++) {
	 for (j=0; j<Nb_2; j++) {
	    int indx = Nb_2*(Nb_2*Zu + i) + j;
	    d->Hzb.Hu[indx] = d->Ham[atom][nbr][i][j].r;
	 }}
	 
	 d->Hzb.indxHU_to_row[Zu] = I;
	 d->Hzb.indxHU_to_col[Zu] = J;
	 Zu++;
      }
   }}
}
#endif


/* Computes y=(H-Is)x where H is the Hamiltonian, s is input shift,
   and x is input vector. */
int Hmult ( rvectr y, qd_struct d, real s, rvectr x )
{
    int i, j, k, l, idum;
    real mxp, mxm, myp, mym, mzp, mzm, Tf;

    Tf = HBAR_MASSFACTOR/d->geo.lattice_x/d->geo.lattice_x;

    for ( l=0; l < d->geo.N_Cell; l++ ){
        i = d->geo.cell__ijk[l][0];
        j = d->geo.cell__ijk[l][1];
        k = d->geo.cell__ijk[l][2];

        /* Average masses in x-dimension */
        mxm = (d->mstar[d->mat[i-1][j][k]] + d->mstar[d->mat[i][j][k]] )/2.0;
	mxp = (d->mstar[d->mat[i+1][j][k]] + d->mstar[d->mat[i][j][k]] )/2.0;

        /* Average masses in y-dimension */ 
        mym = (d->mstar[d->mat[i][j-1][k]] + d->mstar[d->mat[i][j][k]] )/2.0;
	myp = (d->mstar[d->mat[i][j+1][k]] + d->mstar[d->mat[i][j][k]] )/2.0;
      
        /* Average masses in z-dimension */ 
        mzm = (d->mstar[d->mat[i][j][k-1]] + d->mstar[d->mat[i][j][k]] )/2.0;
 	mzp = (d->mstar[d->mat[i][j][k+1]] + d->mstar[d->mat[i][j][k]] )/2.0;
            
        y[l] = 0;

        /* Hamiltonian diagonal */
        y[l] += (Tf*(1/mxm+1/mxp+1/mym+1/myp+1/mzm+1/mzp) - s) * x[l];

        /* Off-diagonal terms */
        idum = d->geo.ijk__cell[i-1][j][k];
        if ( idum >= 0 )
            y[l] += -Tf/mxm * x[idum];

        idum = d->geo.ijk__cell[i+1][j][k];
	if ( idum >= 0 )
            y[l] += -Tf/mxp * x[idum];

        idum = d->geo.ijk__cell[i][j-1][k];
        if ( idum >= 0 )
             y[l] += -Tf/mym * x[idum];

        idum = d->geo.ijk__cell[i][j+1][k]; 
	if ( idum >= 0 )
            y[l] += -Tf/myp * x[idum];

        idum = d->geo.ijk__cell[i][j][k-1];
	if ( idum >= 0 )
            y[l] += -Tf/mzm * x[idum];

        idum = d->geo.ijk__cell[i][j][k+1];
	if ( idum >= 0 )
            y[l] += -Tf/mzp * x[idum];
    }
    return(0);
}


/* Computes y=(H-Is)x where H is the Hamiltonian, s is input shift,
   and x is input vector. */
int Hmult_c ( cvectr y, qd_struct d, real s, cvectr x )
{
    int i, j, k, l, idum;
    real mxp, mxm, myp, mym, mzp, mzm, Tf;

    Tf = HBAR_MASSFACTOR/d->geo.lattice_x/d->geo.lattice_x;

    for ( l=0; l < d->geo.N_Cell; l++ ){
        i = d->geo.cell__ijk[l][0];
        j = d->geo.cell__ijk[l][1];
        k = d->geo.cell__ijk[l][2];

        /* Average masses in x-dimension */


        mxm = (d->mstar[d->mat[i-1][j][k]] + d->mstar[d->mat[i][j][k]] )/2.0;
	mxp = (d->mstar[d->mat[i+1][j][k]] + d->mstar[d->mat[i][j][k]] )/2.0;

        /* Average masses in y-dimension */ 
        mym = (d->mstar[d->mat[i][j-1][k]] + d->mstar[d->mat[i][j][k]] )/2.0;
	myp = (d->mstar[d->mat[i][j+1][k]] + d->mstar[d->mat[i][j][k]] )/2.0;
      
        /* Average masses in z-dimension */ 
        mzm = (d->mstar[d->mat[i][j][k-1]] + d->mstar[d->mat[i][j][k]] )/2.0;
 	mzp = (d->mstar[d->mat[i][j][k+1]] + d->mstar[d->mat[i][j][k]] )/2.0;
            
        y[l].r = y[l].i = 0;

        /* Hamiltonian diagonal */
        y[l].r += (Tf*(1/mxm+1/mxp+1/mym+1/myp+1/mzm+1/mzp) - s) * x[l].r
                  + d->eta*x[l].i;
        y[l].i += -d->eta*x[l].r;

        /* Off-diagonal terms */
        idum = d->geo.ijk__cell[i-1][j][k];
        if ( idum >= 0 )
	{
            y[l].r += -Tf/mxm * x[idum].r;
            y[l].i += -Tf/mxm * x[idum].i;
        } 

        idum = d->geo.ijk__cell[i+1][j][k];
	if ( idum >= 0 )
	{
            y[l].r += -Tf/mxp * x[idum].r;
            y[l].i += -Tf/mxp * x[idum].i;
	}

        idum = d->geo.ijk__cell[i][j-1][k];
        if ( idum >= 0 )
	{
            y[l].r += -Tf/mym * x[idum].r;
	    y[l].i += -Tf/mym * x[idum].i;
        }

        idum = d->geo.ijk__cell[i][j+1][k]; 
	if ( idum >= 0 )
	{  
            y[l].r += -Tf/myp * x[idum].r;
            y[l].i += -Tf/myp * x[idum].i;
        }

        idum = d->geo.ijk__cell[i][j][k-1];
	if ( idum >= 0 )
	{  
            y[l].r += -Tf/mzm * x[idum].r;
	    y[l].i += -Tf/mzm * x[idum].i;
        }

        idum = d->geo.ijk__cell[i][j][k+1];
	if ( idum >= 0 )
	{  
            y[l].r += -Tf/mzp * x[idum].r;
            y[l].i += -Tf/mzp * x[idum].i;
        }	    
    }

    return(0);
}



/* Randomize the on-site energies of the sored Hamiltonian in the single band model. */
int Randomize_H ( qd_struct d, int proc)
{
   int atom;
   real step;
   real disorder=d->opt.ExecParam.Phys.Disorder1;
   

   if (d->NBasisStates!=1){
       die("ERROR in Randomize_H:\n"
	     "Cannot run this function in any other than the single band model\n");
   }
   
   
   for (atom=0;atom<vdim1(d->Ham);atom++) {
       step = disorder*rand()/RAND_MAX - disorder/2.0;
       d->Ham[atom][0][0][0].r += step;
   }
   
   return d->nvmap[proc];
}


// count the number of neighbors on this processor and outside this proc.
void count_nbrs_for_sparse_zincblende_ham(int* Nout, int* Nin, qd_struct d) 
{
   bool has_periodicity = d->geo.hasPeriodicity();

   int Zu=0, Zo=0;
   int imin = d->seg_s[mpi_n3d_id];
   int imax = d->seg_s[mpi_n3d_id] + d->seg_ln[mpi_n3d_id];
   int isx = d->cell_s[mpi_n3d_id];
   int iex = isx + d->cell_ln[mpi_n3d_id];

   /* Allocate space for mapping of upper triangular and off-processor 
      part of H.  This algorithm should remain identical in structure 
      to the subsequent algorithm in Hmult_spds_col() that fills the 
      hamiltonian. */
   for ( int l=isx; l < iex; l++ ) {
      /* These are real space indexes for the unit cell */
      int i = d->geo.cell__ijk[l][0];
      int j = d->geo.cell__ijk[l][1];
      int k = d->geo.cell__ijk[l][2];
      
      /* Loop through all the atoms in the unit cell */
      for ( int m=0; m < d->geo.AtomsPerCellMax(); m++ ) {
         int AtomType_l_m = (int) d->geo.AtomType[l][m];
         if ( !AtomType_l_m ) continue;
         
         /* now, take care of off-diagonal guys */
         for ( int n = 0; n < d->geo.Neighbors(m); n++ ) {
            /* Hamiltonian off-diagonal.  Here we loop through all this 
               atom's nearest neighbors get the coupling hamiltonians and 
               do the matrix vector multiplications. */
            /* Get relative position of unit cell where this neighbor lives. */ 
            int* nnmap_m_n = d->geo.NbrCell(m,n);
            int id = *(nnmap_m_n++);
            int jd = *(nnmap_m_n++);
            int kd = *(nnmap_m_n++);
            int aindx = *(nnmap_m_n);
            
            int cindx = d->geo.ijk__cell[i+id][j+jd][k+kd];

#ifdef ELIMINATE_SSMAP
	    if (cindx < 0 && has_periodicity && ( cindx + d->geo.N_Cell + 1) >= 0) {
	      cindx = cindx + d->geo.N_Cell + 1 ;
	    }
#else /* ELIMINATE_SSMAP */
	    if (cindx < 0 && has_periodicity && d->geo.ssmap[l][m][n]!=-1) {
	      cindx = d->geo.ssmap[l][m][n];
	    }
#endif /* ELIMINATE_SSMAP */
            

            /* Check that neighbor's unit cell is within simulation domain 
               AND that atom is within the simulation domain */
            if ( cindx >= 0 && d->geo.AtomType[cindx][aindx] ) {
               /* Get starting point on vector for this neigbhor.  */
               int isn = d->geo.l_m_2_offset_ham_tot[cindx][aindx];
               
               /* Check if neighbor is outside the x vector range */
               if ( isn >= imax || isn < imin )
                  Zo++;
               else
                  Zu++;
            }
         }
      }
   }
   
   Zu /= 2; /* correct our double count */

   *Nout = Zo;
   *Nin = Zu;
}


