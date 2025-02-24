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
$Header: /repo/nemo3d/src/base/rayleigh_ritz.c,v 1.22 2005/12/24 03:37:14 gekco Exp $ 
*****************************************************************************/

#include "rayleigh_ritz.h"

#include "macopt.h"
#include "macopt_par.h"
#include "cmacopt.h"
#include "ritz.h"
#include "MPI_Timing.h"

#undef COMPLEX_CG	

/* solve (H-shift)^2 * x = r using CG on normal equations */
void linSolv_CGNE(cvectr x, qd_struct d, real shift, cvectr r, real tol, 
                  real* xx, real* xHx) 
{
   //void **argList = (void**) new void* [4];
   //argList[0] = (void*) &d;
   //argList[1] = (void*) &(d->opt.ExecParam.ElCalc.ElAlg.RayleighRitz.EigValGuess[indx]);
   //double shift1 = 0.0;
   //argList[1] = shift1;
   //argList[2] = (void*) &indx;
   //argList[3] = (void*) &eval_scoperti;

   int sz = d->seg_ln[mpi_n3d_id];
   real alpha, beta, rr_loc, rr, rr_old, pHp, pHp_loc, rcv[3], snd[3];

   real ritz=0.0, ritz_old=1e10;

   cvectr p = Cvectr(sz);
   cvectr Hp = Cvectr(sz);
   cvectr HHp = Cvectr(sz);

   FILE* fp = fopen("cg_conv","a");

   masterPrint("entering CGNE (shift=%e)\n", shift);
   
   for (int j=0; j<sz; j++) {
      p[j].r = r[j].r;
      p[j].i = r[j].i;
      x[j].r = x[j].i = 0.0;
      Hp[j].r = Hp[j].i = 0.0;
      HHp[j].r = HHp[j].i = 0.0;
   }

#define MAX_ITER_LINSOLV 100000
   for (int i=0; i < MAX_ITER_LINSOLV; i++) {
      
      //Hmult_spds_complete(Hp, d, shift, p);
#undef CGNE
#ifdef CGNE
      Hmult_spds_complete(HHp, d, shift, p);
      Hmult_spds_complete(Hp, d, shift, HHp);
#else
   //double shift = d->eigv[0].r;
      //matvect_mult((const int**) argList, Hp, p);
      Hmult_spds_complete(Hp, d, shift, p);
#endif

      rr_loc = 0.0;
      pHp_loc = 0.0;
      for (int j=0; j<sz; j++) {
	 rr_loc += r[j].r*r[j].r + r[j].i*r[j].i;
	 pHp_loc += Hp[j].r*p[j].r + Hp[j].i*p[j].i;
      }


      snd[0] = rr_loc;
      snd[1] = pHp_loc;
#if (defined MPI3d && !defined FAKE_MPI)
      MPI_Allreduce(&snd[0], &rcv[0], 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else /* MPI3d && !FAKE_MPI */
      rcv[0]=snd[0];
      rcv[1]=snd[1];
#endif /* MPI3d && !FAKE_MPI */
      rr = rcv[0];  // residual = [(H-shift)^2 * x - b]
      pHp = rcv[1];


      // debuggin infocmp
#ifdef DEBUG_HAMILTONIAN
      char dum[100];
      masterPrint("%d\n", i);
      sprintf(dum, "DB/p%d:%d:%d", mpi_n3d_numprocs, mpi_n3d_id, i);
      FILE *fp_p = fopen(dum,"w");
      sprintf(dum, "DB/Hp%d:%d:%d", mpi_n3d_numprocs, mpi_n3d_id, i);
      FILE *fp_Hp = fopen(dum,"w");
      for (int j=0; j<sz; j++) {
         fprintf(fp_p, "%e %e\n", p[j].r, p[j].i);
         fprintf(fp_Hp, "%e %e\n", Hp[j].r, Hp[j].i);
      }
      fclose(fp_p);
      fclose(fp_Hp);
      if (i==0) exit(1);
#endif      
      


      alpha = rr / pHp; // step length

      for (int j=0; j<sz; j++) {
         // approximate solution
	 x[j].r += alpha*p[j].r;
	 x[j].i += alpha*p[j].i;

         // residual
	 r[j].r -= alpha*Hp[j].r;
	 r[j].i -= alpha*Hp[j].i;
      }
      
      rr_old = rr;
      rr_loc = 0.0;
      for (int j=0; j<sz; j++) {
	 rr_loc += r[j].r*r[j].r + r[j].i*r[j].i;
      }
      
      snd[0] = rr_loc;
#if (defined MPI3d && !defined FAKE_MPI)
      MPI_Allreduce(&snd[0], &rcv[0], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else /* MPI3d && !FAKE_MPI */
      rcv[0]=snd[0];
#endif /* MPI3d && !FAKE_MPI */

      rr = rcv[0];

      beta = rr / rr_old; // improvement this step
      
      for (int j=0; j<sz; j++) {
         // search direction
	 p[j].r = r[j].r + beta*p[j].r;
	 p[j].i = r[j].i + beta*p[j].i;
      }

      /* print out some info every 100 iterations */
      if (i%100==0) {

         if (fabs((ritz-ritz_old)/ritz_old) < tol) {
            masterPrint("%e %e\n", ritz, ritz_old);
            break;
         }
         
	 Hmult_spds_complete(Hp, d, 0.0, x);
	 pHp_loc = 0.0;
	 rr_loc = 0.0;
	 for (int j=0; j<sz; j++) {
	    pHp_loc += Hp[j].r*x[j].r + Hp[j].i*x[j].i;
	    rr_loc  +=  x[j].r*x[j].r +  x[j].i*x[j].i;
	 }
	 snd[0] = rr_loc;
	 snd[1] = pHp_loc;
#if (defined MPI3d && !defined FAKE_MPI)
	 MPI_Allreduce(&snd[0], &rcv[0], 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else /* MPI3d && !FAKE_MPI */
	 rcv[0]=snd[0];
	 rcv[1]=snd[1];
#endif /* MPI3d && !FAKE_MPI */
	 rr = rcv[0];
	 pHp = rcv[1];

         ritz_old = ritz;
         ritz = pHp/rr;

	 if (mpi_n3d_id == mpi_n3d_masterid) {
	    printf("  %23.15e %23.15e %23.15e %23.15e %23.15e\n", 
                   pHp/rr, pHp, rr_old, alpha, beta);
	    fprintf(fp, "  %23.15e %23.15e %23.15e\n", pHp/rr, pHp, rr_old);
	    fflush(stdout);
	 }
      }
   }
   
   {
      Hmult_spds_complete(Hp, d, 0.0, x);
      pHp_loc = 0.0;
      rr_loc = 0.0;
      for (int j=0; j<sz; j++) {
	 pHp_loc += Hp[j].r*x[j].r + Hp[j].i*x[j].i;
	 rr_loc  +=  x[j].r*x[j].r +  x[j].i*x[j].i;
      }
      snd[0] = rr_loc;
      snd[1] = pHp_loc;
#if (defined MPI3d && !defined FAKE_MPI)
      MPI_Allreduce(&snd[0], &rcv[0], 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else /* MPI3d && !FAKE_MPI */
      rcv[0]=snd[0];
      rcv[1]=snd[1];
#endif /* MPI3d && !FAKE_MPI */
      *xx = rcv[0];
      *xHx = rcv[1];
   }

   fclose(fp);
   rm_cvectr(&p);
   rm_cvectr(&Hp);
}



void rayleighQuotientIter(qd_struct d)
{
   extern void verify_ham ( qd_struct d);
   
   d->rc = Cvectr(d->nvmap[mpi_n3d_id]*d->NBasisStates); 
   d->drq = Cvectr(d->n_ham_tot);
   d->wfc = Cmatrix(d->WF.size(), d->seg_ln[mpi_n3d_id]);
   //verify_ham(d);exit(1);

   for (int indx=0; indx < (signed int) d->WF.size(); indx++) {
      double eval_guess = d->WF[indx].E;

      masterPrint("----------------------------------------------------------");
      masterPrint("Rayleigh Quotient Iteration:  eigenvalue guess[%d] = %e",
                  indx, eval_guess);
      masterPrint("----------------------------------------------------------");

      complex* y     = new complex[d->seg_ln[mpi_n3d_id]];
      complex* y_old = new complex[d->seg_ln[mpi_n3d_id]];
      
      real dummy=1.0/d->seg_ln[mpi_n3d_id] / 30.0;
      for (int i=0; i<d->seg_ln[mpi_n3d_id]; i++) {
         // initialize with nonzero 
         y_old[i].r = dummy;
         y_old[i].i = 0.0;

         if (i%10==1) {
            y_old[i].r += 0.3*(i%7)*dummy;
            y_old[i].i += 0.5*dummy;
         }
         

      }
      
      real dE=1e20;
      for (int j=0; j<40; j++) {
         real yy, yHy;
         // solve (H-eval_guess)*y = y_old

         // assume we know eigenvalues exactly
         linSolv_CGNE(y, d, d->WF[indx].E, y_old, 1e-8, &yy, &yHy);
         //linSolv_CGNE(y, d, eval_guess, y_old, 1e-11, &yy, &yHy);

         real sqrt_yy_m1 = 1.0 / sqrt(yy);
         for (int i=0; i<d->seg_ln[mpi_n3d_id]; i++) {
            y_old[i].r = sqrt_yy_m1 * y[i].r;
            y_old[i].i = sqrt_yy_m1 * y[i].i;
         }
         dE = fabs(yHy/yy - eval_guess);
         eval_guess = yHy / yy;

         if (dE < 1.0e-10) {
            masterPrint("Rayleigh quotient iteration done.  (<y|H|y>/<y|y>=%e)\n", eval_guess);
            break;
         }
      }
      for(int jj=0; jj<d->seg_ln[mpi_n3d_id]; jj++) {
        d->wfc[indx][jj].r = y_old[jj].r;
        d->wfc[indx][jj].i = y_old[jj].i;
      }
      writeEigenvector(d, y_old, indx, "-cgne");
      printPsisq(d, y_old, indx);
/*      
      complex* p_new = new complex[d->seg_ln[mpi_n3d_id]];
      Hmult_spds_complete(p_new, d, 0.0, y_old);
      writeEigenvector(d, p_new, indx, "-debug");
      delete [] p_new;
*/      
      delete [] y_old;
      delete [] y;
   }
   
   compute_ritz_residual_norm_deviation(d);
#if 0
   compute_ritz_value(d);
   compute_residual(d);
   compute_eigenvector_norm(d);   
#endif
   rm_cvectr(&d->drq);
   rm_cvectr(&d->rc);
}


real calc_rayleigh_sq( rvectr p, void *data )
{
   int i;
   real xx;
   complex xa2x, rq, cdum, cdumr, cdumi;
   cvectr x, y, z;
   qd_struct d;
   static int itercount=0;

   itercount++;
   d = (qd_struct)(data);
   x = Cvectr(d->n_ham_tot);
   y = Cvectr(d->n_ham_tot);
   z = Cvectr(d->n_ham_tot);

   for ( i=0; i < d->n_ham_tot; i++ ) {
      x[i].r = p[i];
      x[i].i = p[i+d->n_ham_tot];
   } 
            
   Hmult_spds_complete ( y, d, d->opt.ExecParam.ElCalc.ElAlg.RayleighRitz.EigValGuess[0], x );
   Hmult_spds_complete ( z, d, d->opt.ExecParam.ElCalc.ElAlg.RayleighRitz.EigValGuess[0], y ); 

   xx = 0;
   xa2x.r = xa2x.i = 0;
   for ( i=0; i < d->n_ham_tot; i++ ) {
      xx += x[i].r*x[i].r + x[i].i*x[i].i;
      xa2x.r += x[i].r*z[i].r + x[i].i*z[i].i;
      xa2x.i += x[i].r*z[i].i - x[i].i*z[i].r;
      /* xa2x.r += y[i].r*y[i].r + y[i].i*y[i].i; */
   }

   rq.r = xa2x.r/xx;
   rq.i = xa2x.i/xx;

   /* Need complex version of conjugate gradient */
   for ( i=0; i < d->n_ham_tot; i++ ) {
      cdum = csub( z[i], cmul(rq, x[i] ) );
      cdumr = csub( z[i], cmul(rq,Complex(x[i].r,0.0)) );
      cdumi = csub( z[i], cmul(rq,Complex(0.0,x[i].i)) );
      d->drq[i] = cdum;
   }

#ifdef LOW_LEVEL_PRINT_RAYLEIGH
   printf("\nIter = %d Rayleigh = %e %e", itercount, rq.r, rq.i);
#endif /* LOW_LEVEL_PRINT_RAYLEIGH */

   rm_cvectr(&x);
   rm_cvectr(&y);
   rm_cvectr(&z);
   return rq.r; 
}


/*  Computes Rayleigh quotients = <x|H|x> / <x|x>
    Since H is hermitian, this value should be mostly real.
      */

complex calc_rayleigh_compl( cvectr x, int sz, qd_struct d )
{
  
   int i;
   real xx;
   complex xhx, rq;
   cvectr y = Cvectr(sz);
   
#ifdef MPI3d
   real rcv[3], snd[3];
#endif /* MPI3d */

   Hmult_spds_complete(y, d, 0.0, x);
   
   xx = xhx.r = xhx.i = 0;
   for ( i=0; i < sz; i++ ) {
      xx += x[i].r*x[i].r + x[i].i*x[i].i;
      xhx = cadd( xhx, cmul( cconj(x[i]), y[i]));
   }
   
#if (defined MPI3d && !defined FAKE_MPI)
   snd[0] = xx;
   snd[1] = xhx.r;
   snd[2] = xhx.i;
   MPI_Allreduce(&snd[0], &rcv[0], 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   rq.r = rcv[1] / rcv[0];
   rq.i = rcv[2] / rcv[0];
#else /* MPI3d && !FAKE_MPI */
   rq.r = xhx.r/xx;
   rq.i = xhx.i/xx;
#endif /* MPI3d && !FAKE_MPI */
   
   rm_cvectr(&y);
   
   return rq;
}


#ifdef MPI3d
void calc_rayleigh_sq_macopt(rvectr p, int sz, rvectr xi, int _reset_shift, 
			     void *data, int itercount, real *onevalue)
#else /* MPI3d */
void calc_rayleigh_sq_macopt(rvectr p, rvectr xi, void *data, int itercount, 
			     real *onevalue)
#endif /* MPI3d */
{
   /* IMPORTANT NOTE:  The parallel version of this function needs another
      argument that specifies the size of p.  We cannot simply use vdim(p) 
      because it is unknown whether p has been recast from a cvectr.  
      If p has been recast, vdim looks 2*sizeof(double) bytes before p to 
      determine size, instead of sizeof(double) if p truly is an rvector.  
      In general, the current implementation of vectors in NEMO makes using 
      vdim quite dangerous.  Thus, we implement the following kludge:  We add 
      an extra argument to this function and to calc_rayleigh_sq_cmaopt() for 
      the parallel version and promise not to use vdim with these 
      pointers p, xi */

#ifndef MPI3d
   //int sz = vdim(p);
#endif /* MPI3d */

   cvectr cp = (cvectr) p;
   cvectr cxi = (cvectr) xi;
#ifdef MPI3d
   calc_rayleigh_sq_cmacopt(cp, sz/2, cxi, _reset_shift, data, itercount, 
			    onevalue);
#else /* MPI3d */
   calc_rayleigh_sq_cmacopt(cp, cxi, data, itercount, onevalue);
#endif /* MPI3d */
   return;
}


#define H2_min_start 1.0e30

#ifdef MPI3d
void calc_rayleigh_sq_cmacopt( cvectr p, int sz, cvectr xi, int _reset_shift, 
			       void *data, int itercount, real *onevalue )
#else /* MPI3d */
void calc_rayleigh_sq_cmacopt( cvectr p, cvectr xi, void *data, int itercount, 
			       real *onevalue )
#endif /* MPI3d */
{
   int i;
   real xx=0.0;
   complex xa2x=czero, rq;
   cvectr x, y, z;
   register real rx, ix, rz, iz, ry, iy;
    
   static real H2_min=H2_min_start;

   qd_struct d = (qd_struct) data;
   real shift = d->opt.ExecParam.ElCalc.ElAlg.RayleighRitz.EigValGuess[d->opt.ExecParam.ElCalc.ElAlg.RayleighRitz.RunningIndex];

#define GET_RAYLEIGH_AS_WELL
#ifdef GET_RAYLEIGH_AS_WELL
   complex r, xax=czero;
#endif /* GET_RAYLEIGH_AS_WELL */

#ifdef MPI3d
   real snd[5], rcv[5];
#else /* MPI3d */
   int sz = d->n_ham_tot;
#endif /* MPI3d */

   x = Cvectr(sz); /* copy input array p */
   y = Cvectr(sz); /* set to (H-shift) * x */
   z = Cvectr(sz); /* set to (H-shift)^2 * x */

   for ( i=0; i < sz; i++ )
      x[i] = p[i];

   Hmult_spds_complete( y, d, shift, x );
   Hmult_spds_complete( z, d, shift, y ); 

   for ( i=0; i < sz; i++ ){
      rx = x[i].r, ix = x[i].i;
      rz = z[i].r, iz = z[i].i;
      /* xx += x[i].r*x[i].r + x[i].i*x[i].i; */
      xx += rx*rx + ix*ix;
	
      /* xa2x = cadd( xa2x, cmul( cconj(x[i]), z[i]));  */
      xa2x.r += rx*rz + ix*iz;
      xa2x.i += rx*iz - ix*rz; 
      /* xa2x.r += x[i].r*z[i].r + x[i].i*z[i].i;
	 xa2x.i += x[i].r*z[i].i - x[i].i*z[i].r; */
	
#ifdef GET_RAYLEIGH_AS_WELL
      /* xax  = cadd( xax, cmul( cconj(x[i]), y[i]));  */
      ry = y[i].r, iy = y[i].i;
      xax.r += rx*ry + ix*iy;
      xax.i += rx*iy - ix*ry; 
#endif /* GET_RAYLEIGH_AS_WELL */
   }

#if (defined MPI3d && !defined FAKE_MPI)
    snd[0] = xx;
    snd[1] = xa2x.r;
    snd[2] = xa2x.i;
#ifdef GET_RAYLEIGH_AS_WELL
    snd[3] = xax.r;
    snd[4] = xax.i;
#endif /* GET_RAYLEIGH_AS_WELL */
    MPI_Allreduce(&snd[0], &rcv[0], 5, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    xx = rcv[0];
    xa2x.r = rcv[1];
    xa2x.i = rcv[2];
#ifdef GET_RAYLEIGH_AS_WELL
    xax.r = rcv[3];
    xax.i = rcv[4];
#endif /* GET_RAYLEIGH_AS_WELL */
#endif /* MPI3d && !FAKE_MPI */

   rq.r = xa2x.r/xx;
   rq.i = xa2x.i/xx;
#ifdef GET_RAYLEIGH_AS_WELL
   r.r = xax.r/xx;
   r.i = xax.i/xx;
                                                 #ifdef MPI3d
                                                   if (mpi_n3d_id==mpi_n3d_masterid)
                                                 #endif
   {
      if (d->opt.ExecParam.ElCalc.ElOut.CGiter_log2 && 
          mpi_n3d_id==mpi_n3d_masterid) {
	 FILE *fp;
	 char *filename = n3d_strdup_n(d->inputfile);
	 n3d_FileTypeSet(&filename, "nd_eval_iter2", TRUE);
	 fp = fopen( filename, "a" );
	 str_free(filename);
	 fprintf(fp,"%d %e %e\n",itercount, r.r+shift, rq.r);
	 fclose(fp);
      }
      *onevalue = r.r+shift;

      if ( fabs(r.i/r.r) > 1.0e-10 || fabs(rq.i/rq.r) > 1.0e-10 )
	 printf("ACHTUNG:  Rayleigh quotient has a large imaginary value:  %e %e\n",
		r.i, rq.i);
      printf("\n  <x|(H-s)^2|x> = %12.8g,  ", rq.r);
      printf("  <x|H|x>=%12.8g  target=%12.8g\n",r.r+shift, shift), fflush(stdout);;

#ifdef MPI3d
      if (_reset_shift) {
	 /*
	 if (mpi_n3d_id==mpi_n3d_masterid)
	    printf("    >>> resetting shift to %f\n", r.r + shift);

	 d->opt.ExecParam.ElCalc.ElAlg.RayleighRitz.EigValGuess[d->opt.ExecParam.ElCalc.ElAlg.RayleighRitz.RunningIndex] =
	    r.r+shift;
	 */
      }
#endif /* MPI3d */
   }
#endif /* GET_RAYLEIGH_AS_WELL */	

#define WRITE_VEC
#ifdef WRITE_VEC
   {
      int j;
      real dummy;

      dummy=sqrt(rq.r*rq.r+rq.i+rq.i);
      /* printf("sqrt(rq.r*rq.r+rq.i+rq.i)=%g  H2_min=%g\n",dummy,H2_min); */

      if (dummy<H2_min){
                                                 #ifdef MPI3d
                                                   if (mpi_n3d_id==mpi_n3d_masterid)
                                                 #endif
	 {
	    printf("Found a new minimum from the previous value ");
	    printf("%g to the new value %g\n",H2_min,dummy);
	    fflush(stdout);
	 }
	 
	 if (H2_min==H2_min_start){
	    /* Do not write out the very first vector */
	    H2_min = dummy;
	 } else {
	    H2_min = dummy;
#if 0
                                                 #ifdef MPI3d
                                                   if (mpi_n3d_id==mpi_n3d_masterid)
                                                 #endif
	    printf("Prep Write intermediate file\n");
		
#ifdef MPI3d
	    writeSiteInfo(d, "basis", "Data/PAR_rayleigh_vec","{real x.r, real x.i}", x);
	    MPI_Barrier(MPI_COMM_WORLD);
#else /* MPI3d */
	    {
	       FILE *fvec;
	       if ((fvec = fopen("Data/rayleigh_vec","w"))!=NULL){
		  for ( j=0; j < sz; j++ )
		     fprintf(fvec,"%1.15e %1.15e\n",x[j].r, x[j].i);
		  fclose(fvec);
	       } else {
		  die("Could not open the file Data/rayleigh_vec\n");
	       }
	    }
#endif /* MPI3d */
                                                 #ifdef MPI3d
                                                   if (mpi_n3d_id==mpi_n3d_masterid)
                                                 #endif
	    printf("Done	 Write intermediate file\n");
#endif /* 0 */

#ifdef RESTART_GEKCO	
	    {
#ifdef COMPLEX_CG
	       /* NOTE: parallelize complex version later */
	       cmacopt_args *cmac_arg = get_cmacopt_args();
	       cvectr x_restart = cmac_arg->x_restart;
	       if (!x_restart){
		  x_restart=Cvectr(cmac_arg->n);
		  cmac_arg->x_restart=x_restart;
	       }
		    
	       for ( j=0; j < sz; j++ )
		  x_restart[j] = x[j];
#else /* COMPLEX_CG */
#ifdef MPI3d
	       macopt_args *mac_arg = get_macopt_args_par();
#else /* MPI3d */
	       macopt_args *mac_arg = get_macopt_args();
#endif /* MPI3d */
	       rvectr x_restart = mac_arg->x_restart;
	       cvectr cx_restart=NULL;
	       
	       if (!x_restart){
		  x_restart=Rvectr(2*mac_arg->n);
		  mac_arg->x_restart=x_restart;
	       }
	       cx_restart = (cvectr) x_restart;
		    		    
	       for ( j=0; j < sz; j++ )
		  cx_restart[j] = x[j];
#endif /* COMPLEX_CG */

	    }
#endif /* RESTART_GEKCO */	    	
	 }
      }
   }
#endif /* WRITE_VEC */


#define PROPER_NORMALIZATION
#ifdef PROPER_NORMALIZATION
   /* Routine returns grad{ r = <x|H^2|x> / <x|x> } = (H^2 x - rx) * 2/<x|x> */
    
   rz = 2.0/xx;
    
   ry = rq.r,   iy = rq.i;	
   for ( i=0; i < sz; i++ ){
      /* xi[i] = Rcmul(2.0/xx,csub( z[i], cmul(rq, x[i] ) )); */
      rx = x[i].r, ix = x[i].i;
      xi[i].r = rz * (z[i].r - (ry*rx - iy*ix));
      xi[i].i = rz * (z[i].i - (ry*ix + iy*rx));	
   }
#else /* PROPER_NORMALIZATION */
   ry = rq.r,   iy = rq.i;	
   for ( i=0; i < sz; i++ ){
      /* xi[i] = csub( z[i], cmul(rq, x[i] ) ); */
      rx = x[i].r, ix = x[i].i;
      xi[i].r = z[i].r - (ry*rx - iy*ix);
      xi[i].i = z[i].i - (ry*ix + iy*rx);	
   }
#endif /* PROPER_NORMALIZATION */
   

   rm_cvectr(&x);
   rm_cvectr(&y);
   rm_cvectr(&z);
}


real calc_rayleigh( rvectr p, void *data )
{
   int i;
   real xx;
   complex xa2x, rq, cdum;
   cvectr x, y, z;
   qd_struct d;

   d = (qd_struct)(data);
   x = Cvectr(d->n_ham_tot);
   y = Cvectr(d->n_ham_tot);
   z = Cvectr(d->n_ham_tot);

   for ( i=0; i < d->n_ham_tot; i++ ) {
      x[i].r = p[i];
      x[i].i = p[i+d->n_ham_tot];
   }  

   Hmult_spds_complete ( y, d, d->opt.ExecParam.ElCalc.ElAlg.RayleighRitz.EigValGuess[0], x );

   xx = 0;
   xa2x.r = xa2x.i = 0;
   for ( i=0; i < d->n_ham_tot; i++ ) {
      xx += x[i].r*x[i].r + x[i].i*x[i].i;
      xa2x.r += x[i].r*z[i].r + x[i].i*z[i].i;
      xa2x.i += x[i].r*z[i].i - x[i].i*z[i].r;
   }

   rq.r = xa2x.r/xx;
   rq.i = xa2x.i/xx;

   for ( i=0; i < d->n_ham_tot; i++ ) {
      cdum = csub( z[i], cmul(rq, x[i] ) ); 
      d->drq[i] = cdum;
   }

#ifdef LOW_LEVEL_PRINT_RAYLEIGH
   printf("\nRayleigh = %e %e", rq.r, rq.i);
#endif /* LOW_LEVEL_PRINT_RAYLEIGH */

   rm_cvectr(&x);
   rm_cvectr(&y);
   rm_cvectr(&z);
   return rq.r; 
}


void calc_rayleigh_d( rvectr p, rvectr xi, void *data )
{
   int i;
   qd_struct d;

   d = (qd_struct)(data);

   for ( i=0; i < d->n_ham_tot*2; i++ )
      xi[i] = d->drq[i].r;

   return;
}


int min_rayleigh(qd_struct d)
{
   const int& Nguess = d->opt.ExecParam.ElCalc.ElAlg.RayleighRitz.EigValGuess.size();

   cmacopt_args camac;
   macopt_args amac;
   char filename[100];

   int n_tot = d->n_ham_tot;
   char *filename_dyn=NULL;

#ifdef MPI3d
   char outfile[100];
   int i, proc;
   cvectr p = Cvectr(d->seg_ln[mpi_n3d_id]);
#else /* MPI3d */
   cvectr p = Cvectr(d->n_ham_tot);
#endif /* MPI3d */

   d->rc = Cvectr(d->nvmap[mpi_n3d_id]*d->NBasisStates); 
   d->drq = Cvectr(d->n_ham_tot);

   cmacopt_defaults(&camac,
		    d->opt.ExecParam.ElCalc.ElAlg.RayleighRitz.MaxIter,
                    d->opt.ExecParam.ElCalc.ElAlg.RayleighRitz.verbose,
		    d->opt.ExecParam.ElCalc.ElAlg.RayleighRitz.tol, 
                    d->opt.ExecParam.ElCalc.ElAlg.RayleighRitz.ConvCrit);
   macopt_defaults(&amac,
		   d->opt.ExecParam.ElCalc.ElAlg.RayleighRitz.MaxIter,
                   d->opt.ExecParam.ElCalc.ElAlg.RayleighRitz.verbose,
		   d->opt.ExecParam.ElCalc.ElAlg.RayleighRitz.tol, 
                   d->opt.ExecParam.ElCalc.ElAlg.RayleighRitz.ConvCrit);

   amac.grad_tol_tiny = 0.0;
   camac.grad_tol_tiny = 0.0;

   // list of discovered eigenvalues
   int* eval_scoperti = new int[Nguess];

   if (d->opt.ExecParam.ElCalc.ElAlg.ResFind==ElAlg_struct::Rayleigh) {
      d->WF = vector<waveFunction>(Nguess);
      for (int i=0; i<Nguess; i++)
         d->WF[i].psi = vector<complex>(d->seg_ln[mpi_n3d_id]);
   }

   for (int indx=0; indx<Nguess; indx++){
      d->opt.ExecParam.ElCalc.ElAlg.RayleighRitz.RunningIndex = indx;
      
      if (d->opt.ExecParam.ElCalc.ElAlg.RayleighRitz.StartFileRead != "NULL") {
         if (Nguess==1)
            sprintf(filename,"%s",
                    d->opt.ExecParam.ElCalc.ElAlg.RayleighRitz.StartFileRead.c_str());
         else
            sprintf(filename,"%s_%d",
                    d->opt.ExecParam.ElCalc.ElAlg.RayleighRitz.StartFileRead.c_str(), 
                    indx+1);

         readWaveFunction(filename, d, p);

         complex cdummy = calc_rayleigh_compl( p, nml_dcv_extent(p), d );
         masterPrint("Rayleigh(INITIAL) = (%g,%g)\n",cdummy.r,cdummy.i);
         fflush(stdout);
      } else {
	 /* at least as a fall-back position provide a non-zero
	    guess to the eigenvector:  */
	 p[d->seg_ln[mpi_n3d_id]/2].r = 1.0;
      }

      masterPrint("Start the conjugate gradient work now.\n");
      
      if (mpi_n3d_id==mpi_n3d_masterid){
	 if (d->opt.ExecParam.ElCalc.ElOut.CGiter_log) {
	     filename_dyn = n3d_strdup_n(d->inputfile);
	     n3d_FileTypeSet(&filename_dyn, "nd_eval_iter", TRUE);
	     amac.filename_its = filename_dyn;
	  }
      }
#ifdef COMPLEX_CG	
      cmacoptII(p, n, (calc_rayleigh_sq_cmacopt), d, &camac); 
#else /* COMPLEX_CG */

#ifdef MPI3d

      {
         int verbosity = 1;

         char* fil_trace=NULL;
         if (d->opt.ExecParam.ElCalc.ElOut.TraceEigval ){
            fil_trace = n3d_strdup_n(d->inputfile);
            n3d_FileTypeSet(&fil_trace, "nd_eval_trace", TRUE);
         }

         void **opArgList = (void**) new void* [4];
         opArgList[0] = (void*) &d;
         opArgList[1] = (void*) &(d->opt.ExecParam.ElCalc.ElAlg.RayleighRitz.EigValGuess[indx]);
         opArgList[2] = (void*) &indx;
         opArgList[3] = (void*) &eval_scoperti;

         int **argList = (int**) new int*[5];
         argList[0] = (int*) &matvect_mult;
         argList[1] = (int*) opArgList;
         argList[2] = (int*) 
            &(d->opt.ExecParam.ElCalc.ElAlg.RayleighRitz.EigValGuess[indx]);
         argList[3] = (int*) &verbosity;
         argList[4] = (int*) &fil_trace;

         macoptII_par((rvectr) p, 2*d->seg_ln[mpi_n3d_id], 2*n_tot,
                      (gradRayleighQuotient), 
                      argList, (macopt_args*) &amac);


         if (d->opt.ExecParam.ElCalc.ElOut.TraceEigval )   str_free(fil_trace);

         delete [] opArgList;
         delete [] argList;
      }


      writeEigenvector(d, p, indx, "-rayleigh");

#else /* MPI3d */
      macoptII((rvectr) p, 2*n_tot, (calc_rayleigh_sq_macopt), d, &amac);
#endif /* MPI3d */
#endif /* COMPLEX_CG */

      if (amac.filename_its)
	 str_free(amac.filename_its), amac.filename_its=NULL;
      if (amac.filename_last)
	 str_free(amac.filename_last), amac.filename_last=NULL;

      // find a match to the newly discovered eigenvector
      // if we've gotten to this point after having done Lanczos,
      // then use this info
      if (d->opt.ExecParam.ElCalc.ElAlg.ResFind==ElAlg_struct::LanczAndRayleigh) {
         int closest_indx=0;
         real closest_diff = 1e20;
         for (int i=0; i < (signed int) d->WF.size(); i++) {
            real this_diff = fabs(d->WF[i].E - amac.track_one_value);
            if ( this_diff < closest_diff ) {
               closest_indx = i;
               closest_diff = this_diff;
            }
            cout << i << "  " << closest_indx << "   " << closest_diff << endl;
         }
         masterPrint("### found a match:  indx=%d E=(%e,%e)\n", 
                     closest_indx,
                     amac.track_one_value,
                     d->WF[closest_indx].E
                     );
         eval_scoperti[indx] = closest_indx;
         // copy eigenvector to waveFunction data structure
         d->WF[ eval_scoperti[indx] ].copy_wavefunc(p);
      }
      else {
         eval_scoperti[indx] = indx;
         d->WF[indx].indx = indx;
         d->WF[indx].E = amac.track_one_value;
         d->WF[indx].copy_wavefunc(p);
      }

      printPsisq(d, p, indx);
      writeEigenvector(d, p, indx, "-rayleigh");
   }

   if (d->opt.ExecParam.ElCalc.ElAlg.ResFind==ElAlg_struct::Rayleigh) {
      writeEigenvalues(d, 'a');
   }

   delete [] eval_scoperti;

   rm_cvectr(&d->rc);
   rm_cvectr(&d->drq);
   
   rm_cvectr(&p);

   return 0;
}
