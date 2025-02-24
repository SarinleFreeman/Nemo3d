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
$Header: /repo/nemo3d/src/base/eigsys3d.c,v 1.11 2006/09/04 15:58:37 gekco Exp $ 
*****************************************************************************/

#include "eigsys3d.h"

int sym_lanc_it_c( qd_struct d, real shift, real *alpha, real *beta, real beta_m1,
		   cvectr r, cvectr q, cvectr q_m1, cvectr q_m2 )
{
   int i;
   complex alc;
   real bl, al;

   for ( i=0; i < d->n_ham_tot; i++ ) {
      q[i].r = r[i].r/beta_m1;
      q[i].i = r[i].i/beta_m1;
   }
   
   Hmult_spds_complete( r, d, shift, q); 

   alc.r = alc.i = 0;
   for ( i=0; i < d->n_ham_tot; i++ ) {
      alc = cadd( alc, cmul(cconj(q[i]),r[i]));
   } 
   al = alc.r;
   bl = 0;
   
   if (q_m2 && q_m1){
      for ( i=0; i < d->n_ham_tot; i++ ){
	 r[i].r -= al*q[i].r + beta_m1*q_m1[i].r;
	 r[i].i -= al*q[i].i + beta_m1*q_m1[i].i;
	 q_m2[i] = q_m1[i];
	 q_m1[i] = q[i];
	 bl += r[i].r*r[i].r + r[i].i*r[i].i;
      }
   } else{
      for ( i=0; i < d->n_ham_tot; i++ ){
	 r[i].r -= al*q[i].r + beta_m1*q_m1[i].r;
	 r[i].i -= al*q[i].i + beta_m1*q_m1[i].i;
	 bl += r[i].r*r[i].r + r[i].i*r[i].i;
      }
   }
   
   bl = sqrt(bl);
   
   *alpha = al;
   *beta = bl;
   
   return(0);
}


int h_ceigval ( qd_struct d, real s )
{
    int i, j, k, l, idum;
    real mxp, mxm, myp, mym, mzp, mzm, Tf;
    cvectr x;
    cmatrix h;
    FILE *fp;

    fp = fopen( "ceigval.out", "w");

    h = Coperator(d->geo.N_Cell);
    x = Cvectr(d->geo.N_Cell);

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

        /* Hamiltonian diagonal */
        h[l][l].r = Tf*(1/mxm+1/mxp+1/mym+1/myp+1/mzm+1/mzp) - s;
        h[l][l].i = -d->eta;

        /* Off-diagonal terms */
        idum = d->geo.ijk__cell[i-1][j][k];
        if ( idum >= 0 )
            h[l][idum].r = -Tf/mxm;

        idum = d->geo.ijk__cell[i+1][j][k];
	if ( idum >= 0 )
            h[l][idum].r = -Tf/mxp;

        idum = d->geo.ijk__cell[i][j-1][k];
        if ( idum >= 0 )
            h[l][idum].r = -Tf/mym;

        idum = d->geo.ijk__cell[i][j+1][k]; 
	if ( idum >= 0 )
            h[l][idum].r = -Tf/myp;

        idum = d->geo.ijk__cell[i][j][k-1];
	if ( idum >= 0 )
            h[l][idum].r = -Tf/mzm;

        idum = d->geo.ijk__cell[i][j][k+1];
	if ( idum >= 0 )
            h[l][idum].r = -Tf/mzp;

    }

    ceigval_ns_full( x, h );

    for ( i=0; i < d->geo.N_Cell; i++ )
        fprintf(fp, "%e %e\n", x[i].r, x[i].i);

    rm_cvectr(&x);
    rm_cmatrix(&h);
    fclose(fp);
    return(0);
}



/* Calculates neigv eigenvalues and eigenvectors of the hamiltonian in the
energy range [emin:emax].  max_iter is the maximum number of iterations.
Algorithms stops if 100 iterations occur with no new converged eigenvalues */
int readstates_driver_c( qd_struct d )
{
    int j, k, n_converged;
     char line[100], item[30], *token, *tail;
    rvectr wfc_norm;
    cvectr eig_val_empty=NULL;
    
    char filename[100];
    FILE *fvec;
    
#ifdef MPI3d
    die("Parallel version of  readstates_driver_c() isn't working.\n");
#endif

    n_converged = d->opt.ExecParam.ElCalc.ElAlg.ReadInEigen.StateFileNumb;
    d->wfc = Cmatrix(n_converged,d->n_ham_tot);
    wfc_norm = Rvectr(n_converged);
    eig_val_empty = Cvectr(n_converged);

    for (k=0;k<n_converged;k++){
       sprintf(filename,"%s_%d",d->opt.ExecParam.ElCalc.ElAlg.ReadInEigen.StateFileRead.c_str(),k+1);
       if (mpi_n3d_id==mpi_n3d_masterid) printf("Read file %s number %d/%d\n",filename,k+1,n_converged);
       if ( (fvec = fopen(filename,"r")) == 0 ){
          die("File \"%s\" not found in working directory",filename);
       }
      
       for ( j=0; j < d->n_ham_tot; j++ ){
          fgets(line,100,fvec);
          token = strtok(line,DELIMITERS);
          strcpy(item,token);
          d->wfc[k][j].r = strtod(item,&tail);
          token = strtok(NULL,DELIMITERS);
          strcpy(item,token);
          d->wfc[k][j].i = strtod(item,&tail);	  
          wfc_norm[k] += 
             (d->wfc[k][j].r*d->wfc[k][j].r + d->wfc[k][j].i*d->wfc[k][j].i);
       }
       fclose(fvec);
       wfc_norm[k] = sqrt(wfc_norm[k]);

#ifdef FIX_LATER
       if (d->opt.ExecParam.ElCalc.ElOut.Eigval){
          if (mpi_n3d_id==mpi_n3d_masterid) printf("Compute Rayleigh Quotient\n");
          complex cdummy = calc_rayleigh_compl( d->wfc[k], d->n_ham_tot, d );
          eig_val_empty[k]=cdummy;
	
          if (mpi_n3d_id==mpi_n3d_masterid){
             FILE *f_eval;
             printf("Eigenvalue (%g,%g)\n",cdummy.r,cdummy.i);
             if (d->opt.ExecParam.ElCalc.ElOut.Eigval ){
                char label[]="nd_eval.m";
                char *filename=NULL;
                filename = n3d_strdup_n(d->inputfile);
                n3d_FileTypeSet(&filename, label, TRUE);
	      
                if ((f_eval = fopen(filename, "a"))==NULL) {
                   die("Could not open the file \"%s\"",filename);
                }
                str_free(filename);
                fprintf( f_eval, "Eigval(%d) = %e;\n",k+1, 
                         cdummy.r*HAMILTONIAN_SCALE_VALUE);
                fclose(f_eval);
             }
          }

          {
             cvectr q = Cvectr(d->n_ham_tot);
             cvectr q_m1 = Cvectr(d->n_ham_tot);
             real alpha=0.0, beta=0.0, beta_m1=1.0;
             for ( j=0; j < d->n_ham_tot; j++ ){
                d->wfc[k][j].r /= wfc_norm[k];
                d->wfc[k][j].i /= wfc_norm[k];
             }
             sym_lanc_it_c( d, 0.0, &alpha, &beta, beta_m1, d->wfc[k], 
                            q, q_m1, NULL );
	  
             printf("alpha=%g   beta=%g\n",alpha,beta);
             cp_cvectr_dim(d->wfc[k],q,d->n_ham_tot);

             /* print kth eigvect to file now in case other eigvects are faulty  */
             print_psisq_gnu_pdb_k(d->wfc, k, d,0);
#ifdef NEW_CORRUPT
             print_psisq_gnu_pdb(d->wfc, k, 1, d,0,eig_val_empty);
#endif
             rm_cvectr(&q);
             rm_cvectr(&q_m1);
          }
       }
#endif // FIX_LATER      
    }

    if (!d->opt.ExecParam.ElCalc.ElOut.PsiSqr     &&
	!d->opt.ExecParam.ElCalc.ElOut.PsiSqr_cell &&
	!d->opt.ExecParam.ElCalc.ElOut.Psi4_IPR          &&
	!d->opt.ExecParam.ElCalc.ElOut.PsiSqrPdb         ){
	rm_cmatrix(&d->wfc);
	return 0;    
    }

    print_psisq_gnu_pdb(d->wfc,   n_converged,d,0);
#ifdef NEW_CORRUPT
    get_print_psisq(d->wfc,0, n_converged,d,0,eig_val_empty);
#endif

    rm_cvectr(&eig_val_empty);
    
    rm_cmatrix(&d->wfc);
    rm_rvectr(&wfc_norm);

    return(0);
}


/* Calculates neigv eigenvalues and eigenvectors of the hamiltonian in the
energy range [emin:emax].  max_iter is the maximum number of iterations.
Algorithms stops if 100 iterations occur with no new converged eigenvalues */
int ns_lanc_driver(qd_struct d) {
    return(0);
}


int sym_lanc_it( qd_struct d, real shift, real *alpha, real *beta, real beta_m1, rvectr r, rvectr q, rvectr q_m1, rvectr q_m2 ) {
    int i, n;
    real bl, al;

    n = d->geo.N_Cell;

#define LOCAL_REORTH
#ifdef LOCAL_REORTH
    al = 0;
    bl = 0;
    for ( i=0; i < n; i++ )
    {
        al += r[i]*q_m1[i];
        bl += r[i]*q_m2[i];
    }
    for ( i=0; i < n; i++ )
        r[i] -= al*q_m1[i] + bl*q_m2[i];
#endif

    for ( i=0; i < n; i++ )
        q[i] = r[i]/beta_m1;

    Hmult( r, d, shift, q ); 

    al = 0;
    for ( i=0; i < n; i++ )
        al += r[i]*q[i];

    bl = 0;
    for ( i=0; i < n; i++ )
    {
        r[i] -= al*q[i] + beta_m1*q_m1[i];
         q_m2[i] = q_m1[i];
        q_m1[i] = q[i];
        bl += r[i]*r[i];
    }

    bl = sqrt(bl);

    *alpha = al;
    *beta = bl;

    return(0);
}


int nonsym_lanc_it( qd_struct d, real shift, complex *alpha, complex *beta, complex *Gamma, cvectr r, cvectr s, cvectr q, cvectr q_m1, cvectr p, cvectr p_m1 )
{
    int i;
    complex bl, al, gl;

    bl.r = bl.i = 0;
    for ( i=0; i < d->geo.N_Cell; i++ )
        bl = cadd( bl, cmul(r[i],r[i])); 
    bl = cxsqrt(bl);

    gl.r = gl.i = 0;
    for ( i=0; i < d->geo.N_Cell; i++ )
        gl = cadd( gl, cmul(s[i],r[i])); 

    gl = cxdiv( gl, bl );

    for ( i=0; i < d->geo.N_Cell; i++ )
    {
        q[i] = cxdiv( r[i], bl );
        p[i] = cxdiv( s[i], gl );
    }

    Hmult_c( r, d, shift, q ); 
    Hmult_c( s, d, shift, p );

    al.r = al.i = 0;
    for ( i=0; i < d->geo.N_Cell; i++ )
        al = cadd( al, cmul( p[i], r[i]));

    for ( i=0; i < d->geo.N_Cell; i++ )
    {
        r[i] = csub( r[i], cadd( cmul(al,q[i]), cmul(gl,q_m1[i]) ) );
        q_m1[i] = q[i];
        s[i] = csub( s[i], cadd( cmul(al,p[i]), cmul(bl,p_m1[i]) ) );
        p_m1[i] = p[i];
    }

    *alpha = al;
    *beta = bl;
    *Gamma = gl;

    return(0);
}


#ifdef KEEP_AROUND
int bulk_bandstructure( qd_struct d )
{
    int i, j, nk=80, fold=2;
    ivectr BasisAll, spin;
    rvectr k;
    rmatrix nnv;
    cvectr eigv;
    cmatrix h;
    FILE *f_bulk;

    f_bulk = fopen("bulk.out", "w");

    BasisAll = Ivectr(2*d->N_BasisStates);
    spin = Ivectr(2*d->N_BasisStates);
    k = Rvectr(3);
    nnv = Rmatrix(d->geo.NeighborsMax(),3);
    eigv = Cvectr(2*d->N_BasisStates);
    h = Coperator(2*d->N_BasisStates);

    for ( i=0; i < d->N_BasisStates; i++ )
    {
        BasisAll[i] = d->BasisAnion[i];
	BasisAll[i+d->N_BasisStates] = d->BasisCation[i];
	spin[i] = spin[i+d->N_BasisStates] = d->Basis_Spin[i];
    }

    for ( i=0; i < d->geo.NeighborsMax(); i++ )
        for ( j=0; j < 3; j++ )
            nnv[i][j] = d->a * 
               ( d->geo.PositionInCell(d->geo.NbrCell(0,i,3),j) - 
                 d->geo.PositionInCell(0,j) );
    
    for ( i=0; i < nk; i++ )
    {
        k[0] = 2*Pi/d->a * (real)(i)/(real)(nk-1);
        tbh_const_bulk( h, nnv, BasisAll, spin, d->param, k );
        ceigval_ns_full( eigv, h );
        if ( i < nk/2 )
	    fprintf(f_bulk, "%e ", k[0]*d->a/Pi);
	else
	    fprintf(f_bulk, "%e ", 2.0 - k[0]*d->a/Pi);
	for ( j=0; j < 2*d->N_BasisStates; j++ )
            fprintf(f_bulk, " %e", eigv[j].r*HAMILTONIAN_SCALE_VALUE);
	fprintf(f_bulk, "\n");
    }

    rm_ivectr(&BasisAll);
    rm_ivectr(&spin);
    rm_rvectr(&k);
    rm_rmatrix(&nnv);
    rm_cvectr(&eigv);
    rm_cmatrix(&h);

    return(0);
}
#endif /* KEEP_AROUND */




