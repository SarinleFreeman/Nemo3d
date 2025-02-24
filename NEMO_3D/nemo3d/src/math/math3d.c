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
$Header: /repo/nemo3d/src/math/math3d.c,v 1.3 2004/12/02 12:09:56 marek Exp $ 
*****************************************************************************/

#include "math3d.h"



/* Assorts the first n entries of rvector ar[] in ascending order */
void eigvec_srt( cvectr vec, cmatrix cmat, int n )
{
    int i, j;
    complex a;
    int nmax=cdim(cmat);
    cvectr vec_dum=Cvectr(nmax);

    for ( j=1; j < n; j++ )
    {
        a = vec[j];
	cp_carray2carray(cmat[j],vec_dum,nmax);
        i = j-1;
        while ( i >= 0 && vec[i].r > a.r )
	{
            vec[i+1] = vec[i];
	    cp_carray2carray(cmat[i],cmat[i+1],nmax);
            i--;
        }
        vec[i+1] = a;
	cp_carray2carray(vec_dum,cmat[i+1],nmax);
    }

    rm_cvectr(&vec_dum);
}

real invit ( rtridiag h, rvectr x, rvectr y, real lambda )
{
    int i, ret;
    int dim = nml_dv_extent(x);
    real max, nrm2, lambda0, de;
    int iter=0, itermax=30;
    
    rvectr y_old = Rvectr(dim);

    /* copy:  hs <-- h */
    rtridiag hs = Rtridiag(dim);
    for ( i=0; i < dim; i++ ) {
       hs[0][i] = h[0][i];
       hs[1][i] = h[1][i];
       hs[-1][i] = h[-1][i];
    }

    de = lambda;
    lambda0 = lambda;
    max = 1e10;
    

    do {
       iter++;
	
       for ( i=0; i < dim; i++ )
	  hs[0][i] -= de;

       ret = solve_rtridiag_ii(x, hs, y );

       for ( i=0; i < dim; i++ )
	  y_old[i] = y[i];

       /* y <-- x / |x| */
       nrm2 = 0.0;
       for ( i=0; i < dim; i++ )
	  nrm2 += x[i]*x[i];
       for ( i=0; i < dim; i++ )
	  y[i] = x[i] / sqrt(nrm2);
       
       /* max <-- infinity norm of |y| - |y_old| */
       max=0.0;
       for ( i=0; i < dim; i++ ) {
	  if ( fabs(fabs(y[i])-fabs(y_old[i])) > max) {
	     max = fabs(fabs(y[i])-fabs(y_old[i]));
	  }
       }
       
       lambda0 = lambda;
       lambda = 0.0;
       for ( i=0; i < dim; i++ )
	  lambda += y[i]*y[i]*h[0][i];
       for ( i=1; i < dim; i++ )
	  lambda += 2*y[i-1]*y[i]*h[-1][i]; /* assume symmetric matrix */
       
       de = lambda - lambda0;
       
       printf("%25.18e   %25.18e  %e\n", lambda, de, max);
    }
    while(max > 1e-4 && iter<itermax);
    
    if (iter>=itermax)
	printf("Bailed out on inverse iteration\n");
    
    rm_rtridiag(&hs);
    rm_rvectr(&y_old);

    return lambda;
}


/* Called by eigval_s_tridi_QD() below, which is nolonger needed. */
int td_ql_im_QD (rvectr d, rvectr e, roperator z, int n) {
    return(0);
}


/*
 * Called from eigsys3d_par.[ch]::srch_converged_eigval_par(),
 * which is now made obsolete by LAPACK.
 */
rvectr eigval_s_tridi_QD (rvectr d, rtridiag t, int n) {
    int i;
    rvectr e = Rvectr(n);

    for (i = 0; i < n; i++)
    {
	d[i] = t[ 0][i];
	e[i] = t[-1][i];
    }

    td_ql_im_QD (d, e, NULL, n);
    rm_rvectr(&e);

    return(0);
}


