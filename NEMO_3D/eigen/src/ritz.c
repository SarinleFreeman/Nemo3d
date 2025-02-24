/*****************************************************************************
The JPL Parallel Eigensolvers package.
Copyright (C) 2002 California Institute of Technology (Caltech)

This file is part of
The Jet Propulsion Laboratory (JPL) Parallel Eigensolvers package.

This library is free software, which you can redistribute and/or modify
under the terms of the GNU Lesser General Public License as published by the
Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this library; see the file COPYING. If not, write to the
Free Software Foundation, Inc.,
59 Temple Place, Suite 330,
Boston, MA  02111-1307  USA

For additional information, please contact
  Gerhard Klimeck (gekco@jpl.nasa.gov)
  Fabiano Oyafuso (fabiano@jpl.nasa.gov)
  E. Robert Tisdale (E.Robert.Tisdale@jpl.nasa.gov)

Written by: Chris Bowen
            Gerhard Klimeck
            Fabiano Oyafuso
            E. Robert Tisdale

*****************************************************************************
$Source: /repo/eigen/src/ritz.c,v $
*****************************************************************************/

#include "ritz.h"

#undef COMPLEX_CG	

#define A_LARGE_DOUBLE 1.0e30		/* PLEASE EXPLAIN ME		*/

/* The following two functions set xi = d/dp (<p|(op-shift)^2|p>/<p|p>)	*/

void gradRayleighQuotient(
    nml_dvector	*p,
    int		sz,
    nml_dvector	*xi,
    void        *argList,
    int		itercount,
    double	*onevalue
    ) {
  nml_dcvector	*cp  = (nml_dcvector*)p;
  nml_dcvector	*cxi = (nml_dcvector*)xi;

  gradRayleighQuotient_cmplx(cp, sz/2, cxi, argList, itercount, onevalue);
  }

void gradRayleighQuotient_cmplx(
    nml_dcvector *p,
    int		 sz,
    nml_dcvector *xi,
    void        *argument,
    int		 itercount,
    double	 *onevalue
    ) {

  const int **argList = (const int**)argument;
  /* void (*op)(const int**, nml_dcscalar*, nml_dcscalar*) = argList[0];*/
  typedef void (*multiplier)(const int**, nml_dcscalar*, nml_dcscalar*);
  multiplier	op = (multiplier)(argList[0]);
  const int **opArgs = ((const int***) argument)[1];
  double	shift = *((double*) argList[2]);
  int		verbosity = *((int*) argList[3]);
  const char	*fil_trace = *((char**) argList[4]);

  int		i, j;
  double	nrm2_rq;
  double	xx = 0.0;
  const
  nml_dcscalar	czero = nml_dcmplx(0.0, 0.0);
  nml_dcscalar	xa2x = czero, rq;
  nml_dcvector	*x;
  nml_dcvector	*y;
  nml_dcvector	*z;

  static double	H2_min = A_LARGE_DOUBLE;

  nml_dcscalar	r, xax = czero;	/* used if GET_RAYLH_AS_WELL is defined	*/

  double	snd[5], rcv[5];
  double	ry, iy;
  register
  double	rx, ix, rz, iz;

  x = nml_dcv_new(sz);		/* copy input array p			*/
  y = nml_dcv_new(sz);		/* set to (H-shift)*x			*/
  z = nml_dcv_new(sz);		/* set to (H-shift)^2*x			*/

  for (i = 0; i < sz; ++i)
    x[i] = p[i];

  (*op)(opArgs, y, x);
  (*op)(opArgs, z, y);

  for (i = 0; i < sz; ++i) {
    rx = x[i].r, ix = x[i].i;
    rz = z[i].r, iz = z[i].i;

    xx += rx*rx + ix*ix;
	
    xa2x.r += rx*rz + ix*iz;
    xa2x.i += rx*iz - ix*rz;
	
#define	GET_RAYLH_AS_WELL
#ifdef	GET_RAYLH_AS_WELL
    ry = y[i].r, iy = y[i].i;
    xax.r += rx*ry + ix*iy;
    xax.i += rx*iy - ix*ry;
#endif/*GET_RAYLH_AS_WELL		*/
    }

#if	(defined MPI3d && !defined FAKE_MPI)
  snd[0] = xx;
  snd[1] = xa2x.r;
  snd[2] = xa2x.i;
  snd[3] = xax.r;
  snd[4] = xax.i;
  MPI_Allreduce(&snd[0], &rcv[0], 5, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  xx = rcv[0];
  xa2x.r = rcv[1];
  xa2x.i = rcv[2];
  xax.r = rcv[3];
  xax.i = rcv[4];
#endif/*(defined MPI3d && !defined FAKE_MPI)	*/

  rq.r = xa2x.r/xx;
  rq.i = xa2x.i/xx;
#ifdef	GET_RAYLH_AS_WELL
  r.r = xax.r/xx;
  r.i = xax.i/xx;

  *onevalue = r.r+shift;

  /* take care of some I/O						*/
  if (mpi_n3d_masterid == mpi_n3d_id) {
    FILE *f_trace;

    if (fil_trace) {
      if (NULL == (f_trace = fopen(fil_trace, "a"))) {
        error("Could not open the file \"%s\"", fil_trace);
        }
      else {
        fprintf(f_trace, "%d %g %e %e ", itercount, r.r+shift, r.r, rq.r);
        fclose(f_trace);
        }
      }

					/* PLEASE EXPLAIN ME		*/
    if ((fabs(r.i/r.r) > 1.0e-10) || (fabs(rq.i/rq.r) > 1.0e-10))
      printf(
	  "ACHTUNG:  Rayleigh quotient has a large imaginary value:  %e %e\n",
	  r.i, rq.i);

    if (verbosity) {
      printf("\n  <x|(A-s)^2|x> = %g,  ", rq.r);
      printf("  <x|A|x> = %g  target = %g\n", r.r+shift, shift);
      fflush(stdout);
      }
    }
#endif/*GET_RAYLH_AS_WELL		*/

  nrm2_rq = sqrt(rq.r*rq.r + rq.i+rq.i);
  /*
  printf("sqrt(rq.r*rq.r + rq.i+rq.i) = %g  H2_min = %g\n", nrm2_rq, H2_min);
  */

  if (nrm2_rq<H2_min) {

    if ((mpi_n3d_masterid == mpi_n3d_id) && verbosity) {
      printf("Found a new minimum from the previous value ");
      printf("%g to the new value %g\n", H2_min, nrm2_rq);
      fflush(stdout);
      }

    if (A_LARGE_DOUBLE == H2_min) {	/* PLEASE EXPLAIN ME		*/
      /* Do not write out the very first vector.			*/
      H2_min = nrm2_rq;
      }
    else {
#ifdef	RESTART_GEKCO	
#ifdef	COMPLEX_CG
      /* NOTE: parallelize complex version later.			*/
      cmacopt_args*	cmac_arg = get_cmacopt_args();
      nml_dcvector*		x_restart = cmac_arg->x_restart;
      if (!x_restart) {
        x_restart = nml_dcv_new(cmac_arg->n);
        cmac_arg->x_restart = x_restart;
        }

      for (j = 0; j < sz; ++j)
	x_restart[j] = x[j];
#else /*COMPLEX_CG			*/
#ifdef	MPI3d
      macopt_args*	mac_arg = get_macopt_args_par();
#else /*MPI3d				*/
      macopt_args*	mac_arg = get_macopt_args();
#endif/*MPI3d				*/
      nml_dvector*		x_restart = mac_arg->x_restart;
      nml_dcvector*		cx_restart = NULL;

      if (!x_restart) {
        x_restart = nml_dv_new(2*mac_arg->n);
        mac_arg->x_restart = x_restart;
        }
      cx_restart = (nml_dcvector*) x_restart;

      for (j = 0; j < sz; ++j)
        cx_restart[j] = x[j];
#endif/*COMPLEX_CG			*/
#endif/*RESTART_GEKCO			*/	    	
      H2_min = nrm2_rq;
      }
    }

#define	PROPER_NORMALIZATION 1
#ifdef	PROPER_NORMALIZATION
  /* Routine returns grad{ r = <x|H^2|x>/<x|x> } = (H^2 x - rx)*2/<x|x> */

  rz = 2.0/xx;

  ry = rq.r, iy = rq.i;	
  for (i = 0; i < sz; ++i) {
    /* xi[i] = Rcmul(2.0/xx, csub(z[i], cmul(rq, x[i]))); */
    rx = x[i].r, ix = x[i].i;
    xi[i].r = rz*(z[i].r - (ry*rx - iy*ix));
    xi[i].i = rz*(z[i].i - (ry*ix + iy*rx));	
    }
#else /*PROPER_NORMALIZATION		*/
  ry = rq.r, iy = rq.i;	
  for (i = 0; i < sz; ++i) {
    /* xi[i] = csub(z[i], cmul(rq, x[i]));				*/
    rx = x[i].r, ix = x[i].i;
    xi[i].r = z[i].r - (ry*rx - iy*ix);
    xi[i].i = z[i].i - (ry*ix + iy*rx);	
    }
#endif/*PROPER_NORMALIZATION		*/


#define NORMALIZE_P
#ifdef NORMALIZE_P
  double nrm = 1.0 / sqrt(xx);
  for (i = 0; i < sz; ++i) {
     p[i].r *= nrm;
     p[i].i *= nrm;
     xi[i].r *= nrm;
     xi[i].i *= nrm;
  }
#endif /* NORMALIZE_P */

#ifdef	LOW_LEVEL_PRINT_RAYLH
#ifdef	MPI3d
  if (mpi_n3d_masterid == mpi_n3d_id)		/* PLEASE EXPLAIN ME		*/
#endif/*MPI3d				*/
    printf("\nRayleigh = %e %e", rq.r, rq.i);
#endif/*LOW_LEVEL_PRINT_RAYLH		*/

  nml_dcv_delete(x);
  nml_dcv_delete(y);
  nml_dcv_delete(z);
  }

/* -------------------------------------------------------------------- */
#ifdef	_FIX_THIS_STUFF
/* -------------------------------------------------------------------- */

/* solve (H-shift)^2*x = r using CG on normal equations			*/
static void dum_linSolv(
    nml_dcvector*	x,
    qd_struct	d,
    double	shift,
    nml_dcvector*	r,
    double	tol,
    double*	xx,
    double*	xHx
    ) {
#define	MAX_ITER_LINSOLV 150000

  int		i, j;
  int		sz = d->seg_ln[mpi_n3d_id];
  double	alpha, beta, r2_loc, r2, r2_old, pHp, pHp_loc, rcv[3], snd[3];

  nml_dcvector*	p   = nml_dcv_new(sz);
  nml_dcvector*	Hp  = nml_dcv_new(sz);
  nml_dcvector*	HHp = nml_dcv_new(sz);

  FILE*		fp = fopen("cg_conv", "a");

  if (mpi_n3d_masterid == mpi_n3d_id)
    printf("entering dum_linSolv(shift = %e)\n", shift);

  for (j = 0; j < sz; ++j) {
    p[j].r = r[j].r;
    p[j].i = r[j].i;
    x[j].r = x[j].i = 0.0;
    Hp[j].r = Hp[j].i = 0.0;
    HHp[j].r = HHp[j].i = 0.0;
    }

  for (i = 0; i < MAX_ITER_LINSOLV; ++i) {
 /* Hmult_spds_complete(Hp, d, shift, p);	*/
    Hmult_spds_complete(HHp, d, shift, p);
    Hmult_spds_complete(Hp, d, shift, HHp);

    r2_loc = 0.0;
    pHp_loc = 0.0;
    for (j = 0; j < sz; ++j) {
      r2_loc  +=  r[j].r*r[j].r +  r[j].i*r[j].i;
      pHp_loc += Hp[j].r*p[j].r + Hp[j].i*p[j].i;
      }

    snd[0] =  r2_loc;
    snd[1] = pHp_loc;

    Allreduce_MPI_sp(&snd[0], &rcv[0], 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

     r2 = rcv[0];
    pHp = rcv[1];

    alpha = r2/pHp;

    for (j = 0; j < sz; ++j) {
      x[j].r += alpha*p[j].r;
      x[j].i += alpha*p[j].i;

      r[j].r -= alpha*Hp[j].r;
      r[j].i -= alpha*Hp[j].i;
      }

    r2_old = r2;
    r2_loc = 0.0;
    for (j = 0; j < sz; ++j) {
      r2_loc += r[j].r*r[j].r + r[j].i*r[j].i;
      }

    snd[0] = r2_loc;

    Allreduce_MPI_sp(&snd[0], &rcv[0], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    r2 = rcv[0];

    beta = r2/r2_old;

    for (j = 0; j < sz; ++j) {
      p[j].r = r[j].r + beta*p[j].r;
      p[j].i = r[j].i + beta*p[j].i;
      }

    if (r2 < tol)
      break;

    /* Print out some info every 100 iterations.			*/
    if (0 == i%100) {
      Hmult_spds_complete(Hp, d, 0.0, x);
      pHp_loc = 0.0;
       r2_loc = 0.0;
      for (j = 0; j < sz; ++j) {
	pHp_loc += Hp[j].r*x[j].r + Hp[j].i*x[j].i;
	 r2_loc +=  x[j].r*x[j].r +  x[j].i*x[j].i;
	}
      snd[0] = r2_loc;
      snd[1] = pHp_loc;

      Allreduce_MPI_sp(&snd[0], &rcv[0], 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

       r2 = rcv[0];
      pHp = rcv[1];
      if (mpi_n3d_masterid == mpi_n3d_id) {
	printf("  %23.15e %23.15e %23.15e\n", pHp/r2, pHp, r2_old);
	fprintf(fp, "  %23.15e %23.15e %23.15e\n", pHp/r2, pHp, r2_old);
	fflush(stdout);
	}
      }
    }

  {
    Hmult_spds_complete(Hp, d, 0.0, x);
    pHp_loc = 0.0;
     r2_loc = 0.0;
    for (j = 0; j < sz; ++j) {
      pHp_loc += Hp[j].r*x[j].r + Hp[j].i*x[j].i;
       r2_loc +=  x[j].r*x[j].r +  x[j].i*x[j].i;
      }
    snd[0] =  r2_loc;
    snd[1] = pHp_loc;

    Allreduce_MPI_sp(&snd[0], &rcv[0], 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

     *xx = rcv[0];
    *xHx = rcv[1];
    }

  fclose(fp);
  nml_dcv_delete(p);
  nml_dcv_delete(Hp);
  }

/* (p, eval_loc) is initial eigenpair guess				*/
double dum_invIter(
    qd_struct	d,
    nml_dcvector*	p,
    double	eval_loc,
    double	tol,
    int		itnum
    ) {
  double	eval, yy, yHy, deval = 1e20;	/* PLEASE EXPLAIN ME	*/
  int		i, j;

  nml_dcvector*	y     = nml_dcv_new(d->seg_ln[mpi_n3d_id]);
  nml_dcvector*	p_old = nml_dcv_new(d->seg_ln[mpi_n3d_id]);
  nml_dcvector*	p_new = nml_dcv_new(d->seg_ln[mpi_n3d_id]);

  for (i = 0; i < d->seg_ln[mpi_n3d_id]; ++i) {
    p_old[i] = p[i];
    }

  writeSiteInfo(d, "basis", "old.jnk", "{real old.r, real old.i}", (void*) p_old);

  for (j = 0; j < 10; ++j) {
    dum_linSolv(y, d, eval_loc, p_old, tol, &yy, &yHy);
    for (i = 0; i < d->seg_ln[mpi_n3d_id]; ++i) {
      p_old[i].r = y[i].r/yy;
      p_old[i].i = y[i].i/yy;
      }
    deval = yHy/yy - eval_loc;
    eval_loc = yHy/yy;

    if (fabs(deval) < 1.0e-6) {		/* PLEASE EXPLAIN ME		*/
      if (mpi_n3d_masterid == mpi_n3d_id)
	printf("Terminating Rayliegh iteration.\n");
      break;
      }
    }

  writeSiteInfo(d, "basis", "test.jnk", "{real test.r, real test.i}", (void*) p_old);
  Hmult_spds_complete(p_new, d, eval_loc, y);
  writeSiteInfo(d, "basis", "new.jnk", "{real new.r, real new.i}", (void*) p_new);

  printPsisq(d, p_old, itnum);

  nml_dcv_delete(y);
  nml_dcv_delete(p_old);
  nml_dcv_delete(p_new);

  eval = eval_loc;
  return eval;
  }

double calc_rayleigh_sq(nml_dvector* p, void* data) {
  int		i;
  double	xx;
  nml_dcscalar	xa2x, rq, cdum, cdumr, cdumi;
  nml_dcvector*	x, y, z;
  qd_struct	d;
  static int	itercount = 0;

  ++itercount;
  d = (qd_struct)(data);
  x = nml_dcv_new(d->n_ham_tot);
  y = nml_dcv_new(d->n_ham_tot);
  z = nml_dcv_new(d->n_ham_tot);

  for (i = 0; i < d->n_ham_tot; ++i) {
    x[i].r = p[i];
    x[i].i = p[i+d->n_ham_tot];
    }

  Hmult_spds_complete(y, d, d->I_N3D->Rayleigh->EigValGuess[0], x);
  Hmult_spds_complete(z, d, d->I_N3D->Rayleigh->EigValGuess[0], y);

  xx = 0;
  xa2x.r = xa2x.i = 0;
  for (i = 0; i < d->n_ham_tot; ++i) {
    xx     += x[i].r*x[i].r + x[i].i*x[i].i;
    xa2x.r += x[i].r*z[i].r + x[i].i*z[i].i;
    xa2x.i += x[i].r*z[i].i - x[i].i*z[i].r;
 /* xa2x.r += y[i].r*y[i].r + y[i].i*y[i].i; */
  }

  rq.r = xa2x.r/xx;
  rq.i = xa2x.i/xx;

  /* Need complex version of conjugate gradient */
  for (i = 0; i < d->n_ham_tot; ++i) {
    cdum  = csub(z[i], cmul(rq, x[i]));
    cdumr = csub(z[i], cmul(rq, Complex(x[i].r, 0.0)));
    cdumi = csub(z[i], cmul(rq, Complex(0.0, x[i].i)));
    d->drq[i] = cdum;
    }

#ifdef	LOW_LEVEL_PRINT_RAYLH
  printf("\nIter = %d Rayleigh = %e %e", itercount, rq.r, rq.i);
#endif/*LOW_LEVEL_PRINT_RAYLH	*/

  nml_dcv_delete(x);
  nml_dcv_delete(y);
  nml_dcv_delete(z);
  return rq.r;
  }

/* Computes Rayleigh quotients = <x|H|x>/<x|x>
   Since H is hermitian, this value should be mostly real.		*/

nml_dcscalar calc_rayleigh_compl(nml_dcvector* x, int sz, qd_struct d) {
  int		i;
  double	xx;
  nml_dcscalar	xhx, rq;
  nml_dcvector*	y = nml_dcv_new(sz);

#ifdef	MPI3d
  double rcv[3], snd[3];
#endif/*MPI3d		*/

  Hmult_spds_complete(y, d, 0.0, x);

  xx = xhx.r = xhx.i = 0;
  for (i = 0; i < sz; ++i) {
    xx += x[i].r*x[i].r + x[i].i*x[i].i;
    xhx = cadd(xhx, cmul(cconj(x[i]), y[i]));
    }

#if	(defined MPI3d && !defined FAKE_MPI)
  snd[0] = xx;
  snd[1] = xhx.r;
  snd[2] = xhx.i;
  MPI_Allreduce(&snd[0], &rcv[0], 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  rq.r = rcv[1]/rcv[0];
  rq.i = rcv[2]/rcv[0];
#else /*(defined MPI3d && !defined FAKE_MPI)	*/
  rq.r = xhx.r/xx;
  rq.i = xhx.i/xx;
#endif/*(defined MPI3d && !defined FAKE_MPI)	*/

  nml_dcv_delete(y);

  return rq;
  }

double calc_rayleigh(nml_dvector* p, void *data) {
  int		i;
  double	xx;
  nml_dcscalar	xa2x, rq, cdum;
  nml_dcvector*	x, y, z;
  qd_struct	d;

  d = (qd_struct)(data);
  x = nml_dcv_new(d->n_ham_tot);
  y = nml_dcv_new(d->n_ham_tot);
  z = nml_dcv_new(d->n_ham_tot);

  for (i = 0; i < d->n_ham_tot; ++i) {
    x[i].r = p[i];
    x[i].i = p[i+d->n_ham_tot];
    }

  Hmult_spds_complete (y, d, d->I_N3D->Rayleigh->EigValGuess[0], x);

  xx = 0;
  xa2x.r = xa2x.i = 0;
  for (i = 0; i < d->n_ham_tot; ++i) {
    xx     += x[i].r*x[i].r + x[i].i*x[i].i;
    xa2x.r += x[i].r*z[i].r + x[i].i*z[i].i;
    xa2x.i += x[i].r*z[i].i - x[i].i*z[i].r;
    }

  rq.r = xa2x.r/xx;
  rq.i = xa2x.i/xx;

  for (i = 0; i < d->n_ham_tot; ++i) {
    cdum = csub(z[i], cmul(rq, x[i]));
    d->drq[i] = cdum;
    }

#ifdef	LOW_LEVEL_PRINT_RAYLH
   printf("\nRayleigh = %e %e", rq.r, rq.i);
#endif/*LOW_LEVEL_PRINT_RAYLH	*/

  nml_dcv_delete(x);
  nml_dcv_delete(y);
  nml_dcv_delete(z);
  return rq.r;
  }

int min_rayleigh(qd_struct d) {
  int		states;
  cmacopt_args	camac;
  macopt_args	amac;
  char		filename[100];

  int		n_tot = d->n_ham_tot;
  char*		filename_dyn = NULL;

#ifdef	MPI3d
  char		outfile[100];
  int		i, proc;
  nml_dcvector*	p = nml_dcv_new(d->seg_ln[mpi_n3d_id]);
  nml_dcvector*	p_full = nml_dcv_new(d->n_ham_tot);
#else /*MPI3d			*/
  nml_dcvector*	p = nml_dcv_new(d->n_ham_tot);
#endif/*MPI3d			*/

  d->rc  = nml_dcv_new(d->nvmap[mpi_n3d_id]*d->N_BasisStates);
  d->drq = nml_dcv_new(d->n_ham_tot);

  cmacopt_defaults(&camac,
      d->I_N3D->Rayleigh->MaxIter, d->I_N3D->Rayleigh->verbose,
      d->I_N3D->Rayleigh->tol, d->I_N3D->Rayleigh->ConvCrit);
   macopt_defaults(&amac,
      d->I_N3D->Rayleigh->MaxIter, d->I_N3D->Rayleigh->verbose,
      d->I_N3D->Rayleigh->tol, d->I_N3D->Rayleigh->ConvCrit);

  for (states = 0; states < nml_dcextent(d->I_N3D->Rayleigh->EigValGuess); ++states) {
    d->I_N3D->Rayleigh->RunningIndex = states;

    if (d->I_N3D->Rayleigh->StartFileRead) {

#ifdef	MPI3d
      if (mpi_n3d_masterid == mpi_n3d_id) {
	if (1 == nml_dcextent(d->I_N3D->Rayleigh->EigValGuess)) {
	  sprintf(filename, "%s", d->I_N3D->Rayleigh->StartFileRead);
	  }
	else {
	  sprintf(filename, "%s_%d", d->I_N3D->Rayleigh->StartFileRead, states+1);
	  }

	printf("Read the complete initial guess\n"), fflush(stdout);;
	initStartVect(filename, d, p_full);
	printf("Send initial guess to CPU"), fflush(stdout);;
	for  (proc = 1; proc < mpi_n3d_numprocs; ++proc) {
	  /* might use non-blocking send?				*/
	  printf(" %d", proc), fflush(stdout);
	  MPI_Send(&p_full[d->seg_s[proc]], 2*d->seg_ln[proc], MPI_DOUBLE,
	      proc, mpi_n3d_masterid, MPI_COMM_WORLD);
	  }
	for (i = 0; i < d->seg_ln[mpi_n3d_id]; ++i) {
	  p[i] = p_full[i];
	  }
	printf(" Done\n"), fflush(stdout);
	}
      else {
	MPI_Recv(p, 2*d->seg_ln[mpi_n3d_id], MPI_DOUBLE, mpi_n3d_masterid,
	    MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_n3d_status);
	}
      MPI_Barrier(MPI_COMM_WORLD);
#else /*MPI3d				*/
      if (1 == nml_dcextent(d->I_N3D->Rayleigh->EigValGuess)) {
	sprintf(filename, "%s", d->I_N3D->Rayleigh->StartFileRead);
	}
      else {
	sprintf(filename, "%s_%d", d->I_N3D->Rayleigh->StartFileRead, states+1);
	}

      initStartVect(filename, d, p);
#endif/*MPI3d				*/
      /* Print rayleigh coefficient starting from the read in guess.	*/
      {
	nml_dcscalar	cdummy = calc_rayleigh_compl(p, nml_dcv_extent(p), d);
	if (mpi_n3d_masterid == mpi_n3d_id)
	  printf("Rayleigh(INITIAL) = (%g, %g)\n", cdummy.r, cdummy.i);
	fflush(stdout);
	}

      }
    else {
      /* At least, as a fall-back position,
	 provide a non-zero guess to the eigenvector:			*/
      p[d->seg_ln[mpi_n3d_id]/2].r = 1.0;
      }

#ifdef	_TEST_CG
    dum_invIter(d, p, 0.856572, 1e-4, states);/* PLEASE EXPLAIN ME	*/
    return 0;
#endif/*_TEST_CG			*/

#ifdef	MPI3d
    if (mpi_n3d_masterid == mpi_n3d_id)
      printf("Start the conjugate gradient work now.\n");
#endif/*MPI3d				*/

    if (mpi_n3d_masterid == mpi_n3d_id) {
      if (d->I_N3D->Out1->EigenEnergyCGiter) {
	filename_dyn = strdup_n(d->inputfile);
	FileTypeSet(&filename_dyn, "nd_eval_iter", TRUE);
	amac.filename_its = filename_dyn;
	}
      if (d->I_N3D->Out1->Eigval) {
	filename_dyn = strdup_n(d->inputfile);
	FileTypeSet(&filename_dyn, "nd_eval", TRUE);
	amac.filename_last = filename_dyn;
	}
      }
#ifdef	COMPLEX_CG	
      cmacoptII(p, n, (gradRayleighQuotient_cmplx), d, &camac);
#else /*COMPLEX_CG			*/
#ifdef	MPI3d
    macoptII_par((nml_dvector*) p, 2*d->seg_ln[mpi_n3d_id], 2*n_tot,
	(gradRayleighQuotient), d, &amac);

    if (d->I_N3D->Out1->FullEigvec) {
      sprintf(outfile, "nd_evec_%d", states+1);
      filename_dyn = strdup_n(d->inputfile);
      FileTypeSet(&filename_dyn, outfile, TRUE);
      writeSiteInfo(d, "basis", filename_dyn, "{real wf.r, real wf.i}", (void*) p);
      mfree(filename_dyn);
      }
#else /*MPI3d				*/
      macoptII((nml_dvector*) p, 2*n_tot, (gradRayleighQuotient), d, &amac);
#endif/*MPI3d				*/
#endif/*COMPLEX_CG			*/

    if (amac.filename_its)
      mfree(amac.filename_its),  amac.filename_its  = NULL;
    if (amac.filename_last)
      mfree(amac.filename_last), amac.filename_last = NULL;

#if	0
    {
      nml_dcscalar	cdummy = calc_rayleigh_compl(p, nml_dcextent(p), d);
      if (mpi_n3d_masterid == mpi_n3d_id) {
	if (d->I_N3D->Out1->Eigval) {
	  FILE*	fp_eval = NULL;
	  filename_dyn = strdup_n(d->inputfile);
	  FileTypeSet(&filename_dyn, "nd_eval", TRUE);
	  fp_eval = fopen(filename_dyn, "a");
	  fprintf(fp_eval, "%g\n", cdummy.r);
	  fclose(fp_eval);
	  }
	printf("Rayleigh(FINAL) = (%g, %g)\n", cdummy.r, cdummy.i);
	printf("Printing out the %d th eigenstate.\n", states);
	}
      }
#endif/*0				*/

    printPsisq(d, p, states);
    }

  nml_dcv_delete(d->rc);
  nml_dcv_delete(d->drq);

  nml_dcv_delete(p);
#ifdef	MPI3d
   nml_dcv_delete(p_full);
#endif/*MPI3d				*/

  return 0;
  }

#endif /* FIX_THIS_STUFF */


