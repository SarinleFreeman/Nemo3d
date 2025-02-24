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
$Source: /repo/eigen/src/macopt.c,v $
*****************************************************************************/

/*
This file is part of
The Jet Propulsion Laboratory (JPL) Parallel Eigensolvers package.

macopt library source file	release 1.1	gradient-based optimizer

Copyright	(c) 2002	David J.C. MacKay

This library is free software; you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License
as published by the Free Software Foundation;
either version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this library; if not, write to the Free Software Foundation, Inc.,
59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

GNU licenses are here:

    http://www.gnu.org/licenses/licenses.html

Author contact details are here:

    http://www.inference.phy.cam.ac.uk/mackay/c/macopt.html
    mackay@mrao.cam.ac.uk

If you find macopt useful,
please feel free to make a donation to support David MacKay's research group.

Written by David J.C. MacKay
Modified by Chris Bowen, Gerhard Klimeck, Fabiano A. Oyafuso and
E. Robert Tisdale
*/

/*
   Please do not use macopt without understanding a little about how it works;
   there are some control parameters which the user MUST set!

   David MacKay's optimizer, based on conjugate gradient ideas, 
   but using bracketing of the zero of the inner product 

             (gradient).(line_search_direction)

   to do the line minimization. Only derivative calculations are required.
   The length of the first step in the line search (often set to "1.0"
   in other code)is adapted here so that, if 0.00001 is a better step size, 
   it soon cottons on to that and saves ~log(10000)bracketing operations.
   The result is that (with rich set to 0)the program can use 
   as few as 2 derivatives per line search. (If rich is set to 1, it does 
   an extra derivative calculation at the beginning of each line search 
   making a minimum of 3 per line search. Set rich = 0 if you think 
   that the surface is locally quite quadratic.)If the program does average 
   2 derivatives per line search then it must be superior to most cg methods 
   including use of Rbackprop (which costs 2 derivatives straight off)

   A possible modification: where the function can be returned at same 
   time as the dfunction --- there is nothing clever to do with the value
   but it could be used as a sanity check and a convergence criterion. 

   See

	http://www.inference.phy.cam.ac.uk/mackay/c/macopt.html 

   for further discussion.

   NB: The value of "tol" is totally arbitrary
   and must be set by you to a value that works well for your problem. 
   It depends completely on the typical value of the gradient/step size. 

   Tol specifies a magnitude of gradient at which a halt is called. 
   or a step size.

   This program MINIMIZES a function.
*/

#include "macopt.h"

void macoptII(
  nml_dvector*	p,	/* starting vector				*/
  int		n,	/* number of dimensions				*/
  void	(*dfunc)(nml_dvector *, nml_dvector *, void *, int, double *), 
			/* evaluates the gradient of the optimized function */
  void*	dfunc_arg,	/* arguments that get passed to dfunc		    */
  macopt_args*	a	/* structure in which optimizer arguments stored    */
 ) {			/* Note, (*func)(nml_dvector*, void*) is not used   */
  int j;
  double gg, gam, dgg;
  nml_dvector *g, *h, *xi;
  int end_if_small_grad = 1 - a->end_if_small_step;
  double step = 0.0, tmpd;

  /* A total of 7 nml_dvector* 0..n-1 are used by this optimizer. 
     p		is provided when the optimizer is called 
     pt		is used by the line minimizer as the temporary vector. 
		this could be cut out with minor rewriting, using p alone
		g, h and xi are used by the cg method as in NR -
	       	could one of these be cut out?
     the line minimizer uses an extra gx and gy to evaluate two gradients. 
  */
  
  a->n  = n; 
  a->g  = nml_dv_new(n); g = a->g;	/* vectors as in NR code */
  a->h  = nml_dv_new(n); h = a->h;    /*                       */
  a->xi = nml_dv_new(n); xi = a->xi;   /*                       */
  a->pt = nml_dv_new(n); /* scratch vector for sole use of macprod */
  a->gx = nml_dv_new(n); /* scratch gradients             */
  a->gy = nml_dv_new(n); /* used by maclinmin and macprod */

#ifdef RESTART_GEKCO
  a->x_restart = NULL;
  set_macopt_args(a);
#endif	
  a->its = 0;
  a->track_one_value = 0.0;

  printf("First Gradient Evaluation before the start of the algorithm.\n");
  fflush(stdout);
  (*dfunc)(p, xi, dfunc_arg, a->its, &a->track_one_value);
  macopt_restart (a, 1); 
  for (a->its = 1; a->its <= a->itmax; a->its ++) {

    for (gg = 0.0, j = 0; j < n; ++j) 
      gg += g[j]*g[j];          /* find the magnitude of the old gradient */
    a->gtyp = sqrt (gg / (double)(n)); 

    if (a->verbose > 0) 
      if (a->end_if_small_step)
	printf ("mac_it %d of %d : converge on (step = %6.3g with tol = %6.3g), monitor gradient = %6.3g, grad_tiny = %6.3g: ", a->its, a->itmax, step, a->tol,  gg, a->grad_tol_tiny), fflush(stdout);
      else
	printf ("mac_it %d of %d : converge on (gradient = %6.3g with tol = %6.3g grad_tiny = %6.3g) monitor step = %6.3g step_tiny = %6.3g ", a->its, a->itmax, gg, a->tol, a->grad_tol_tiny, step, a->step_tol_tiny), fflush(stdout);
    

    if ((end_if_small_grad && (gg <= a->tol)) 
	|| (gg <= a->grad_tol_tiny)) {
      if (gg <= a->grad_tol_tiny){
	printf("\nAchieved convergence though gradient!! Gradient is less than grad_tol_tiny\n"), fflush(stdout);
	}
      else {
	printf("\nAchieved convergence though gradient!!\n"), fflush(stdout);
	}
      printf ("mac_it %d of %d : converge on (gradient = %6.3g with tol = %6.3g grad_tiny = %6.3g) monitor step = %6.3g step_tiny = %6.3g ", a->its, a->itmax, gg, a->tol, a->grad_tol_tiny, step, a->step_tol_tiny), fflush(stdout);

      macopt_free (a);
      if (a->verbose > 0) printf ("\n");
      return;
      }

    step = maclinminII (p, dfunc, dfunc_arg, a); 

    if (a->filename_its) {
      FILE *fp = NULL;
      fp = fopen(a->filename_its, "a");
      fprintf(fp, "%d %g %g %g\n", a->its, a->track_one_value, step, gg);
      fclose(fp);
      }

    if (a->restart == 0) {
      if (a->verbose > 1) printf (" (step %9.5g)", step);
      if (a->verbose > 0) printf ("\n");
      if ((a->end_if_small_step  && (step <= a->tol)) 
	  || (step <= a->step_tol_tiny)) {
	if (step <= a->step_tol_tiny) {
	  printf("\nAchieved convergence though step size!! Gradient is less than step_tol_tiny\n"), fflush(stdout);
	  }
	else {
	  printf("\nAchieved convergence though step size!!\n"), fflush(stdout);
	  }
	printf ("mac_it %d of %d : converge on (step = %6.3g with tol = %6.3g), monitor gradient = %6.3g, grad_tiny = %6.3g: ", a->its, a->itmax, step, a->tol,  gg, a->grad_tol_tiny), fflush(stdout);
	macopt_free (a);
	return;
	}
      }

    /* if we are feeling rich, evaluate the gradient at the new `minimum'.
     * alternatively, linmin has already estimated this gradient
     * by linear combination of the last two evaluations and left it in xi */
    if (a->rich || a->restart) { 
      (*dfunc)(p, xi, dfunc_arg, a->its, &a->track_one_value); 
      }
    if (a->restart) {
      fprintf(stderr, "Restarting macopt within the main iteration loop\n"), fflush(stderr); 
      macopt_restart (a, 0);
      /* This is not quite right.  It should distinguish whether
       * there was an overrun indicating that the value of lastx
       * needs to be bigger / smaller; in which case
       * resetting lastx to default value may be a bad idea, 
       * giving an endless loop of resets 
       */
      }
    else {
      dgg = 0.0;
      for (j = 0; j < n; ++j) {
	dgg += (xi[j] + g[j]) * xi[j];
	}
      gam = dgg/gg;
      for (tmpd = 0.0, j = 0; j < n; ++j) {
	g[j] = -xi[j];		/* g stores (-) the most recent gradient */
	xi[j] = h[j] = g[j] + gam*h[j];
	/* h stores xi, the current line direction			*/
	/* check that the inner product of gradient and line search is < 0 */
	tmpd -= xi[j] * g[j]; 
	}

      if (tmpd > 0.0  || a->verbose > 2) {
	fprintf(stderr, "new line search has inner prod %9.4g, this will result in a restart\n", tmpd), fflush(stderr); 
	}
      if (tmpd > 0.0) { 
	if (a->rich == 0) {
	  fprintf (stderr, "Setting rich to 1; "); 
	  a->rich = 1; 
	  }
	a->restart = 2; /* signifies that g[j] = -xi[j] is already done */
	fprintf(stderr, "Restarting macopt (2) due to positive inner product within the line search\n\\n"), fflush(stderr); 
	macopt_restart (a, 0);
	}
      }
    }
  fprintf(stderr, "Reached iteration limit in macopt; continuing.\n"), fflush(stderr); 
  macopt_free (a);	
  return;
  } /* NB this leaves the best value of p in the p vector,
       but the function has not been evaluated there if rich = 0	*/

/* MaclinminII: evaluate gradient at a sequence of points
 * and calculate the inner product with the line search direction.
 * Continue until a bracketing is achieved (i.e a change in sign).
 */
double maclinminII(
    nml_dvector*	p, 
    void	(*dfunc)(nml_dvector *, nml_dvector *, void *, int, double *)
    /* evaluates the gradient */,
    void*	arg,
    macopt_args*	a
    ) {

  int n = a->n; 

  double x, y;
  double s, t, m;
  int    its = 1, i;
  double step, tmpd; 
  double  *gx = a->gx, *gy = a->gy;

  /* at x = 0, the gradient (uphill) satisfies s < 0			*/
  if (a->verbose > 2) { /* check this is true:				*/
    /* (no need to do this really as it is already checked
     * at the end of the main loop of macopt)				*/
  /*
#define TESTS 5
    x = a->lastx / a->gtyp;
    fprintf (stderr, "inner product at:\n"); 
    for (i = -TESTS; i <= TESTS; i += 2) {
      step = x * 2.0 * (double) i / (double) TESTS; 
      fprintf (stderr, "%9.5g %9.5g\n", step,
	  tmpd = macprodII (p, gy, step, dfunc, arg, a)); 
      }
  */
    fprintf (stderr, "inner product at zero = %9.4g\n",
	tmpd = macprodII (p, gy, 0.0, dfunc, arg, a)), fflush(stderr); 
    if (tmpd > 0.0) { 
      a->restart = 1; 
      return 0.0; 
      }
    }

  x = a->lastx / a->gtyp;
  s = macprodII (p, gx, x, dfunc, arg, a); 
  
  if (s < 0)  {			/* we need to go further		*/
    do {
      y = x * a->linmin_g1;
      t = macprodII (p, gy, y, dfunc, arg, a); 
      if (a->verbose > 1) 
	printf ("s = %6.3g: t = %6.3g; x = %6.3g y = %6.3g\n", s, t, x, y), fflush(stdout);
      if (t >= 0.0) break;
      x = y; s = t; a->gunused = gx; gx = gy; gy = a->gunused; 
      ++its;
      /* replaces: for (i = 0; i < n; ++i) gx[i] = gy[i];		*/
      } while (its <= a->linmin_maxits);
    }
  else
  if (s > 0) {			/* need to step back inside interval	*/
    do {
      y = x * a->linmin_g3;
      t = macprodII (p, gy, y, dfunc, arg, a); 
      if (a->verbose > 1) 
	printf ("s = %6.3g: t = %6.3g; x = %6.3g y = %6.3g\n", s, t, x, y), fflush(stdout);
      if (t <= 0.0) break;
      x = y; s = t; a->gunused = gx; gx = gy; gy = a->gunused; 
      ++its;
      } while (its <= a->linmin_maxits);
    }
  else {			/* hole in one s = 0.0			*/
    t = 1.0; y = x;
    }

  if (its > a->linmin_maxits) {
    fprintf (stderr, "\n\nWarning! maclinmin overran the maximum number of linmin iterations: %d\n", a->linmin_maxits);
    /* this can happen where the function goes \_ and doesn't buck up again;
     * it also happens if the initial `gradient' does not satisfy gradient.
     * `gradient' > 0, so that there is no minimum
     * in the supposed downhill direction.
     * I don't know if this actually happens...
     * If it does then I guess a->rich should be 1.
     *
     * If the overrun is because too big a step was taken
     * then the interpolation should be made
     * between zero and the most recent measurement. 
     * 
     * If the overrun is because too small a step was taken
     * then the best place to go is the most distant point. 
     * I will assume that this doesn't happen for the moment.
     * 
     * Also need to check up what happens to t and s in the case of overrun.
     *  And gx and gy. 
     * 
     * Maybe sort this out when writing a macopt
     * that makes use of the gradient at zero? 
     */
    fprintf (stderr, "- inner product at zero = %9.4g\n",
	tmpd = macprodII (p, gy, 0.0, dfunc, arg, a)); 
    if (tmpd > 0 && a->rich == 0) {
      fprintf (stderr, "setting rich to 1\n");       a->rich = 1; 
      }
    if (tmpd > 0) a->restart = 1; 
    }

  /*  Linear interpolate between the last two. 
      This assumes that x and y do bracket.				*/
  if (s < 0.0) s = - s;
  if (t < 0.0) t = - t;
  m = (s + t);
  s /= m; t /= m;
  
  m =  s * y + t * x; 
  /* evaluate the step length, not that it necessarily means anything	*/
  for (step = 0.0, i = 0; i < n; ++i) {
    tmpd = m * a->xi[i];
    p[i] += tmpd; /* this is the point where the parameter vector steps	*/
    step += fabs (tmpd); 
    a->xi[i] = s * gy[i] + t * gx[i];
    /* send back the estimated gradient in xi (NB not like linmin)	*/
    }
  a->lastx = m * a->linmin_g2 *  a->gtyp;
  
  return (step / (double) (n)); 
  }

double macprodII(
    nml_dvector * p,
    nml_dvector * gy,
    double y, 
    void (*dfunc)(nml_dvector *, nml_dvector *, void *, int,  double *), 
    void *arg, 
    macopt_args *a
   ) {
  nml_dvector * pt = a->pt; 
  nml_dvector * xi = a->xi; 
  /* finds pt = p + y xi and gets gy there, returning gy . xi		*/
  int n = a->n; 

  int i;
  double s = 0.0;

  for (i = 0; i < n; ++i) {
    pt[i] = p[i] + y * xi[i];
    }
  dfunc(pt, gy, arg, a->its, &a->track_one_value);

  for (i = 0; i < n; ++i) {
    s += gy[i] * xi[i];
    }
  return s;
  }

void macopt_defaults(
    macopt_args*	a,
    int		itermax,
    int		verbose,
    double	tol,
    int		end_if_small_step
    ) {
  a->verbose = verbose;
  /*	If verbose = 1, then there is one report for each line minimization.
	If verbose = 2, then there is an additional
	  report for each step of the line minimization.
	If verbose = 3, then extra debugging routines kick in.		*/

  a->tol = tol;			/* Do fiddle with this			*/
  a->end_if_small_step = end_if_small_step;
				/* Change this to 0/1 if you prefer	*/
  a->itmax = itermax;		/* You may wish to change this.		*/
  a->rich = 1;	/* If this is 1, then the program runs a bit slower.	*/
  a->stepmax = 0.5; 

  a->grad_tol_tiny = 1e-28;	/* Probably not worth fiddling with	*/
  a->step_tol_tiny = 0.0;	/* Probably not worth fiddling with	*/
  a->linmin_maxits = 20;	/* Probably not worth fiddling with	*/
  a->lastx = 0.01;		/* only has a transient effect		*/
  a->lastx_default = 0.01;	/* -- This defines typical distance
    in parameter space at which the line minimum is expected;
    Both these should be set.  The default is consulted
    if something goes badly wrong and a reset is demanded.		*/

  /* Don't fiddle with the following, unless you really mean it.	*/
  a->linmin_g1 = 2.0; 
  a->linmin_g2 = 1.25; 
  a->linmin_g3 = 0.5; 
  a->restart = 0; 
  a->filename_its = NULL;
  a->filename_last = NULL;
  }

void macopt_free (macopt_args *a) {
  if (a->filename_last && (mpi_n3d_masterid == mpi_n3d_id)) {
    FILE*	fp = NULL;
    fp = fopen(a->filename_last, "a");
    fprintf(fp, "%d %g\n", a->its, a->track_one_value);
    fclose(fp);
    }

  nml_dv_delete(a->xi);
  nml_dv_delete(a->h);
  nml_dv_delete(a->g);
  nml_dv_delete(a->pt);
  nml_dv_delete(a->gx);
  nml_dv_delete(a->gy); 
  nml_dv_delete(a->x_restart); 
  }

void macopt_restart (macopt_args *a, int start) {
  /* if start == 1, then this is the start of a fresh macopt, not a restart */
  int j, n = a->n; 
  nml_dvector	*g;
  nml_dvector	*h;
  nml_dvector	*xi;
#ifdef	RESTART_GEKCO
  nml_dvector* x_restart = a->x_restart;
#endif/*RESTART_GEKCO	*/

  g = a->g;  h = a->h;  xi = a->xi; 

  if (0 == start)
    a->lastx = a->lastx_default; 
  /* It is assumed that (*dfunc)(p, xi, dfunc_arg); has happened.	*/
#ifdef RESTART_GEKCO
#define RESTART2sign -	/* This sign is - in the original code.....	*/
  if (a->x_restart) {
    if (2 != a->restart) {
      for (j = 0; j < n; ++j) {
	g[j] = RESTART2sign x_restart[j];
	xi[j] = h[j] = g[j];
	}
      }
    else {
      for (j = 0; j < n; ++j) {
	xi[j] = h[j] = g[j] = x_restart[j];
	}
      }
    }
  else {
    if (2 != a->restart) {
      for (j = 0; j < n; ++j) {
	g[j] = RESTART2sign xi[j];
	xi[j] = h[j] = g[j];
	}
      }
    else {
      for (j = 0; j < n; ++j) {
	xi[j] = h[j] = g[j] = xi[j];
	}
      }
    }
#else /*RESTART_GEKCO	*/
  for (j = 0; j < n; ++j) {
    if (2 != a->restart)
      g[j] = -xi[j];
    xi[j] = h[j] = g[j];
    }
#endif/*RESTART_GEKCO	*/

  a->restart = 0; 
  }

#ifdef	RESTART_GEKCO
static macopt_args *macopt_global = NULL;

void set_macopt_args(macopt_args *a) {
  macopt_global = a;
  return;
  }

macopt_args *get_macopt_args(void) {
  return macopt_global;
  }

#endif/*RESTART_GEKCO	*/

