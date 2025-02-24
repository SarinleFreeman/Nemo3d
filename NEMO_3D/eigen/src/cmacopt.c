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
$Header: /repo/eigen/src/cmacopt.c,v 1.2 2003/10/08 16:15:51 hook Exp $
*****************************************************************************/

#include "cmacopt.h"

/*
   http://131.111.48.24/mackay/c/macopt.html       mackay@mrao.cam.ac.uk

   Please do not use macopt without understanding a little about how it works;
   there are some control parameters which the user MUST set!

   David MacKay's optimizer, based on conjugate gradient ideas,
   but using bracketing of the zero of the inner product

             (gradient).(line_search_direction)

   to do the line minimization. Only derivative calculations are required.
   The length of the first step in the line search (often set to "1.0"
   in other code) is adapted here so that, if 0.00001 is a better step size,
   it soon cottons on to that and saves ~log(10000) bracketing operations.
   The result is that (with rich set to 0) the program can use
   as few as 2 derivatives per line search. (If rich is set to 1, it does
   an extra derivative calculation at the beginning of each line search
   making a minimum of 3 per line search. Set rich = 0 if you think
   that the surface is locally quite quadratic.) If the program does average
   2 derivatives per line search then it must be superior to most cg methods
   including use of Rbackprop (which costs 2 derivatives straight off)

   A possible modification: where the function can be returned at same
   time as the dfunction --- there is nothing clever to do with the
   value, but it could be used as a sanity check and a convergence criterion.

   See http://131.111.48.24/mackay/c/macopt.html for further discussion.

   NB: The value of "tol" is totally arbitrary and must be set by
   you to a value that works well for your problem.
   It depends completely on the typical value of the gradient/step size.

   Tol specifies a magnitude of gradient at which a halt is called.
   or a step size.

   This program MINIMIZES a function.
*/

void cmacoptII(
  nml_dcvector*	p,	/* starting vector				*/
  int		n,	/* number of dimensions				*/
  void		(*dfunc)(nml_dcvector*, nml_dcvector*, void*),
			/* evaluates the gradient of the optimized function */
  void*	dfunc_arg,	/* arguments that get passed to dfunc		*/
  cmacopt_args*	a	/* structure in which optimizer arguments stored*/
  ) {			/* Note: (*func)(cvectr, void *) is not used	*/
  int j;
  complex gg = czero, gam = czero, dgg = czero;
  cvectr g, h, xi;
  int end_if_small_grad = 1 - a->end_if_small_step;
  complex step = czero, tmpd = czero;

  /* A total of 7 cvectrs 0..n-1 are used by this optimizer.
     p		is provided when the optimizer is called
     pt		is used by the line minimizer as the temporary vector.
		this could be cut out with minor rewriting, using p alone
		g, h and xi are used by the cg method as in NR -
	       	could one of these be cut out?
     the line minimizer uses an extra gx and gy to evaluate two gradients.
  */

  a->n	= n;
  a->g  = Cvectr (n); g  = a->g;	/* vectors as in NR code	*/
  a->h  = Cvectr (n); h  = a->h;	/*				*/
  a->xi = Cvectr (n); xi = a->xi;	/*				*/
  a->pt = Cvectr (n);	/* scratch vector for sole use of macprod	*/
  a->gx = Cvectr (n);	/* scratch gradients				*/
  a->gy = Cvectr (n);	/* used by maclinmin and macprod		*/

#ifdef RESTART_GEKCO
  a->x_restart = NULL;
  set_cmacopt_args(a);
#endif	

  printf("First Gradient Evaluation before the start of the algorithm.\n");
  (*dfunc)(p, xi, dfunc_arg);
  cmacopt_restart (a, 1);

  for (a->its = 1; a->its <= a->itmax; a->its ++) {

    gg = czero;
    for (j = 0; j < n; j++) {
      /* gg = cadd(gg, cmul(cconj(g[j]), g[j]));			*/
      /* find the magnitude of the old gradient				*/
      gg.r += (g[j].r*g[j].r + g[j].i*g[j].i);
      }

    /* a->gtyp = cxsqrt (cxdiv(gg, Complex(n, 0)));			*/
    a->gtyp.r = sqrt(gg.r/n);
    a->gtyp.i = 0.0;

    if (a->verbose > 0)
      if (a->end_if_small_step)
	printf("cmac_it %d of %d : converge on (step = (%6.3g, %6.3g) with tol = %6.3g), monitor: (gradient = %6.3g, grad_tiny = %6.3g:) ",
	a->its, a->itmax, step.r, step.i, a->tol,  c_abs(gg), a->grad_tol_tiny);
      else
	/* print out gradient oriented convergence			*/
	printf("cmac_it %d of %d : converge on (gradient = %6.3g tol = %6.3g grad_tiny = %6.3g), monitor (step = (%6.3g, %6.3g) step_tiny = %6.3g) ",
	a->its, a->itmax, c_abs(gg), a->tol, a->grad_tol_tiny, step.r, step.i, a->step_tol_tiny);

    if ((end_if_small_grad && (c_abs(gg) <= a->tol))
	|| (c_abs(gg) <= a->grad_tol_tiny)) {
      printf("\nAchieved convergence though gradient!!\n");
      printf("cmac_it %d of %d : converge on (gradient = %6.3g tol = %6.3g grad_tiny = %6.3g), monitor (step = (%6.3g, %6.3g)  step_tiny = %6.3g)\n ",
	  a->its, a->itmax, c_abs(gg), a->tol, a->grad_tol_tiny, step.r, step.i, a->step_tol_tiny);
      cmacopt_free (a);
      if (a->verbose > 0) printf("\n");
      return;
      }

    step = cmaclinminII (p, dfunc, dfunc_arg, a);

    if (a->restart == 0) {
      if (a->verbose > 1) printf(" (step (%9.5g, %9.5g))", step.r, step.i);
      if (a->verbose > 0) printf("\n");
      if ((a->end_if_small_step  && (c_abs(step) <= a->tol))
	  || (c_abs(step) <= a->step_tol_tiny)) {
	cmacopt_free (a);
	return;
	}
      }

    /* if we are feeling rich, evaluate the gradient at the new `minimum'.
     * alternatively, linmin has already estimated this gradient
     * by linear combination of the last two evaluations and left it in xi */
    if (a->rich || a->restart) {
      (*dfunc)(p, xi, dfunc_arg);
      }
    if (a->restart) {
      fprintf(stderr, "Restarting cmacopt within the main iteration loop\n");
      cmacopt_restart (a, 0);
      /* This is not quite right.  It should distinguish
       * whether there was an overrun indicating that the value of lastx
       * needs to be bigger/smaller;
       * in which case resetting lastx to default value may be a bad idea,
       * giving an endless loop of resets				*/
      }
    else {
      dgg = czero;
      for (j = 0; j < n; ++j) {
	dgg = cadd(dgg, cmul(cadd(xi[j], g[j]), xi[j]));
	}
      gam = cxdiv(dgg, gg);
      tmpd.r = tmpd.i = 0.0;
      for (j = 0; j < n; ++j) {
	g[j] = RCmul(-1.0, xi[j]);
			/* g stores (-) the most recent gradient	*/
	h[j] = cadd(g[j], cmul(gam, h[j]));
	xi[j] = h[j];
	/* h stores xi, the current line direction */
	/* check that the inner product of gradient and line search is < 0 */
	tmpd = csub(tmpd, cmul(xi[j], g[j]));
	}
      if (tmpd.r > 0.0  || a->verbose > 2) {
	fprintf(stderr, "cmacopI error: new line search has positive inner prod (%9.4g, %9.4g)\n", tmpd.r, tmpd.i);
	}
      if (tmpd.r > 0.0) {
	if (a->rich == 0) {
	  fprintf(stderr, "cmacopII: Setting rich to 1; ");
	  a->rich = 1;
	  }
	a->restart = 2; /* signifies that g[j] = -xi[j] is already done */
	fprintf(stderr, "\n\ncmacopII: Restarting cmacopt due to positive inner product within the line search\n\n");
	cmacopt_restart (a, 0);
	}
      }
    }
  fprintf(stderr, "Reached iteration limit in cmacopt; continuing.\n");
  cmacopt_free (a);	
  return;
  } /* NB this leaves the best value of p in the p vector,
       but the function has not been evaluated there if rich = 0	*/

/* Examines objective function and d_objective function
 * to see if they agree for a step of size epsilon			*/
void cmaccheckgrad(
  nml_dcvector*	p,
  int		n,
  nml_dcscalar	epsilon,
  nml_dcscalar	(*func)(nml_dcvector*, void *),
  void*		func_arg,
  void		(*dfunc)(nml_dcvector*, nml_dcvector*, void *),
  void*		dfunc_arg,
  int		stopat	/* stop at this component. If 0, do the lot.	*/
   ) {

  int j;
  complex f1;
  cvectr g, h;
  complex tmpp;

  h = Cvectr(n);
  g = Cvectr(n);
  f1 = (*func)(p, func_arg);
  (*dfunc)(p, g, dfunc_arg);
  if (stopat <= 0 || stopat > n) stopat = n;

  printf("Testing gradient evaluation\n");
  printf("      analytic     1st_diffs    difference\n");
  for (j = 1; j <= stopat; ++j) {
    tmpp = p[j];
    p[j] = cadd(p[j], epsilon);
    h[j] = csub((*func)(p, func_arg), f1);
    p[j] = tmpp;

    /*
    printf("%2d %9.5g %9.5g %9.5g\n", j, c_abs(g[j]),
	c_abs(cxdiv(h[j], epsilon)), c_abs(csub(g[j], cxdiv(h[j], epsilon))));
    */
    fflush(stdout);
    }
  rm_cvectr(&h);
  rm_cvectr(&g);
  printf("      --------     ---------\n");
  }

void cmacopt_defaults (cmacopt_args *a, int itermax, int verbose,
    double tol, int end_if_small_step) {

  /* if verbose = 1 then there is one report for each line minimization.
     if verbose = 2 then there is an additional report
		    for each step of the line minimization.
     if verbose = 3 then extra debugging routines kick in.		*/
  a->verbose = verbose;

  a->tol = tol;		/* Do fiddle with this				*/
  a->end_if_small_step = end_if_small_step;
			/* Change this to 0/1 if you prefer		*/
  a->itmax = itermax;	/* You may wish to change this			*/
  a->rich = 1;		/* if this is 1, then the program runs a bit slower */
  a->stepmax = 0.5;

  a->grad_tol_tiny = 1e-16;	/* Probably not worth fiddling with	*/
  a->step_tol_tiny = 0.0;	/* Probably not worth fiddling with	*/
  a->linmin_maxits = 20;	/* Probably not worth fiddling with	*/
  a->lastx = Complex(0.01, 0.0);	/* only has a transient effect	*/
  a->lastx = Complex(0.0001, 0.0);	/* only has a transient effect	*/
  a->lastx_default = Complex(0.01, 0.0);
  a->lastx_default = Complex(0.0001, 0.0);
				/* -- defines typical distance in parameter
				space at which the line minimum is expected;
				both these should be set. the default is
				consulted if something goes badly wrong and
				a reset is demanded. */

  /* don't fiddle with the following, unless you really mean it		*/
  a->linmin_g1 = Complex(2.0, 0.0);
  a->linmin_g2 = Complex(1.25, 0.0);
  a->linmin_g3 = Complex(0.5, 0.0);
  a->restart = 0;
  }

void cmacopt_free (cmacopt_args *a) {
  rm_cvectr(&a->xi);
  rm_cvectr(&a->h);
  rm_cvectr(&a->g);
  rm_cvectr(&a->pt);
  rm_cvectr(&a->gx);
  rm_cvectr(&a->gy);
#ifdef RESTART_GEKCO
  rm_cvectr(&a->x_restart);
#endif	
  }

/* maclinmin.
   Method:
       evaluate gradient at a sequence of points and calculate the inner
       product with the line search direction. Continue until a
       bracketing is achieved (i.e a change in sign). */
complex cmaclinminII(
  cvectr	p,
  void (*dfunc)(cvectr, cvectr, void *)	/* evaluates the gradient	*/,
  void*	arg,
  cmacopt_args*	a
  ) {
  int n = a->n;

  complex x, y, cret;
  complex s, t, m;
  int its = 1, i;
  complex step, tmpd;
  cvectr  gx = a->gx;
  cvectr  gy = a->gy;

  /* at x = 0, the gradient (uphill) satisfies s < 0 */
  if (a->verbose > 2) { /* check this is true:				*/
    /* (no need to do this really as it is already checked
     * at the end of the main loop of macopt)				*/
    /*
#define TESTS 5
    x = a->lastx/a->gtyp;
    fprintf(stderr, "inner product at:\n");
    for (i = -TESTS; i <= TESTS; i += 2) {
      step = x*2.0*(real) i/(real) TESTS;
      fprintf(stderr, "%9.5g %9.5g\n", step,
	  tmpd = cmacprodII (p, gy, step, dfunc, arg, a));
      }
    */
    tmpd = cmacprodII (p, gy, Complex(0.0, 0.0), dfunc, arg, a);
    fprintf(stderr, "inner product at zero = %9.4g %9.4g\n", tmpd.r, tmpd.i);
    if (tmpd.r > 0.0) {
      a->restart = 1;
      return Complex(0.0, 0.0);
      }
    }
  x = cxdiv(a->lastx, a->gtyp);
  s = cmacprodII (p, gx, x, dfunc, arg, a);

#define OLD_BRACKETING
#ifdef OLD_BRACKETING
  if (s.r < 0) {	/* we need to go further			*/
    do {
      y = cmul(x, a->linmin_g1);
      t = cmacprodII (p, gy, y, dfunc, arg, a);
      if (a->verbose > 1)
	printf(" cminII s = (%1.2g %1.2g); t = %6.3g; x = %6.3g y = %6.3g\n",
	    s.r, s.i, t.r, x.r, y.r);
      if (t.r >= 0.0) break;
      x = y; s = t; a->gunused = gx; gx = gy; gy = a->gunused;
      ++its;
      /* replaces: for (i = 0; i < n; ++i) gx[i] = gy[i];		*/
      } while (its <= a->linmin_maxits);
    }
  else
  if (s.r > 0) {	/* need to step back inside interval		*/
    do {
      y = cmul(x, a->linmin_g3);
      t = cmacprodII (p, gy, y, dfunc, arg, a);
      if (a->verbose > 1)
	printf(" cmacII s = (%1.2g %1.2g) t = %6.3g; x = %6.3g y = %6.3g\n",
	    s.r, s.i, t.r, x.r, y.r);
      if (t.r <= 0.0) break;
      x = y; s = t; a->gunused = gx; gx = gy; gy = a->gunused;
      ++its;
      } while (its <= a->linmin_maxits);
    }
  else {		/* hole in one s = 0.0				*/
    t = Complex(1.0, 0.0); y = x;
    }
#else
#undef RUSH2112
#ifdef RUSH2112
  if (s.r < 0)  {	/* we need to go further			*/
    do {
      y.r = x.r*a->linmin_g1.r;
      t = cmacprodII (p, gy, y, dfunc, arg, a);
      if (a->verbose > 1)
	printf("s = %6.3g: t = %6.3g; x = %6.3g y = %6.3g\n", s, t, x, y);
      if (t.r >= 0.0) break;
      x = y; s = t; a->gunused = gx; gx = gy; gy = a->gunused;
      ++its;
      /* replaces: for (i = 0; i < n; ++i) gx[i] = gy[i]; */
      } while (its <= a->linmin_maxits);
    }
  else
  if (s.r > 0) {	/* need to step back inside interval		*/
    do {
      y.r = x.r*a->linmin_g3.r;
      t = cmacprodII (p, gy, y, dfunc, arg, a);
      if (a->verbose > 1)
	printf("s = %6.3g: t = %6.3g; x = %6.3g y = %6.3g\n", s, t, x, y);
      if (t.r <= 0.0) break;
      x = y; s = t; a->gunused = gx; gx = gy; gy = a->gunused;
      ++its;
      } while (its <= a->linmin_maxits);
    }
  else {			/* hole in one s = 0.0			*/
    t = Complex(1.0, 0.0); y = x;
    }
#endif
  if (s.r < 0 && s.i < 0)  {	/* we need to go further */
    t.r = t.i = -1;
    do {
      y.r = x.r*a->linmin_g1.r;
      y.i = x.i*a->linmin_g1.r;
      t = cmacprodII (p, gy, y, dfunc, arg, a);
      if (a->verbose > 1)
	printf(" minII s = %1.1e %1.1e: t = %1.1e %1.1e; x = %1.1e %1.1e y = %1.1e %1.1e\n", s.r, s.i, t.r, t.i, x.r, x.i, y.r, y.i);
      if (t.r >= 0.0 || t.i >= 0.0) break;
      x = y; s = t; a->gunused = gx; gx = gy; gy = a->gunused;
      ++its;
      } while (its <= a->linmin_maxits);
    }
  else
  if (s.r > 0 && s.i > 0) {	/* need to step back inside interval	*/
    t.r = t.i = 1;
    do {
      y.r = x.r*a->linmin_g3.r;
      y.i = x.i*a->linmin_g3.r;
      t = cmacprodII (p, gy, y, dfunc, arg, a);
      if (a->verbose > 1)
	printf("s = %1.1e %1.1e: t = %1.1e %1.1e; x = %1.1e %1.1e y = %1.1e %1.1e\n", s.r, s.i, t.r, t.i, x.r, x.i, y.r, y.i);
      if (t.r <= 0.0 || t.i <= 0.0) break;
      x = y; s = t; a->gunused = gx; gx = gy; gy = a->gunused;
      ++its;
      } while (its <= a->linmin_maxits);
    }
  else
    if (s.r < 0 && s.i > 0) {	/* need to step back inside interval	*/
      t.r = -1;
      t.i = 1;
      do {
	y.r = x.r*a->linmin_g1.r;
	y.i = x.i*a->linmin_g3.r;
	t = cmacprodII (p, gy, y, dfunc, arg, a);
	if (a->verbose > 1)
	  printf("s = %1.1e %1.1e: t = %1.1e %1.1e; x = %1.1e %1.1e y = %1.1e %1.1e\n", s.r, s.i, t.r, t.i, x.r, x.i, y.r, y.i);
	if (t.r >= 0.0 || t.i <= 0.0) break;
	x = y; s = t; a->gunused = gx; gx = gy; gy = a->gunused;
	++its;
	} while (its <= a->linmin_maxits);
      }
    else
    if (s.r > 0 && s.i < 0) {	/* need to step back inside interval	*/
      t.r = 1;
      t.i = -1;
      do {
	y.r = x.r*a->linmin_g3.r;
	y.i = x.i*a->linmin_g1.r;
	t = cmacprodII (p, gy, y, dfunc, arg, a);
	if (a->verbose > 1)
	printf("s = %1.1e %1.1e: t = %1.1e %1.1e; x = %1.1e %1.1e y = %1.1e %1.1e\n", s.r, s.i, t.r, t.i, x.r, x.i, y.r, y.i);
	if (t.r <= 0.0 || t.i >= 0.0) break;
	x = y; s = t; a->gunused = gx; gx = gy; gy = a->gunused;
	++its;
	} while (its <= a->linmin_maxits);
      }
    else {			/* hole in one s = 0.0			*/
      t = Complex(1.0, 0.0); y = x;
      }
#endif

  if (its > a->linmin_maxits) {
    fprintf(stderr, "\n\nWarning! cmaclinmin overran the maximum number of linmin iterations: %d\n", a->linmin_maxits);
    /* this can happen where the function goes \_ and doesn't buck up again;
     * it also happens if the initial `gradient' does not satisfy gradient.
     * `gradient' > 0, so that there is no minimum
     * in the supposed downhill direction.
     * I don't know if this actually happens...
     * If it does then I guess a->rich should be 1.
     *
     * If the overrun is because too big a step was taken
     * then the interpolation should be made between zero
     * and the most recent measurement.
     *
     * If the overrun is because too small a step was taken,
     * then the best place to go is the most distant point.
     * I will assume that this doesn't happen for the moment.
     *
     * Also need to check up what happens to t and s in the case of overrun.
     * And gx and gy.
     *
     * Maybe sort this out when writing a macopt
     * that makes use of the gradient at zero?
     */
    tmpd = cmacprodII (p, gy, Complex(0.0, 0.0), dfunc, arg, a);
    fprintf(stderr, "- inner product at zero = (%9.4g, %9.4g)\n",
	tmpd.r, tmpd.i);

    if (tmpd.r > 0 && a->rich == 0) {
      fprintf(stderr, "setting rich to 1\n");       a->rich = 1;
      }
    if (tmpd.r > 0) a->restart = 1;
    }

  /*  Linear interpolate between the last two.
      This assumes that x and y do bracket.				*/
  if (s.r < 0.0) s.r = -s.r;
  if (s.i < 0.0) s.i = -s.i;
  if (t.r < 0.0) t.r = -t.r;
  if (t.i < 0.0) t.i = -t.i;
  m = cadd(s, t);
  s = cxdiv(s, m);
  t = cxdiv(t, m);

  m = cadd(cmul(s, y), cmul(t, x));
  /* evaluate the step length, not that it necessarily means anything	*/
  step = Complex(0.0, 0.0);
  for (i = 0; i < n; ++i) {
    tmpd = cmul(m, a->xi[i]);
    p[i] = cadd(p[i], tmpd); /* this is the point
				where the parameter vector steps	*/
    step = cadd(step, Complex(c_abs(tmpd), 0.0));
    a->xi[i] = cadd(cmul(s, gy[i]), cmul(t, gx[i]));
    /* send back the estimated gradient in xi (NB not like linmin)	*/
    }
  a->lastx = cmul(m, cmul(a->linmin_g2, a->gtyp));

  cret = cxdiv(step, Complex(n, 0));
  return (cret);
  }

complex cmacprodII(
  cvectr	p,
  cvectr	gy,
  complex	y,
  void		(*dfunc)(cvectr, cvectr, void *),
  void*		arg,
  cmacopt_args*	a
  ) {
  cvectr pt = a->pt;
  cvectr xi = a->xi;
  /* finds pt = p + y xi and gets gy there, returning gy . xi		*/
  int n = a->n;

  int i;
  complex s = Complex (0.0, 0.0);

  for (i = 0; i < n; i++) {
    pt[i] = cadd(p[i], cmul(y, xi[i]));
    }

  dfunc(pt, gy, arg);

  for (i = 0; i < n; i++) {
    s = cadd(s, cmul(gy[i], xi[i]));
    }

  return s;
  }

/* if start == 1 then this is the start of a fresh macopt, not a restart */
void cmacopt_restart (cmacopt_args *a, int start) {
  int j, n = a->n;
  cvectr g, h, xi;
#ifdef RESTART_GEKCO
  cvectr x_restart;
  x_restart = a->x_restart;
#endif /* RESTART_GEKCO */
  g = a->g;  h = a->h;  xi = a->xi;

  if (start == 0) a->lastx = a->lastx_default;
  /* it is assumed that (*dfunc)(p, xi, dfunc_arg);  has happened	*/

#ifdef RESTART_GEKCO
#define RESTART2sign +	/* This sign is - in the original code.		*/
  if (a->x_restart){
    if (a->restart != 2) {
      for (j = 0; j < n; ++j) {
	g[j].r = RESTART2sign x_restart[j].r;
	g[j].i = RESTART2sign x_restart[j].i;
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
    if (a->restart != 2) {
      for (j = 0; j < n; ++j) {
	g[j].r = RESTART2sign xi[j].r;
	g[j].i = RESTART2sign xi[j].i;
	xi[j] = h[j] = g[j];
	}
      }
    else {
      for (j = 0; j < n; ++j) {
	xi[j] = h[j] = g[j] = xi[j];
	}
      }
    }
#else /* RESTART_GEKCO */	
  if (start == 0) a->lastx = a->lastx_default;
  /* it is assumed that (*dfunc)(p, xi, dfunc_arg);  has happened	*/
  for (j = 0; j < n; ++j) {
    if (a->restart != 2) g[j] = RCmul(-1.0, xi[j]);
    xi[j] = h[j] = g[j];
    }
#endif /* RESTART_GEKCO */

  a->restart = 0;
  }

#ifdef RESTART_GEKCO
static cmacopt_args *a_global = NULL;

void set_cmacopt_args(cmacopt_args *a) {
  a_global = a;
  return;
  }

cmacopt_args *get_cmacopt_args(void) {
  return a_global;
  }

#endif /* RESTART_GEKCO */
