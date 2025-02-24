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
$Source: /repo/eigen/src/macopt_par.c,v $
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

#include "macopt_par.h"
#include "MPI_Timing.h"

void macoptII_par(
    nml_dvector	*p,	/* starting vector				*/
    int		n_loc,	/* size of p (this proc)			*/
    int		n,	/* number of dimensions				*/
			/* evaluates the gradient of the optimized function */
    void	(*dfunc)(nml_dvector*, int, nml_dvector*, void*, int, double*),
    void         *dfunc_arg,	/* arguments that get passed to dfunc	*/
    macopt_args	*a	/* structure in which optimizer arguments stored*/
    ) {		/* Note, (*func)(nml_dvector*, void *) is not used.	*/
  int j, flag_done;
  double gg_loc, gg, gam, dgg_loc, dgg;
  nml_dvector *g;
  nml_dvector *h;
  nml_dvector *xi;
  int end_if_small_grad = 1 - a->end_if_small_step;
  double step = 0.0, tmpd, tmpd_loc;

  /* A total of 7 nml_dvector*s 0..n-1 are used by this optimizer.
     p	   is provided when the optimizer is called
     pt	   is used by the line minimizer as the temporary vector.
	   this could be cut out with minor rewriting, using p alone
	   g, h and xi are used by the cg method as in NR -
	   could one of these be cut out?
     the line minimizer uses an extra gx and gy to evaluate two gradients.
  */

  int		i_am_root = (mpi_n3d_id == mpi_n3d_masterid);

  nml_memory_report("macoptII_par  beginning");
  a->n  = n_loc;
  a->g  = nml_dv_new(n_loc); g  = a->g;	/* vectors as in NR code*/
  a->h  = nml_dv_new(n_loc); h  = a->h;	/*			*/
  a->xi = nml_dv_new(n_loc); xi = a->xi;	/*			*/
  nml_memory_report("macoptII_par  after three basic vector assignments");
  a->pt = nml_dv_new(n_loc); /* scratch vector for sole use of macprod	*/
  a->gx = nml_dv_new(n_loc); /* scratch gradients			*/
  a->gy = nml_dv_new(n_loc); /* used by maclinmin and macprod		*/

  nml_memory_report("macoptII_par  after initial memory creation");

#ifdef	RESTART_GEKCO
  a->x_restart = NULL;
  set_macopt_args_par(a);
#endif/*RESTART_GEKCO	*/
  a->its = 0;
  a->track_one_value = 0.0;

  if (i_am_root)
    printf("First Gradient Evaluation before the start of the algorithm.\n");

  (*dfunc)(p, n_loc, xi, dfunc_arg, a->its, &a->track_one_value);
  macopt_restart(a, 1);
  for (a->its = 1; a->its <= a->itmax; ++(a->its)) {

    for (gg_loc = 0.0, j = 0; j < n_loc; ++j)
      gg_loc += g[j]*g[j];	/* find the magnitude of the old gradient */

    Allreduce_MPI_sp(&gg_loc, &gg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    a->gtyp = sqrt(gg/(double)(n));

    flag_done = (end_if_small_grad && (gg <= a->tol)) || (gg <= a->grad_tol_tiny);

    if (i_am_root) {
       if (a->verbose > 0 && i_am_root) {
          if (a->end_if_small_step) {
             printf("\nmac_it = %d / %d : ", a->its, a->itmax);
             printf("step = %6.3e < %6.3e   grad = %6.3e < %6.3e", 
                    step, a->tol, gg, a->grad_tol_tiny);
	  }
          else {
             printf("\nmac_it = %d / %d : ", a->its, a->itmax);
             printf("grad = %6.3e < %6.3e   step = %6.3e < %6.3e", 
                    gg, a->grad_tol_tiny, step, a->tol);
	  }
       }
       if (flag_done) {
          if (gg <= a->grad_tol_tiny) {
             printf("\nAchieved convergence though gradient.");
          }
          else {
             printf("\nAchieved convergence but grad > tolerance_grad!!\n");
	  }

          printf("\nmac_it = %d / %d : ", a->its, a->itmax);
          printf("grad = %6.3e < %6.3e   step = %6.3e < %6.3e", 
                 gg, a->grad_tol_tiny, step, a->tol);

          if (a->verbose > 0) printf("\n");
       }
    }

    if (flag_done) {
      macopt_free(a);
      return;
    }

    step = maclinminII_par(p, dfunc, dfunc_arg, a, n);

    if (a->filename_its && mpi_n3d_id == mpi_n3d_masterid) {
      FILE *fp = NULL;
      fp = fopen(a->filename_its, "a");
      fprintf(fp, "%d %g %g %g\n", a->its, a->track_one_value, step, gg);
      fclose(fp);
      }

    flag_done = (a->end_if_small_step && step <= a->tol) || (step <= a->step_tol_tiny);

    if (a->restart == 0 && i_am_root) {
       if (a->verbose > 1) printf(" (step %9.5g)", step);
       if (a->verbose > 0) printf("\n");
       if (flag_done) {
          if (step <= a->step_tol_tiny) {
             printf("\nAchieved convergence though step size.");
	  }
          else {
             printf("\nAchieved convergence but step size > tolerance!!\n");
	  }

          printf("\nmac_it = %d / %d : ", a->its, a->itmax);
          printf("step = %6.3e < %6.3e   grad = %6.3e < %6.3e", 
                 step, a->tol, gg, a->grad_tol_tiny);
       }
    }

    if (flag_done && a->restart == 0) {
      macopt_free(a);
      return;
      }

    /* if we are feeling rich, evaluate the gradient at the new `minimum'.
       Alternatively, linmin has already estimated this gradient
       by linear combination of the last two evaluations and left it in xi.
    */
    if (a->rich || a->restart) {
      (*dfunc)(p, n_loc, xi, dfunc_arg, a->its, &a->track_one_value);
      }

    if (a->restart) {
      if (i_am_root)
	fprintf(stderr, "Restarting macopt within the main iteration loop\n");
      macopt_restart(a, 0);
      /* This is not quite right.  It should distinguish
       * whether there was an overrun indicating that
       * the value of lastx needs to be bigger/smaller;
       * in which case resetting lastx to default value may be a bad idea,
       * giving an endless loop of resets				*/
      }
    else {

      dgg_loc = 0.0;
      for (j = 0; j < n_loc; ++j) {
	dgg_loc += (xi[j] + g[j])*xi[j];
	}

      Allreduce_MPI_sp(&dgg_loc, &dgg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      gam = dgg/gg;

      /* GLOBINFO: tmpd	*/
      for (tmpd_loc = 0.0, j = 0; j < n_loc; ++j) {
	g[j] = -xi[j];	/* g stores (-) the most recent gradient	*/
	xi[j] = h[j] = g[j] + gam*h[j];
	/* h stores xi, the current line direction			*/
	/* check that the inner prod of gradient and line search is <0	*/
	tmpd_loc -= xi[j]*g[j];
	}

      Allreduce_MPI_sp(&tmpd_loc, &tmpd, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      if (i_am_root && (tmpd > 0.0 || a->verbose > 2)) {
	fprintf(stderr, "new line search has inner prod %9.4g, ", tmpd);
	fprintf(stderr, "this will result in a restart\n");
	}
      if (tmpd > 0.0) {
	if (a->rich == 0) {
	  if (i_am_root)  fprintf(stderr, "Setting rich to 1; ");
	  a->rich = 1;
	  }
	a->restart = 2; /* signifies that g[j] = -xi[j] is already done */
	if (i_am_root) {
	  fprintf(stderr,
	      "Restarting macopt (2) due to positive inner product ");
	  fprintf(stderr, "within the line search\n\\n");
	  }
	macopt_restart(a, 0);
	}
      }
    }
  if (i_am_root)
    fprintf(stderr, "Reached iteration limit in macopt; continuing.\n");
  macopt_free(a);	

  fflush(stdout);

  return;
  }	/* macoptII_par
	   NB this leaves the best value of p in the p vector
	   but the function has not been evaluated there if rich = 0	*/

/* maclinmin.
   Method:
	evaluate gradient at a sequence of points and calculate
	the inner product with the line search direction.
	Continue until a bracketing is achieved (i.e a change in sign).	*/
double maclinminII_par(
 nml_dvector *p,
 void (*dfunc)(nml_dvector*, int, nml_dvector*, void*, int, double*),
			/* evaluates the gradient */
 void *arg,
 macopt_args *a,
 int n_tot
 ) {
  int		n = a->n;
  double	x, y;
  double	s, t, m;
  int		its = 1, i;
  double	step_tot, step, tmpd;
  double	*gx = a->gx, *gy = a->gy;

  /* at x = 0, the gradient (uphill) satisfies s < 0			*/
  if (a->verbose > 2) {		/* check that this is true:
    (no need to do this really as it is already checked
    at the end of the main loop of macopt)				*/
    tmpd = macprodII_par(p, gy, 0.0, dfunc, arg, a);
    if (mpi_n3d_id == mpi_n3d_masterid)
      fprintf(stderr, "inner product at zero = %9.4g\n", tmpd);
    if (tmpd > 0.0) {
      a->restart = 1;
      return 0.0;
      }
    }

  x = a->lastx/a->gtyp;
  s = macprodII_par(p, gx, x, dfunc, arg, a);

  if (s < 0) {	/* we need to go further				*/
    do {
      y = x*a->linmin_g1;
      t = macprodII_par(p, gy, y, dfunc, arg, a);
      if (a->verbose > 1 && mpi_n3d_id == mpi_n3d_masterid)
	printf("     s = %6.3g: t = %6.3g; x = %6.3g y = %6.3g\n", s, t, x, y);
      if (t >= 0.0) break;
	x = y; s = t; a->gunused = gx; gx = gy; gy = a->gunused;
	++its;
	/* replaces: for (i = 0; i < n; ++i) gx[i] = gy[i];		*/
      } while (its <= a->linmin_maxits);
    }
  else
  if (s > 0) {	/* need to step back inside interval			*/
    do {
      y = x*a->linmin_g3;
      t = macprodII_par(p, gy, y, dfunc, arg, a);
      if (a->verbose > 1 && mpi_n3d_id == mpi_n3d_masterid)
	printf("     s = %6.3g: t = %6.3g; x = %6.3g y = %6.3g\n", s, t, x, y);
      if (t <= 0.0) break;
      x = y; s = t; a->gunused = gx; gx = gy; gy = a->gunused;
      ++its;
      } while (its <= a->linmin_maxits);
    }
  else {	/* hole in one s = 0.0					*/
    t = 1.0; y = x;
    }

  if (its > a->linmin_maxits) {
    if (mpi_n3d_id == mpi_n3d_masterid) {
      fprintf(stderr, "\n\nWarning! maclinmin overran the maximum number ");
      fprintf(stderr, "of linmin iterations: %d\n", a->linmin_maxits);
      }
    /* this can happen where the function goes \_ and doesn't buck up again;
       it also happens if the initial `gradient' does not satisfy gradient.
       `gradient' > 0, so that
       there is no minimum in the supposed downhill direction.
       I don't know if this actually happens...
       If it does then I guess a->rich should be 1.

       If the overrun is because too big a step was taken,
       then the interpolation should be made between zero
       and the most recent measurement.

       If the overrun is because too small a step was taken,
       then the best place to go is the most distant point.
       I will assume that this doesn't happen for the moment.

       Also need to check up what happens to t and s in the case of overrun.
       And gx and gy.

       Maybe sort this out when writing a macopt
       that makes use of the gradient at zero?
    */

    tmpd = macprodII_par(p, gy, 0.0, dfunc, arg, a);
    if (mpi_n3d_id == mpi_n3d_masterid)
      fprintf(stderr, "- inner product at zero = %9.4g\n", tmpd);
    if (tmpd > 0 && a->rich == 0) {
      if (mpi_n3d_id == mpi_n3d_masterid)
	fprintf(stderr, "setting rich to 1\n");
      a->rich = 1;
      }
    if (tmpd > 0) a->restart = 1;
    }

  /* Linearly interpolate between the last two.
     This assumes that x and y do bracket.				*/
  if (s < 0.0) s = - s;
  if (t < 0.0) t = - t;
  m = (s + t);
  s /= m; t /= m;

  /* GLOBINFO: step	*/
  m =  s*y + t*x;
  /* evaluate the step length, not that it necessarily means anything	*/
  for (step = 0.0, i = 0; i < n; ++i) {
    tmpd = m*a->xi[i];
    p[i] += tmpd; /* this is the point where the parameter vector steps */
    step += fabs(tmpd);
    a->xi[i] = s*gy[i] + t*gx[i];
    /* send back the estimated gradient in xi (NB not like linmin)	*/
    }
  a->lastx = m*a->linmin_g2*a->gtyp;

  Allreduce_MPI_sp(&step, &step_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return (step_tot/(double)(n_tot));
  }	/* maclinminII_par	*/

double macprodII_par(
    nml_dvector* p, nml_dvector* gy, double y,
    void (*dfunc)(nml_dvector*, int, nml_dvector*, void*, int, double *),
    void *arg,
    macopt_args *a
    ) {
  nml_dvector* pt = a->pt;
  nml_dvector* xi = a->xi;
  /* finds pt = p + y xi and gets gy there, returning gy . xi		*/
  int n = a->n;

  int i;
  double sum_local = 0.0, sum_tot;

  for (i = 0; i < n; ++i) {
    pt[i] = p[i] + y*xi[i];
  }
  dfunc(pt, n, gy, arg, a->its, &a->track_one_value);

  for (i = 0; i < n; ++i) {
    sum_local += gy[i]*xi[i];
  }
  Allreduce_MPI_sp(&sum_local, &sum_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return sum_tot;
  } /* macprodII_par */


#ifdef	RESTART_GEKCO

void set_macopt_args_par(macopt_args *a) {
  int i;

  if (mpi_n3d_id == mpi_n3d_masterid) {
    for (i = 0; i<MAX_NPROC_macopt; ++i) {
      macopt_global[i] = NULL;
      }
    }
  MPI_Barrier(MPI_COMM_WORLD);

  macopt_global[mpi_n3d_id] = a;
  }

macopt_args *get_macopt_args_par(void) {
  return macopt_global[mpi_n3d_id];
  }

#endif /*RESTART_GEKCO	*/

