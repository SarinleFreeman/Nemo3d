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
$Header: /repo/eigen/include/macopt.h,v 1.3 2003/10/08 16:15:44 hook Exp $
*****************************************************************************/

#ifndef	MACOPT_H
#define	MACOPT_H 1

/* http://131.111.48.24/mackay/c/macopt.html       mackay@mrao.cam.ac.uk

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
   making a minimum of 3 per line search. Set rich=0 if you think
   that the surface is locally quite quadratic.) If the program does average
   2 derivatives per line search then it must be superior to most cg methods
   including use of Rbackprop (which costs 2 derivatives straight off)

   A possible modification: where the function can be returned at the same
   time as the dfunction --- there is nothing clever to do with the value
   but it could be used as a sanity check and a convergence criterion.

   See http://131.111.48.24/mackay/c/macopt.html for further discussion.

   NB: The value of "tol" is totally arbitrary
   and must be set by you to a value that works well for your problem.
   It depends completely on the typical value of the gradient/step size.

   Tol specifies a magnitude of gradient at which a halt is called.
   or a step size.

   This program MINIMIZES a function.
*/

#include "mpi_local.h"
#include "nml_dvector.h"

#define RESTART_GEKCO

/* structure for macopt */
typedef struct {
  double	tol;		/* Convergence is declared
    when the gradient vector is smaller in magnitude than this
    or when the mean absolute step is less than this (see above).	*/
  double	grad_tol_tiny;	/* If the gradient is less than this,
    we definitely stop even if we are not using a gradient tolerance.	*/
  double	step_tol_tiny;	/* If the step is less than this,
    we stop even if we are not using a step tolerance			*/
  int	end_if_small_step;	/* This defines the role of tol --
    the alternative is end_on_small_grad.				*/
  int	its;			/* number of its			*/
  int	itmax;			/* max					*/
  int	rich;			/* This determines whether to do
    the extra gradient evaluation at the beginning of each new line min.*/
  int	verbose;
  double stepmax;	/* largest step permitted (not used in macopt)	*/

  int	linmin_maxits;	/* in maclinmin					*/
  double linmin_g1;     /* factors for growing and shrinking the interval */
  double linmin_g2;
  double linmin_g3;
  double lastx;		/* keeps track of typical step length		*/
  double lastx_default;	/* If maclinmin is reset, lastx is set to this.	*/

  /* These should not be touched by the user.
     They are handy pointers for macopt to use.				*/
  double gtyp;		/* stores the rms gradient for linmin		*/
  nml_dvector	*pt;
  nml_dvector	*gx;
  nml_dvector	*gy;
  nml_dvector	*gunused;
  nml_dvector	*xi;
  nml_dvector	*g;
  nml_dvector	*h;
#ifdef RESTART_GEKCO
  nml_dvector	*x_restart;
#endif
  int n;       /* dimension of parameter space */
  int restart; /* whether to restart macopt - fresh cg directions */

  double track_one_value;	/* One double value
    to be tracked and returned in each gradient evaluation,		*/
  char *filename_its;
  char *filename_last;
  } macopt_args;

/* lastx :--- 1.0 might make general sense, (cf N.R.)
   but the best setting of all is to have a prior idea of the eigenvalues.
   If the objective function is equal to sum of N terms then set this to 1/N,
   for example Err on the small side to be conservative.		*/
void macoptII(
   nml_dvector*,	/* starting vector				*/
   int,			/* number of dimensions				*/
   void (*dfunc)(nml_dvector*, nml_dvector*, void *, int, double *),
			/* evaluates the gradient			*/
   void*,
   macopt_args *
   );

double maclinminII(
   nml_dvector*,
   void (*dfunc)(nml_dvector*, nml_dvector*, void *, int, double *),
			/* evaluates the gradient			*/
   void*,
   macopt_args*
   );

double macprodII(
   nml_dvector*,
   nml_dvector*,
   double,
   void (*dfunc)(nml_dvector*, nml_dvector*, void *, int, double *),
   void*,
   macopt_args*
   );

void macopt_defaults(
   macopt_args*,
   int itermax,
   int verbose,
   double tol,
   int end_if_small_step
   );

void macopt_free(macopt_args*);

/* The following functions could be declared static within macopt.c	*/

void macopt_restart(macopt_args*, int);

#ifdef RESTART_GEKCO
void set_macopt_args(macopt_args *a);
macopt_args *get_macopt_args(void);
#endif /* RESTART_GEKCO */

#endif /* MACOPT_H */

