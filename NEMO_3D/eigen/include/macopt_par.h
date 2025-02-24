#ifndef	MACOPT_PAR_H
#define	MACOPT_PAR_H 1

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
$Header: /repo/eigen/include/macopt_par.h,v 1.3 2003/10/08 16:15:45 hook Exp $
*****************************************************************************/


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


#include "macopt.h"
#include "mpi_local.h"
#include "mpi_fake.h"
#include "nml_dvector.h"

#define RESTART_GEKCO

void macoptII_par(
   nml_dvector*,	/* starting vector				*/
   int,			/* number of dimensions, this proc		*/
   int,			/* number of dimensions				*/
   void (*dfunc)(nml_dvector*, int, nml_dvector*, void*, int, double*),
			/* evaluates the gradient			*/
   void*,
   macopt_args*
   );

/* The following functions could be declared static within macopt.c	*/

double maclinminII_par(
   nml_dvector*,
   void (*dfunc)(nml_dvector*, int, nml_dvector*, void*, int, double*),
			/* evaluates the gradient			*/
   void*,
   macopt_args *,
   int
   );


double macprodII_par(
   nml_dvector*,
   nml_dvector*,
   double,
   void (*dfunc)(nml_dvector*, int, nml_dvector*, void*, int, double *),
   void*,
   macopt_args*
   );

void macopt_restart_par(macopt_args*, int);

#ifdef RESTART_GEKCO
#define MAX_NPROC_macopt 512  /* ASSUME we have less than 512 processors */
static macopt_args *macopt_global[MAX_NPROC_macopt];
void set_macopt_args_par(macopt_args *);
macopt_args *get_macopt_args_par(void);
#endif /* RESTART_GEKCO */


#endif /* MACOPT_PAR_H */

