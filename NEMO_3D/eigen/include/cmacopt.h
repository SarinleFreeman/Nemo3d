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
$Header: /repo/eigen/include/cmacopt.h,v 1.3 2003/10/08 16:15:41 hook Exp $
*****************************************************************************/

#ifndef	CMACOPT_H
#define	CMACOPT_H 1

/*
This file is part of
The Jet Propulsion Laboratory (JPL) Parallel Eigensolvers package.

macopt library header file	release 1.1	gradient-based optimizer

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
Modified by Chris Bowen, Gerhard Klimeck, Fabiano A. Oyafuso and E. Robert Tisdale
*/

#include "nml_dcvector.h"
#include "rvector.h"
#include "cvector.h"
#include "complex.h"


#define RESTART_GEKCO

typedef struct {	/* structure for macopt				*/
  double tol;		/* Convergence is declared
    when the gradient vector is smaller in magnitude than this
    or when the mean absolute step is less than this (see above)	*/
  double grad_tol_tiny;	/* If gradient is less than this,
    we definitely stop, even if we are not using a gradient tolerance.	*/
  double step_tol_tiny;	/* If step is less than this,
    we stop, even if we are not using a step tolerance.			*/
  int end_if_small_step;/* This defines the role of tol --
    the alternative is end_on_small_grad.				*/
  int its;		/* number of its				*/
  int itmax;		/* max						*/
  int rich;		/* This determines whether to do
    the extra gradient evaluation at the beginning of each new line min.*/
  int verbose;
  double stepmax;	/* largest step permitted (not used in macopt).	*/

  int linmin_maxits;	/* in maclinmin					*/
  nml_dcscalar linmin_g1;
			/* factors for growing and shrinking the interval */
  nml_dcscalar linmin_g2;
  nml_dcscalar linmin_g3;
  nml_dcscalar lastx;	/* keeps track of typical step length.		*/
  nml_dcscalar lastx_default;
			/* If maclinmin is reset, lastx is set to this.	*/

/* These should not be touched by the user.
   They are handy pointers for macopt to use.				*/
  nml_dcscalar gtyp;	/* stores the rms gradient for linmin		*/
  nml_dcvector* pt;
  nml_dcvector* gx;
  nml_dcvector* gy;
  nml_dcvector* gunused;
  nml_dcvector* xi;
  nml_dcvector* g;
  nml_dcvector* h;
#ifdef	RESTART_GEKCO
  nml_dcvector* x_restart;
#endif/*RESTART_GEKCO	*/
  int n;		/* dimension of parameter space			*/
  int restart;		/* whether to restart macopt - fresh cg directions */
  } cmacopt_args;

/* lastx :--- 1.0 might make general sense, (cf N.R.)
   but the best setting of all is to have a prior idea of the eigenvalues.
   If the objective function is equal to sum of N terms then set this to 1/N,
   for example Err on the small side to be conservative.		*/
void cmacoptII(
   nml_dcvector*,	/* starting vector				*/
   int,			/* number of dimensions				*/
   void (*dfunc)(nml_dcvector*, nml_dcvector*, void *),
			/* evaluates the gradient			*/
   void*,
   cmacopt_args*
   );

void cmaccheckgrad(
   nml_dcvector*,
   int,
   nml_dcscalar,
   nml_dcscalar (*f)(nml_dcvector*, void *),
   void *,
   void (*g)(nml_dcvector*, nml_dcvector*, void *),
   void *,
   int
   );	

void cmacopt_defaults(
   cmacopt_args*,
   int,
   int,
   double,
   int
   );

void cmacopt_free(cmacopt_args*);

/* The following functions could be declared static within macopt.c	*/

nml_dcscalar cmaclinminII(
   nml_dcvector*,
   void (*dfunc)(nml_dcvector*, nml_dcvector*, void *),
			/* evaluates the gradient			*/
   void*,
   cmacopt_args*
   );

nml_dcscalar cmacprodII(
   nml_dcvector*,
   nml_dcvector*,
   nml_dcscalar,
   void (*dfunc)(nml_dcvector*, nml_dcvector*, void *),
   void*,
   cmacopt_args*
   );

void cmacopt_restart(cmacopt_args *, int);

#ifdef	RESTART_GEKCO
void set_cmacopt_args(cmacopt_args *a);
cmacopt_args *get_cmacopt_args(void);
#endif/*RESTART_GEKCO	*/

#endif/*CMACOPT_H	*/

