#ifndef	RITZ_H
#define	RITZ_H 1

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
$Header: /repo/eigen/include/ritz.h,v 1.3 2003/10/08 16:15:48 hook Exp $
*****************************************************************************/


#include "io.h"
#include "mpi_local.h"
#include "macopt.h"
#include "macopt_par.h"

#include "nml_dcvector.h"

#define	RESTART_GEKCO	
#ifdef	RESTART_GEKCO	
#include "cmacopt.h"
#endif /* RESTART_GEKCO */

#ifdef	_GET_RID_OF_THIS
int	min_rayleigh(qd_struct d);
double	calc_rayleigh(nml_dvector* p, void* data);
double	calc_rayleigh_sq(nml_dvector* p, void* data);
void	calc_rayleigh_d(nml_dvector* p, nml_dvector* xi, void* data);
complex	calc_rayleigh_compl(nml_dcvector* x, int sz, qd_struct d);
double dum_invIter(
    qd_struct	 d,
    nml_dcvector *p,
    double	 eval_loc,
    double	 tol,
    int		 itnum
    );
#endif /* _GET_RID_OF_THIS */

void	gradRayleighQuotient(
    nml_dvector	*p,
    int		sz,
    nml_dvector	*xi, 
    void        *argList,
    int		itercount,
    double	*onevalue
    );

void	gradRayleighQuotient_cmplx(
    nml_dcvector *p,
    int		sz,
    nml_dcvector *xi, 
    void        *argList, 
    int		itercount,
    double      *onevalue
    );

#endif /* RITZ_H */

