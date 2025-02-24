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
$Header: /repo/nemo3d/src/ham/mb_ham_spds_diag.h,v 1.2 2003/10/08 16:17:24 hook Exp $
*****************************************************************************/

/*  diag(  vpar_diag, int basis_diag )  */

static int diag_map[20];
diag_map[STa]       = 4;
diag_map[Sa]        = 0;
diag_map[Pxa]       = 1;
diag_map[Pya]       = 1;
diag_map[Pza]       = 1;
diag_map[Dxya]      = 6;
diag_map[Dyza]      = 6;
diag_map[Dzxa]      = 6;
diag_map[Dx2my2a]   = 6;
diag_map[Dz2a]      = 6;
diag_map[STc]       = 5;
diag_map[Sc]        = 2;
diag_map[Pxc]       = 3;
diag_map[Pyc]       = 3;
diag_map[Pzc]       = 3;
diag_map[Dxyc]      = 7;
diag_map[Dyzc]      = 7;
diag_map[Dzxc]      = 7;
diag_map[Dx2my2c]   = 7;
diag_map[Dz2c]      = 7;


