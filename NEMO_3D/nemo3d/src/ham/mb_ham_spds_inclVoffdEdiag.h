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
$Header: /repo/nemo3d/src/ham/mb_ham_spds_inclVoffdEdiag.h,v 1.4 2007/03/23 16:12:15 bhaley Exp $
*****************************************************************************/

#include "mb_ham_spds_VoffdEdiag_constructors.h"

/*
 * Inline code that returns the off-diagonal matrix elements
 * as well as the additional corrections to the diagonal elements
 * due to strain.
 */

    switch (basis1_econst){
    case STa:
        switch(basis2_econst){
           case STc: get_STa_STc(Out_econst, Out_diag_corr); break;
           case Sc: get_STa_Sc(Out_econst, Out_diag_corr); break;
           case Pxc: get_STa_Pxc(Out_econst, Out_diag_corr); break;
           case Pyc: get_STa_Pyc(Out_econst, Out_diag_corr); break;
           case Pzc: get_STa_Pzc(Out_econst, Out_diag_corr); break;
           case Dxyc: get_STa_Dxyc(Out_econst, Out_diag_corr); break;
           case Dyzc: get_STa_Dyzc(Out_econst, Out_diag_corr); break;
           case Dzxc: get_STa_Dzxc(Out_econst, Out_diag_corr); break;
           case Dx2my2c: get_STa_Dx2my2c(Out_econst, Out_diag_corr); break;
           case Dz2c: get_STa_Dz2c(Out_econst, Out_diag_corr); break;
        }
        break;
    case Sa:
        switch(basis2_econst){
           case STc: get_Sa_STc(Out_econst, Out_diag_corr); break;
           case Sc: get_Sa_Sc(Out_econst, Out_diag_corr); break;
           case Pxc: get_Sa_Pxc(Out_econst, Out_diag_corr); break;
           case Pyc: get_Sa_Pyc(Out_econst, Out_diag_corr); break;
           case Pzc: get_Sa_Pzc(Out_econst, Out_diag_corr); break;
           case Dxyc: get_Sa_Dxyc(Out_econst, Out_diag_corr); break;
           case Dyzc: get_Sa_Dyzc(Out_econst, Out_diag_corr); break;
           case Dzxc: get_Sa_Dzxc(Out_econst, Out_diag_corr); break;
           case Dx2my2c: get_Sa_Dx2my2c(Out_econst, Out_diag_corr); break;
           case Dz2c: get_Sa_Dz2c(Out_econst, Out_diag_corr); break;
        }
        break;
    case Pxa:
        switch(basis2_econst){
           case STc: get_Pxa_STc(Out_econst, Out_diag_corr); break;
           case Sc: get_Pxa_Sc(Out_econst, Out_diag_corr); break;
           case Pxc: get_Pxa_Pxc(Out_econst, Out_diag_corr); break;
           case Pyc: get_Pxa_Pyc(Out_econst, Out_diag_corr); break;
           case Pzc: get_Pxa_Pzc(Out_econst, Out_diag_corr); break;
           case Dxyc: get_Pxa_Dxyc(Out_econst, Out_diag_corr); break;
           case Dyzc: get_Pxa_Dyzc(Out_econst, Out_diag_corr); break;
           case Dzxc: get_Pxa_Dzxc(Out_econst, Out_diag_corr); break;
           case Dx2my2c: get_Pxa_Dx2my2c(Out_econst, Out_diag_corr); break;
           case Dz2c: get_Pxa_Dz2c(Out_econst, Out_diag_corr); break;
        }
        break;
    case Pya:
        switch(basis2_econst){
           case STc: get_Pya_STc(Out_econst, Out_diag_corr); break;
           case Sc: get_Pya_Sc(Out_econst, Out_diag_corr); break;
           case Pxc: get_Pya_Pxc(Out_econst, Out_diag_corr); break;
           case Pyc: get_Pya_Pyc(Out_econst, Out_diag_corr); break;
           case Pzc: get_Pya_Pzc(Out_econst, Out_diag_corr); break;
           case Dxyc: get_Pya_Dxyc(Out_econst, Out_diag_corr); break;
           case Dyzc: get_Pya_Dyzc(Out_econst, Out_diag_corr); break;
           case Dzxc: get_Pya_Dzxc(Out_econst, Out_diag_corr); break;
           case Dx2my2c: get_Pya_Dx2my2c(Out_econst, Out_diag_corr); break;
           case Dz2c: get_Pya_Dz2c(Out_econst, Out_diag_corr); break;
        }
        break;
    case Pza:
        switch(basis2_econst){
           case STc: get_Pza_STc(Out_econst, Out_diag_corr); break;
           case Sc: get_Pza_Sc(Out_econst, Out_diag_corr); break;
           case Pxc: get_Pza_Pxc(Out_econst, Out_diag_corr); break;
           case Pyc: get_Pza_Pyc(Out_econst, Out_diag_corr); break;
           case Pzc: get_Pza_Pzc(Out_econst, Out_diag_corr); break;
           case Dxyc: get_Pza_Dxyc(Out_econst, Out_diag_corr); break;
           case Dyzc: get_Pza_Dyzc(Out_econst, Out_diag_corr); break;
           case Dzxc: get_Pza_Dzxc(Out_econst, Out_diag_corr); break;
           case Dx2my2c: get_Pza_Dx2my2c(Out_econst, Out_diag_corr); break;
           case Dz2c: get_Pza_Dz2c(Out_econst, Out_diag_corr); break;
        }
        break;
    case Dxya:
        switch(basis2_econst){
           case STc: get_Dxya_STc(Out_econst, Out_diag_corr); break;
           case Sc: get_Dxya_Sc(Out_econst, Out_diag_corr); break;
           case Pxc: get_Dxya_Pxc(Out_econst, Out_diag_corr); break;
           case Pyc: get_Dxya_Pyc(Out_econst, Out_diag_corr); break;
           case Pzc: get_Dxya_Pzc(Out_econst, Out_diag_corr); break;
           case Dxyc: get_Dxya_Dxyc(Out_econst, Out_diag_corr); break;
           case Dyzc: get_Dxya_Dyzc(Out_econst, Out_diag_corr); break;
           case Dzxc: get_Dxya_Dzxc(Out_econst, Out_diag_corr); break;
           case Dx2my2c: get_Dxya_Dx2my2c(Out_econst, Out_diag_corr); break;
           case Dz2c: get_Dxya_Dz2c(Out_econst, Out_diag_corr); break;
        }
        break;
    case Dyza:
        switch(basis2_econst){
           case STc: get_Dyza_STc(Out_econst, Out_diag_corr); break;
           case Sc: get_Dyza_Sc(Out_econst, Out_diag_corr); break;
           case Pxc: get_Dyza_Pxc(Out_econst, Out_diag_corr); break;
           case Pyc: get_Dyza_Pyc(Out_econst, Out_diag_corr); break;
           case Pzc: get_Dyza_Pzc(Out_econst, Out_diag_corr); break;
           case Dxyc: get_Dyza_Dxyc(Out_econst, Out_diag_corr); break;
           case Dyzc: get_Dyza_Dyzc(Out_econst, Out_diag_corr); break;
           case Dzxc: get_Dyza_Dzxc(Out_econst, Out_diag_corr); break;
           case Dx2my2c: get_Dyza_Dx2my2c(Out_econst, Out_diag_corr); break;
           case Dz2c: get_Dyza_Dz2c(Out_econst, Out_diag_corr); break;
        }
        break;
    case Dzxa:
        switch(basis2_econst){
           case STc: get_Dzxa_STc(Out_econst, Out_diag_corr); break;
           case Sc: get_Dzxa_Sc(Out_econst, Out_diag_corr); break;
           case Pxc: get_Dzxa_Pxc(Out_econst, Out_diag_corr); break;
           case Pyc: get_Dzxa_Pyc(Out_econst, Out_diag_corr); break;
           case Pzc: get_Dzxa_Pzc(Out_econst, Out_diag_corr); break;
           case Dxyc: get_Dzxa_Dxyc(Out_econst, Out_diag_corr); break;
           case Dyzc: get_Dzxa_Dyzc(Out_econst, Out_diag_corr); break;
           case Dzxc: get_Dzxa_Dzxc(Out_econst, Out_diag_corr); break;
           case Dx2my2c: get_Dzxa_Dx2my2c(Out_econst, Out_diag_corr); break;
           case Dz2c: get_Dzxa_Dz2c(Out_econst, Out_diag_corr); break;
        }
        break;
    case Dx2my2a:
        switch(basis2_econst){
           case STc: get_Dx2my2a_STc(Out_econst, Out_diag_corr); break;
           case Sc: get_Dx2my2a_Sc(Out_econst, Out_diag_corr); break;
           case Pxc: get_Dx2my2a_Pxc(Out_econst, Out_diag_corr); break;
           case Pyc: get_Dx2my2a_Pyc(Out_econst, Out_diag_corr); break;
           case Pzc: get_Dx2my2a_Pzc(Out_econst, Out_diag_corr); break;
           case Dxyc: get_Dx2my2a_Dxyc(Out_econst, Out_diag_corr); break;
           case Dyzc: get_Dx2my2a_Dyzc(Out_econst, Out_diag_corr); break;
           case Dzxc: get_Dx2my2a_Dzxc(Out_econst, Out_diag_corr); break;
           case Dx2my2c: get_Dx2my2a_Dx2my2c(Out_econst, Out_diag_corr); break;
           case Dz2c: get_Dx2my2a_Dz2c(Out_econst, Out_diag_corr); break;
        }
        break;
    case Dz2a:
        switch(basis2_econst){
           case STc: get_Dz2a_STc(Out_econst, Out_diag_corr); break;
           case Sc: get_Dz2a_Sc(Out_econst, Out_diag_corr); break;
           case Pxc: get_Dz2a_Pxc(Out_econst, Out_diag_corr); break;
           case Pyc: get_Dz2a_Pyc(Out_econst, Out_diag_corr); break;
           case Pzc: get_Dz2a_Pzc(Out_econst, Out_diag_corr); break;
           case Dxyc: get_Dz2a_Dxyc(Out_econst, Out_diag_corr); break;
           case Dyzc: get_Dz2a_Dyzc(Out_econst, Out_diag_corr); break;
           case Dzxc: get_Dz2a_Dzxc(Out_econst, Out_diag_corr); break;
           case Dx2my2c: get_Dz2a_Dx2my2c(Out_econst, Out_diag_corr); break;
           case Dz2c: get_Dz2a_Dz2c(Out_econst, Out_diag_corr); break;
        }
        break;
    case STc:
        switch(basis2_econst){
           case STa: get_STc_STa(Out_econst, Out_diag_corr); break;
           case Sa: get_STc_Sa(Out_econst, Out_diag_corr); break;
           case Pxa: get_STc_Pxa(Out_econst, Out_diag_corr); break;
           case Pya: get_STc_Pya(Out_econst, Out_diag_corr); break;
           case Pza: get_STc_Pza(Out_econst, Out_diag_corr); break;
           case Dxya: get_STc_Dxya(Out_econst, Out_diag_corr); break;
           case Dyza: get_STc_Dyza(Out_econst, Out_diag_corr); break;
           case Dzxa: get_STc_Dzxa(Out_econst, Out_diag_corr); break;
           case Dx2my2a: get_STc_Dx2my2a(Out_econst, Out_diag_corr); break;
           case Dz2a: get_STc_Dz2a(Out_econst, Out_diag_corr); break;
        }
        break;
    case Sc:
        switch(basis2_econst){
           case STa: get_Sc_STa(Out_econst, Out_diag_corr); break;
           case Sa: get_Sc_Sa(Out_econst, Out_diag_corr); break;
           case Pxa: get_Sc_Pxa(Out_econst, Out_diag_corr); break;
           case Pya: get_Sc_Pya(Out_econst, Out_diag_corr); break;
           case Pza: get_Sc_Pza(Out_econst, Out_diag_corr); break;
           case Dxya: get_Sc_Dxya(Out_econst, Out_diag_corr); break;
           case Dyza: get_Sc_Dyza(Out_econst, Out_diag_corr); break;
           case Dzxa: get_Sc_Dzxa(Out_econst, Out_diag_corr); break;
           case Dx2my2a: get_Sc_Dx2my2a(Out_econst, Out_diag_corr); break;
           case Dz2a: get_Sc_Dz2a(Out_econst, Out_diag_corr); break;
        }
        break;
    case Pxc:
        switch(basis2_econst){
           case STa: get_Pxc_STa(Out_econst, Out_diag_corr); break;
           case Sa: get_Pxc_Sa(Out_econst, Out_diag_corr); break;
           case Pxa: get_Pxc_Pxa(Out_econst, Out_diag_corr); break;
           case Pya: get_Pxc_Pya(Out_econst, Out_diag_corr); break;
           case Pza: get_Pxc_Pza(Out_econst, Out_diag_corr); break;
           case Dxya: get_Pxc_Dxya(Out_econst, Out_diag_corr); break;
           case Dyza: get_Pxc_Dyza(Out_econst, Out_diag_corr); break;
           case Dzxa: get_Pxc_Dzxa(Out_econst, Out_diag_corr); break;
           case Dx2my2a: get_Pxc_Dx2my2a(Out_econst, Out_diag_corr); break;
           case Dz2a: get_Pxc_Dz2a(Out_econst, Out_diag_corr); break;
        }
        break;
    case Pyc:
        switch(basis2_econst){
           case STa: get_Pyc_STa(Out_econst, Out_diag_corr); break;
           case Sa: get_Pyc_Sa(Out_econst, Out_diag_corr); break;
           case Pxa: get_Pyc_Pxa(Out_econst, Out_diag_corr); break;
           case Pya: get_Pyc_Pya(Out_econst, Out_diag_corr); break;
           case Pza: get_Pyc_Pza(Out_econst, Out_diag_corr); break;
           case Dxya: get_Pyc_Dxya(Out_econst, Out_diag_corr); break;
           case Dyza: get_Pyc_Dyza(Out_econst, Out_diag_corr); break;
           case Dzxa: get_Pyc_Dzxa(Out_econst, Out_diag_corr); break;
           case Dx2my2a: get_Pyc_Dx2my2a(Out_econst, Out_diag_corr); break;
           case Dz2a: get_Pyc_Dz2a(Out_econst, Out_diag_corr); break;
        }
        break;
    case Pzc:
        switch(basis2_econst){
           case STa: get_Pzc_STa(Out_econst, Out_diag_corr); break;
           case Sa: get_Pzc_Sa(Out_econst, Out_diag_corr); break;
           case Pxa: get_Pzc_Pxa(Out_econst, Out_diag_corr); break;
           case Pya: get_Pzc_Pya(Out_econst, Out_diag_corr); break;
           case Pza: get_Pzc_Pza(Out_econst, Out_diag_corr); break;
           case Dxya: get_Pzc_Dxya(Out_econst, Out_diag_corr); break;
           case Dyza: get_Pzc_Dyza(Out_econst, Out_diag_corr); break;
           case Dzxa: get_Pzc_Dzxa(Out_econst, Out_diag_corr); break;
           case Dx2my2a: get_Pzc_Dx2my2a(Out_econst, Out_diag_corr); break;
           case Dz2a: get_Pzc_Dz2a(Out_econst, Out_diag_corr); break;
        }
        break;
    case Dxyc:
        switch(basis2_econst){
           case STa: get_Dxyc_STa(Out_econst, Out_diag_corr); break;
           case Sa: get_Dxyc_Sa(Out_econst, Out_diag_corr); break;
           case Pxa: get_Dxyc_Pxa(Out_econst, Out_diag_corr); break;
           case Pya: get_Dxyc_Pya(Out_econst, Out_diag_corr); break;
           case Pza: get_Dxyc_Pza(Out_econst, Out_diag_corr); break;
           case Dxya: get_Dxyc_Dxya(Out_econst, Out_diag_corr); break;
           case Dyza: get_Dxyc_Dyza(Out_econst, Out_diag_corr); break;
           case Dzxa: get_Dxyc_Dzxa(Out_econst, Out_diag_corr); break;
           case Dx2my2a: get_Dxyc_Dx2my2a(Out_econst, Out_diag_corr); break;
           case Dz2a: get_Dxyc_Dz2a(Out_econst, Out_diag_corr); break;
        }
        break;
    case Dyzc:
        switch(basis2_econst){
           case STa: get_Dyzc_STa(Out_econst, Out_diag_corr); break;
           case Sa: get_Dyzc_Sa(Out_econst, Out_diag_corr); break;
           case Pxa: get_Dyzc_Pxa(Out_econst, Out_diag_corr); break;
           case Pya: get_Dyzc_Pya(Out_econst, Out_diag_corr); break;
           case Pza: get_Dyzc_Pza(Out_econst, Out_diag_corr); break;
           case Dxya: get_Dyzc_Dxya(Out_econst, Out_diag_corr); break;
           case Dyza: get_Dyzc_Dyza(Out_econst, Out_diag_corr); break;
           case Dzxa: get_Dyzc_Dzxa(Out_econst, Out_diag_corr); break;
           case Dx2my2a: get_Dyzc_Dx2my2a(Out_econst, Out_diag_corr); break;
           case Dz2a: get_Dyzc_Dz2a(Out_econst, Out_diag_corr); break;
        }
        break;
    case Dzxc:
        switch(basis2_econst){
           case STa: get_Dzxc_STa(Out_econst, Out_diag_corr); break;
           case Sa: get_Dzxc_Sa(Out_econst, Out_diag_corr); break;
           case Pxa: get_Dzxc_Pxa(Out_econst, Out_diag_corr); break;
           case Pya: get_Dzxc_Pya(Out_econst, Out_diag_corr); break;
           case Pza: get_Dzxc_Pza(Out_econst, Out_diag_corr); break;
           case Dxya: get_Dzxc_Dxya(Out_econst, Out_diag_corr); break;
           case Dyza: get_Dzxc_Dyza(Out_econst, Out_diag_corr); break;
           case Dzxa: get_Dzxc_Dzxa(Out_econst, Out_diag_corr); break;
           case Dx2my2a: get_Dzxc_Dx2my2a(Out_econst, Out_diag_corr); break;
           case Dz2a: get_Dzxc_Dz2a(Out_econst, Out_diag_corr); break;
        }
        break;
    case Dx2my2c:
        switch(basis2_econst){
           case STa: get_Dx2my2c_STa(Out_econst, Out_diag_corr); break;
           case Sa: get_Dx2my2c_Sa(Out_econst, Out_diag_corr); break;
           case Pxa: get_Dx2my2c_Pxa(Out_econst, Out_diag_corr); break;
           case Pya: get_Dx2my2c_Pya(Out_econst, Out_diag_corr); break;
           case Pza: get_Dx2my2c_Pza(Out_econst, Out_diag_corr); break;
           case Dxya: get_Dx2my2c_Dxya(Out_econst, Out_diag_corr); break;
           case Dyza: get_Dx2my2c_Dyza(Out_econst, Out_diag_corr); break;
           case Dzxa: get_Dx2my2c_Dzxa(Out_econst, Out_diag_corr); break;
           case Dx2my2a: get_Dx2my2c_Dx2my2a(Out_econst, Out_diag_corr); break;
           case Dz2a: get_Dx2my2c_Dz2a(Out_econst, Out_diag_corr); break;
        }
        break;
    case Dz2c:
        switch(basis2_econst){
           case STa: get_Dz2c_STa(Out_econst, Out_diag_corr); break;
           case Sa: get_Dz2c_Sa(Out_econst, Out_diag_corr); break;
           case Pxa: get_Dz2c_Pxa(Out_econst, Out_diag_corr); break;
           case Pya: get_Dz2c_Pya(Out_econst, Out_diag_corr); break;
           case Pza: get_Dz2c_Pza(Out_econst, Out_diag_corr); break;
           case Dxya: get_Dz2c_Dxya(Out_econst, Out_diag_corr); break;
           case Dyza: get_Dz2c_Dyza(Out_econst, Out_diag_corr); break;
           case Dzxa: get_Dz2c_Dzxa(Out_econst, Out_diag_corr); break;
           case Dx2my2a: get_Dz2c_Dx2my2a(Out_econst, Out_diag_corr); break;
           case Dz2a: get_Dz2c_Dz2a(Out_econst, Out_diag_corr); break;
        }
        break;
    }
