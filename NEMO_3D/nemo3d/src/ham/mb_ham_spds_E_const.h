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
$Header: /repo/nemo3d/src/ham/mb_ham_spds_E_const.h,v 1.2 2003/10/08 16:17:23
hook Exp $
*****************************************************************************/

/*
 * Inline code that returns the off-diagonal matrix elements
 */

/* E_const_get2( Out_econst , vpar_e_const, basis1_econst, basis2_econst, l, m,
 * n, l2, m2, n2)      */
switch (basis1_econst) {
case STa:
  switch (basis2_econst) {
  case STc:
    Out_econst = ststs;
    break;
  case Sc:
    Out_econst = stass;
    break;
  case Pxc:
    Out_econst = l * stapcs;
    break;
  case Pyc:
    Out_econst = m * stapcs;
    break;
  case Pzc:
    Out_econst = n * stapcs;
    break;
  case Dxyc:
    Out_econst = sqr3 * l * m * stadcs;
    break;
  case Dyzc:
    Out_econst = sqr3 * m * n * stadcs;
    break;
  case Dzxc:
    Out_econst = sqr3 * n * l * stadcs;
    break;
  case Dx2my2c:
    Out_econst = 0.5 * sqr3 * (l2 - m2) * stadcs;
    break;
  case Dz2c:
    Out_econst = (n2 - 0.5 * (l2 + m2)) * stadcs;
  }
  break;
case Sa:
  switch (basis2_econst) {
  /* case Sa is for single band */
  /* case Sa:       Out_econst = sss;
     break; */
  case STc:
    Out_econst = sasts;
    break;
  case Sc:
    Out_econst = sss;
    break;
  case Pxc:
    Out_econst = l * sapcs;
    break;
  case Pyc:
    Out_econst = m * sapcs;
    break;
  case Pzc:
    Out_econst = n * sapcs;
    break;
  case Dxyc:
    Out_econst = sqr3 * l * m * sadcs;
    break;
  case Dyzc:
    Out_econst = sqr3 * m * n * sadcs;
    break;
  case Dzxc:
    Out_econst = sqr3 * l * n * sadcs;
    break;
  case Dx2my2c:
    Out_econst = 0.5 * sqr3 * (l2 - m2) * sadcs;
    break;
  case Dz2c:
    Out_econst = (n2 - 0.5 * (l2 + m2)) * sadcs;
  }
  break;
case Pxa:
  switch (basis2_econst) {
  case STc:
    Out_econst = l * stcpas;
    break;
  case Sc:
    Out_econst = l * scpas;
    break;
  case Pxc:
    Out_econst = l2 * pps + (1 - l2) * ppp;
    break;
  case Pyc:
    Out_econst = l * m * pps - l * m * ppp;
    break;
  case Pzc:
    Out_econst = l * n * pps - l * n * ppp;
    break;
  case Dxyc:
    Out_econst = sqr3 * l2 * m * padcs + m * (1 - 2 * l2) * padcp;
    break;
  case Dyzc:
    Out_econst = sqr3 * l * m * n * padcs - 2 * l * m * n * padcp;
    break;
  case Dzxc:
    Out_econst = sqr3 * l2 * n * padcs + n * (1 - 2 * l2) * padcp;
    break;
  case Dx2my2c:
    Out_econst = 0.5 * sqr3 * l * (l2 - m2) * padcs + l * (1 - l2 + m2) * padcp;
    break;
  case Dz2c:
    Out_econst = l * (n2 - 0.5 * (l2 + m2)) * padcs - sqr3 * l * n2 * padcp;
  }
  break;
case Pya:
  switch (basis2_econst) {
  case STc:
    Out_econst = m * stcpas;
    break;
  case Sc:
    Out_econst = m * scpas;
    break;
  case Pxc:
    Out_econst = l * m * pps - l * m * ppp;
    break;
  case Pyc:
    Out_econst = m2 * pps + (1 - m2) * ppp;
    break;
  case Pzc:
    Out_econst = m * n * pps - m * n * ppp;
    break;
  case Dxyc:
    Out_econst = sqr3 * l * m2 * padcs + l * (1 - 2 * m2) * padcp;
    break;
  case Dyzc:
    Out_econst = sqr3 * n * m2 * padcs + n * (1 - 2 * m2) * padcp;
    break;
  case Dzxc:
    Out_econst = sqr3 * l * m * n * padcs - 2 * l * m * n * padcp;
    break;
  case Dx2my2c:
    Out_econst = 0.5 * sqr3 * m * (l2 - m2) * padcs - m * (1 + l2 - m2) * padcp;
    break;
  case Dz2c:
    Out_econst = m * (n2 - 0.5 * (l2 + m2)) * padcs - sqr3 * m * n2 * padcp;
  }
  break;
case Pza:
  switch (basis2_econst) {
  case STc:
    Out_econst = n * stcpas;
    break;
  case Sc:
    Out_econst = n * scpas;
    break;
  case Pxc:
    Out_econst = l * n * pps - l * n * ppp;
    break;
  case Pyc:
    Out_econst = m * n * pps - m * n * ppp;
    break;
  case Pzc:
    Out_econst = n2 * pps + (1 - n2) * ppp;
    break;
  case Dxyc:
    Out_econst = sqr3 * l * m * n * padcs - 2 * l * m * n * padcp;
    break;
  case Dyzc:
    Out_econst = sqr3 * m * n2 * padcs + m * (1 - 2 * n2) * padcp;
    break;
  case Dzxc:
    Out_econst = sqr3 * l * n2 * padcs + l * (1 - 2 * n2) * padcp;
    break;
  case Dx2my2c:
    Out_econst = 0.5 * sqr3 * n * (l2 - m2) * padcs - n * (l2 - m2) * padcp;
    break;
  case Dz2c:
    Out_econst =
        n * (n2 - 0.5 * (l2 + m2)) * padcs + sqr3 * n * (l2 + m2) * padcp;
  }
  break;
case Dxya:
  switch (basis2_econst) {
  case STc:
    Out_econst = sqr3 * l * m * stcdas;
    break;
  case Sc:
    Out_econst = sqr3 * l * m * scdas;
    break;
  case Pxc:
    Out_econst = sqr3 * l2 * m * pcdas + m * (1 - 2 * l2) * pcdap;
    break;
  case Pyc:
    Out_econst = sqr3 * l * m2 * pcdas + l * (1 - 2 * m2) * pcdap;
    break;
  case Pzc:
    Out_econst = sqr3 * l * m * n * pcdas - 2 * l * m * n * pcdap;
    break;
  case Dxyc:
    Out_econst = 3 * l2 * m2 * dds + (l2 + m2 - 4 * l2 * m2) * ddp +
                 (n2 + l2 * m2) * ddd;
    break;
  case Dyzc:
    Out_econst = 3 * l * m2 * n * dds + l * n * (1 - 4 * m2) * ddp +
                 l * n * (m2 - 1) * ddd;
    break;
  case Dzxc:
    Out_econst = 3 * l2 * m * n * dds + m * n * (1 - 4 * l2) * ddp +
                 m * n * (l2 - 1) * ddd;
    break;
  case Dx2my2c:
    Out_econst = 1.5 * l * m * (l2 - m2) * dds + 2 * l * m * (m2 - l2) * ddp +
                 0.5 * l * m * (l2 - m2) * ddd;
    break;
  case Dz2c:
    Out_econst = sqr3 * l * m * (n2 - 0.5 * (l2 + m2)) * dds -
                 sqr3 * 2 * l * m * n2 * ddp +
                 0.5 * sqr3 * l * m * (1 + n2) * ddd;
  }
  break;
case Dyza:
  switch (basis2_econst) {
  case STc:
    Out_econst = sqr3 * m * n * stcdas;
    break;
  case Sc:
    Out_econst = sqr3 * m * n * scdas;
    break;
  case Pxc:
    Out_econst = sqr3 * l * m * n * pcdas - 2 * l * m * n * pcdap;
    break;
  case Pyc:
    Out_econst = sqr3 * n * m2 * pcdas + n * (1 - 2 * m2) * pcdap;
    break;
  case Pzc:
    Out_econst = sqr3 * m * n2 * pcdas + m * (1 - 2 * n2) * pcdap;
    break;
  case Dxyc:
    Out_econst = 3 * l * m2 * n * dds + l * n * (1 - 4 * m2) * ddp +
                 l * n * (m2 - 1) * ddd;
    break;
  case Dyzc:
    Out_econst = 3 * m2 * n2 * dds + (m2 + n2 - 4 * m2 * n2) * ddp +
                 (l2 + m2 * n2) * ddd;
    break;
  case Dzxc:
    Out_econst = 3 * m * n2 * l * dds + m * l * (1 - 4 * n2) * ddp +
                 m * l * (n2 - 1) * ddd;
    break;
  case Dx2my2c:
    Out_econst = 1.5 * m * n * (l2 - m2) * dds -
                 m * n * (1 + 2 * (l2 - m2)) * ddp +
                 m * n * (1 + 0.5 * (l2 - m2)) * ddd;
    break;
  case Dz2c:
    Out_econst = sqr3 * m * n * (n2 - 0.5 * (l2 + m2)) * dds +
                 sqr3 * m * n * (l2 + m2 - n2) * ddp -
                 0.5 * sqr3 * m * n * (l2 + m2) * ddd;
  }
  break;
case Dzxa:
  switch (basis2_econst) {
  case STc:
    Out_econst = sqr3 * n * l * stcdas;
    break;
  case Sc:
    Out_econst = sqr3 * n * l * scdas;
    break;
  case Pxc:
    Out_econst = sqr3 * l2 * n * pcdas + n * (1 - 2 * l2) * pcdap;
    break;
  case Pyc:
    Out_econst = sqr3 * l * m * n * pcdas - 2 * l * m * n * pcdap;
    break;
  case Pzc:
    Out_econst = sqr3 * l * n2 * pcdas + l * (1 - 2 * n2) * pcdap;
    break;
  case Dxyc:
    Out_econst = 3 * l2 * m * n * dds + m * n * (1 - 4 * l2) * ddp +
                 m * n * (l2 - 1) * ddd;
    break;
  case Dyzc:
    Out_econst = 3 * n2 * l * m * dds + l * m * (1 - 4 * n2) * ddp +
                 l * m * (n2 - 1) * ddd;
    break;
  case Dzxc:
    Out_econst = 3 * n2 * l2 * dds + (n2 + l2 - 4 * n2 * l2) * ddp +
                 (m2 + n2 * l2) * ddd;
    break;
  case Dx2my2c:
    Out_econst = 1.5 * n * l * (l2 - m2) * dds +
                 n * l * (1 - 2 * (l2 - m2)) * ddp -
                 n * l * (1 - 0.5 * (l2 - m2)) * ddd;
    break;
  case Dz2c:
    Out_econst = sqr3 * l * n * (n2 - 0.5 * (l2 + m2)) * dds +
                 sqr3 * l * n * (l2 + m2 - n2) * ddp -
                 0.5 * sqr3 * l * n * (l2 + m2) * ddd;
  }
  break;
case Dx2my2a:
  switch (basis2_econst) {
  case STc:
    Out_econst = 0.5 * sqr3 * (l2 - m2) * stcdas;
    break;
  case Sc:
    Out_econst = 0.5 * sqr3 * (l2 - m2) * scdas;
    break;
  case Pxc:
    Out_econst = 0.5 * sqr3 * l * (l2 - m2) * pcdas + l * (1 - l2 + m2) * pcdap;
    break;
  case Pyc:
    Out_econst = 0.5 * sqr3 * m * (l2 - m2) * pcdas - m * (1 + l2 - m2) * pcdap;
    break;
  case Pzc:
    Out_econst = 0.5 * sqr3 * n * (l2 - m2) * pcdas - n * (l2 - m2) * pcdap;
    break;
  case Dxyc:
    Out_econst = 1.5 * l * m * (l2 - m2) * dds + 2 * l * m * (m2 - l2) * ddp +
                 0.5 * l * m * (l2 - m2) * ddd;
    break;
  case Dyzc:
    Out_econst = 1.5 * m * n * (l2 - m2) * dds -
                 m * n * (1 + 2 * (l2 - m2)) * ddp +
                 m * n * (1 + 0.5 * (l2 - m2)) * ddd;
    break;
  case Dzxc:
    Out_econst = 1.5 * n * l * (l2 - m2) * dds +
                 n * l * (1 - 2 * (l2 - m2)) * ddp -
                 n * l * (1 - 0.5 * (l2 - m2)) * ddd;
    break;
  case Dx2my2c:
    Out_econst = 0.75 * (l2 - m2) * (l2 - m2) * dds +
                 (l2 + m2 - (l2 - m2) * (l2 - m2)) * ddp +
                 (n2 + 0.25 * (l2 - m2) * (l2 - m2)) * ddd;
    break;
  case Dz2c:
    Out_econst = 0.5 * sqr3 * (l2 - m2) * (n2 - 0.5 * (l2 + m2)) * dds +
                 sqr3 * n2 * (m2 - l2) * ddp +
                 0.25 * sqr3 * (1 + n2) * (l2 - m2) * ddd;
  }
  break;
case Dz2a:
  switch (basis2_econst) {
  case STc:
    Out_econst = (n2 - 0.5 * (l2 + m2)) * stcdas;
    break;
  case Sc:
    Out_econst = (n2 - 0.5 * (l2 + m2)) * scdas;
    break;
  case Pxc:
    Out_econst = l * (n2 - 0.5 * (l2 + m2)) * pcdas - sqr3 * l * n2 * pcdap;
    break;
  case Pyc:
    Out_econst = m * (n2 - 0.5 * (l2 + m2)) * pcdas - sqr3 * m * n2 * pcdap;
    break;
  case Pzc:
    Out_econst =
        n * (n2 - 0.5 * (l2 + m2)) * pcdas + sqr3 * n * (l2 + m2) * pcdap;
    break;
  case Dxyc:
    Out_econst = sqr3 * l * m * (n2 - 0.5 * (l2 + m2)) * dds -
                 sqr3 * 2 * l * m * n2 * ddp +
                 0.5 * sqr3 * l * m * (1 + n2) * ddd;
    break;
  case Dyzc:
    Out_econst = sqr3 * m * n * (n2 - 0.5 * (l2 + m2)) * dds +
                 sqr3 * m * n * (l2 + m2 - n2) * ddp -
                 0.5 * sqr3 * m * n * (l2 + m2) * ddd;
    break;
  case Dzxc:
    Out_econst = sqr3 * l * n * (n2 - 0.5 * (l2 + m2)) * dds +
                 sqr3 * l * n * (l2 + m2 - n2) * ddp -
                 0.5 * sqr3 * l * n * (l2 + m2) * ddd;
    break;
  case Dx2my2c:
    Out_econst = 0.5 * sqr3 * (l2 - m2) * (n2 - 0.5 * (l2 + m2)) * dds +
                 sqr3 * n2 * (m2 - l2) * ddp +
                 0.25 * sqr3 * (1 + n2) * (l2 - m2) * ddd;
    break;
  case Dz2c:
    Out_econst = (n2 - 0.5 * (l2 + m2)) * (n2 - 0.5 * (l2 + m2)) * dds +
                 3 * n2 * (l2 + m2) * ddp + 0.75 * (l2 + m2) * (l2 + m2) * ddd;
  }
  break;
case STc:
  switch (basis2_econst) {
  case STa:
    Out_econst = ststs;
    break;
  case Sa:
    Out_econst = sasts;
    break;
  case Pxa:
    Out_econst = l * stcpas;
    break;
  case Pya:
    Out_econst = m * stcpas;
    break;
  case Pza:
    Out_econst = n * stcpas;
    break;
  case Dxya:
    Out_econst = sqr3 * l * m * stcdas;
    break;
  case Dyza:
    Out_econst = sqr3 * m * n * stcdas;
    break;
  case Dzxa:
    Out_econst = sqr3 * n * l * stcdas;
    break;
  case Dx2my2a:
    Out_econst = 0.5 * sqr3 * (l2 - m2) * stcdas;
    break;
  case Dz2a:
    Out_econst = (n2 - 0.5 * (l2 + m2)) * stcdas;
  }
  break;
case Sc:
  switch (basis2_econst) {
  case STa:
    Out_econst = stass;
    break;
  case Sa:
    Out_econst = sss;
    break;
  case Pxa:
    Out_econst = l * scpas;
    break;
  case Pya:
    Out_econst = m * scpas;
    break;
  case Pza:
    Out_econst = n * scpas;
    break;
  case Dxya:
    Out_econst = sqr3 * l * m * scdas;
    break;
  case Dyza:
    Out_econst = sqr3 * m * n * scdas;
    break;
  case Dzxa:
    Out_econst = sqr3 * l * n * scdas;
    break;
  case Dx2my2a:
    Out_econst = 0.5 * sqr3 * (l2 - m2) * scdas;
    break;
  case Dz2a:
    Out_econst = (n2 - 0.5 * (l2 + m2)) * scdas;
  }
  break;
case Pxc:
  switch (basis2_econst) {
  case STa:
    Out_econst = l * stapcs;
    break;
  case Sa:
    Out_econst = l * sapcs;
    break;
  case Pxa:
    Out_econst = l2 * pps + (1 - l2) * ppp;
    break;
  case Pya:
    Out_econst = l * m * pps - l * m * ppp;
    break;
  case Pza:
    Out_econst = l * n * pps - l * n * ppp;
    break;
  case Dxya:
    Out_econst = sqr3 * l2 * m * pcdas + m * (1 - 2 * l2) * pcdap;
    break;
  case Dyza:
    Out_econst = sqr3 * l * m * n * pcdas - 2 * l * m * n * pcdap;
    break;
  case Dzxa:
    Out_econst = sqr3 * l2 * n * pcdas + n * (1 - 2 * l2) * pcdap;
    break;
  case Dx2my2a:
    Out_econst = 0.5 * sqr3 * l * (l2 - m2) * pcdas + l * (1 - l2 + m2) * pcdap;
    break;
  case Dz2a:
    Out_econst = l * (n2 - 0.5 * (l2 + m2)) * pcdas - sqr3 * l * n2 * pcdap;
  }
  break;
case Pyc:
  switch (basis2_econst) {
  case STa:
    Out_econst = m * stapcs;
    break;
  case Sa:
    Out_econst = m * sapcs;
    break;
  case Pxa:
    Out_econst = l * m * pps - l * m * ppp;
    break;
  case Pya:
    Out_econst = m2 * pps + (1 - m2) * ppp;
    break;
  case Pza:
    Out_econst = m * n * pps - m * n * ppp;
    break;
  case Dxya:
    Out_econst = sqr3 * l * m2 * pcdas + l * (1 - 2 * m2) * pcdap;
    break;
  case Dyza:
    Out_econst = sqr3 * n * m2 * pcdas + n * (1 - 2 * m2) * pcdap;
    break;
  case Dzxa:
    Out_econst = sqr3 * l * m * n * pcdas - 2 * l * m * n * pcdap;
    break;
  case Dx2my2a:
    Out_econst = 0.5 * sqr3 * m * (l2 - m2) * pcdas - m * (1 + l2 - m2) * pcdap;
    break;
  case Dz2a:
    Out_econst = m * (n2 - 0.5 * (l2 + m2)) * pcdas - sqr3 * m * n2 * pcdap;
  }
  break;
case Pzc:
  switch (basis2_econst) {
  case STa:
    Out_econst = n * stapcs;
    break;
  case Sa:
    Out_econst = n * sapcs;
    break;
  case Pxa:
    Out_econst = l * n * pps - l * n * ppp;
    break;
  case Pya:
    Out_econst = m * n * pps - m * n * ppp;
    break;
  case Pza:
    Out_econst = n2 * pps + (1 - n2) * ppp;
    break;
  case Dxya:
    Out_econst = sqr3 * l * m * n * pcdas - 2 * l * m * n * pcdap;
    break;
  case Dyza:
    Out_econst = sqr3 * m * n2 * pcdas + m * (1 - 2 * n2) * pcdap;
    break;
  case Dzxa:
    Out_econst = sqr3 * l * n2 * pcdas + l * (1 - 2 * n2) * pcdap;
    break;
  case Dx2my2a:
    Out_econst = 0.5 * sqr3 * n * (l2 - m2) * pcdas - n * (l2 - m2) * pcdap;
    break;
  case Dz2a:
    Out_econst =
        n * (n2 - 0.5 * (l2 + m2)) * pcdas + sqr3 * n * (l2 + m2) * pcdap;
  }
  break;
case Dxyc:
  switch (basis2_econst) {
  case STa:
    Out_econst = sqr3 * l * m * stadcs;
    break;
  case Sa:
    Out_econst = sqr3 * l * m * sadcs;
    break;
  case Pxa:
    Out_econst = sqr3 * l2 * m * padcs + m * (1 - 2 * l2) * padcp;
    break;
  case Pya:
    Out_econst = sqr3 * l * m2 * padcs + l * (1 - 2 * m2) * padcp;
    break;
  case Pza:
    Out_econst = sqr3 * l * m * n * padcs - 2 * l * m * n * padcp;
    break;
  case Dxya:
    Out_econst = 3 * l2 * m2 * dds + (l2 + m2 - 4 * l2 * m2) * ddp +
                 (n2 + l2 * m2) * ddd;
    break;
  case Dyza:
    Out_econst = 3 * l * m2 * n * dds + l * n * (1 - 4 * m2) * ddp +
                 l * n * (m2 - 1) * ddd;
    break;
  case Dzxa:
    Out_econst = 3 * l2 * m * n * dds + m * n * (1 - 4 * l2) * ddp +
                 m * n * (l2 - 1) * ddd;
    break;
  case Dx2my2a:
    Out_econst = 1.5 * l * m * (l2 - m2) * dds + 2 * l * m * (m2 - l2) * ddp +
                 0.5 * l * m * (l2 - m2) * ddd;
    break;
  case Dz2a:
    Out_econst = sqr3 * l * m * (n2 - 0.5 * (l2 + m2)) * dds -
                 sqr3 * 2 * l * m * n2 * ddp +
                 0.5 * sqr3 * l * m * (1 + n2) * ddd;
  }
  break;
case Dyzc:
  switch (basis2_econst) {
  case STa:
    Out_econst = sqr3 * m * n * stadcs;
    break;
  case Sa:
    Out_econst = sqr3 * m * n * sadcs;
    break;
  case Pxa:
    Out_econst = sqr3 * l * m * n * padcs - 2 * l * m * n * padcp;
    break;
  case Pya:
    Out_econst = sqr3 * n * m2 * padcs + n * (1 - 2 * m2) * padcp;
    break;
  case Pza:
    Out_econst = sqr3 * m * n2 * padcs + m * (1 - 2 * n2) * padcp;
    break;
  case Dxya:
    Out_econst = 3 * l * m2 * n * dds + l * n * (1 - 4 * m2) * ddp +
                 l * n * (m2 - 1) * ddd;
    break;
  case Dyza:
    Out_econst = 3 * m2 * n2 * dds + (m2 + n2 - 4 * m2 * n2) * ddp +
                 (l2 + m2 * n2) * ddd;
    break;
  case Dzxa:
    Out_econst = 3 * m * n2 * l * dds + m * l * (1 - 4 * n2) * ddp +
                 m * l * (n2 - 1) * ddd;
    break;
  case Dx2my2a:
    Out_econst = 1.5 * m * n * (l2 - m2) * dds -
                 m * n * (1 + 2 * (l2 - m2)) * ddp +
                 m * n * (1 + 0.5 * (l2 - m2)) * ddd;
    break;
  case Dz2a:
    Out_econst = sqr3 * m * n * (n2 - 0.5 * (l2 + m2)) * dds +
                 sqr3 * m * n * (l2 + m2 - n2) * ddp -
                 0.5 * sqr3 * m * n * (l2 + m2) * ddd;
  }
  break;
case Dzxc:
  switch (basis2_econst) {
  case STa:
    Out_econst = sqr3 * n * l * stadcs;
    break;
  case Sa:
    Out_econst = sqr3 * n * l * sadcs;
    break;
  case Pxa:
    Out_econst = sqr3 * l2 * n * padcs + n * (1 - 2 * l2) * padcp;
    break;
  case Pya:
    Out_econst = sqr3 * l * m * n * padcs - 2 * l * m * n * padcp;
    break;
  case Pza:
    Out_econst = sqr3 * l * n2 * padcs + l * (1 - 2 * n2) * padcp;
    break;
  case Dxya:
    Out_econst = 3 * l2 * m * n * dds + m * n * (1 - 4 * l2) * ddp +
                 m * n * (l2 - 1) * ddd;
    break;
  case Dyza:
    Out_econst = 3 * n2 * l * m * dds + l * m * (1 - 4 * n2) * ddp +
                 l * m * (n2 - 1) * ddd;
    break;
  case Dzxa:
    Out_econst = 3 * n2 * l2 * dds + (n2 + l2 - 4 * n2 * l2) * ddp +
                 (m2 + n2 * l2) * ddd;
    break;
  case Dx2my2a:
    Out_econst = 1.5 * n * l * (l2 - m2) * dds +
                 n * l * (1 - 2 * (l2 - m2)) * ddp -
                 n * l * (1 - 0.5 * (l2 - m2)) * ddd;
    break;
  case Dz2a:
    Out_econst = sqr3 * l * n * (n2 - 0.5 * (l2 + m2)) * dds +
                 sqr3 * l * n * (l2 + m2 - n2) * ddp -
                 0.5 * sqr3 * l * n * (l2 + m2) * ddd;
  }
  break;
case Dx2my2c:
  switch (basis2_econst) {
  case STa:
    Out_econst = 0.5 * sqr3 * (l2 - m2) * stadcs;
    break;
  case Sa:
    Out_econst = 0.5 * sqr3 * (l2 - m2) * sadcs;
    break;
  case Pxa:
    Out_econst = 0.5 * sqr3 * l * (l2 - m2) * padcs + l * (1 - l2 + m2) * padcp;
    break;
  case Pya:
    Out_econst = 0.5 * sqr3 * m * (l2 - m2) * padcs - m * (1 + l2 - m2) * padcp;
    break;
  case Pza:
    Out_econst = 0.5 * sqr3 * n * (l2 - m2) * padcs - n * (l2 - m2) * padcp;
    break;
  case Dxya:
    Out_econst = 1.5 * l * m * (l2 - m2) * dds + 2 * l * m * (m2 - l2) * ddp +
                 0.5 * l * m * (l2 - m2) * ddd;
    break;
  case Dyza:
    Out_econst = 1.5 * m * n * (l2 - m2) * dds -
                 m * n * (1 + 2 * (l2 - m2)) * ddp +
                 m * n * (1 + 0.5 * (l2 - m2)) * ddd;
    break;
  case Dzxa:
    Out_econst = 1.5 * n * l * (l2 - m2) * dds +
                 n * l * (1 - 2 * (l2 - m2)) * ddp -
                 n * l * (1 - 0.5 * (l2 - m2)) * ddd;
    break;
  case Dx2my2a:
    Out_econst = 0.75 * (l2 - m2) * (l2 - m2) * dds +
                 (l2 + m2 - (l2 - m2) * (l2 - m2)) * ddp +
                 (n2 + 0.25 * (l2 - m2) * (l2 - m2)) * ddd;
    break;
  case Dz2a:
    Out_econst = 0.5 * sqr3 * (l2 - m2) * (n2 - 0.5 * (l2 + m2)) * dds +
                 sqr3 * n2 * (m2 - l2) * ddp +
                 0.25 * sqr3 * (1 + n2) * (l2 - m2) * ddd;
  }
  break;
case Dz2c:
  switch (basis2_econst) {
  case STa:
    Out_econst = (n2 - 0.5 * (l2 + m2)) * stadcs;
    break;
  case Sa:
    Out_econst = (n2 - 0.5 * (l2 + m2)) * sadcs;
    break;
  case Pxa:
    Out_econst = l * (n2 - 0.5 * (l2 + m2)) * padcs - sqr3 * l * n2 * padcp;
    break;
  case Pya:
    Out_econst = m * (n2 - 0.5 * (l2 + m2)) * padcs - sqr3 * m * n2 * padcp;
    break;
  case Pza:
    Out_econst =
        n * (n2 - 0.5 * (l2 + m2)) * padcs + sqr3 * n * (l2 + m2) * padcp;
    break;
  case Dxya:
    Out_econst = sqr3 * l * m * (n2 - 0.5 * (l2 + m2)) * dds -
                 sqr3 * 2 * l * m * n2 * ddp +
                 0.5 * sqr3 * l * m * (1 + n2) * ddd;
    break;
  case Dyza:
    Out_econst = sqr3 * m * n * (n2 - 0.5 * (l2 + m2)) * dds +
                 sqr3 * m * n * (l2 + m2 - n2) * ddp -
                 0.5 * sqr3 * m * n * (l2 + m2) * ddd;
    break;
  case Dzxa:
    Out_econst = sqr3 * l * n * (n2 - 0.5 * (l2 + m2)) * dds +
                 sqr3 * l * n * (l2 + m2 - n2) * ddp -
                 0.5 * sqr3 * l * n * (l2 + m2) * ddd;
    break;
  case Dx2my2a:
    Out_econst = 0.5 * sqr3 * (l2 - m2) * (n2 - 0.5 * (l2 + m2)) * dds +
                 sqr3 * n2 * (m2 - l2) * ddp +
                 0.25 * sqr3 * (1 + n2) * (l2 - m2) * ddd;
    break;
  case Dz2a:
    Out_econst = (n2 - 0.5 * (l2 + m2)) * (n2 - 0.5 * (l2 + m2)) * dds +
                 3 * n2 * (l2 + m2) * ddp + 0.75 * (l2 + m2) * (l2 + m2) * ddd;
    break;
  }
}
