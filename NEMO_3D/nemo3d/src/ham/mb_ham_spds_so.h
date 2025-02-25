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
$Header: /repo/nemo3d/src/ham/mb_ham_spds_so.h,v 1.2 2003/10/08 16:17:28 hook
Exp $
*****************************************************************************/

/*
 * spin_orbit_spds replacement for inlining
 */

/* complex spin_orbit_spds2( Out_e_so, e_so_basis1,  e_so_basis2,  e_so_spin1,
 * e_so_spin2,  e_so_delta_a,  e_so_delta_c ) */
{

  Out_e_so.r = Out_e_so.i = 0.0;

  if (e_so_spin1 == 1 && e_so_spin2 == 0) {
    switch (e_so_basis1) {
    case Pxa:
      switch (e_so_basis2) {
      case Pza:
        Out_e_so.r = e_so_delta_a;
        break;
      }
      break;
    case Pya:
      switch (e_so_basis2) {
      case Pza:
        Out_e_so.i = -e_so_delta_a;
        break;
      }
      break;
    case Pza:
      switch (e_so_basis2) {
      case Pxa:
        Out_e_so.r = -e_so_delta_a;
        break;
      case Pya:
        Out_e_so.i = e_so_delta_a;
        break;
      }
      break;
    case Pxc:
      switch (e_so_basis2) {
      case Pzc:
        Out_e_so.r = e_so_delta_c;
        break;
      }
      break;
    case Pyc:
      switch (e_so_basis2) {
      case Pzc:
        Out_e_so.i = -e_so_delta_c;
        break;
      }
      break;
    case Pzc:
      switch (e_so_basis2) {
      case Pxc:
        Out_e_so.r = -e_so_delta_c;
        break;
      case Pyc:
        Out_e_so.i = e_so_delta_c;
        break;
      }
    }
  }
  if (e_so_spin1 == 0 && e_so_spin2 == 1) {
    switch (e_so_basis1) {
    case Pxa:
      switch (e_so_basis2) {
      case Pza:
        Out_e_so.r = -e_so_delta_a;
        break;
      }
      break;
    case Pya:
      switch (e_so_basis2) {
      case Pza:
        Out_e_so.i = -e_so_delta_a;
        break;
      }
      break;
    case Pza:
      switch (e_so_basis2) {
      case Pxa:
        Out_e_so.r = e_so_delta_a;
        break;
      case Pya:
        Out_e_so.i = e_so_delta_a;
        break;
      }
      break;
    case Pxc:
      switch (e_so_basis2) {
      case Pzc:
        Out_e_so.r = -e_so_delta_c;
        break;
      }
      break;
    case Pyc:
      switch (e_so_basis2) {
      case Pzc:
        Out_e_so.i = -e_so_delta_c;
        break;
      }
      break;
    case Pzc:
      switch (e_so_basis2) {
      case Pxc:
        Out_e_so.r = e_so_delta_c;
        break;
      case Pyc:
        Out_e_so.i = e_so_delta_c;
        break;
      }
    }
  }
  if (e_so_spin1 == 1 && e_so_spin2 == 1) {
    switch (e_so_basis1) {
    case Pxa:
      switch (e_so_basis2) {
      case Pya:
        Out_e_so.i = -e_so_delta_a;
        break;
      }
      break;
    case Pya:
      switch (e_so_basis2) {
      case Pxa:
        Out_e_so.i = e_so_delta_a;
        break;
      }
      break;
    case Pxc:
      switch (e_so_basis2) {
      case Pyc:
        Out_e_so.i = -e_so_delta_c;
        break;
      }
      break;
    case Pyc:
      switch (e_so_basis2) {
      case Pxc:
        Out_e_so.i = e_so_delta_c;
        break;
      }
    }
  }
  if (e_so_spin1 == 0 && e_so_spin2 == 0) {
    switch (e_so_basis1) {
    case Pxa:
      switch (e_so_basis2) {
      case Pya:
        Out_e_so.i = e_so_delta_a;
        break;
      }
      break;
    case Pya:
      switch (e_so_basis2) {
      case Pxa:
        Out_e_so.i = -e_so_delta_a;
        break;
      }
      break;
    case Pxc:
      switch (e_so_basis2) {
      case Pyc:
        Out_e_so.i = e_so_delta_c;
        break;
      }
      break;
    case Pyc:
      switch (e_so_basis2) {
      case Pxc:
        Out_e_so.i = -e_so_delta_c;
        break;
      }
    }
  }
}
