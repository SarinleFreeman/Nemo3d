#!/usr/bin/env perl

# *****************************************************************************
# The Jet Propulsion Laboratory (JPL) NanoElectronicMOdeling-3D package.
# Copyright (C) 2002 California Institute of Technology (Caltech)
#
# This application is free software, which you can redistribute and/or modify
# under the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation; either version 2.1 of the License, or (at your
# option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library; see the file COPYING. If not, write to the
# Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA  02111-1307  USA
#
# For additional information, please contact
#   Gerhard Klimeck (gekco@jpl.nasa.gov)
#   Fabiano Oyafuso (fabiano@jpl.nasa.gov)
#
# Written by:  Chris Bowen
#              Gerhard Klimeck
#              Fabiano Oyafuso
#              Seungwon Lee
#              Olga Lazarenkova
#              Hook Hua
#
# This product includes software developed by the Apache Software Foundation
# (http://www.apache.org/).
#
# *****************************************************************************
# $Header: /repo/nemo3d/misc/xtr.pl,v 1.2 2003/10/08 16:16:07 hook Exp $
# *****************************************************************************

# -------------------------------------------------------
#    xtract.pm
# -------------------------------------------------------
#   FORMAT:  xtract.pm 
#
#   OPTIONS:    +/-g  print out energy gap
#               +/-l  print out lattice constant
# -------------------------------------------------------

sub suffix {
  my ($fl) = @_;

  $prfx = $fl . "-";

  @sfx = ();
  while (defined($line = <$prfx*>)) {
    ($x) = $line =~ m/$prfx(\S+)/;
    push @sfx, $x;
  }
  @sfx;
}

sub suffixIntersection {
  my ($f1) = @_;
  shift;
  my ($f2) = @_;

  @x1=suffix($f1);
  @x2=suffix($f2);

  %x1=map{$_ =>1} @x1;
  %x2=map{$_ =>1} @x2;
  @x1_and_x2 = grep( $x1{$_}, @x2 );

  print $f1 . "-$_\n" foreach (@x1);
  print $f2 . "-$_\n" foreach (@x2);

  @x1_and_x2;
}

@y=suffixIntersection("crud", "crap");

print "intersection:\t$_\n" foreach (@y);

