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
# $Header: /repo/nemo3d/misc/stdout2species.pl,v 1.2 2003/10/08 16:16:06 hook Exp $
# *****************************************************************************

# -------------------------------------------------------
#    stdout2species.pm
# -------------------------------------------------------
#   FORMAT:  stdout2species.pm  [filename] 
#

die("must supply exactly one argument") if ($#ARGV != 0);

$flag=0;
$fil_in = $ARGV[$1];

open(FP, $fil_in) || die("Could not open file:  ", $fil_in);
while (!$flag & defined($line = <FP>)) {
  ($flag) = $line =~ m/^Start processing an inner shell/;
  
}

# found string; now search for atom counts
while ($flag & defined($line = <FP>)) {
  ($N, $atom) = $line =~ m/^(\S+) atoms of type (\S+)/;
  print STDOUT $atom, "  ", $N , "\n" if $N;
}

close(FP);

