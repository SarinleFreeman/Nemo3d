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
# $Header: /repo/nemo3d/misc/reth.pl,v 1.2 2003/10/08 16:16:05 hook Exp $
# *****************************************************************************

# -------------------------------------------------------
#    run.pm
# -------------------------------------------------------
#   FORMAT:  run.pm [-r <template>.in | -np INT ] 
#
#   OPTIONS:    -r     default run is "alloy"
#               -np    number of processors
#               -x0    start concentration
#               -xf    end concentration
#               -dx    conentration increment
#               +/- bs do bisection (-, default)
#               +/- bl do bondlength (-, default)
# -------------------------------------------------------

use Sys::Hostname;

sub printHelp {
  my ($arg) = @_;

  system("head -n14 `which run.pm` | tail -n13");
  die("") if ($arg eq "error");
  exit 0 if ($arg eq "help");
}

# default values
$exec = "./run3d";
$fil_in = "rnd.in";
$fil_out = "tmp.in";
$fmt_ext = "%4.3f";
$Nproc = 8;
#$Nproc = 1;
$eps = 1e-6;
$x0 = 0.0;
$xf = 1.0 + $eps;
$dx = .05;
$vary_conc = 0;
$vary_seed = 1;
$do_bondlength = 0;
$do_bisect = 0;

# for random alloy distribution calculation
($host) = hostname() =~ m/pluto0*(\S+).jpl.nasa.gov/;
$seed0 = $host*5000 + 1000;
$seedf = ($host+1)*5000 - 3001;

# overide defaults from command line
for ($i=0; $i<=$#ARGV; $i++) {
  if ($ARGV[$i] eq "-h") {
    printHelp("help");
  }
  elsif ($ARGV[$i] eq "-r") {
    $fil_in = $ARGV[++$i] . ".in";
  }
  elsif ($ARGV[$i] eq "-np") {
    $Nproc = $ARGV[++$i];
  }
  elsif ($ARGV[$i] eq "-x0") {
    $x0 = $ARGV[++$i];
  }
  elsif ($ARGV[$i] eq "-xf") {
    $xf = $ARGV[++$i] + $eps;
  }
  elsif ($ARGV[$i] eq "-dx") {
    $dx = $ARGV[++$i];
  }
  elsif ($ARGV[$i] eq "-e") {
    $cmd = $ARGV[++$i];
    print STDOUT "Evaluating command:  " . $cmd . "\n";
    eval $cmd;
  }
  # boolean options
  elsif ($ARGV[$i] eq "+bs") {
    $do_bisect = 1;
  }
  elsif ($ARGV[$i] eq "-bs") {
    $do_bisect = 0;
  }
  elsif ($ARGV[$i] eq "+bl") {
    $do_bondlength = 1;
  }
  elsif ($ARGV[$i] eq "-bl") {
    $do_bondlength = 0;
  }
  else {
    printHelp("error");
  }
}

# fill concentration array
@conc=();
if ($do_bisect) {
  push @conc, $x0;
  push @conc, $xf unless ($x0==$xf);
  for ($n=1, $delx=0.5*($xf-$x0); $delx > 0.5*$dx; $n++, $delx*=.5) {
    for ($j=1; $j<2**$n; $j+=2) {
      push @conc, $x0 + $j*$delx;
    }
  } 
} else {
  for ($x=$x0; $x <= $xf; $x += $dx) {
    push @conc, $x;
  }
}

open(FP, $fil_in) || die("Could not open file:  ", $fil_in);
while (defined($line = <FP>)) {
  ($y) = $line =~ m/^ResFind\s+=\s+(\S+)/;
  $lanczos = (($y eq "Iterative") ? 1 : 0) if ($y ne "");
  ($y) = $line =~ m/^AlloyTreatment\s+=\s+(\S+)/;
  $is_vca_alloy = (($y eq "VCA") ? 1 : 0) if ($y ne "");
}
close(FP);


if ($vary_seed) {
  for ($i=$seed0; $i<=$seedf; $i++) {

    $x_fmt = sprintf("%d",$i);

    open(FP, $fil_in) || die("Could not open file:  ", $fil_in);
    open(FP_OUT, ">".$fil_out);
    while (defined($line = <FP>)) {
      ($y) = $line =~ m/^RandomSeed/;
      if ($y!="") {
	print FP_OUT sprintf("RandomSeed = %d \n", $i);
      }
      else {
	print FP_OUT $line;
      }
    }
    close(FP_OUT);
    close(FP);
    
    system("mpirun -nolocal -machinefile mach -np " . $Nproc . " " . $exec . " " . $fil_out);
#    system("mpirun.ch_gm --gm-f gmfile -np " . $Nproc . " " . $exec . " " . $fil_out);
#    system("./run.py");
    
#    while(<tmp.nd_evec_*>) {
#      system("scp",$_,"alhena:/silo/fabiano/" . $_ . "-" . $x_fmt);
#    }

    while(<tmp.nd_evec_*>) {
      system("/bin/rm -f", $_);
    }

#    system("scp", "species_count.out", "alhena:/silo/fabiano/species_count-" . $x_fmt);
#    system("scp", "tmp.nd_Ek", "alhena:/silo/fabiano/Egamma-" . $x_fmt);

  }
}
