# *****************************************************************************
# The JPL Parallel Eigensolvers package.
# Copyright (C) 2002 California Institute of Technology (Caltech)
#
# This file is part of
# The Jet Propulsion Laboratory (JPL) Parallel Eigensolvers package.
#
# This library is free software, which you can redistribute and/or modify
# under the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation; either version 2.1 of the License, or (at your
# option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
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
#   E. Robert Tisdale (E.Robert.Tisdale@jpl.nasa.gov)

# Written by: Chris Bowen
#             Gerhard Klimeck
#             Fabiano Oyafuso
#             E. Robert Tisdale
#
# *****************************************************************************
# $Header: /repo/eigen/Makefile,v 1.6 2006/10/05 17:58:50 stevenclark Exp $
# *****************************************************************************


#================================================================
# Makefile to generate the eigsolver library
#================================================================

#----------------------------------------------------------------
# System configuration: compiler and libraries
#----------------------------------------------------------------

  MAKE	= make
  SHELL	= /bin/sh

  VERSION   =
  PLATFORM  = $(HOSTTYPE)

#----------------------------------------------------------------
#  Platform specific options
#----------------------------------------------------------------

# NML options
	AR         = /usr/bin/ar
	AR_OPTIONS = rcvs
# LAPACK options
#	TCS
 	FORTRAN	= f90
# eigen options
	CCMKD = g++ -MM
#	CCOPT = -Wall -O3 -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING
	CCOPT = -Wall -O3 -DFORTRAN_UNDERSCORE -DMPI3d
#	TCS C++
	CC    = g++
	LD    = g++

#----------------------------------------------------------------
# User rules
#----------------------------------------------------------------

all:
	(mkdir -p lib/$(PLATFORM);)
	(mkdir -p obj/$(PLATFORM);)
	(cd src; $(MAKE) PLATFORM='$(PLATFORM)' CCMKD='$(CCMKD)' \
	CC='$(CC)' CCOPT='$(CCOPT)' LD='$(LD)' AR='$(AR) -$(AR_OPTIONS)')

clean:
	(cd src; $(MAKE) PLATFORM='$(PLATFORM)' clean)
	(rm -rf obj/$(PLATFORM);)

cleanall: clean
	(rm -rf lib/$(PLATFORM);)

tar:
	tar cvf ../eigen.tar ../eigen/

