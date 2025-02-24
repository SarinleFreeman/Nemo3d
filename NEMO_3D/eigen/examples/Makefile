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
# $Header: /repo/eigen/examples/Makefile,v 1.2 2003/10/08 16:15:38 hook Exp $
# *****************************************************************************


  NML_HOME=../nml
  LAPACK_HOME=../LAPACK
# Parallel C or C++ compiler.
# PCC=$(MPI_HOME)/bin/mpicc
  PCC=$(MPI_HOME)/bin/mpiCC
# COPTIONS=-O2 -ansi -pedantic -Wall
  COPTIONS=-Wall -O3
# C preprocessor macro definitions
# LAPACK compiled with ABSOFT f77 -f -B108
#		  or with GNU f77 -Wall -funroll-all-loops -fno-f2c -O3
# DEFINES=-DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING
  DEFINES=-DFORTRAN_UNDERSCORE -DMPI3d
# LAPACK compiled with ABSOFT f77 -N109
# DEFINES=-DFORTRAN_ALLCAPS -DMPI3d -DMPI_TIMING
# DEFINES=-DFORTRAN_ALLCAPS -DMPI3d
# 
# Define the NML_INLINE macro if your C compiler uses inline functions.
# DEFINES+=-DNML_INLINE
# Define the NML_MACROS macro if your C compiler can't accept inline functions.
# DEFINES+=-DNML_MACROS

  INCLUDE=-I../include -I$(NML_HOME)/include
  PLATFORM=$(HOSTTYPE)
# LAPACK compiled with GNU    f77 compiler
  LIBRARY=-L../lib/$(PLATFORM) -L$(NML_HOME)/lib/$(PLATFORM)
# LAPACK compiled with ABSOFT f77 compiler
# LIBRARY=-L../lib/$(PLATFORM) -L$(NML_HOME)/lib/$(PLATFORM) -L$(ABSOFT)/lib
  COMPILE=$(PCC) $(DEFINES) $(INCLUDE) $(LIBRARY) $(COPTIONS)

# Parallel Fortran compiler.
  PFC=$(MPI_HOME)/bin/mpif77
# Use the same Fortran options used to compile LAPACK
# LAPACK compiled with GNU    f77 compiler
  FOPTIONS=-Wall -funroll-all-loops -fno-f2c -O3
# LAPACK compiled with ABSOFT f77 compiler
# FOPTIONS=-O -N109	# upper case
# FOPTIONS=-O -f -B108	# lower case underscore
  FORTRAN=$(PFC) $(FOPTIONS)
# LAPACK compiled with GNU    f77 compiler
  library=-leigen -lnml -llapack -lblas -lg2c
# LAPACK compiled with ABSOFT f77 compiler
# library=-leigen -lnml -lm -llapack -lblas -lfio -lf77math -lf90math
  SOURCES=pexample.c
  HEADERS=
  OBJECTS=
  ARCHIVE=

pexample:	pexample.c pzsmvmul.o
	$(COMPILE) -o pexample pexample.c pzsmvmul.o $(library)

pzsmvmul.o:	pzsmvmul.f
	$(FORTRAN) -c pzsmvmul.f

clean: 
	/bin/rm -f pexample pexample.o pzsmvmul.o

