# *****************************************************************************
# The NEMO Math Library.
# Copyright (C) 2002-2004 California Institute of Technology (Caltech)
#
# This file is part of
# The NanoElectrictronic MOdeling (NEMO) Math Library.
#
# This library is free software which you can redistribute and/or modify
# under the terms of the GNU Library General Public License
# as published by the Free Software Foundation;
# either version 2, or (at your option) any later version.
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
# Written by E. Robert Tisdale
# Modified by Marek J. Korkusinski, October, 2004
# *****************************************************************************
# $Header: /repo/nml/examples/Makefile,v 1.6 2004/10/21 21:15:34 marek Exp $
# *****************************************************************************


# Set the default C or C++ compiler.  Note that
# the C++ compiler "mangles" function names to create the symbols
# that it leaves behind in object files for the link editor.
# Sequential version:
# CC=gcc
# CC=g++
# Parallel version:
# CC=$(MPI_HOME)/bin/mpicc
#  CC=$(MPI_HOME)/bin/mpiCC
ifneq ($(MPI_HOME),)
	CC = $(MPI_HOME)/bin/mpiCC
	LD = $(MPI_HOME)/bin/mpiCC
else
	CC = g++
	LD = g++
endif
INCLUDE=-I../include
#DEFINES=-DNML_MEM_EXT
# DEFINES=-DNML_INLINE
# DEFINES=-DNML_MACROS
OPTIONS=-O2 -Wall -ansi -pedantic
#OPTIONS= -O2 -Wall -felide-constructors
PLATFORM=$(HOSTTYPE)
LIBRARY=-L../lib/$(PLATFORM)
libraries= -lnml

COMPILE=$(CC) $(INCLUDE) $(LIBRARY) $(DEFINES) $(OPTIONS)

 example:	example.c
	$(COMPILE) -o example example.c $(libraries)

r2example:	r2example.c
	$(COMPILE) -o r2example r2example.c $(libraries)

s2example:	s2example.c
	$(COMPILE) -o s2example s2example.c $(libraries)

f2example:	f2example.c
	$(COMPILE) -o f2example f2example.c $(libraries)

b2example:	b2example.c
	$(COMPILE) -o b2example b2example.c $(libraries)

rmexample:	rmexample.c
	$(COMPILE) -o rmexample rmexample.c $(libraries)

cmexample:	cmexample.c
	$(COMPILE) -o cmexample cmexample.c $(libraries)

kmexample:	kmexample.c
	$(COMPILE) -o kmexample kmexample.c $(libraries)

   clean:
	rm -f example example.o r2example r2example.o \
	rmexample rmexample.o cmexample cmexample.o s2example.o s2example \
	b2example.o b2example f2example.o f2example kmexample kmexample.o

