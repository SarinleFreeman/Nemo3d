# *****************************************************************************
# The NEMO Math Library.
# Copyright (C) 2002 California Institute of Technology (Caltech)
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
#
# *****************************************************************************
# $Header: /repo/nml/src/matrix/dmatrix/Makefile,v 1.3 2004/10/14 16:02:41 marek Exp $
# *****************************************************************************

MAKE = make
INCPATH=../../../include
LIBPATH=../../../lib
BINPATH=../../../bin
CC=gcc
INCLUDE=-I$(INCPATH)
DEFINES=-DNCL_REVEAL -DNML_MEM_EXT
OPTIONS=-O2 -Wall -ansi -pedantic -felide-constructors
AR=ar
AR_OPTIONS=rcvs
PLATFORM=$(HOSTTYPE)
CATALOG=$(LIBPATH)/$(PLATFORM)
ARCHIVE=$(CATALOG)/libnml_dmatrix.a
HEADER0=$(INCPATH)/nml_dmatrix.h
HEADER1=$(INCPATH)/nml_d2DBase.h
HEADER2=$(INCPATH)/nml_d3bands.h
HEADERS=$(HEADER0) $(HEADER1) $(HEADER2)
COMPILE=$(CC) $(INCLUDE) $(DEFINES) $(OPTIONS)
GENCLASS=$(BINPATH)/genclass

SOURCES=nml_dmatrix.c nml_d2DBase.c nml_d3bands.c

OBJECTS=$(PLATFORM)/nml_dmatrix.o $(PLATFORM)/rmatrix.o \
	$(PLATFORM)/nml_d2DBase.o $(PLATFORM)/r2tensor.o \
	$(PLATFORM)/nml_d3bands.o $(PLATFORM)/rtridiag.o

all:
	(cd ../../vector/dvector; $(MAKE) \
		CC='$(CC)' DEFINES='$(DEFINES)' OPTIONS='$(OPTIONS)' \
		AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)' PLATFORM='$(PLATFORM)' MAKE='$(MAKE)')
	($(MAKE) $(ARCHIVE) \
		CC='$(CC)' DEFINES='$(DEFINES)' OPTIONS='$(OPTIONS)' \
		AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)' PLATFORM='$(PLATFORM)' MAKE='$(MAKE)')

$(ARCHIVE):		$(CATALOG) $(OBJECTS)
	$(AR) $(AR_OPTIONS) $(ARCHIVE) $(OBJECTS)

$(CATALOG):
	mkdir -p $(CATALOG)

$(PLATFORM)/nml_dmatrix.o:	$(PLATFORM) $(HEADER0) nml_dmatrix.c
	$(COMPILE) -c nml_dmatrix.c
	mv nml_dmatrix.o $(PLATFORM)

$(PLATFORM)/rmatrix.o:		$(PLATFORM) $(INCPATH)/rmatrix.h rmatrix.c
	$(COMPILE) -c rmatrix.c
	mv rmatrix.o $(PLATFORM)

$(HEADER0):		../matrix.hP
	$(GENCLASS) double d matrix
	mv nml_dmatrix.h $(HEADER0)

nml_dmatrix.c:		../matrix.cP
	$(GENCLASS) double d matrix
	mv nml_dmatrix.h $(HEADER0)

$(PLATFORM)/nml_d2DBase.o:	$(PLATFORM) $(HEADER1) nml_d2DBase.c
	$(COMPILE) -c nml_d2DBase.c
	mv nml_d2DBase.o $(PLATFORM)

$(PLATFORM)/r2tensor.o:		$(PLATFORM) $(INCPATH)/r2tensor.h r2tensor.c
	$(COMPILE) -c r2tensor.c
	mv r2tensor.o $(PLATFORM)

$(HEADER1):		../2DBase.hP
	$(GENCLASS) double d 2DBase
	mv nml_d2DBase.h $(HEADER1)

nml_d2DBase.c:		../2DBase.cP
	$(GENCLASS) double d 2DBase
	mv nml_d2DBase.h $(HEADER1)

$(PLATFORM)/nml_d3bands.o:	$(PLATFORM) $(HEADER2) nml_d3bands.c
	$(COMPILE) -c nml_d3bands.c
	mv nml_d3bands.o $(PLATFORM)

$(PLATFORM)/rtridiag.o:		$(PLATFORM) $(INCPATH)/rtridiag.h rtridiag.c
	$(COMPILE) -c rtridiag.c
	mv rtridiag.o $(PLATFORM)

$(HEADER2):		../3bands.hP
	$(GENCLASS) double d 3bands
	mv nml_d3bands.h $(HEADER2)

nml_d3bands.c:		../3bands.cP
	$(GENCLASS) double d 3bands
	mv nml_d3bands.h $(HEADER2)

$(PLATFORM):
	mkdir -p $(PLATFORM)

clean:
	rm -f $(OBJECTS)

sterile:
	rm -f $(ARCHIVE) $(OBJECTS) $(SOURCES) $(HEADERS)
	rm -rf $(PLATFORM)

