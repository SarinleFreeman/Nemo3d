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
# $Header: /repo/nml/src/tensor/dtensor/Makefile,v 1.4 2005/04/09 18:04:55 marek Exp $
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
ARCHIVE=$(CATALOG)/libnml_dtensor.a
HEADER0=$(INCPATH)/nml_dtensor.h
HEADER1=$(INCPATH)/nml_d3DBase.h
HEADERS=$(HEADER0) $(HEADER1)
COMPILE=$(CC) $(INCLUDE) $(DEFINES) $(OPTIONS)
GENCLASS=$(BINPATH)/genclass

SOURCES=nml_dtensor.c nml_d3DBase.c

OBJECTS=$(PLATFORM)/nml_dtensor.o \
	$(PLATFORM)/nml_d3DBase.o $(PLATFORM)/r3tensor.o

all:
	(cd ../../matrix/dmatrix; $(MAKE) \
		CC='$(CC)' DEFINES='$(DEFINES)' OPTIONS='$(OPTIONS)' \
		AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)' PLATFORM='$(PLATFORM)' MAKE='$(MAKE)')
	($(MAKE) $(ARCHIVE) \
		CC='$(CC)' DEFINES='$(DEFINES)' OPTIONS='$(OPTIONS)' \
		AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)' PLATFORM='$(PLATFORM)' MAKE='$(MAKE)')

$(ARCHIVE):		$(CATALOG) $(OBJECTS)
	$(AR) $(AR_OPTIONS) $(ARCHIVE) $(OBJECTS)

$(CATALOG):
	mkdir -p $(CATALOG)

$(PLATFORM)/nml_dtensor.o:	$(PLATFORM) $(HEADER0) nml_dtensor.c
	$(COMPILE) -c nml_dtensor.c
	mv nml_dtensor.o $(PLATFORM)

$(HEADER0):		../tensor.hP
	$(GENCLASS) double d tensor
	mv nml_dtensor.h $(HEADER0)

nml_dtensor.c:		../tensor.cP
	$(GENCLASS) double d tensor
	mv nml_dtensor.h $(HEADER0)

$(PLATFORM)/nml_d3DBase.o:	$(PLATFORM) $(HEADER1) nml_d3DBase.c
	$(COMPILE) -c nml_d3DBase.c
	mv nml_d3DBase.o $(PLATFORM)

$(PLATFORM)/r3tensor.o:		$(PLATFORM) $(INCPATH)/r3tensor.h r3tensor.c
	$(COMPILE) -c r3tensor.c
	mv r3tensor.o $(PLATFORM)

$(HEADER1):		../3DBase.hP
	$(GENCLASS) double d 3DBase
	mv nml_d3DBase.h $(HEADER1)

nml_d3DBase.c:		../3DBase.cP
	$(GENCLASS) double d 3DBase
	mv nml_d3DBase.h $(HEADER1)

$(PLATFORM):
	mkdir -p $(PLATFORM)

clean:
	rm -f $(OBJECTS)

sterile:
	rm -f $(ARCHIVE) $(OBJECTS) $(SOURCES) $(HEADERS)
	rm -rf $(PLATFORM)

