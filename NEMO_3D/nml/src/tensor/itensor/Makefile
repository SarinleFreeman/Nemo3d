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
# $Header: /repo/nml/src/tensor/itensor/Makefile,v 1.3 2004/10/14 16:02:41 marek Exp $
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
ARCHIVE=$(CATALOG)/libnml_itensor.a
HEADER0=$(INCPATH)/nml_itensor.h
HEADER1=$(INCPATH)/nml_i3DBase.h
HEADERS=$(HEADER0) $(HEADER1)
COMPILE=$(CC) $(INCLUDE) $(DEFINES) $(OPTIONS)
GENCLASS=$(BINPATH)/genclass

SOURCES=nml_itensor.c nml_i3DBase.c

OBJECTS=$(PLATFORM)/nml_itensor.o \
	$(PLATFORM)/nml_i3DBase.o $(PLATFORM)/i3tensor.o

all:
	(cd ../../matrix/imatrix; $(MAKE) \
		CC='$(CC)' DEFINES='$(DEFINES)' OPTIONS='$(OPTIONS)' \
		AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)' PLATFORM='$(PLATFORM)' MAKE='$(MAKE)')
	($(MAKE) $(ARCHIVE) \
		CC='$(CC)' DEFINES='$(DEFINES)' OPTIONS='$(OPTIONS)' \
		AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)' PLATFORM='$(PLATFORM)' MAKE='$(MAKE)')

$(ARCHIVE):		$(CATALOG) $(OBJECTS)
	$(AR) $(AR_OPTIONS) $(ARCHIVE) $(OBJECTS)

$(CATALOG):
	mkdir -p $(CATALOG)

$(PLATFORM)/nml_itensor.o:	$(PLATFORM) $(HEADER0) nml_itensor.c
	$(COMPILE) -c nml_itensor.c
	mv nml_itensor.o $(PLATFORM)

$(HEADER0):		../tensor.hP
	$(GENCLASS) integer i tensor
	mv nml_itensor.h $(HEADER0)

nml_itensor.c:		../tensor.cP
	$(GENCLASS) integer i tensor
	mv nml_itensor.h $(HEADER0)

$(PLATFORM)/nml_i3DBase.o:	$(PLATFORM) $(HEADER1) nml_i3DBase.c
	$(COMPILE) -c nml_i3DBase.c
	mv nml_i3DBase.o $(PLATFORM)

$(PLATFORM)/i3tensor.o:		$(PLATFORM) $(INCPATH)/i3tensor.h i3tensor.c
	$(COMPILE) -c i3tensor.c
	mv i3tensor.o $(PLATFORM)

$(HEADER1):		../3DBase.hP
	$(GENCLASS) integer i 3DBase
	mv nml_i3DBase.h $(HEADER1)

nml_i3DBase.c:		../3DBase.cP
	$(GENCLASS) integer i 3DBase
	mv nml_i3DBase.h $(HEADER1)

$(PLATFORM):
	mkdir -p $(PLATFORM)

clean:
	rm -f $(OBJECTS)

sterile:
	rm -f $(ARCHIVE) $(OBJECTS) $(SOURCES) $(HEADERS)
	rm -rf $(PLATFORM)
