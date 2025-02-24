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
# $Header: /repo/nml/src/vector/ivector/Makefile,v 1.3 2004/10/14 16:02:41 marek Exp $
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
ARCHIVE=$(CATALOG)/libnml_ivector.a
HEADER0=$(INCPATH)/nml_ivector.h
HEADER1=$(INCPATH)/nml_i1DBase.h
HEADERS=$(HEADER0) $(HEADER1)
COMPILE=$(CC) $(INCLUDE) $(DEFINES) $(OPTIONS)
GENCLASS=$(BINPATH)/genclass

SOURCES=nml_ivector.c nml_i1DBase.c

OBJECTS=$(PLATFORM)/nml_ivector.o $(PLATFORM)/ivector.o \
	$(PLATFORM)/nml_i1DBase.o

all:
	(cd ../../scalar/iscalar; $(MAKE) \
		CC='$(CC)' DEFINES='$(DEFINES)' OPTIONS='$(OPTIONS)' \
		AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)' PLATFORM='$(PLATFORM)' MAKE='$(MAKE)')
	($(MAKE) $(ARCHIVE) \
		CC='$(CC)' DEFINES='$(DEFINES)' OPTIONS='$(OPTIONS)' \
		AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)' PLATFORM='$(PLATFORM)' MAKE='$(MAKE)')

$(ARCHIVE):		$(CATALOG) $(OBJECTS)
	$(AR) $(AR_OPTIONS) $(ARCHIVE) $(OBJECTS)

$(CATALOG):
	mkdir -p $(CATALOG)

$(PLATFORM)/nml_ivector.o:	$(PLATFORM) $(HEADER0) nml_ivector.c
	$(COMPILE) -c nml_ivector.c
	mv nml_ivector.o $(PLATFORM)

$(PLATFORM)/ivector.o:		$(PLATFORM) $(INCPATH)/ivector.h ivector.c
	$(COMPILE) -c ivector.c
	mv ivector.o $(PLATFORM)

$(HEADER0):		../vector.hP
	$(GENCLASS) integer i vector
	mv nml_ivector.h $(HEADER0)

nml_ivector.c:		../vector.cP
	$(GENCLASS) integer i vector
	mv nml_ivector.h $(HEADER0)

$(PLATFORM)/nml_i1DBase.o:	$(PLATFORM) $(HEADER1) nml_i1DBase.c
	$(COMPILE) -c nml_i1DBase.c
	mv nml_i1DBase.o $(PLATFORM)

$(PLATFORM)/i1tensor.o:		$(PLATFORM) $(INCPATH)/i1tensor.h i1tensor.c
	$(COMPILE) -c i1tensor.c
	mv i1tensor.o $(PLATFORM)

$(HEADER1):		../1DBase.hP
	$(GENCLASS) integer i 1DBase
	mv nml_i1DBase.h $(HEADER1)

nml_i1DBase.c:		../1DBase.cP
	$(GENCLASS) integer i 1DBase
	mv nml_i1DBase.h $(HEADER1)

$(PLATFORM):
	mkdir -p $(PLATFORM)

clean:
	rm -f $(OBJECTS)

sterile:
	rm -f $(ARCHIVE) $(OBJECTS) $(SOURCES) $(HEADERS)
	rm -rf $(PLATFORM)

