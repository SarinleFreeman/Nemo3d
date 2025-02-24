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
# $Header: /repo/nml/src/vector/dvector/Makefile,v 1.3 2004/10/14 16:02:41 marek Exp $
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
ARCHIVE=$(CATALOG)/libnml_dvector.a
HEADER0=$(INCPATH)/nml_dvector.h
HEADER1=$(INCPATH)/nml_d1DBase.h
HEADERS=$(HEADER0) $(HEADER1)
COMPILE=$(CC) $(INCLUDE) $(DEFINES) $(OPTIONS)
GENCLASS=$(BINPATH)/genclass

SOURCES=nml_dvector.c nml_d1DBase.c

OBJECTS=$(PLATFORM)/nml_dvector.o $(PLATFORM)/rvector.o \
	$(PLATFORM)/nml_d1DBase.o

all:
	(cd ../../scalar/dscalar; $(MAKE) \
		CC='$(CC)' DEFINES='$(DEFINES)' OPTIONS='$(OPTIONS)' \
		AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)' PLATFORM='$(PLATFORM)' MAKE='$(MAKE)')
	($(MAKE) $(ARCHIVE) \
		CC='$(CC)' DEFINES='$(DEFINES)' OPTIONS='$(OPTIONS)' \
		AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)' PLATFORM='$(PLATFORM)' MAKE='$(MAKE)')

$(ARCHIVE):		$(CATALOG) $(OBJECTS)
	$(AR) $(AR_OPTIONS) $(ARCHIVE) $(OBJECTS)

$(CATALOG):
	mkdir -p $(CATALOG)

$(PLATFORM)/nml_dvector.o:	$(PLATFORM) $(HEADER0) nml_dvector.c
	$(COMPILE) -c nml_dvector.c
	mv nml_dvector.o $(PLATFORM)

$(PLATFORM)/rvector.o:		$(PLATFORM) $(INCPATH)/rvector.h rvector.c
	$(COMPILE) -c rvector.c
	mv rvector.o $(PLATFORM)

$(HEADER0):		../vector.hP
	$(GENCLASS) double d vector
	mv nml_dvector.h $(HEADER0)

nml_dvector.c:		../vector.cP
	$(GENCLASS) double d vector
	mv nml_dvector.h $(HEADER0)

$(PLATFORM)/nml_d1DBase.o:	$(PLATFORM) $(HEADER1) nml_d1DBase.c
	$(COMPILE) -c nml_d1DBase.c
	mv nml_d1DBase.o $(PLATFORM)

$(PLATFORM)/r2tensor.o:		$(PLATFORM) $(INCPATH)/r2tensor.h r2tensor.c
	$(COMPILE) -c r2tensor.c
	mv r2tensor.o $(PLATFORM)

$(HEADER1):		../1DBase.hP
	$(GENCLASS) double d 1DBase
	mv nml_d1DBase.h $(HEADER1)

nml_d1DBase.c:		../1DBase.cP
	$(GENCLASS) double d 1DBase
	mv nml_d1DBase.h $(HEADER1)

$(PLATFORM):
	mkdir -p $(PLATFORM)

clean:
	rm -f $(OBJECTS)

sterile:
	rm -f $(ARCHIVE) $(OBJECTS) $(SOURCES) $(HEADERS)
	rm -rf $(PLATFORM)

