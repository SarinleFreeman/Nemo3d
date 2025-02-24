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
# Written by Marek J. Korkusinski
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
ARCHIVE=$(CATALOG)/libnml_svector.a
HEADER0=$(INCPATH)/nml_svector.h
HEADER1=$(INCPATH)/nml_s1DBase.h
HEADERS=$(HEADER0) $(HEADER1)
COMPILE=$(CC) $(INCLUDE) $(DEFINES) $(OPTIONS)
GENCLASS=$(BINPATH)/genclass

SOURCES=nml_svector.c nml_s1DBase.c

OBJECTS=$(PLATFORM)/nml_svector.o $(PLATFORM)/svector.o \
	$(PLATFORM)/nml_s1DBase.o

all:
	(cd ../../scalar/sscalar; $(MAKE) \
		CC='$(CC)' DEFINES='$(DEFINES)' OPTIONS='$(OPTIONS)' \
		AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)' PLATFORM='$(PLATFORM)' MAKE='$(MAKE)')
	($(MAKE) $(ARCHIVE) \
		CC='$(CC)' DEFINES='$(DEFINES)' OPTIONS='$(OPTIONS)' \
		AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)' PLATFORM='$(PLATFORM)' MAKE='$(MAKE)')

$(ARCHIVE):		$(CATALOG) $(OBJECTS)
	$(AR) $(AR_OPTIONS) $(ARCHIVE) $(OBJECTS)

$(CATALOG):
	mkdir -p $(CATALOG)

$(PLATFORM)/nml_svector.o:	$(PLATFORM) $(HEADER0) nml_svector.c
	$(COMPILE) -c nml_svector.c
	mv nml_svector.o $(PLATFORM)

$(PLATFORM)/svector.o:		$(PLATFORM) $(INCPATH)/svector.h svector.c
	$(COMPILE) -c svector.c
	mv svector.o $(PLATFORM)

$(HEADER0):		../vector.hP
	$(GENCLASS) short s vector
	mv nml_svector.h $(HEADER0)

nml_svector.c:		../vector.cP
	$(GENCLASS) short s vector
	mv nml_svector.h $(HEADER0)

$(PLATFORM)/nml_s1DBase.o:	$(PLATFORM) $(HEADER1) nml_s1DBase.c
	$(COMPILE) -c nml_s1DBase.c
	mv nml_s1DBase.o $(PLATFORM)

$(PLATFORM)/s2tensor.o:		$(PLATFORM) $(INCPATH)/s2tensor.h s2tensor.c
	$(COMPILE) -c s2tensor.c
	mv s2tensor.o $(PLATFORM)

$(HEADER1):		../1DBase.hP
	$(GENCLASS) short s 1DBase
	mv nml_s1DBase.h $(HEADER1)

nml_s1DBase.c:		../1DBase.cP
	$(GENCLASS) short s 1DBase
	mv nml_s1DBase.h $(HEADER1)

$(PLATFORM):
	mkdir -p $(PLATFORM)

clean:
	rm -f $(OBJECTS)

sterile:
	rm -f $(ARCHIVE) $(OBJECTS) $(SOURCES) $(HEADERS)
	rm -rf $(PLATFORM)

