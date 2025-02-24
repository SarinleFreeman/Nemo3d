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
ARCHIVE=$(CATALOG)/libnml_bvector.a
HEADER0=$(INCPATH)/nml_bvector.h
HEADER1=$(INCPATH)/nml_b1DBase.h
HEADERS=$(HEADER0) $(HEADER1)
COMPILE=$(CC) $(INCLUDE) $(DEFINES) $(OPTIONS)
GENCLASS=$(BINPATH)/genclass

SOURCES=nml_bvector.c nml_b1DBase.c

OBJECTS=$(PLATFORM)/nml_bvector.o $(PLATFORM)/bvector.o \
	$(PLATFORM)/nml_b1DBase.o

all:
	(cd ../../scalar/bscalar; $(MAKE) \
		CC='$(CC)' DEFINES='$(DEFINES)' OPTIONS='$(OPTIONS)' \
		AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)' PLATFORM='$(PLATFORM)' MAKE='$(MAKE)')
	($(MAKE) $(ARCHIVE) \
		CC='$(CC)' DEFINES='$(DEFINES)' OPTIONS='$(OPTIONS)' \
		AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)' PLATFORM='$(PLATFORM)' MAKE='$(MAKE)')

$(ARCHIVE):		$(CATALOG) $(OBJECTS)
	$(AR) $(AR_OPTIONS) $(ARCHIVE) $(OBJECTS)

$(CATALOG):
	mkdir -p $(CATALOG)

$(PLATFORM)/nml_bvector.o:	$(PLATFORM) $(HEADER0) nml_bvector.c
	$(COMPILE) -c nml_bvector.c
	mv nml_bvector.o $(PLATFORM)

$(PLATFORM)/bvector.o:		$(PLATFORM) $(INCPATH)/bvector.h bvector.c
	$(COMPILE) -c bvector.c
	mv bvector.o $(PLATFORM)

$(HEADER0):		../vector.hP
	$(GENCLASS) char b vector
	mv nml_bvector.h $(HEADER0)

nml_bvector.c:		../vector.cP
	$(GENCLASS) char b vector
	mv nml_bvector.h $(HEADER0)

$(PLATFORM)/nml_b1DBase.o:	$(PLATFORM) $(HEADER1) nml_b1DBase.c
	$(COMPILE) -c nml_b1DBase.c
	mv nml_b1DBase.o $(PLATFORM)

$(PLATFORM)/b2tensor.o:		$(PLATFORM) $(INCPATH)/b2tensor.h b2tensor.c
	$(COMPILE) -c b2tensor.c
	mv b2tensor.o $(PLATFORM)

$(HEADER1):		../1DBase.hP
	$(GENCLASS) char b 1DBase
	mv nml_b1DBase.h $(HEADER1)

nml_b1DBase.c:		../1DBase.cP
	$(GENCLASS) char b 1DBase
	mv nml_b1DBase.h $(HEADER1)

$(PLATFORM):
	mkdir -p $(PLATFORM)

clean:
	rm -f $(OBJECTS)

sterile:
	rm -f $(ARCHIVE) $(OBJECTS) $(SOURCES) $(HEADERS)
	rm -rf $(PLATFORM)

