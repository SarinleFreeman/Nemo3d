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
ARCHIVE=$(CATALOG)/libnml_bmatrix.a
HEADER0=$(INCPATH)/nml_bmatrix.h
HEADER1=$(INCPATH)/nml_b2DBase.h
HEADERS=$(HEADER0) $(HEADER1)
COMPILE=$(CC) $(INCLUDE) $(DEFINES) $(OPTIONS)
GENCLASS=$(BINPATH)/genclass

SOURCES=nml_bmatrix.c nml_b2DBase.c

OBJECTS=$(PLATFORM)/nml_bmatrix.o $(PLATFORM)/bmatrix.o \
	$(PLATFORM)/nml_b2DBase.o $(PLATFORM)/b2tensor.o

all:
	(cd ../../vector/bvector; $(MAKE) \
		CC='$(CC)' DEFINES='$(DEFINES)' OPTIONS='$(OPTIONS)' \
		AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)' PLATFORM='$(PLATFORM)' MAKE='$(MAKE)')
	($(MAKE) $(ARCHIVE) \
		CC='$(CC)' DEFINES='$(DEFINES)' OPTIONS='$(OPTIONS)' \
		AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)' PLATFORM='$(PLATFORM)' MAKE='$(MAKE)')

$(ARCHIVE):		$(CATALOG) $(OBJECTS)
	$(AR) $(AR_OPTIONS) $(ARCHIVE) $(OBJECTS)

$(CATALOG):
	mkdir -p $(CATALOG)

$(PLATFORM)/nml_bmatrix.o:	$(PLATFORM) $(HEADER0) nml_bmatrix.c
	$(COMPILE) -c nml_bmatrix.c
	mv nml_bmatrix.o $(PLATFORM)

$(PLATFORM)/bmatrix.o:		$(PLATFORM) $(INCPATH)/bmatrix.h bmatrix.c
	$(COMPILE) -c bmatrix.c
	mv bmatrix.o $(PLATFORM)

$(HEADER0):		../matrix.hP
	$(GENCLASS) char b matrix
	mv nml_bmatrix.h $(HEADER0)

nml_bmatrix.c:		../matrix.cP
	$(GENCLASS) char b matrix
	mv nml_bmatrix.h $(HEADER0)

$(PLATFORM)/nml_b2DBase.o:	$(PLATFORM) $(HEADER1) nml_b2DBase.c
	$(COMPILE) -c nml_b2DBase.c
	mv nml_b2DBase.o $(PLATFORM)

$(PLATFORM)/b2tensor.o:		$(PLATFORM) $(INCPATH)/b2tensor.h b2tensor.c
	$(COMPILE) -c b2tensor.c
	mv b2tensor.o $(PLATFORM)

$(HEADER1):		../2DBase.hP
	$(GENCLASS) char b 2DBase
	mv nml_b2DBase.h $(HEADER1)

nml_b2DBase.c:		../2DBase.cP
	$(GENCLASS) char b 2DBase
	mv nml_b2DBase.h $(HEADER1)

$(PLATFORM):
	mkdir -p $(PLATFORM)

clean:
	rm -f $(OBJECTS)

sterile:
	rm -f $(ARCHIVE) $(OBJECTS) $(SOURCES) $(HEADERS)
	rm -rf $(PLATFORM)

