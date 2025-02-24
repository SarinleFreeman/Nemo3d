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
#
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
ARCHIVE=$(CATALOG)/libnml_fmatrix.a
HEADER0=$(INCPATH)/nml_fmatrix.h
HEADER1=$(INCPATH)/nml_f2DBase.h
HEADER2=$(INCPATH)/nml_f3bands.h
HEADERS=$(HEADER0) $(HEADER1) $(HEADER2)
COMPILE=$(CC) $(INCLUDE) $(DEFINES) $(OPTIONS)
GENCLASS=$(BINPATH)/genclass

SOURCES=nml_fmatrix.c nml_f2DBase.c nml_f3bands.c

OBJECTS=$(PLATFORM)/nml_fmatrix.o $(PLATFORM)/fmatrix.o \
	$(PLATFORM)/nml_f2DBase.o $(PLATFORM)/f2tensor.o \
	$(PLATFORM)/nml_f3bands.o $(PLATFORM)/ftridiag.o

all:
	(cd ../../vector/fvector; $(MAKE) \
		CC='$(CC)' DEFINES='$(DEFINES)' OPTIONS='$(OPTIONS)' \
		AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)' PLATFORM='$(PLATFORM)' MAKE='$(MAKE)')
	($(MAKE) $(ARCHIVE) \
		CC='$(CC)' DEFINES='$(DEFINES)' OPTIONS='$(OPTIONS)' \
		AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)' PLATFORM='$(PLATFORM)' MAKE='$(MAKE)')

$(ARCHIVE):		$(CATALOG) $(OBJECTS)
	$(AR) $(AR_OPTIONS) $(ARCHIVE) $(OBJECTS)

$(CATALOG):
	mkdir -p $(CATALOG)

$(PLATFORM)/nml_fmatrix.o:	$(PLATFORM) $(HEADER0) nml_fmatrix.c
	$(COMPILE) -c nml_fmatrix.c
	mv nml_fmatrix.o $(PLATFORM)

$(PLATFORM)/fmatrix.o:		$(PLATFORM) $(INCPATH)/fmatrix.h fmatrix.c
	$(COMPILE) -c fmatrix.c
	mv fmatrix.o $(PLATFORM)

$(HEADER0):		../matrix.hP
	$(GENCLASS) float f matrix
	mv nml_fmatrix.h $(HEADER0)

nml_fmatrix.c:		../matrix.cP
	$(GENCLASS) float f matrix
	mv nml_dmatrix.h $(HEADER0)

$(PLATFORM)/nml_f2DBase.o:	$(PLATFORM) $(HEADER1) nml_f2DBase.c
	$(COMPILE) -c nml_f2DBase.c
	mv nml_f2DBase.o $(PLATFORM)

$(PLATFORM)/f2tensor.o:		$(PLATFORM) $(INCPATH)/f2tensor.h f2tensor.c
	$(COMPILE) -c f2tensor.c
	mv f2tensor.o $(PLATFORM)

$(HEADER1):		../2DBase.hP
	$(GENCLASS) float f 2DBase
	mv nml_f2DBase.h $(HEADER1)

nml_f2DBase.c:		../2DBase.cP
	$(GENCLASS) float f 2DBase
	mv nml_f2DBase.h $(HEADER1)

$(PLATFORM)/nml_f3bands.o:	$(PLATFORM) $(HEADER2) nml_f3bands.c
	$(COMPILE) -c nml_f3bands.c
	mv nml_f3bands.o $(PLATFORM)

$(PLATFORM)/ftridiag.o:		$(PLATFORM) $(INCPATH)/ftridiag.h ftridiag.c
	$(COMPILE) -c ftridiag.c
	mv ftridiag.o $(PLATFORM)

$(HEADER2):		../3bands.hP
	$(GENCLASS) float f 3bands
	mv nml_f3bands.h $(HEADER2)

nml_f3bands.c:		../3bands.cP
	$(GENCLASS) float f 3bands
	mv nml_f3bands.h $(HEADER2)

$(PLATFORM):
	mkdir -p $(PLATFORM)

clean:
	rm -f $(OBJECTS)

sterile:
	rm -f $(ARCHIVE) $(OBJECTS) $(SOURCES) $(HEADERS)
	rm -rf $(PLATFORM)

