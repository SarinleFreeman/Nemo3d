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
# $Header: /repo/nml/src/matrix/dcmatrix/Makefile,v 1.3 2004/10/14 16:02:41 marek Exp $
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
ARCHIVE=$(CATALOG)/libnml_dcmatrix.a
HEADER0=$(INCPATH)/nml_dcmatrix.h
HEADER2=$(INCPATH)/nml_dc3bands.h
HEADERS=$(HEADER0) $(HEADER2)
COMPILE=$(CC) $(INCLUDE) $(DEFINES) $(OPTIONS)
GENCLASS=$(BINPATH)/genclass

SOURCES=nml_dcmatrix.c nml_dc3bands.c

OBJECTS=$(PLATFORM)/nml_dcmatrix.o $(PLATFORM)/cmatrix.o \
	$(PLATFORM)/nml_dc3bands.o $(PLATFORM)/ctridiag.o

all:
	(cd ../../vector/dcvector; $(MAKE) \
		CC='$(CC)' DEFINES='$(DEFINES)' OPTIONS='$(OPTIONS)' \
		AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)' PLATFORM='$(PLATFORM)' MAKE='$(MAKE)')
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

$(PLATFORM)/nml_dcmatrix.o:	$(PLATFORM) $(HEADER0) nml_dcmatrix.c
	$(COMPILE) -c nml_dcmatrix.c
	mv nml_dcmatrix.o $(PLATFORM)

$(PLATFORM)/cmatrix.o:		$(PLATFORM) $(INCPATH)/cmatrix.h cmatrix.c
	$(COMPILE) -c cmatrix.c
	mv cmatrix.o $(PLATFORM)

$(HEADER0):		../matrix.hP
	$(GENCLASS) double d complex c matrix
	mv nml_dcmatrix.h $(HEADER0)

nml_dcmatrix.c: 	../matrix.cP
	$(GENCLASS) double d complex c matrix
	mv nml_dcmatrix.h $(HEADER0)

$(PLATFORM)/nml_dc3bands.o:	$(PLATFORM) $(HEADER2) nml_dc3bands.c
	$(COMPILE) -c nml_dc3bands.c
	mv nml_dc3bands.o $(PLATFORM)

$(PLATFORM)/ctridiag.o:		$(PLATFORM) $(INCPATH)/ctridiag.h ctridiag.c
	$(COMPILE) -c ctridiag.c
	mv ctridiag.o $(PLATFORM)

$(HEADER2):		../3bands.hP
	$(GENCLASS) double d complex c 3bands
	mv nml_dc3bands.h $(HEADER2)

nml_dc3bands.c: 	../3bands.cP
	$(GENCLASS) double d complex c 3bands
	mv nml_dc3bands.h $(HEADER2)

$(PLATFORM):
	mkdir -p $(PLATFORM)

clean:
	rm -f $(OBJECTS)

sterile:
	rm -f $(ARCHIVE) $(OBJECTS) $(SOURCES) $(HEADERS)
	rm -rf $(PLATFORM)

