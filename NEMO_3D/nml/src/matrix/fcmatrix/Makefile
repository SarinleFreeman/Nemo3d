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
ARCHIVE=$(CATALOG)/libnml_fcmatrix.a
HEADER0=$(INCPATH)/nml_fcmatrix.h
HEADER2=$(INCPATH)/nml_fc3bands.h
HEADERS=$(HEADER0) $(HEADER2)
COMPILE=$(CC) $(INCLUDE) $(DEFINES) $(OPTIONS)
GENCLASS=$(BINPATH)/genclass

SOURCES=nml_fcmatrix.c nml_fc3bands.c

OBJECTS=$(PLATFORM)/nml_fcmatrix.o $(PLATFORM)/kmatrix.o \
	$(PLATFORM)/nml_fc3bands.o $(PLATFORM)/ktridiag.o

all:
	(cd ../../vector/fcvector; $(MAKE) \
		CC='$(CC)' DEFINES='$(DEFINES)' OPTIONS='$(OPTIONS)' \
		AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)' PLATFORM='$(PLATFORM)' MAKE='$(MAKE)')
	(cd ../../matrix/fmatrix; $(MAKE) \
		CC='$(CC)' DEFINES='$(DEFINES)' OPTIONS='$(OPTIONS)' \
		AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)' PLATFORM='$(PLATFORM)' MAKE='$(MAKE)')
	($(MAKE) $(ARCHIVE) \
		CC='$(CC)' DEFINES='$(DEFINES)' OPTIONS='$(OPTIONS)' \
		AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)' PLATFORM='$(PLATFORM)' MAKE='$(MAKE)')

$(ARCHIVE):		$(CATALOG) $(OBJECTS)
	$(AR) $(AR_OPTIONS) $(ARCHIVE) $(OBJECTS)

$(CATALOG):
	mkdir -p $(CATALOG)

$(PLATFORM)/nml_fcmatrix.o:	$(PLATFORM) $(HEADER0) nml_fcmatrix.c
	$(COMPILE) -c nml_fcmatrix.c
	mv nml_fcmatrix.o $(PLATFORM)

$(PLATFORM)/kmatrix.o:		$(PLATFORM) $(INCPATH)/kmatrix.h kmatrix.c
	$(COMPILE) -c kmatrix.c
	mv kmatrix.o $(PLATFORM)

$(HEADER0):		../matrix.hP
	$(GENCLASS) float f complex c matrix
	mv nml_fcmatrix.h $(HEADER0)

nml_fcmatrix.c: 	../matrix.cP
	$(GENCLASS) float f complex c matrix
	mv nml_fcmatrix.h $(HEADER0)

$(PLATFORM)/nml_fc3bands.o:	$(PLATFORM) $(HEADER2) nml_fc3bands.c
	$(COMPILE) -c nml_fc3bands.c
	mv nml_fc3bands.o $(PLATFORM)

$(PLATFORM)/ktridiag.o:		$(PLATFORM) $(INCPATH)/ktridiag.h ktridiag.c
	$(COMPILE) -c ktridiag.c
	mv ktridiag.o $(PLATFORM)

$(HEADER2):		../3bands.hP
	$(GENCLASS) float f complex c 3bands
	mv nml_fc3bands.h $(HEADER2)

nml_fc3bands.c: 	../3bands.cP
	$(GENCLASS) float f complex c 3bands
	mv nml_fc3bands.h $(HEADER2)

$(PLATFORM):
	mkdir -p $(PLATFORM)

clean:
	rm -f $(OBJECTS)

sterile:
	rm -f $(ARCHIVE) $(OBJECTS) $(SOURCES) $(HEADERS)
	rm -rf $(PLATFORM)

