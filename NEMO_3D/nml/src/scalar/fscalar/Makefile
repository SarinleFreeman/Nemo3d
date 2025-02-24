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
DEFINES= -DNML_MEM_EXT
OPTIONS=-O2 -Wall -ansi -pedantic
AR=ar
AR_OPTIONS=rcvs
PLATFORM=$(HOSTTYPE)
CATALOG=$(LIBPATH)/$(PLATFORM)
ARCHIVE=$(CATALOG)/libnml_fscalar.a
HEADER=$(INCPATH)/nml_fscalar.h
COMPILE=$(CC) $(INCLUDE) $(DEFINES) $(OPTIONS)
GENCLASS=$(BINPATH)/genclass

SOURCES=nml_fscalar.c

OBJECTS=$(PLATFORM)/nml_fscalar.o

all:
	(cd ../../global/global; $(MAKE) \
		CC='$(CC)' DEFINES='$(DEFINES)' OPTIONS='$(OPTIONS)' \
		AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)' PLATFORM='$(PLATFORM)' MAKE='$(MAKE)')
	($(MAKE) $(ARCHIVE) \
		CC='$(CC)' DEFINES='$(DEFINES)' OPTIONS='$(OPTIONS)' \
		AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)' PLATFORM='$(PLATFORM)' MAKE='$(MAKE)')

$(ARCHIVE):		$(CATALOG) $(OBJECTS)
	$(AR) $(AR_OPTIONS) $(ARCHIVE) $(OBJECTS)

$(CATALOG):
	mkdir -p $(CATALOG)

$(PLATFORM)/nml_fscalar.o:	$(PLATFORM) $(HEADER) nml_fscalar.c
	$(COMPILE) -c nml_fscalar.c
	mv nml_fscalar.o $(PLATFORM)

$(PLATFORM):
	mkdir -p $(PLATFORM)

$(HEADER):		../scalar.hP
	$(GENCLASS) float f scalar
	mv nml_fscalar.h $(HEADER)

nml_fscalar.c:		../scalar.cP
	$(GENCLASS) float f scalar
	mv nml_fscalar.h $(HEADER)

clean:
	rm -f $(OBJECTS)

sterile:
	rm -f $(ARCHIVE) $(OBJECTS) $(SOURCES) $(HEADER)
	rm -rf $(PLATFORM)

