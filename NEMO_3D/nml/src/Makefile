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
# Written by E. Robert Tisdale
# Modified by Marek J. Korkusinski, October 2004
# *****************************************************************************


# Set the default C or C++ compiler.  Note that
# the C++ compiler "mangles" function names to create the symbols
# that it leaves behind in object files for the link editor.
# Sequential version:
# CC=gcc
# CC=g++
# Parallel version:
# CC=$(MPI_HOME)/bin/mpicc

MAKE = make

ifneq ($(MPI_HOME),)
	CC = $(MPI_HOME)/bin/mpiCC
	LD = $(MPI_HOME)/bin/mpiCC
else
	CC = g++
	LD = g++
endif

#DEFINES = -DNML_ERRANT
DEFINES=-DNML_ERRANT -DNML_MEM_EXT
# DEFINES=-DNML_INLINE
# DEFINES=-DNML_MACROS
OPTIONS=-O2 -Wall -ansi -pedantic
AR=ar
AR_OPTIONS=rcvs
PLATFORM=$(HOSTTYPE)
CATALOG=../lib/$(PLATFORM)
ARCHIVE=$(CATALOG)/libnml.a

all:
	(cd tensor/itensor; $(MAKE) \
		CC='$(CC)' DEFINES='$(DEFINES)' OPTIONS='$(OPTIONS)' \
		AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)' PLATFORM='$(PLATFORM)' MAKE='$(MAKE)')
	(cd tensor/stensor; $(MAKE) \
		CC='$(CC)' DEFINES='$(DEFINES)' OPTIONS='$(OPTIONS)' \
		AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)' PLATFORM='$(PLATFORM)' MAKE='$(MAKE)')
	(cd tensor/btensor; $(MAKE) \
		CC='$(CC)' DEFINES='$(DEFINES)' OPTIONS='$(OPTIONS)' \
		AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)' PLATFORM='$(PLATFORM)' MAKE='$(MAKE)')
	(cd tensor/ftensor; $(MAKE) \
		CC='$(CC)' DEFINES='$(DEFINES)' OPTIONS='$(OPTIONS)' \
		AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)' PLATFORM='$(PLATFORM)' MAKE='$(MAKE)')
	(cd tensor/dtensor; $(MAKE) \
		CC='$(CC)' DEFINES='$(DEFINES)' OPTIONS='$(OPTIONS)' \
		AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)' PLATFORM='$(PLATFORM)' MAKE='$(MAKE)')
	(cd tensor/dctensor; $(MAKE) \
		CC='$(CC)' DEFINES='$(DEFINES)' OPTIONS='$(OPTIONS)' \
		AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)' PLATFORM='$(PLATFORM)' MAKE='$(MAKE)')
	(cd tensor/fctensor; $(MAKE) \
		CC='$(CC)' DEFINES='$(DEFINES)' OPTIONS='$(OPTIONS)' \
		AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)' PLATFORM='$(PLATFORM)' MAKE='$(MAKE)')
	($(AR) $(AR_OPTIONS) $(ARCHIVE) */*/$(PLATFORM)/*.o */$(PLATFORM)/*.o )

clean:
	(cd tensor/dctensor;	$(MAKE) clean)
	(cd tensor/fctensor;	$(MAKE) clean)
	(cd tensor/dtensor;	$(MAKE) clean)
	(cd tensor/itensor;	$(MAKE) clean)
	(cd tensor/btensor;	$(MAKE) clean)
	(cd tensor/stensor;	$(MAKE) clean)
	(cd tensor/ftensor;	$(MAKE) clean)
	(cd matrix/dcmatrix;	$(MAKE) clean)
	(cd matrix/fcmatrix;	$(MAKE) clean)
	(cd matrix/dmatrix;	$(MAKE) clean)
	(cd matrix/imatrix;	$(MAKE) clean)
	(cd matrix/bmatrix;	$(MAKE) clean)
	(cd matrix/smatrix;	$(MAKE) clean)
	(cd matrix/fmatrix;	$(MAKE) clean)
	(cd vector/dcvector;	$(MAKE) clean)
	(cd vector/fcvector;	$(MAKE) clean)
	(cd vector/dvector;	$(MAKE) clean)
	(cd vector/ivector;	$(MAKE) clean)
	(cd vector/bvector;	$(MAKE) clean)
	(cd vector/svector;	$(MAKE) clean)
	(cd vector/fvector;	$(MAKE) clean)
	(cd scalar/dcscalar;	$(MAKE) clean)
	(cd scalar/fcscalar;	$(MAKE) clean)
	(cd scalar/dscalar;	$(MAKE) clean)
	(cd scalar/iscalar;	$(MAKE) clean)
	(cd scalar/bscalar;	$(MAKE) clean)
	(cd scalar/sscalar;	$(MAKE) clean)
	(cd scalar/fscalar;	$(MAKE) clean)
	(cd global/global;	$(MAKE) clean)

cleanall: clean
	(rm -f $(ARCHIVE); rm -rf $(CATALOG))

sterile:
	(cd tensor/dctensor;	$(MAKE) sterile)
	(cd tensor/fctensor;	$(MAKE) sterile)
	(cd tensor/dtensor;	$(MAKE) sterile)
	(cd tensor/itensor;	$(MAKE) sterile)
	(cd tensor/btensor;	$(MAKE) sterile)
	(cd tensor/stensor;	$(MAKE) sterile)
	(cd tensor/ftensor;	$(MAKE) sterile)
	(cd matrix/dcmatrix;	$(MAKE) sterile)
	(cd matrix/fcmatrix;	$(MAKE) sterile)
	(cd matrix/dmatrix;	$(MAKE) sterile)
	(cd matrix/imatrix;	$(MAKE) sterile)
	(cd matrix/bmatrix;	$(MAKE) sterile)
	(cd matrix/smatrix;	$(MAKE) sterile)
	(cd matrix/fmatrix;	$(MAKE) sterile)
	(cd vector/dcvector;	$(MAKE) sterile)
	(cd vector/fcvector;	$(MAKE) sterile)
	(cd vector/dvector;	$(MAKE) sterile)
	(cd vector/ivector;	$(MAKE) sterile)
	(cd vector/bvector;	$(MAKE) sterile)
	(cd vector/svector;	$(MAKE) sterile)
	(cd vector/fvector;	$(MAKE) sterile)
	(cd scalar/dcscalar;	$(MAKE) sterile)
	(cd scalar/fcscalar;	$(MAKE) sterile)
	(cd scalar/dscalar;	$(MAKE) sterile)
	(cd scalar/iscalar;	$(MAKE) sterile)
	(cd scalar/bscalar;	$(MAKE) sterile)
	(cd scalar/sscalar;	$(MAKE) sterile)
	(cd scalar/fscalar;	$(MAKE) sterile)
	(cd global/global;	$(MAKE) sterile)
	(rm -f $(ARCHIVE); rm -rf $(CATALOG))

