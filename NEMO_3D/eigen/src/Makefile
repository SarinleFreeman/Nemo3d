# *****************************************************************************
# The JPL Parallel Eigensolvers package.
# Copyright (C) 2002 California Institute of Technology (Caltech)
#
# This file is part of
# The Jet Propulsion Laboratory (JPL) Parallel Eigensolvers package.
#
# This library is free software, which you can redistribute and/or modify
# under the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation; either version 2.1 of the License, or (at your
# option) any later version.
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
# For additional information, please contact
#   Gerhard Klimeck (gekco@jpl.nasa.gov)
#   Fabiano Oyafuso (fabiano@jpl.nasa.gov)
#   E. Robert Tisdale (E.Robert.Tisdale@jpl.nasa.gov)

# Written by: Chris Bowen
#             Gerhard Klimeck
#             Fabiano Oyafuso
#             E. Robert Tisdale
#
# *****************************************************************************
# $Header: /repo/eigen/src/Makefile,v 1.5 2007/08/10 15:00:13 baeh Exp $
# *****************************************************************************


#================================================================
# Makefile to generate the eigsolver library
#================================================================

#----------------------------------------------------------------
# System configuration: compiler and libraries
#----------------------------------------------------------------


  NML_HOME= ../../nml
  DEPENDS = Makefile.depend

  VERSION = 
  PLATFORM= $(HOSTTYPE)
  OBJ_DIR = ../obj/$(PLATFORM)
  INC_DIR = ../include
  SRC_DIR = ../src
  LIB_DIR = ../lib

ifeq ($(MPI_HOME),)
  INCLUDE = -I$(INC_DIR) -I$(NML_HOME)/include
else
  INCLUDE = -I$(INC_DIR) -I$(NML_HOME)/include -I$(MPI_HOME)/include
endif
	INCLUDE += -I../../nemo3d/src/base -I../../nemo3d/src/math

  CATALOG = $(LIB_DIR)/$(PLATFORM)
  ARCHIVE = $(CATALOG)/libeigen.a
  POINTER = $(LIB_DIR)/libeigen.a

#----------------------------------------------------------------
#  Platform specific options
#----------------------------------------------------------------

# Command to generate Makefile dependencies
  CCMKD = g++ -MM
# LAPACK compiled with ABSOFT f77 -f -B108
# 		  or with GNU f77 -Wall -funroll-all-loops -fno-f2c -O3
# CCOPT = -c -Wall -O3 -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING
  CCOPT = -c -Wall -O3 -DFORTRAN_UNDERSCORE -DMPI3d
# LAPACK compiled with ABSOFT f77 -N109
# CCOPT = -c -Wall -O3 -DFORTRAN_ALLCAPS -DMPI3d -DMPI_TIMING
# CCOPT = -c -Wall -O3 -DFORTRAN_ALLCAPS -DMPI3d
# MPI C
# CC    = $(MPI_HOME)/bin/mpicc
# LD    = $(MPI_HOME)/bin/mpicc
# MPI C++
# CC    = $(MPI_HOME)/bin/mpiCC
# LD    = $(MPI_HOME)/bin/mpiCC
  AR    = /usr/bin/ar -r

#----------------------------------------------------------------
# Source and object files + program name
#----------------------------------------------------------------

  SOURCES = \
	$(SRC_DIR)/cmacopt.c \
	$(SRC_DIR)/io.c \
	$(SRC_DIR)/lanczos.c \
	$(SRC_DIR)/lanczos_blk.c \
	$(SRC_DIR)/macopt.c \
	$(SRC_DIR)/macopt_par.c \
	$(SRC_DIR)/mpi_fake.c \
	$(SRC_DIR)/ritz.c

  HEADERS = \
	$(INC_DIR)/cmacopt.h \
	$(INC_DIR)/io.h \
	$(INC_DIR)/lanczos.h \
	$(INC_DIR)/lanczos_blk.h \
	$(INC_DIR)/macopt.h \
	$(INC_DIR)/macopt_par.h \
	$(INC_DIR)/mpi_fake.h \
	$(INC_DIR)/ritz.h

  OBJECTS = \
	$(OBJ_DIR)/cmacopt.o \
	$(OBJ_DIR)/io.o \
	$(OBJ_DIR)/lanczos.o \
	$(OBJ_DIR)/lanczos_blk.o \
	$(OBJ_DIR)/macopt.o \
	$(OBJ_DIR)/macopt_par.o \
	$(OBJ_DIR)/mpi_fake.o \
	$(OBJ_DIR)/ritz.o

#----------------------------------------------------------------
# User rules
#----------------------------------------------------------------

lib:	$(ARCHIVE)
	@ echo
	@ echo "Library '$(ARCHIVE)' is up-to-date."
	@ echo
#	/bin/rm -rf $(POINTER)
#	ln -f -s $(ARCHIVE) $(POINTER)

depend:	$(DEPENDS)
	@ echo "   " $(DEPENDS) "is now up-to-date."

clean: 
	/bin/rm -f $(OBJ_DIR)/*.o $(DEPENDS)

tar:
	tar cvf $(HOME)/ml-src.tar ../src

all: lib

#----------------------------------------------------------------
# Local rules
#----------------------------------------------------------------

$(DEPENDS):	$(HEADERS)
	@ /bin/rm -f $(DEPENDS)
	@ $(CCMKD) $(INCLUDE) $(SOURCES) > $(DEPENDS)-temp
	@ sed "s?^.*\.o:?$(OBJ_DIR)/&?" $(DEPENDS)-temp > $(DEPENDS)
	@ /bin/rm $(DEPENDS)-temp

$(ARCHIVE):	$(OBJECTS)
	$(AR) $(@) $(OBJECTS)

$(OBJ_DIR)/%.o:	$(SRC_DIR)/%.c	
	$(CC) $(CCOPT) $(INCLUDE) $(<) -o $(@)

$(OBJ_DIR)/%.o:	$(SRC_DIR)/%.f
	$(F77) $(<) -o $(@)

#----------------------------------------------------------------
# Dependencies description
#----------------------------------------------------------------

-include $(DEPENDS)

# DO NOT DELETE THIS LINE -- make depend depends on it.
