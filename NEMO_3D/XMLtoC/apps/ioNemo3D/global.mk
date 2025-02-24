# *****************************************************************************
# The Jet Propulsion Laboratory (JPL) XML-to-C++ package.
# Copyright (C) 2002 California Institute of Technology (Caltech)
#
# This library is free software, which you can redistribute and/or modify
# under the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation; either version 2.1 of the License, or (at your
# option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library; see the file COPYING. If not, write to the
# Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA  02111-1307  USA
#
# For additional information, please contact
#   Fabiano Oyafuso (fabiano@jpl.nasa.gov)
#
# Written by:  Fabiano Oyafuso
#              Hook Hua
#              Ed Vinyard
#
# This product includes software developed by the Apache Software Foundation
# (http://www.apache.org/).
#
# *****************************************************************************
# $Header: /repo/XMLtoC/global.mk,v 1.3 2003/10/08 16:15:19 hook Exp $
# *****************************************************************************
 

#================================================================
# This file contains the default directory structure and build
# configurations for a number of applications developed at JPL.
# Typically this file will be included in the top level makefile
# for a particular application.
#================================================================

MAKE = make

SHELL = /bin/sh
OS        := $(shell uname)
PLATFORM = $(HOSTTYPE)

ifeq ($(PJCT),)
	$(error "Undefined environment variable, $PJCT, needed by global.mk")
endif

ifeq ($($(PJCT)_ROOT),)
	$(error "Undefined environment variable, $(PJCT)_ROOT, needed by global.mk")
endif

ifeq ($($(PJCT)_BLD),)
	$(error "Undefined environment variable, $(PJCT)_BLD, needed by global.mk")
endif

# location of makefile includes
$(PJCT)_TOPMK   = $($(PJCT)_ROOT)/Makefile
$(PJCT)_DEPMK   = $($(PJCT)_ROOT)/$(PJCT)_dependencies.mk
$(PJCT)_FILMK   = $($(PJCT)_ROOT)/$(PJCT)_files.mk
$(PJCT)_CMPMK   = $($(PJCT)_ROOT)/$(PJCT)_cmplopts.mk

# directory structure
$(PJCT)_SRCD   = $($(PJCT)_ROOT)/src
$(PJCT)_INCLD  = $($(PJCT)_ROOT)/src
$(PJCT)_OBJD   = $($(PJCT)_ROOT)/obj/$(PLATFORM)
$(PJCT)_BIND   = $($(PJCT)_ROOT)/bin
$(PJCT)_LIBD   = $($(PJCT)_ROOT)/lib/$(PLATFORM)



# Absoft Fortran libraries
ifneq ($(ABSOFT),)
	LIB_F77MATH = $(ABSOFT)/lib/libf77math.a
	LIB_F90MATH = $(ABSOFT)/lib/libf90math.a
	LIB_FIO     = $(ABSOFT)/lib/libfio.a
endif


#----------------------------------------------------------------
#  Default platform specific options
#----------------------------------------------------------------

ifeq ($($(PJCT)_BLD),gnuLinux)
	CCMKD    = g++ -MM
	CC       = g++ -c
	CCOPT   = -Wall -O2
	LD       = g++
	AR       = ar
	AR_OPTIONS = r
	MISCLIB = -ldl -lpthread -lm
endif


ifeq ($($(PJCT)_BLD),absoft_mpi_fast)
	CCMKD    = g++ -MM
	CCOPT   = -Wall -O2
	ifneq ($(MPI_HOME),)
		CC = $(MPI_HOME)/bin/mpiCC -c
		LD = $(MPI_HOME)/bin/mpiCC
	else
		CC = g++ -c
		LD = g++
	endif
	AR         = ar
	AR_OPTIONS = r
	MISCLIB = $(LIB_F77MATH) $(LIB_F90MATH) $(LIB_FIO) -ldl -lpthread -lm
endif


# PYTHON
PYTHON_ROOT = /usr/lib/python1.5
PYTHON_EXPORT_INCL = -I/usr/include/python1.5
PYTHON_EXPORT_LIB = $(PYTHON_ROOT)/config/libpython1.5.a


# XERCES-C++
XERCESCROOT = ../xerces-c-src2_1_0
XERCESC_EXPORT_INCL = \
	-I$(XERCESCROOT)/include \
	-I$(XERCESCROOT)/include/xercesc \
	-I$(XERCESCROOT)/include/xercesc/dom/deprecated
XERCESC_EXPORT_LIB = $(XERCESCROOT)/lib/libxerces-c.so
