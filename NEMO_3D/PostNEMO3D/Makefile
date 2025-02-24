# *****************************************************************************
# The Jet Propulsion Laboratory (JPL) NanoElectronicMOdeling-3D
# PostProcessing package.
# Copyright (C) 2002 California Institute of Technology (Caltech)
#
# This application is free software, which you can redistribute and/or modify
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
#   Seungwon Lee (Seungwon.Lee@jpl.nasa.gov)
#
# Written by:  Seungwon Lee
#
# *****************************************************************************


#================================================================
# Postprocess Makefile
#================================================================

#----------------------------------------------------------------
# System configuration: compiler and libraries
#----------------------------------------------------------------

include ../build/make.inc

#----------------------------------------------------------------
#  Platform specific options
#----------------------------------------------------------------

RM     = rm

#----------------------------------------------------------------
#  Directory Setup
#----------------------------------------------------------------

SRC_DIR       = $(POSTPROCESS_ROOT)/src
LIB_DIR       = $(POSTPROCESS_ROOT)/lib
OBJ_DIR       = $(POSTPROCESS_ROOT)/obj
BIN_DIR       = $(POSTPROCESS_ROOT)/bin
INC_DIR       = $(POSTPROCESS_ROOT)/include

extra_INCL = -I $(NML_ROOT)/include   -I $(NEMO3D_ROOT)/src/math -I $(NEMO3D_ROOT)/src//util
NEMO3D_LIBD   = $(NEMO3D_ROOT)/lib/$(PLATFORM)
NEMO3D_LIB_MATH = $(NEMO3D_LIBD)/libnemo3d-math-$(PLATFORM).a
NEMO3D_LIB_UTIL = $(NEMO3D_LIBD)/libnemo3d-util-$(PLATFORM).a
N3D_ADDTL_LIB = $(NEMO3D_LIB_MATH) $(NEMO3D_LIB_UTIL) $(LIB_LAPACK) $(LIB_BLAS) $(MISCLIB) 

#----------------------------------------------------------------
#  Library Files and Includes
#----------------------------------------------------------------

DEP_FILE = Makefile.depend

# The executable to build
EXEC     = $(BIN_DIR)/PostProcess-$(PLATFORM).ex
EXEC_LNK = $(BIN_DIR)/PostProcess.ex

# Build Libraries and Includes
INC_SRCH_DIR = -I$(INC_DIR) -I$(MPI_HOME)/include $(extra_INCL)
LIBS  = $(LIB_NML) $(N3D_ADDTL_LIB)  -lm



#----------------------------------------------------------------
# Source and object files + program name
#----------------------------------------------------------------

SRC = \
	$(SRC_DIR)/PostProcess.cpp \
	$(SRC_DIR)/Coulomb.cpp \
	$(SRC_DIR)/Dipole.cpp \
	$(SRC_DIR)/ReadData.cpp \
	$(SRC_DIR)/SimpleParser.cpp \
	$(SRC_DIR)/Projection.cpp \
        $(SRC_DIR)/FCCproj.cpp \
        $(SRC_DIR)/FCCguess.cpp \
        $(SRC_DIR)/CubicProj.cpp \
        $(SRC_DIR)/CubicGuess.cpp
#       $(SRC_DIR)/DispersionUnfold.cpp \ 


HEADERS = \
	$(INC_DIR)/FakeMPI.h \
	$(INC_DIR)/Coulomb.h \
	$(INC_DIR)/Dipole.h \
	$(INC_DIR)/ReadData.h \
	$(INC_DIR)/SimpleParser.h \
	$(INC_DIR)/Projection.h \
        $(INC_DIR)/FCCproj.h \
        $(INC_DIR)/FCCguess.h \
        $(INC_DIR)/CubicProj.h \
        $(INC_DIR)/CubicGuess.h
#       $(INC_DIR)/DispersionUnfold.h \


OBJS = \
	$(OBJ_DIR)/PostProcess.o \
	$(OBJ_DIR)/Coulomb.o \
	$(OBJ_DIR)/Dipole.o \
	$(OBJ_DIR)/ReadData.o \
	$(OBJ_DIR)/SimpleParser.o \
	$(OBJ_DIR)/Projection.o \
        $(OBJ_DIR)/FCCproj.o \
        $(OBJ_DIR)/FCCguess.o \
        $(OBJ_DIR)/CubicProj.o \
        $(OBJ_DIR)/CubicGuess.o
#       $(OBJ_DIR)/DispersionUnfold.o \



#----------------------------------------------------------------
# User rules
#----------------------------------------------------------------

all : $(EXEC)
	@ echo "Executable $(EXEC_LNK) is up-to-date"

depend : $(DEP_FILE)
	@ echo "   " $(DEP_FILE) "is now up-to-date"

clean :
	$(RM) -f $(OBJ_DIR)/*.o $(DEP_FILE) core.*

cleanall : clean
	$(RM) -f $(EXEC)
	$(RM) -f $(EXEC-LNK)


#----------------------------------------------------------------
# Local rules
#----------------------------------------------------------------

$(DEP_FILE) : $(HEADERS)
	@ $(RM) -f $(DEP_FILE)
	@ $(CCMKD) $(INC_SRCH_DIR) $(SRC) > $(DEP_FILE)-temp
	@ sed "s?^.*\.o:?$(OBJ_DIR)/&?" $(DEP_FILE)-temp > $(DEP_FILE)
	@ $(RM) $(DEP_FILE)-temp

$(EXEC) : $(OBJS)
	$(LD) $(LDOPT) $(OBJS) $(LIBS) -o $(@)
	ln -s -f $(EXEC) $(EXEC_LNK)
	@ echo "==== CREATED postprocess EXECUTABLE ===="

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp
	$(CC) $(CCOPT) $(INC_SRCH_DIR) $(<) -o $(@)

#----------------------------------------------------------------
# Dependencies description
#----------------------------------------------------------------

-include $(DEP_FILE)

# DO NOT DELETE THIS LINE -- make depend depends on it.

