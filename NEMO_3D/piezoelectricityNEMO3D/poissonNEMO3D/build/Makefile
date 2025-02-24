############################################################################
# Copyright (C) 2005 Shaikh Shahid Ahmed
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
# Shaikh S. Ahmed 
# Department of Electrical and Computer Engineering
# Purdue University
# West Lafayette, IN, USA.
# Email: ssahmed@purdue.edu
########################################################################


SHELL = /bin/sh
OS    = $(shell uname)
HOSTTYPE = $(shell echo $${HOSTTYPE})

OSTYPE =
ifeq ($(OS),Linux)
        OSTYPE = LINUX
endif
ifeq ($(OS),Darwin)
        OSTYPE = MACOSX
endif
ifeq ($(OS),IRIX64)
        OSTYPE = IRIX
endif
ifeq ($(OS),AIX)
        OSTYPE = AIX
endif

ifeq ($(BUILD_TYPE),)
#       BUILD_TYPE = intel8_p4_mpi_static
#       BUILD_TYPE = intel9_32_mpi_intelfast
#       BUILD_TYPE = intel9_64_mpi_intelfast
#       BUILD_TYPE = intel8_p4_openmp_static
#       BUILD_TYPE = aix_static64
#       BUILD_TYPE = sgi_intel
#       BUILD_TYPE = sgi_f90
#	BUILD_TYPE = gnu_g95
#	BUILD_TYPE = gnu_gfortran
       BUILD_TYPE = intel_ifort
endif

PLATFORM = $(HOSTTYPE)_$(BUILD_TYPE)

RM     = /bin/rm
MKDIR  = /bin/mkdir
AR     = /usr/bin/ar
AROPT  =  rcvs
MAKE   = make
RANLIB = ranlib

ROOT             = $(shell pwd)/..
SRC_DIR          = $(ROOT)/src
LIB_DIR          = $(ROOT)/lib/$(PLATFORM)
OBJ_DIR          = $(ROOT)/obj/$(PLATFORM)
BIN_DIR          = $(ROOT)/bin
INC_DIR          = $(ROOT)/include
EXA_DIR          = $(ROOT)/example

#DEPMK            = $(ROOT)/src/Makefile.depend

POISSON3D_EXEC     = $(BIN_DIR)/poisson3d-$(PLATFORM).ex
POISSON3D_EXEC_LNK = $(BIN_DIR)/poisson3d.ex


LIB_POISSON3D  = $(LIB_DIR)/libpoisson3d.a  
LIB_RAPPTURE = 
LIB_MISC    = -ldl -lpthread -lm

ifeq ($(OS),Darwin)
	LIB_PATHS =
else
	LIB_PATHS =
endif


# Platform project build profiles (default set above) 

# Intel Pentium4 optimized MPI Static
ifeq ($(BUILD_TYPE),intel8_p4_mpi_static)
	AR = xiar
         ifeq ($(MPI_HOME),)
           @ echo "Error - Undefined environment variable,MPI_HOME"
         endif
	F90 = $(MPI_HOME)/bin/mpif90
	F90MKD = $(MPI_HOME)/bin/mpif90  
	F90OPT = -c -g  -axW -tpp7 -xW -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING  -I/opt/intel_cc_80/include/c++ -FI -D$(OSTYPE)
	LD = $(MPI_HOME)/bin/mpif90
	LDOPT  = -g -Wl,-static  -axW -tpp7 -xW -FI
	MISCLIB0 = -static -lpthread -lm /opt/intel_fc_80/lib/libifcore.a /opt/intel_fc_80/lib/libifport.a -ldl -lsvml
	MISCINCL0  =
endif

# Intel Pentium4 optimized MPI Static
ifeq ($(BUILD_TYPE),intel9_32_mpi_intelfast)
	AR = xiar
         ifeq ($(MPI_HOME),)
           @ echo "Error - Undefined environment variable,MPI_HOME"
         endif
	F90 = $(MPI_HOME)/bin/mpif90
	F90MKD = $(MPI_HOME)/bin/mpif90  
	F90OPT = -c -g  -axW -tpp7 -xW -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING  -I/opt/cc/9.1.045/include/c++ -FI -D$(OSTYPE)
	LD = $(MPI_HOME)/bin/mpif90
	LDOPT  = -g -Wl,-static  -axW -tpp7 -xW -FI
	MISCLIB0 = -static -lpthread -lm /opt/intel/fc/9.1.040/lib/libifcore.a /opt/intel/fc/9.1.040/lib/libifport.a -ldl -lsvml
	MISCINCL0  =
endif

ifeq ($(BUILD_TYPE),intel9_64_mpi_intelfast)
	AR = xiar
         ifeq ($(MPI_HOME),)
           @ echo "Error - Undefined environment variable,MPI_HOME"
         endif
	F90 = $(MPI_HOME)/bin/mpif90
	F90MKD = $(MPI_HOME)/bin/mpif90  
	F90OPT = -c -g  -axW -tpp7 -xW -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING  -I/opt/cce/9.1.045/include/c++ -FI -D$(OSTYPE)
	LD = $(MPI_HOME)/bin/mpif90
	LDOPT  = -g -Wl,-static  -axW -tpp7 -xW -FI
	MISCLIB0 =  -static -lpthread -lm /opt/intel/fce/9.1.040/lib/libifcore.a /opt/intel/fce/9.1.040/lib/libifport.a -ldl -lsvml
	MISCINCL0  =
endif

# Intel ifort
ifeq ($(BUILD_TYPE),intel_ifort)
        AR = xiar
        F90 = ifort
        F90MKD = $(F90)
        F90OPT = -c -g -axW -xW -FI -D$(OSTYPE)
	LD = $(F90)
        LDOPT  = -g -Wl,-static -axW -xW -FI
        MISCLIB0   = -static
        MISCINCL0  =
endif

# Intel Pentium4 optimized OpenMP Static
ifeq ($(BUILD_TYPE),intel8_p4_openmp_static)
	AR = xiar
	F90 = ifort
	F90MKD = $(F90)
	F90OPT = -c -O3 -axW -tpp7 -xW -openmp -FI -D$(OSTYPE) 
	LD = $(F90)
	LDOPT  = -O3 -Wl,-static  -axW -tpp7 -xW -openmp -FI 
	MISCLIB0   = -static
	MISCINCL0  =
endif

# AIX 64bit 
ifeq ($(BUILD_TYPE),aix_static64)
	HOSTTYPE=AIX
	MAKE=gmake
	AR = ar
	F90 = xlf90_r -qsuffix=f=f90
	F90MKD = $(F90)
	F90OPT = -c -O3 -qsmp=omp  -qfixed -q64 -qarch=auto -qtune=auto -qlanglvl=extended -qsave -D$(OSTYPE) 
	LD = $(F90)
	LDOPT  = -O3 -qsmp=omp -qfixed -q64 -qarch=auto -qtune=auto -qlanglvl=extended -qsave -brename:.flush,.flush_

	MISCLIB0   = 
	MISCINCL0  =
endif

# SGI with Intel and OpenMP Static
ifeq ($(BUILD_TYPE),sgi_intel)
	AR = xiar
	F90 = ifort
	F90MKD = $(F90)
	F90OPT = -c -g -FI -D$(OSTYPE) -cm 
	LD = $(F90)
	LDOPT  = -g -Wl,-static -FI -cm
	MISCLIB0   = -static
	MISCINCL0  =
endif

# SGI with f90 and OpenMP Static
ifeq ($(BUILD_TYPE),sgi_f90)
	AR = xiar
	F90 = f90
	F90MKD = $(F90)
	F90OPT = -c -g -FI
	LD = $(F90)
	LDOPT  = -g -static -FI  
	MISCLIB0   = -static
	MISCINCL0  =
endif

# GNU fortran compiler
ifeq ($(BUILD_TYPE),gnu_g95)
	AR = ar
	F90 = g95
	F90MKD = $(F90)
	F90OPT = -c -O    
	LD = $(F90)
	LDOPT  = -O -static    
	MISCLIB0   = -static
	MISCINCL0  =
endif

# GNU fortran compiler
ifeq ($(BUILD_TYPE),gnu_gfortran)
	AR = ar
	F90 = gfortran
	F90MKD = $(F90)
	#F90OPT = -c -O    
	F90OPT = -c -FI     
	LD = $(F90)
	LDOPT  = -O -static    
	LDOPT  = -FI  
	MISCLIB0   = -lpthread -lm
	MISCINCL0  =
endif

# Misc additional library and include flags
ifeq ($(OS),Linux)
	MISCLIB  = $(MISCLIB0) -lutil $(LIB_PATHS) 
	MISCLIB  = $(MISCLIB0)  $(LIB_PATHS) 
else
	MISCLIB  = $(MISCLIB0) $(LIB_PATHS)   
endif
ifeq ($(OS),Darwin)
	MISCINCL = $(MISCINCL0) -I/usr/include/objc 
else
	MISCINCL = $(MISCINCL0) 
endif

INC_SRCH_DIR = -I$(INC_DIR) $(MISCINCL)
LIBS = $(LIB_MISC) $(MISCLIB) 



# Files  

HEADERS = \
	$(INC_DIR)/params.h \
	$(INC_DIR)/common_par.inc  
SRC =  \
	$(SRC_DIR)/POISSON.f 
OBJS =  \
	$(OBJ_DIR)/POISSON.o 


# Rules

all:    $(POISSON3D_EXEC)
 
#depend : $(DEPMK) 
#	@ echo " " $(DEPMK) " now up-to-date"

#$(DEPMK) : $(HEADERS) 
#	@ echo "===== MAKING DEPENDENCIES ====="
#	@ cd $(SRC_DIR)
#	@ $(RM) -f  $(DEPMK)
#	@ $(F90MKD) $(INC_SRCH_DIR) $(SRC) > $(DEPMK)-temp
#	@ sed "s?^.*\.o:?$(OBJ_DIR)/&?" $(DEPMK)-temp > $(DEPMK)
#	@ rm $(DEPMK)-temp


$(POISSON3D_EXEC): $(LIB_POISSON3D)  
	$(MKDIR) -p $(BIN_DIR)
	$(LD) $(LDOPT) $(INC_SRCH_DIR) -o $(@) $(LIB_POISSON3D) $(LIBS)
	ln -s -f $(POISSON3D_EXEC) $(POISSON3D_EXEC_LNK) 
	@ echo "==== CREATED POISSON3D EXECUTABLE ===="
	@ echo

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.f
	$(MKDIR) -p $(OBJ_DIR)
	$(F90) $(F90OPT) $(INC_SRCH_DIR) $(<) -o $(@)

$(LIB_POISSON3D) : $(OBJS)
	$(MKDIR) -p $(LIB_DIR)
	$(AR) $(AROPT) $(@) $(OBJS)
	@ echo "==== CREATED POISSON3D LIBRARY ===="
	@ echo


# Clean 

clean:
	$(RM) $(OBJ_DIR)/*.o #$(SRC_DIR)/*.depend 
	@ echo "==== CLEANED *.o FILES ===="
	@ echo

cleanall: 
	$(RM) $(BIN_DIR)/*.ex
	$(RM) $(OBJ_DIR)/*.o #$(SRC_DIR)/*.depend 
	$(RM) $(LIB_POISSON3D) 
	@ echo "==== CLEANED OBJ LIB BIN FILES ===="
	@ echo

#-include $(DEPMK)

# DO NOT DELETE THIS LINE: make depend depends on it

