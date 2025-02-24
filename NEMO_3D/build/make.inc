# *****************************************************************************
# The Jet Propulsion Laboratory (JPL) NanoElectronicMOdeling-3D package.
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
#   Gerhard Klimeck (gekco@jpl.nasa.gov)
#   Fabiano Oyafuso (fabiano@jpl.nasa.gov)
#
# Written by:  Chris Bowen
#              Gerhard Klimeck
#              Fabiano Oyafuso
#              Seungwon Lee
#              Olga Lazarenkova
#              Hook Hua
#
# This product includes software developed by the Apache Software Foundation
# (http://www.apache.org/).
#
# *****************************************************************************

#================================================================
# Makefile include to store global details shared by all components.
#================================================================


#----------------------------------------------------------------
# System vars
#----------------------------------------------------------------

SHELL = /bin/sh
OS    = $(shell uname)
HOSTTYPE = $(shell echo $${HOSTTYPE})

ifeq ($(HOSTTYPE),unknown)
	HOSTTYPE=itanium
	HOSTTYPE=EM64T
endif

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
ifeq ($(OS),CYGWIN_NT-5.1)
	OSTYPE = CYGWIN
endif
ifeq ($(OS),OSF1)
	OSTYPE = TRU64
endif

# make sure XERCESROOT is properly set!
# env XERCESCROOT="`pwd`/../xerces-c-src2_1_0"
# setenv XERCESCROOT `pwd`/../xerces-c-src2_1_0
export XERCESCROOT `pwd`/../xerces-c-src2_1_0
#----------------------------------------------------------------
#  Default platform specific project build options
#----------------------------------------------------------------

# -------------- Different platform build types --------------

# Set Default Build Type if none is specified
# ( Build type profiles defined below )
ifeq ($(BUILD_TYPE),)
#	BUILD_TYPE = gnu_debug_static
#	BUILD_TYPE = gnu_debug
#	BUILD_TYPE = gnu_fast
#	BUILD_TYPE = gnu_fast_static
#	BUILD_TYPE = gnu_mpi_fast_static
#	BUILD_TYPE = gnu_mpi_fast
#	BUILD_TYPE = gnu_mpi_debug_static
#	BUILD_TYPE = gnu_mpi_debug
#	BUILD_TYPE = gnu64_fast
#	BUILD_TYPE = tcs_mpi_debug
#	BUILD_TYPE = tcs_mpi_fast
#	BUILD_TYPE = tcs_mpi_fast_debug
#	BUILD_TYPE = absoft_mpi_fast
#	BUILD_TYPE = absoft_mpi_fast_p4
#	BUILD_TYPE = absoft_mpi_fast_static
#	BUILD_TYPE = absoft_gcc296_mpi_fast_static
#	BUILD_TYPE = absoft_fast_static
#	BUILD_TYPE = absoft_debug_static
#	BUILD_TYPE = absoft_mpi_debug
#	BUILD_TYPE = absoft_debug
#	BUILD_TYPE = intel_p4_mpi_static
#	BUILD_TYPE = intel_p4_mpi
#	BUILD_TYPE = intel7_p4
#	BUILD_TYPE = intel8_p4_mpi_static_debug
#	BUILD_TYPE = intel8_p4_mpi_static
#	BUILD_TYPE = intel8_p4_mpi_static_infband
#	BUILD_TYPE = intel8_p4_mpi_intelfast
#	BUILD_TYPE = intel9_64_inteldebug
#	BUILD_TYPE = intel9_64_intelfast
#	BUILD_TYPE = intel9_64_intelfast_91045
#	BUILD_TYPE = intel9_64_mpi_intelfast
# 	BUILD_TYPE = intel9_64_mpi_inteldebug
#	BUILD_TYPE = intel9_64_mpi_intelfast_instr
#	BUILD_TYPE = intel10_64_mpi_intelfast
#       BUILD_TYPE = intel15_64_mpi_intelfast
#       BUILD_TYPE = intel19_64_mpi_intelfast
#       BUILD_TYPE = intel19_64_mpi_intelfast_gadi
	BUILD_TYPE = intel20_64_openmpi_gadi
#	BUILD_TYPE = intel_mpi_static
#	BUILD_TYPE = intel9_32_mpi_intelfast
#	BUILD_TYPE = intel_piii_mpi_static
#	BUILD_TYPE = intel_piii_static
#	BUILD_TYPE = itanium2_gnu_debug
#	BUILD_TYPE = itanium2_gnu_fast
#	BUILD_TYPE = opteron_pgi_fast
#	BUILD_TYPE = crayxt_backend_pgi_fast
#	BUILD_TYPE = crayxt_frontend_pgi_fast
#	BUILD_TYPE = crayxt_backend_pgi_debug
#	BUILD_TYPE = crayxt_frontend_pgi_debug
#	BUILD_TYPE = crayxt_backend_pgi_fast_profile
#	BUILD_TYPE = crayxt_frontend_pgi_fast_profile
#	BUILD_TYPE = itanium2_intel8_static
#	BUILD_TYPE = itanium2_intel7_mpi_static
#	BUILD_TYPE = itanium2_intel8_mpi_static
#	BUILD_TYPE = itanium2_mpichg2_static
#	BUILD_TYPE = macosx_gnu_debug
#	BUILD_TYPE = macosx_intel_g95_debug
#	BUILD_TYPE = macosx_intel_g95_fast
#	BUILD_TYPE = macosx_intel_gfortran_fast
#	BUILD_TYPE = macosx_intel_gfortran_local
#	BUILD_TYPE = macosx_intel_gfortran_mpi_fast
#	BUILD_TYPE = macosx_powerpc_gfortran_mpi_fast
#	BUILD_TYPE = macosx_gnu_debug_64bit
#	BUILD_TYPE = macosx_absoft_mpi_fast
#	BUILD_TYPE = macosx_absoft_fast
#	BUILD_TYPE = macosx_absoft_debug
#	BUILD_TYPE = macosx_absoft_mwerks_debug
#	BUILD_TYPE = macosx_ibm_fast
#	BUILD_TYPE = macosx_ibm_debug
#	BUILD_TYPE = macosx_ibm_fast
#	BUILD_TYPE = aix_mpi_static
#  	BUILD_TYPE = aix_mpi_static64
# 	BUILD_TYPE = pathscale_mpi_fast
#	BUILD_TYPE = cygwin_gnu_debug
#	BUILD_TYPE = ibm64_mpich_instr
#	BUILD_TYPE = ranger_intel10_openmpi
#	BUILD_TYPE = ranger_intel10_mvapich2
#	BUILD_TYPE = ranger_pgi_mvapich2
#       BUILD_TYPE = nanohub_gcc_mpi
#	BUILD_TYPE = nanohub_intel10_mpi
#	BUILD_TYPE = steele_intel10_64_mpi_intelfast
#	BUILD_TYPE = steele_intel11_64_mpi_intelfast
endif

# Set Postprocess Type
ifeq ($(POSTPROCESS),)
#	POSTPROCESS = Off 
	POSTPROCESS = Self_Standing
#	POSTPROCESS = Embedded
endif

# Note: If you are running on a processor that does not support MPI, You should not uncomment USE_TRACEMIN_PAR
# It had to be synchronized with MPI3d options
# Include tracemin / TRACEMIN?
TRACEMIN_MACRO_DEF=
#TRACEMIN_MACRO_DEF=-DUSE_TRACEMIN_PAR
#TRACEMIN_MACRO_DEF=-DUSE_TRACEMIN_SER

# Include python?
PYTHON_MACRO_DEF=-DDO_NOT_USE_PYTHON
#PYTHON_MACRO_DEF=-DUSE_PYTHON

# Have single precision Hamiltonian storage
 HAM_STORE_DEF=-DHamMemorySave
#HAM_STORE_DEF=-UHamMemorySave

# Note: If you want to use MKL library instead of LAPACK/BLAS included in this software, please switch on this option
#MKL_MACRO_DEF=
MKL_MACRO_DEF=-DUSE_MKL

#----------------CHOICE OF COMPILERS------------------------------
#------------THE DECLARATIONS BELOW JUST ESTABLISH MAKEFILE VARIABLES
#------------IF NEEDED, THEY SHOULD BE DEFINED IN THE BUILD OPTIONS
#----- IF SOME OF THEM ARE LEFT UNDEFINED, THEY WILL BE ASSIGNED
#----- A DEFAULT SETTING AT THE END OF THIS FILE 
#-----------(BELOW ALL SPECIFIC BUILD TYPE DEFINITIONS)

#---------------AGAIN I TELL YOU, HANDS OFF!

XERCES_CC=
XERCES_CCMKD=
XERCES_CCOPT=
XERCES_LD=
XERCES_LDOPT=

XMLtoC_CC=
XMLtoC_CCMKD=
XMLtoC_CCOPT=
XMLtoC_LD=
XMLtoC_LDOPT=

FORTRAN90=
SCTRACEMIN_OPTS=
PCTRACEMIN_OPTS=
SQTRACEMIN_OPTS=
PQTRACEMIN_OPTS=
# -------------- Common to all platform build types --------------

MAKE = make
AR = ar
AR_OPTIONS = rcvs
RANLIB = ranlib


# The platform and build type for which lib and obj subdirectories will be organized in
PLATFORM = $(HOSTTYPE)_$(BUILD_TYPE)


#----------------------------------------------------------------
# Application components
#----------------------------------------------------------------

# The application root which contains the various project component libraries.
APP_ROOT = $(shell pwd)/..

# Project components used for this application
NEMO3D_ROOT = $(APP_ROOT)/nemo3d
XMLtoC_ROOT = $(APP_ROOT)/XMLtoC
EIGEN_ROOT  = $(APP_ROOT)/eigen
NML_ROOT    = $(APP_ROOT)/nml
POSTPROCESS_ROOT    = $(APP_ROOT)/PostNEMO3D
TEST_ROOT   = $(APP_ROOT)/nemo3d/example/TestRun
XERCESCROOT = $(APP_ROOT)/xerces-c-src2_1_0
#SVMTL_ROOT    = $(APP_ROOT)/svmtl

# XML file used to generate I/O classes for NEMO3D
# NOTE: The XML file must be relative from "$(XMLtoC_ROOT)/xml"
ioNemo3D_XML0 = ioNemo3D.xml


#----------------------------------------------------------------
# Application dependent libraries
#----------------------------------------------------------------

#LIB_SVMTL = $(SVMTL_ROOT)/lib/$(PLATFORM)/libsvmtl.a
LIB_NML    = $(NML_ROOT)/lib/$(PLATFORM)/libnml.a
LIB_EIGEN  = $(EIGEN_ROOT)/lib/$(PLATFORM)/libeigen.a
LIB_XMLtoC = $(XMLtoC_ROOT)/apps/ioNemo3D/lib/$(PLATFORM)/libioNemo3D-$(PLATFORM).a


#----------------------------------------------------------------
# External library root & include paths
#----------------------------------------------------------------


# LAPACK : Linear Algebra Package Fortran Library
# NOTE: have to define LAPACK_ROOT as relative because the LAPACK/BLAS uses the library definitions as relative paths.
LAPACK_ROOT  = ../LAPACK

# BLAS : Basic Linear Algebra Fortran Library
BLAS_ROOT    = $(LAPACK_ROOT)/BLAS

# ARPACK: ARPACK is a collection of Fortran77 subroutines designed to solve large scale eigenvalue problems.
ARPACK_ROOT  = $(APP_ROOT)/ARPACK

# PARPACK: Parallel ARPACK
PARPACK_ROOT  = $(ARPACK_ROOT)/PARPACK

##################################################################################
# TRACEMIN: Sequential and Parallel Implementation of Trace Minimization Algorithm 
TRACEMIN_ROOT  = $(APP_ROOT)/TRACEMIN
##################################################################################

# XERCES-C++ : XML Parser
XERCESCROOT = $(APP_ROOT)/xerces-c-src2_1_0
XERCESC_EXPORT_INCL = \
	-I$(XERCESCROOT)/include \
	-I$(XERCESCROOT)/include/xercesc

# PYTHON : Scripting Interpreter
ifeq ($(PYTHON_MACRO_DEF),-DUSE_PYTHON)
	PYTHON_ROOT = /usr/lib/python2.2
	PYTHON_EXPORT_INCL = -I/usr/include/python2.2
endif

#----------------------------------------------------------------
# External libraries
#----------------------------------------------------------------

MKL=$(INTEL_MKL_BASE)
ifeq ($(MKL_MACRO_DEF), )
   # BLAS
   LIB_BLAS    = $(LAPACK_ROOT)/libblas_$(PLATFORM).a
   # LAPACK
   LIB_LAPACK  = $(LAPACK_ROOT)/liblapack_$(PLATFORM).a
   MKL_INCLUDE = -mkl 
   MKL_LINK = -mkl
else
   LIB_BLAS = 
   LIB_LAPACK = 
   MKL_INCLUDE = -I${MKL}/include
   MKL_LINK =  -L${MKL}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lifport -lifcore -lirc -lsvml -limf -lintlc -lmpi_mpifh
endif

# ARPACK
LIB_ARPACK  = $(ARPACK_ROOT)/libarpack_$(PLATFORM).a

# PARPACK
LIB_PARPACK  = $(ARPACK_ROOT)/libparpack_$(PLATFORM).a

############################################################
# TRACEMIN
ifeq ($(TRACEMIN_MACRO_DEF),-DUSE_TRACEMIN_PAR)
	LIB_PCTRACEMIN = $(TRACEMIN_ROOT)/Lib/libpctracemin_$(PLATFORM).a
	LIB_PQTRACEMIN = $(TRACEMIN_ROOT)/Lib/libpqtracemin_$(PLATFORM).a
	LIB_TRACEMIN = $(LIB_PCTRACEMIN) $(LIB_PQTRACEMIN)
endif
ifeq ($(TRACEMIN_MACRO_DEF),-DUSE_TRACEMIN_SER)
	LIB_SCTRACEMIN = $(TRACEMIN_ROOT)/Lib/libsctracemin_$(PLATFORM).a
	LIB_SQTRACEMIN = $(TRACEMIN_ROOT)/Lib/libsqtracemin_$(PLATFORM).a
	LIB_TRACEMIN = $(LIB_SCTRACEMIN) $(LIB_SQTRACEMIN)
endif
############################################################

# PYTHON
ifeq ($(PYTHON_MACRO_DEF),-DUSE_PYTHON)
	LIB_PYTHON  = $(PYTHON_ROOT)/config/libpython2.2.a
endif

LIB_XERCESC = $(XERCESCROOT)/lib/libxerces-c_$(PLATFORM).a


#----------------------------------------------------------------
# NEMO3D Complete Library
#----------------------------------------------------------------

# This is the complete library of NEMO3D. It is really a merge
# of all of NEMO3D's dependent libraries
# i.e. blas, lapack, arpack, parpack, nml, eigen, XMLtoC
# Xerces-C++, and the sub-libraries of nemo3d.

LIB_NEMO3D = $(APP_ROOT)/build/libnemo3d_$(PLATFORM).a


#----------------------------------------------------------------
# External libraries search path
#----------------------------------------------------------------

# Search paths to shared objects for compiler and link editor.
# Specifies directories containing dynamic libraries that are NOT
# in the usual searched library paths.
ifeq ($(OS),Darwin)
	LIB_PATHS = -L$(XERCESCROOT)/lib
else
	LIB_PATHS = -L$(XERCESCROOT)/lib -Wl,-rpath,$(XERCESCROOT)/lib
endif


#----------------------------------------------------------------
#  Platform project build profiles ( Default set above )
#----------------------------------------------------------------

# Compaq TCS system debug
ifeq ($(BUILD_TYPE),tcs_mpi_debug)
        MAKE = gmake

        CCMKD = g++ -MM
        # Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
        CC = g++
        CCOPT = -g -c -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) -Dmem_64BIT
        ifeq ($(POSTPROCESS), Embedded)
                CCOPT = -g -c -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) -Dmem_64BIT
        endif
        LD = g++
        LDOPT = 

        MISCLIB0   = -lfmpi -lmpi -lcxml -lUfor -lfor -lrt -lpthread -lm
        MISCINCL0  =

        #XMLtoC specific options
        XMLtoC_CC=g++
        XMLtoC_CCOPT=$(CCOPT)
        XMLtoC_LD=$(XMLtoC_CC)

        # NML-specific options
        NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT -Dmem_64BIT
        NML_OPTIONS = -Wall -ansi -pedantic

        FORTRAN = f90

        # LAPACK-specific options
        LAPACK_NOOPT =
        ifeq ($(OS),Darwin)
                LAPACK_OPTS   = -fast -N11 -g
                LAPACK_LDOPTS = -fast
        else
                LAPACK_OPTS   = -fast -funderscoring -funroll-all-loops -g
                LAPACK_LDOPTS = -fast
        endif

        ifeq ($(OS),Darwin)
                ARPACK_OPTS   = -fast -N11 -g
        else
                ARPACK_OPTS   = -fast -funderscoring -funroll-all-loops -g
        endif

endif

# Compaq TCS system fast
ifeq ($(BUILD_TYPE),tcs_mpi_fast)
        MAKE = gmake

        CCMKD = g++ -MM
        # Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
        CC = g++
        CCOPT = -c -O3 -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) -Dmem_64BIT
        ifeq ($(POSTPROCESS), Embedded)
                CCOPT = -c -O3 -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) -Dmem_64BIT
        endif
        LD = g++
        LDOPT = 

        MISCLIB0   = -lfmpi -lmpi -lcxml -lUfor -lfor -lrt -lpthread -lm
        MISCINCL0  =

        #XMLtoC specific options
        XMLtoC_CC=g++
        XMLtoC_CCOPT=$(CCOPT)
        XMLtoC_LD=$(XMLtoC_CC)

        # NML-specific options
        NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT -Dmem_64BIT
        NML_OPTIONS = -O3 -Wall -ansi -pedantic

        FORTRAN = f90

        # LAPACK-specific options
        LAPACK_NOOPT =
        ifeq ($(OS),Darwin)
                LAPACK_OPTS   = -fast -N11
                LAPACK_LDOPTS = -fast
        else
                LAPACK_OPTS   = -fast -funderscoring -funroll-all-loops
                LAPACK_LDOPTS = -fast
        endif

        ifeq ($(OS),Darwin)
                ARPACK_OPTS   = -fast -N11
        else
                ARPACK_OPTS   = -fast -funderscoring -funroll-all-loops
        endif

endif

# Compaq TCS system fast debug
ifeq ($(BUILD_TYPE),tcs_mpi_fast_debug)
        MAKE = gmake

        CCMKD = g++ -MM
        # Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
        CC = g++
        CCOPT = -g -c -O3 -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) -Dmem_64BIT
        ifeq ($(POSTPROCESS), Embedded)
                CCOPT = -g -c -O3 -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) -Dmem_64BIT
        endif
        LD = g++
        LDOPT = 

        MISCLIB0   = -lfmpi -lmpi -lcxml -lUfor -lfor -lrt -lpthread -lm
        MISCINCL0  =

        #XMLtoC specific options
        XMLtoC_CC=g++
        XMLtoC_CCOPT=$(CCOPT)
        XMLtoC_LD=$(XMLtoC_CC)

        # NML-specific options
        NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT -Dmem_64BIT
        NML_OPTIONS = -O3 -Wall -ansi -pedantic

        FORTRAN = f90

        # LAPACK-specific options
        LAPACK_NOOPT =
        ifeq ($(OS),Darwin)
                LAPACK_OPTS   = -fast -N11 -g
                LAPACK_LDOPTS = -fast
        else
                LAPACK_OPTS   = -fast -funderscoring -funroll-all-loops -g
                LAPACK_LDOPTS = -fast
        endif

        ifeq ($(OS),Darwin)
                ARPACK_OPTS   = -fast -N11 -g
        else
                ARPACK_OPTS   = -fast -funderscoring -funroll-all-loops -g
        endif

endif


# GNU Debug Static
ifeq ($(BUILD_TYPE),gnu_debug_static)

	CCMKD  = g++ -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY -DSANS_PY
	CC     = g++
	CCOPT  = -c -g -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE  -DFAKE_MPI -DMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) 
	LIB_PARPACK =
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT  = -c -g -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE  -DFAKE_MPI -DMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	endif

	LD     = g++
	LDOPT  = -g -Wl,-static

	MISCLIB0   = -ldl -lpthread -lm -lg2c
	MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
	NML_OPTIONS = -g -Wall -ansi -pedantic

	# LAPACK-specific options
	ifeq ($(OS),Darwin)
		FORTRAN = f77
		LAPACK_NOOPT =
		LAPACK_OPTS   = -g
		LAPACK_LDOPTS = -g
		ARPACK_OPTS   = -g
	else
		FORTRAN = g77
		LAPACK_NOOPT =
		LAPACK_OPTS   = -g -funderscoring
		LAPACK_LDOPTS = -g -funderscoring
		ARPACK_OPTS   = -g -funderscoring
	endif

	# If -DNO_PARPACK is used in CCOPT, then we need to make sure that the PARPACK library doesn't get linked in
	LIB_PARPACK =

endif






# GNU Debug Static
ifeq ($(BUILD_TYPE),gnu_debug)

        CCMKD  = g++ -MM
        # Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
        CC     = g++
        CCOPT  = -c -g -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE  -DFAKE_MPI -UMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) 
	LIB_PARPACK =
        ifeq ($(POSTPROCESS), Embedded)
                CCOPT  = -c -g -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE  -DFAKE_MPI -UMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
        endif
        LD     = g++
        LDOPT  = -g -Wl
        LDOPT  = -g 

        MISCLIB0   = -ldl -lpthread -lm -lg2c
        MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

        # NML-specific options
        NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
        NML_OPTIONS = -g -Wall -ansi -pedantic
        NML_OPTIONS = -g -Wall -ansi 

        # LAPACK-specific options
        ifeq ($(OS),Darwin)
                FORTRAN = f77
                LAPACK_NOOPT =
                LAPACK_OPTS   = -g
                LAPACK_LDOPTS = -g
                ARPACK_OPTS   = -g
        else
                FORTRAN = g77
                LAPACK_NOOPT =
                LAPACK_OPTS   = -g -funderscoring
                LAPACK_LDOPTS = -g -funderscoring
                ARPACK_OPTS   = -g -funderscoring
        endif

        # If -DNO_PARPACK is used in CCOPT, then we need to make sure that the PARPACK library doesn't get linked in
        LIB_PARPACK =

endif

# GNU Fast
ifeq ($(BUILD_TYPE),gnu_fast)

	CCMKD  = g++ -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC     = g++
	CCOPT  = -c -O3 -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DFAKE_MPI -DMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) 
	LIB_PARPACK =
        ifeq ($(POSTPROCESS), Embedded)
		CCOPT  = -c -O3 -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DFAKE_MPI -DMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	LIB_PARPACK =
	endif
	LD     = g++
	LDOPT  = 

	MISCLIB0   = -ldl -lpthread -lm -lg2c
	MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
	NML_OPTIONS = -O3 -Wall -ansi -pedantic
	NML_OPTIONS = -O3 -Wall -ansi 

	# LAPACK/ARPARCK options
	ifeq ($(OS),Darwin)
		FORTRAN       = f77
		LAPACK_NOOPT  =
		LAPACK_OPTS   = -O2
		LAPACK_LDOPTS = -O2
		ARPACK_OPTS   = -O2
	else
		FORTRAN       = g77
		LAPACK_NOOPT  =
		LAPACK_OPTS   = -O2 -funderscoring -funroll-all-loops
		LAPACK_LDOPTS = -O2 -funderscoring -funroll-all-loops
		ARPACK_OPTS   = -O2 -funderscoring -funroll-all-loops
	endif

	# If -DNO_PARPACK is used in CCOPT, then we need to make sure that the PARPACK library doesn't get linked in
	LIB_PARPACK =

endif



# GNU Fast
ifeq ($(BUILD_TYPE),gnu_fast_static)

	CCMKD  = g++ -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC     = g++
	CCOPT  = -c -O3 -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DFAKE_MPI -UMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) 
	LIB_PARPACK =
        ifeq ($(POSTPROCESS), Embedded)
		CCOPT  = -c -O3 -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE -UFAKE_MPI -DMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	LIB_PARPACK =
	endif
	LD     = g++
	LDOPT  = -static

	MISCLIB0   = -ldl -lpthread -lm -lg2c
	MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
	NML_OPTIONS = -O3 -Wall -ansi -pedantic
	NML_OPTIONS = -O3 -Wall -ansi 

	# LAPACK/ARPARCK options
	ifeq ($(OS),Darwin)
		FORTRAN       = f77
		LAPACK_NOOPT  =
		LAPACK_OPTS   = -O2
		LAPACK_LDOPTS = -O2
		ARPACK_OPTS   = -O2
	else
		FORTRAN       = g77
		LAPACK_NOOPT  =
		LAPACK_OPTS   = -O2 -funderscoring -funroll-all-loops
		LAPACK_LDOPTS = -O2 -funderscoring -funroll-all-loops
		ARPACK_OPTS   = -O2 -funderscoring -funroll-all-loops
	endif

	# If -DNO_PARPACK is used in CCOPT, then we need to make sure that the PARPACK library doesn't get linked in
	LIB_PARPACK =

endif


# MPICH GNU MPI Fast
ifeq ($(BUILD_TYPE),gnu_mpi_fast)

	ifeq ($(MPI_HOME),)
		@echo "Error - Undefined environment variable, MPI_HOME"
		exit
	endif

	CCMKD = g++ -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC = /home/packages/mpich-1.2.4/gnu/bin/mpiCC
        CC = $(MPI_HOME)/bin/mpiCC
	CCOPT = -c -O3 -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) 
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT = -c -O3 -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	endif
	LD = /home/packages/mpich-1.2.4/gnu/bin/mpiCC
        LD = $(MPI_HOME)/bin/mpiCC
	LDOPT  = 

	MISCLIB0   = -ldl -lpthread -lm -lg2c
	MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
	NML_OPTIONS = -O3 -Wall -ansi -pedantic
	NML_OPTIONS = -O3 -Wall -ansi 

	FORTRAN = /home/packages/mpich-1.2.3/gnu/bin/mpif77
        FORTRAN = $(MPI_HOME)/bin/mpif77

	# LAPACK-specific options
	LAPACK_NOOPT =
	ifeq ($(OS),Darwin)
		LAPACK_OPTS   = -O2 -N11
		LAPACK_LDOPTS = -O2 -N11
	else
		LAPACK_OPTS   = -O2 -funderscoring -funroll-all-loops
		LAPACK_LDOPTS = -O2 -funderscoring -funroll-all-loops
	endif

        ifeq ($(OS),Darwin)
		ARPACK_OPTS   = -O2 -N11
	else
		ARPACK_OPTS   = -O2 -funderscoring -funroll-all-loops
	endif
endif


# MPICH GNU MPI Debug
ifeq ($(BUILD_TYPE),gnu_mpi_debug)

        ifeq ($(MPI_HOME),)
                @echo "Error - Undefined environment variable, MPI_HOME"
                exit
        endif

        CCMKD = g++ -MM
        # Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
        CC = $(MPI_HOME)/bin/mpiCC
        CCOPT = -c -g -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING
        ifeq ($(POSTPROCESS), Embedded)
                CCOPT = -c -g -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING -DEmbedded_Postprocess
        endif
        LD = $(MPI_HOME)/bin/mpiCC
        LDOPT  =

        MISCLIB0   = -ldl -lpthread -lm -lg2c
        MISCINCL0  =

        #XMLtoC specific options
        XMLtoC_CCOPT=$(CCOPT)

        # NML-specific options
        NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
        NML_OPTIONS = -g -Wall -ansi -pedantic

        FORTRAN = $(MPI_HOME)/bin/mpif77

        # LAPACK-specific options
        LAPACK_NOOPT =
        ifeq ($(OS),Darwin)
                LAPACK_OPTS   = -g -N11
                LAPACK_LDOPTS = -g -N11
        else
                LAPACK_OPTS   = -g -funderscoring -funroll-all-loops
                LAPACK_LDOPTS = -g -funderscoring -funroll-all-loops
        endif

        ifeq ($(OS),Darwin)
                ARPACK_OPTS   = -g -N11
        else
                ARPACK_OPTS   = -g -funderscoring -funroll-all-loops
        endif
endif


# MPICH GNU MPI Debug
ifeq ($(BUILD_TYPE),gnu_mpi_debug_static)

	ifeq ($(MPI_HOME),)
		@echo "Error - Undefined environment variable, MPI_HOME"
		exit
	endif

	CCMKD = g++ -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC = /home/packages/mpich-1.2.4/gnu/bin/mpiCC
	CC = mpiCC
        CC = $(MPI_HOME)/bin/mpiCC
	CCOPT = -c -g -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT = -c -g -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING -DEmbedded_Postprocess
	endif
	LD = /home/packages/mpich-1.2.4/gnu/bin/mpiCC
	LD = mpiCC
        LD = $(MPI_HOME)/bin/mpiCC -static
	LDOPT  = 

	MISCLIB0   = -ldl -lpthread -lm -lg2c
	MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
	NML_OPTIONS = -g -Wall -ansi -pedantic

	FORTRAN = /home/packages/mpich-1.2.3/gnu/bin/mpif77
        FORTRAN = $(MPI_HOME)/bin/mpif77

	# LAPACK-specific options
	LAPACK_NOOPT =
	ifeq ($(OS),Darwin)
		LAPACK_OPTS   = -g -N11
		LAPACK_LDOPTS = -g -N11
	else
		LAPACK_OPTS   = -g -funderscoring -funroll-all-loops
		LAPACK_LDOPTS = -g -funderscoring -funroll-all-loops
	endif

        ifeq ($(OS),Darwin)
		ARPACK_OPTS   = -g -N11
	else
		ARPACK_OPTS   = -g -funderscoring -funroll-all-loops
	endif
endif


# MPICH GNU MPI Fast
ifeq ($(BUILD_TYPE),gnu_mpi_fast_static)

        ifeq ($(MPI_HOME),)
                @echo "Error - Undefined environment variable, MPI_HOME"
                exit
        endif

        CCMKD = g++ -MM
        # Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
        CC = /home/packages/mpich-1.2.4/gnu/bin/mpiCC
        CC = $(MPI_HOME)/bin/mpiCC
        CCOPT = -c -O3 -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) 
        ifeq ($(POSTPROCESS), Embedded)
                CCOPT = -c -O3 -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
        endif
        LD = /home/packages/mpich-1.2.4/gnu/bin/mpiCC
        LD = $(MPI_HOME)/bin/mpiCC
        LDOPT  = -g -Wl,-static

        MISCLIB0   = -ldl -lpthread -lm -lg2c
        MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

        # NML-specific options
        NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
        NML_OPTIONS = -O3 -Wall -ansi -pedantic
        NML_OPTIONS = -O3 -Wall -ansi 

        FORTRAN = /home/packages/mpich-1.2.3/gnu/bin/mpif77
        FORTRAN = $(MPI_HOME)/bin/mpif77

        # LAPACK-specific options
        LAPACK_NOOPT =
        ifeq ($(OS),Darwin)
                LAPACK_OPTS   = -O2 -N11
                LAPACK_LDOPTS = -O2 -N11
        else
                LAPACK_OPTS   = -O2 -funderscoring -funroll-all-loops
                LAPACK_LDOPTS = -O2 -funderscoring -funroll-all-loops
        endif

        ifeq ($(OS),Darwin)
                ARPACK_OPTS   = -O2 -N11
        else
                ARPACK_OPTS   = -O2 -funderscoring -funroll-all-loops
        endif
endif


# GNU Fast
ifeq ($(BUILD_TYPE),gnu64_fast)

	CCMKD  = g++-4.1 -m64  -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC     = g++-4.1 -m64 
	CCOPT  = -c -O3 -Wall -D$(OSTYPE) -Dmem_64BIT -DFORTRAN_UNDERSCORE -DFAKE_MPI -DMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) 
	LIB_PARPACK =
        ifeq ($(POSTPROCESS), Embedded)
		CCOPT  = -c -O3 -Wall -D$(OSTYPE) -Dmem_64BIT -DFORTRAN_UNDERSCORE -DFAKE_MPI -DMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	LIB_PARPACK =
	endif
	LD     = g++-4.1 -m64
	LDOPT  = 

	MISCLIB0   = -ldl -lpthread -lm -lg2c
	MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -Dmem_64BIT -DNML_ERRANT  -DNML_MEM_EXT
	NML_OPTIONS = -O3 -Wall -ansi -pedantic
	NML_OPTIONS = -O3 -Wall -ansi 

	# LAPACK/ARPARCK options
	FORTRAN       = g77 -m64
	LAPACK_NOOPT  =
	LAPACK_OPTS   = -O2 -funderscoring -funroll-all-loops
	LAPACK_LDOPTS = -O2 -funderscoring -funroll-all-loops
	ARPACK_OPTS   = -O2 -funderscoring -funroll-all-loops

	# If -DNO_PARPACK is used in CCOPT, then we need to make sure that the PARPACK library doesn't get linked in
	LIB_PARPACK =

endif



# MPICH Absoft MPI Fast
ifeq ($(BUILD_TYPE),absoft_mpi_fast)

	ifeq ($(ABSOFT),)
#		@echo "Error - Undefined environment variable, ABSOFT"
	endif

	ifeq ($(MPI_HOME),)
#		@echo "Error - Undefined environment variable, MPI_HOME"
	endif

	CCMKD = g++ -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC = $(MPI_HOME)/bin/mpiCC
	CCOPT = -c -O3 -Wall -D$(OSTYPE) -DFORTRAN_ALLCAPS -DMPI3d -DMPI_TIMING  $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) 
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT = -c -O3 -Wall -D$(OSTYPE) -DFORTRAN_ALLCAPS -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	endif
	LD = $(MPI_HOME)/bin/mpiCC
	LDOPT  =

	MISCLIB0   = $(ABSOFT)/lib/libU77.a $(ABSOFT)/lib/libf77math.a $(ABSOFT)/lib/libf90math.a $(ABSOFT)/lib/libfio.a -ldl -lpthread -lm -lg2c
	MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
	NML_OPTIONS = -O3 -Wall -ansi -pedantic

	FORTRAN = $(MPI_HOME)/bin/mpif77

	# LAPACK-specific options
	LAPACK_NOOPT =
	ifeq ($(OS),Darwin)
		LAPACK_OPTS = -O2 -N11
	else
		# Use with FORTRAN_UNDERSCORE
		#LAPACK_OPTS = -O2 -cpu:host -f -N22 -B108
		# Use with FORTRAN_ALLCAPS
		LAPACK_OPTS   = -O2 -cpu:host
		LAPACK_LDOPTS = -O2 -cpu:host
	endif

	ifeq ($(OS),Darwin)
		ARPACK_OPTS   = -O2 -N11
	else
		# Use with FORTRAN_UNDERSCORE
		#ARPACK_OPTS  = -O2 -cpu:host -f -N22 -B108
		# Use with FORTRAN_ALLCAPS
		ARPACK_OPTS   = -O2 -cpu:host
	endif

endif



# MPICH Absoft MPI Fast, for GCC 3.x, Optimized for Pentium4
ifeq ($(BUILD_TYPE),absoft_mpi_fast_p4)

	ifeq ($(ABSOFT),)
#		@echo "Error - Undefined environment variable, ABSOFT"
	endif

	ifeq ($(MPI_HOME),)
#		@echo "Error - Undefined environment variable, MPI_HOME"
	endif

	CCMKD = g++ -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC = $(MPI_HOME)/bin/mpiCC
	CCOPT = -c -O2 -march=pentium4 -mcpu=pentium4 -mfpmath=sse -msse2 -mmmx -Wall -D$(OSTYPE) -DFORTRAN_ALLCAPS -DMPI3d  $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DMPI_TIMING $(HAM_STORE_DEF) 
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT = -c -O2 -march=pentium4 -mcpu=pentium4 -mfpmath=sse -msse2 -mmmx -Wall -D$(OSTYPE) -DFORTRAN_ALLCAPS -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	endif
	LD = $(MPI_HOME)/bin/mpiCC
	LDOPT  =

	MISCLIB0   = $(ABSOFT)/lib/libU77.a $(ABSOFT)/lib/libf77math.a $(ABSOFT)/lib/libf90math.a $(ABSOFT)/lib/libfio.a -ldl -lpthread -lm -lg2c
	MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
	NML_OPTIONS = -O2 -march=pentium4 -mcpu=pentium4 -mfpmath=sse -msse2 -mmmx -Wall -ansi -pedantic

	FORTRAN = $(MPI_HOME)/bin/mpif77

	# LAPACK-specific options
	LAPACK_NOOPT =
	ifeq ($(OS),Darwin)
		LAPACK_OPTS = -O2 -N11
	else
		# Use with FORTRAN_UNDERSCORE
		#LAPACK_OPTS = -O2 -cpu:p7 -f -N22 -B108
		# Use with FORTRAN_ALLCAPS
		LAPACK_OPTS   = -O2 -cpu:p7
		LAPACK_LDOPTS = -O2 -cpu:p7
	endif

	ifeq ($(OS),Darwin)
		ARPACK_OPTS   = -O2 -N11
	else
		# Use with FORTRAN_UNDERSCORE
		#ARPACK_OPTS  = -O2 -cpu:p7 -f -N22 -B108
		# Use with FORTRAN_ALLCAPS
		ARPACK_OPTS   = -O2 -cpu:p7
	endif

endif



# MPICH Absoft MPI Debug
ifeq ($(BUILD_TYPE),absoft_mpi_debug)

	ifeq ($(ABSOFT),)
#		$(error "Error - Undefined environment variable, ABSOFT")
	endif

	ifeq ($(MPI_HOME),)
#		$(error "Error - Undefined environment variable, MPI_HOME")
	endif

	CCMKD = g++ -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC = $(MPI_HOME)/bin/mpiCC
	CCOPT = -c -g -Wall -D$(OSTYPE) -DFORTRAN_ALLCAPS -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) 
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT = -c -g -Wall -D$(OSTYPE) -DFORTRAN_ALLCAPS -DMPI3d -DMPI_TIMING -DEmbedded_Postprocess
	endif
	LD = $(MPI_HOME)/bin/mpiCC
	LDOPT  = -g

	MISCLIB0   = $(ABSOFT)/lib/libU77.a $(ABSOFT)/lib/libf77math.a $(ABSOFT)/lib/libf90math.a $(ABSOFT)/lib/libfio.a -ldl -lpthread -lm -lg2c
	MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
	NML_OPTIONS = -g -Wall -ansi -pedantic

	FORTRAN = $(MPI_HOME)/bin/mpif77

	# LAPACK-specific options
	LAPACK_NOOPT =
	ifeq ($(OS),Darwin)
		LAPACK_OPTS   = -g -N11
	else
		# Use with FORTRAN_UNDERSCORE
		#LAPACK_OPTS  = -g -cpu:host -f -N22 -B108
		# Use with FORTRAN_ALLCAPS
		LAPACK_OPTS   = -g -cpu:host
		LAPACK_LDOPTS = -g -cpu:host
	endif

	ifeq ($(OS),Darwin)
		ARPACK_OPTS   = -g -N11
	else
		# Use with FORTRAN_UNDERSCORE
		#ARPACK_OPTS  = -g -cpu:host -f -N22 -B108
		# Use with FORTRAN_ALLCAPS
		ARPACK_OPTS   = -g -cpu:host
	endif

endif



# MPICH Absoft MPI Fast Static
ifeq ($(BUILD_TYPE),absoft_mpi_fast_static)

	ifeq ($(ABSOFT),)
		#@ echo "Error - Undefined environment variable, ABSOFT"
	endif

	ifeq ($(MPI_HOME),)
		#@ echo "Error - Undefined environment variable, MPI_HOME"
	endif

	CCMKD = g++ -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC = $(MPI_HOME)/bin/mpiCC
	CCOPT = -c -O3 -Wall -D$(OSTYPE) -DFORTRAN_ALLCAPS -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF)
	ifeq ($(POSTPROCESS),Embedded)
		CCOPT = -c -O3 -Wall -D$(OSTYPE) -DFORTRAN_ALLCAPS -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	endif
	LD = $(MPI_HOME)/bin/mpiCC
	LDOPT  = -Wl,-static

	MISCLIB0   = $(ABSOFT)/lib/libU77.a $(ABSOFT)/lib/libf77math.a $(ABSOFT)/lib/libf90math.a $(ABSOFT)/lib/libfio.a -ldl -lpthread -lm -lg2c
	MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
	NML_OPTIONS = -O3 -Wall -ansi -pedantic

	FORTRAN = $(MPI_HOME)/bin/mpif77

	# LAPACK-specific options
	LAPACK_NOOPT =
	ifeq ($(OS),Darwin)
		LAPACK_OPTS   = -O2 -N11
		LAPACK_LDOPTS = -O2 -N11
	else
		# Use with FORTRAN_UNDERSCORE
		#LAPACK_OPTS = -O2 -cpu:host -f -N22 -B108
		# Use with FORTRAN_ALLCAPS
		LAPACK_OPTS   = -O2 -cpu:host
		LAPACK_LDOPTS = -O2 -cpu:host
	endif

	ifeq ($(OS),Darwin)
		ARPACK_OPTS   = -O2 -N11
	else
		# Use with FORTRAN_UNDERSCORE
		#ARPACK_OPTS  = -O2 -cpu:host -f -N22 -B108
		# Use with FORTRAN_ALLCAPS
		ARPACK_OPTS   = -O2 -cpu:host
	endif

endif


# MPICH Absoft Fast Static
ifeq ($(BUILD_TYPE),absoft_fast_static)

	ifeq ($(ABSOFT),)
		#@ echo "Error - Undefined environment variable, ABSOFT"
	endif

	CCMKD = g++ -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC = g++
	CCOPT = -c -O3 -Wall -D$(OSTYPE) -DFORTRAN_ALLCAPS -DFAKE_MPI -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) 
	LIB_PARPACK =
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT = -c -O3 -Wall -D$(OSTYPE) -DFORTRAN_ALLCAPS -DFAKE_MPI -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	endif
	LD = g++
	LDOPT  = -Wl,-static

	MISCLIB0   = $(LIB_U77) $(LIB_F77MATH) $(LIB_F90MATH) $(LIB_FIO) -ldl -lpthread -lm -lg2c
	MISCLIB0   = $(ABSOFT)/lib/libU77.a $(ABSOFT)/lib/libf77math.a $(ABSOFT)/lib/libf90math.a $(ABSOFT)/lib/libfio.a  -ldl -lpthread -lm -lg2c
	MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
	NML_OPTIONS = -O3 -Wall -ansi -pedantic

#	FORTRAN = $(MPI_HOME)/bin/mpif77
	FORTRAN = $(ABSOFT)/bin/f90

	# LAPACK-specific options
	LAPACK_NOOPT =
	ifeq ($(OS),Darwin)
		LAPACK_OPTS   = -O2 -N11
		LAPACK_LDOPTS = -O2 -N11
	else
		# Use with FORTRAN_UNDERSCORE
		#LAPACK_OPTS = -O2 -cpu:host -f -N22 -B108
		# Use with FORTRAN_ALLCAPS
		LAPACK_OPTS   = -O2 -cpu:host
		LAPACK_LDOPTS = -O2 -cpu:host
	endif

	ifeq ($(OS),Darwin)
		ARPACK_OPTS   = -O2 -N11
	else
		# Use with FORTRAN_UNDERSCORE
		#ARPACK_OPTS  = -O2 -cpu:host -f -N22 -B108
		# Use with FORTRAN_ALLCAPS
		ARPACK_OPTS   = -O2 -cpu:host
	endif

endif

# MPICH Absoft Debug Static
ifeq ($(BUILD_TYPE),absoft_debug_static)

	ifeq ($(ABSOFT),)
		#@ echo "Error - Undefined environment variable, ABSOFT"
	endif

	CCMKD = g++ -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC = g++
	CCOPT = -c -g -Wall -D$(OSTYPE) -DFORTRAN_ALLCAPS -DFAKE_MPI -DNO_PARPACK
	LIB_PARPACK =
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT = -c -g -Wall -D$(OSTYPE) -DFORTRAN_ALLCAPS -DFAKE_MPI -DNO_PARPACK -DEmbedded_Postprocess
	endif
	LD = g++
	LDOPT  = -Wl,-static

	MISCLIB0   = $(LIB_U77) $(LIB_F77MATH) $(LIB_F90MATH) $(LIB_FIO) -ldl -lpthread -lm -lg2c
	MISCLIB0   = $(ABSOFT)/lib/libU77.a $(ABSOFT)/lib/libf77math.a $(ABSOFT)/lib/libf90math.a $(ABSOFT)/lib/libfio.a  -ldl -lpthread -lm -lg2c
	MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
	NML_OPTIONS = -g -Wall -ansi -pedantic

#	FORTRAN = $(MPI_HOME)/bin/mpif77
	FORTRAN = $(ABSOFT)/bin/f90

	# LAPACK-specific options
	LAPACK_NOOPT =
	ifeq ($(OS),Darwin)
		LAPACK_OPTS   = -g -N11
		LAPACK_LDOPTS = -g -N11
	else
		# Use with FORTRAN_UNDERSCORE
		#LAPACK_OPTS = -g -cpu:host -f -N22 -B108
		# Use with FORTRAN_ALLCAPS
		LAPACK_OPTS   = -g -cpu:host
		LAPACK_LDOPTS = -g -cpu:host
	endif

	ifeq ($(OS),Darwin)
		ARPACK_OPTS   = -g -N11
	else
		# Use with FORTRAN_UNDERSCORE
		#ARPACK_OPTS  = -g -cpu:host -f -N22 -B108
		# Use with FORTRAN_ALLCAPS
		ARPACK_OPTS   = -g -cpu:host
	endif

endif



# this line is a kludge for gcc 2.96
# MPICH Absoft MPI Fast Static
ifeq ($(BUILD_TYPE),absoft_gcc296_mpi_fast_static)

	ifeq ($(ABSOFT),)
		#@ echo "Error - Undefined environment variable, ABSOFT"
	endif

	ifeq ($(MPI_HOME),)
		#@ echo "Error - Undefined environment variable, MPI_HOME"
	endif

	CCMKD = g++ -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC = $(MPI_HOME)/bin/mpiCC
	CCOPT = -c -O3 -Wall -D$(OSTYPE) -DFORTRAN_ALLCAPS -DMPI3d -DMPI_TIMING -DGCC_2_96 $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) 

	ifeq ($(POSTPROCESS),Embedded)
		CCOPT = -c -O3 -Wall -D$(OSTYPE) -DFORTRAN_ALLCAPS -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	endif
	LD = $(MPI_HOME)/bin/mpiCC
	LDOPT  = -Wl,-static

	MISCLIB0   = $(ABSOFT)/lib/libU77.a $(ABSOFT)/lib/libf77math.a $(ABSOFT)/lib/libf90math.a $(ABSOFT)/lib/libfio.a -ldl -lpthread -lm -lg2c
	MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
	NML_OPTIONS = -O3 -Wall -ansi -pedantic

#	FORTRAN = $(MPI_HOME)/bin/mpif77
	FORTRAN = $(ABSOFT)/bin/f90

	# LAPACK-specific options
	LAPACK_NOOPT =
	ifeq ($(OS),Darwin)
		LAPACK_OPTS   = -O2 -N11
		LAPACK_LDOPTS = -O2 -N11
	else
		# Use with FORTRAN_UNDERSCORE
		#LAPACK_OPTS = -O2 -cpu:host -f -N22 -B108
		# Use with FORTRAN_ALLCAPS
		LAPACK_OPTS   = -O2 -cpu:host
		LAPACK_LDOPTS = -O2 -cpu:host
	endif

	ifeq ($(OS),Darwin)
		ARPACK_OPTS   = -O2 -N11
	else
		# Use with FORTRAN_UNDERSCORE
		#ARPACK_OPTS  = -O2 -cpu:host -f -N22 -B108
		# Use with FORTRAN_ALLCAPS
		ARPACK_OPTS   = -O2 -cpu:host
	endif

endif


# MPICH Intel Pentium4 optimized MPI Static
ifeq ($(BUILD_TYPE),intel_p4_mpi_static)
	AR = xiar
	ifeq ($(MPI_HOME),)
		#@ echo "Error - Undefined environment variable, MPI_HOME"
	endif

	CCMKD = $(MPI_HOME)/bin/mpiCC -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC = $(MPI_HOME)/bin/mpiCC
	CCOPT = -c -O3 -Kc++ -march=pentium4 -mcpu=pentium4 -axW -tpp7 -xW -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) 
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT = -c -O3 -Kc++ -march=pentium4 -mcpu=pentium4 -axW -tpp7 -xW -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	endif
	LD = $(MPI_HOME)/bin/mpiCC
	LDOPT  = -O3 -Wl,-static -Wall -march=pentium4 -mcpu=pentium4 -axW -tpp7 -xW

	MISCLIB0   = -lPEPCF90 -lIEPCF90 -lF90 -lcprts -ldl -lpthread -lm
	MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
	NML_OPTIONS = -O3 -Wall -ansi -march=pentium4 -mcpu=pentium4 -axW -tpp7 -xW

	FORTRAN = $(MPI_HOME)/bin/mpif77

	# LAPACK-specific options
	LAPACK_NOOPT  = -O3 -mp -fp_port
	LAPACK_OPTS   = -O3 -mp -fp_port -axW -tpp7 -xW
	LAPACK_LDOPTS = -O3 -mp -fp_port -axW -tpp7 -xW -Vaxlib

	# ARPACK-specific options
	ARPACK_OPTS   = -O3 -mp -axW -tpp7 -xW

endif


# MPICH Intel Pentium4 optimized MPI
ifeq ($(BUILD_TYPE),intel_p4_mpi)
	AR = xiar
	ifeq ($(MPI_HOME),)
		#@ echo "Error - Undefined environment variable, MPI_HOME"
	endif

	CCMKD = $(MPI_HOME)/bin/mpiCC -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC = $(MPI_HOME)/bin/mpiCC
	CCOPT = -c -O3 -Kc++ -march=pentium4 -mcpu=pentium4 -axW -tpp7 -xW -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) 
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT = -c -O3 -Kc++ -march=pentium4 -mcpu=pentium4 -axW -tpp7 -xW -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	endif
	LD = $(MPI_HOME)/bin/mpiCC
	LDOPT  = -O3 -i_dynamic -Wall -march=pentium4 -mcpu=pentium4 -axW -tpp7 -xW

	MISCLIB0   = -ldl -lpthread -lm -lPEPCF90 -lIEPCF90 -lF90 -lintrins
	MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
	NML_OPTIONS = -O3 -Wall -ansi -march=pentium4 -mcpu=pentium4 -axW -tpp7 -xW

	FORTRAN = $(MPI_HOME)/bin/mpif77

	# LAPACK-specific options
	LAPACK_NOOPT  = -O3 -mp -fp_port
	LAPACK_OPTS   = -O3 -mp -fp_port -axW -tpp7 -xW
	LAPACK_LDOPTS = -O3 -mp -fp_port -axW -tpp7 -xW -Vaxlib

	# ARPACK-specific options
	ARPACK_OPTS   = -O3 -mp -axW -tpp7 -xW

endif

# MPICH Intel Pentium4 optimized MPI
ifeq ($(BUILD_TYPE),intel7_p4)
	AR = xiar

	CCMKD = icpc -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY -DNoGnuStatementExpre -DNO_PARPACK
	CC = icpc
	CCOPT = -c -O3 -Kc++ -march=pentium4 -mcpu=pentium4 -axW -tpp7 -xW -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DFAKE_MPI  $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DNoGnuStatementExpre -DNO_PARPACK $(HAM_STORE_DEF) 
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT = -c -O3 -Kc++ -march=pentium4 -mcpu=pentium4 -axW -tpp7 -xW -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	endif
	LD = icpc
	LDOPT  = -O3 -i_dynamic -Wall -march=pentium4 -mcpu=pentium4 -axW -tpp7 -xW

	MISCLIB0   = -ldl -lpthread -lm -lPEPCF90 -lIEPCF90 -lF90 -lintrins
	MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
	NML_OPTIONS = -O3 -Wall -ansi -march=pentium4 -mcpu=pentium4 -axW -tpp7 -xW

	FORTRAN = ifc

	# LAPACK-specific options
	LAPACK_NOOPT  = -O3 -mp -fp_port
	LAPACK_OPTS   = -O3 -mp -fp_port -axW -tpp7 -xW
	LAPACK_LDOPTS = -O3 -mp -fp_port -axW -tpp7 -xW -Vaxlib

	# ARPACK-specific options
	ARPACK_OPTS   = -O3 -mp -axW -tpp7 -xW

endif




# MPICH Intel Pentium4 optimized MPI Static
ifeq ($(BUILD_TYPE),intel8_p4_mpi_static_debug)
	AR = xiar
	ifeq ($(MPI_HOME),)
		#@ echo "Error - Undefined environment variable, MPI_HOME"
	endif

	CCMKD = $(MPI_HOME)/bin/mpiCC -MM   -I/opt/intel_cc_80/include/c++
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC = $(MPI_HOME)/bin/mpiCC
	CCOPT = -c -g -Kc++ -march=pentium4 -mcpu=pentium4 -axW -tpp7 -xW -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF)   -I/opt/intel_cc_80/include/c++ $(HAM_STORE_DEF) 
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT = -c -g -Kc++ -march=pentium4 -mcpu=pentium4 -axW -tpp7 -xW -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	endif
	LD = $(MPI_HOME)/bin/mpiCC
	LDOPT  = -g -Wl,-static -Wall -march=pentium4 -mcpu=pentium4 -axW -tpp7 -xW

	MISCLIB0   = -lpthread -lm /opt/intel_fc_80/lib/libifcore.a /opt/intel_fc_80/lib/libifport.a -ldl  -lsvml

	MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
	NML_OPTIONS = -g -Kc++ -Wall -ansi -march=pentium4 -mcpu=pentium4 -axW -tpp7 -xW

	FORTRAN = $(MPI_HOME)/bin/mpif77

	# LAPACK-specific options
	LAPACK_NOOPT  = -g -mp -fp_port
	LAPACK_OPTS   = -g -mp -fp_port -axW -tpp7 -xW
	LAPACK_LDOPTS = -g -mp -fp_port -axW -tpp7 -xW -Vaxlib

	# ARPACK-specific options
	ARPACK_OPTS   = -g -mp -axW -tpp7 -xW

        ##############################################################
        # TRACEMIN-specific options (these are set here only for choice of BUILD_TYPE equal to intel9_64_mpi_intelfast (lear.rcac.purdue.edu))
        SCTRACEMIN_OPTS = -O3
        PCTRACEMIN_OPTS = -O3
        SQTRACEMIN_OPTS = -O3
        PQTRACEMIN_OPTS = -O3
        FORTRAN90 = $(MPI_HOME)/bin/mpif90  #tracemin fortran compiler
        ##############################################################
endif



# MPICH Intel Pentium4 optimized MPI Static
ifeq ($(BUILD_TYPE),intel8_p4_mpi_static)
	AR = xiar
	ifeq ($(MPI_HOME),)
		#@ echo "Error - Undefined environment variable, MPI_HOME"
	endif

	CCMKD = $(MPI_HOME)/bin/mpiCC -MM   -I/opt/intel_cc_80/include/c++  -I/opt/intel_cc_80/include/c++/backward
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC = $(MPI_HOME)/bin/mpiCC
	CCOPT = -c -O3 -Kc++ -march=pentium4 -mcpu=pentium4 -axW -tpp7 -xW -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF)   -I/opt/intel_cc_80/include/c++ $(HAM_STORE_DEF) 
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT = -c -O3 -Kc++ -march=pentium4 -mcpu=pentium4 -axW -tpp7 -xW -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	endif
	LD = $(MPI_HOME)/bin/mpiCC
	LDOPT  = -O3 -Wl,-static -Wall -march=pentium4 -mcpu=pentium4 -axW -tpp7 -xW

	MISCLIB0   = -lpthread -lm /opt/intel_fc_80/lib/libifcore.a /opt/intel_fc_80/lib/libifport.a -ldl -lsvml
	MISCINCL0  =

	#XERCES specific options
	XERCES_CC =icc
	XERCES_CCMKD =icc -MM  -I/opt/intel_cc_80/include/c++  -I/opt/intel_cc_80/include/c++/backward
	XERCES_CCOPT=-g -Kc++ -c
	XERCES_LD=icc
	XERCES_LDOPT=-g -Kc++

        #XMLtoC specific options
	XMLtoC_CC=icc 
	XMLtoC_CCOPT=-g -Kc++ -c
	XMLtoC_CCMKD=icc -MM  -I/opt/intel_cc_80/include/c++  -I/opt/intel_cc_80/include/c++/backward
	XMLtoC_LD=icc
	XMLtoC_LDOPT=-g -Kc++

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
	NML_OPTIONS = -O3 -Kc++ -Wall -ansi -march=pentium4 -mcpu=pentium4 -axW -tpp7 -xW

	FORTRAN = $(MPI_HOME)/bin/mpif77

	# LAPACK-specific options
	LAPACK_NOOPT  = -g -O0 -mp 
	LAPACK_OPTS   = -g -O0 -mp 
	LAPACK_LDOPTS = -g -O0 -mp 
	LAPACK_NOOPT  =  -mp -fp_port
	LAPACK_OPTS   =  -mp -fp_port -axW -tpp7 
	LAPACK_LDOPTS =  -mp -fp_port -axW -tpp7 -Vaxlib

	# ARPACK-specific options
	ARPACK_OPTS   = -O3 -mp -axW -tpp7 -xW

        ##############################################################
        # TRACEMIN-specific options (these are set here only for choice of BUILD_TYPE equal to intel9_64_mpi_intelfast (lear.rcac.purdue.edu))
        SCTRACEMIN_OPTS = -O3
        PCTRACEMIN_OPTS = -O3
        SQTRACEMIN_OPTS = -O3
        PQTRACEMIN_OPTS = -O3
        FORTRAN90 = $(MPI_HOME)/bin/mpif90  #tracemin fortran compiler
        ##############################################################

endif

# MPICH Intel Pentium4 optimized MPI Static
ifeq ($(BUILD_TYPE),intel8_p4_mpi_static_infband)
	AR = xiar
	ifeq ($(MPI_HOME),)
		#@ echo "Error - Undefined environment variable, MPI_HOME"
	endif

	CCMKD = $(MPI_HOME)/bin/mpiCC -MM   -I/opt/intel_cc_80/include/c++
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC = $(MPI_HOME)/bin/mpiCC
	CCOPT = -c -O3 -Kc++ -march=pentium4 -mcpu=pentium4 -axW -tpp7 -xW -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF)   -I/opt/intel_cc_80/include/c++ $(HAM_STORE_DEF) 
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT = -c -O3 -Kc++ -march=pentium4 -mcpu=pentium4 -axW -tpp7 -xW -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	endif
	LD = $(MPI_HOME)/bin/mpiCC
	LDOPT  = -O3 -Wl,-static -Wall -march=pentium4 -mcpu=pentium4 -axW -tpp7 -xW

	MISCLIB0   = -lpthread -lm /opt/intel_fc_80/lib/libifcore.a /opt/intel_fc_80/lib/libifport.a -ldl -lsvml
	MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
	NML_OPTIONS = -O3 -Kc++ -Wall -ansi -march=pentium4 -mcpu=pentium4 -axW -tpp7 -xW

	FORTRAN = $(MPI_HOME)/bin/mpif77.i

	# LAPACK-specific options
	LAPACK_NOOPT  = -O3 -mp -fp_port
	LAPACK_OPTS   = -O3 -mp -fp_port -axW -tpp7 -xW
	LAPACK_LDOPTS = -O3 -mp -fp_port -axW -tpp7 -xW -Vaxlib

	# ARPACK-specific options
	ARPACK_OPTS   = -O3 -mp -axW -tpp7 -xW
        ##############################################################
        # TRACEMIN-specific options (these are set here only for choice of BUILD_TYPE equal to intel9_64_mpi_intelfast (lear.rcac.purdue.edu))
        SCTRACEMIN_OPTS = -O3
        PCTRACEMIN_OPTS = -O3
        SQTRACEMIN_OPTS = -O3
        PQTRACEMIN_OPTS = -O3
        FORTRAN90 = $(MPI_HOME)/bin/mpif90  #tracemin fortran compiler
        ##############################################################
endif

# MPICH Intel optimized with SSE3
ifeq ($(BUILD_TYPE),intel8_p4_mpi_intelfast)
	AR = xiar
	ifeq ($(MPI_HOME),)
		#@ echo "Error - Undefined environment variable, MPI_HOME"
	endif

	CCMKD = $(MPI_HOME)/bin/mpiCC -MM   -I/opt/intel_cc_80/include/c++
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC = $(MPI_HOME)/bin/mpiCC

	CFLAGS = -O2 -Kc++ -xN

	    CCOPT = -c $(CFLAGS) -UINTEL_MODIFIED -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF)   -I/opt/intel_cc_80/include/c++ $(HAM_STORE_DEF) 
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT = -c -O3 -Kc++ -march=pentium4 -mcpu=pentium4 -axW -tpp7 -xW -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	endif
	LD = $(MPI_HOME)/bin/mpiCC
	    LDOPT  = $(CFLAGS) -Wl,-static -Wall


	    MISCLIB0   = -lpthread -lm /opt/intel_fc_80/lib/libifcore.a /opt/intel_fc_80/lib/libifport.a -ldl -lsvml 
	MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
	    NML_OPTIONS = $(CCOPT)

	FORTRAN = $(MPI_HOME)/bin/mpif77

	# LAPACK-specific options
	LAPACK_NOOPT  = -O3 -mp -fp_port
	LAPACK_OPTS   = -O3 -mp -fp_port -axW -tpp7 -xW
	LAPACK_LDOPTS = -O3 -mp -fp_port -axW -tpp7 -xW -Vaxlib

	# ARPACK-specific options
	ARPACK_OPTS   = -O3 -mp -axW -tpp7 -xW

        ##############################################################
        # TRACEMIN-specific options (these are set here only for choice of BUILD_TYPE equal to intel9_64_mpi_intelfast (lear.rcac.purdue.edu))
        SCTRACEMIN_OPTS = -O3
        PCTRACEMIN_OPTS = -O3
        SQTRACEMIN_OPTS = -O3
        PQTRACEMIN_OPTS = -O3
        FORTRAN90 = $(MPI_HOME)/bin/mpif90  #tracemin fortran compiler
        ##############################################################
endif

# MPICH Intel 9
ifeq ($(BUILD_TYPE),intel9_64_inteldebug)
	AR = xiar

	CCMKD = icc -MM   -I/opt/intel/cc/9.0/include/c++ 
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC = icc 

	CFLAGS = -g -Kc++ 
	CFLAGS = -O3 -xP -Kc++ 
	CFLAGS = -xP -Kc++ 

	    CCOPT = -c $(CFLAGS) -UINTEL_MODIFIED -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DFAKE_MPI  -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF)   -I/opt/intel/cc/9.0/include/c++ $(HAM_STORE_DEF) -Dmem_64BIT
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT = -c -O3 -xP -Kc++ -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DFAKE_MPI  -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) -Dmem_64BIT
	endif
	LIB_PARPACK= 

	LD = icc
	    LDOPT  = $(CFLAGS) -Wl,-static -Wall
	    LDOPT  = $(CFLAGS)  -Wall


	    MISCLIB0   = -lpthread -lm /opt/intel/fc/9.0/lib/libifcore.a /opt/intel/fc/9.0/lib/libifport.a -ldl -lsvml 
	    MISCLIB0   = -lpthread -lm /opt/intel/cc/9.0/lib/libcprts.a /opt/intel/cc/9.0/lib/libcxa.a /opt/intel/cc/9.0/lib/libunwind.a /opt/intel/fc/9.0/lib/libifcore.a /opt/intel/fc/9.0/lib/libifport.a -ldl  /opt/intel/cc/9.0/lib/libsvml.a
	MISCINCL0  =

	#XERCES specific options
	XERCES_CCOPT= -c -g -Kc++ -Dmem_64BIT

        #XMLtoC specific options
        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT  -Dmem_64BIT
	    NML_OPTIONS = $(CCOPT)

	FORTRAN =ifort 

	# LAPACK-specific options
	LAPACK_NOOPT  = -g -mp -fp_port
	LAPACK_OPTS   = -g -mp -fp_port  -tpp7 
	LAPACK_LDOPTS = -g -mp -fp_port  -tpp7 -Vaxlib

	# ARPACK-specific options
	ARPACK_OPTS   = -g -mp -tpp7 

endif




# MPICH Intel optimized with SSE3
ifeq ($(BUILD_TYPE),intel9_64_intelfast)
	AR = xiar

	CCMKD = icc -MM   -I/opt/intel/cc/9.0/include/c++ 
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC = icc 

	CFLAGS = -O2 -Kc++ -axP

	    CCOPT = -c $(CFLAGS) -DINTEL_MODIFIED -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DFAKE_MPI  -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF)   -I/opt/intel/cc/9.0/include/c++ $(HAM_STORE_DEF) 
#	    CCOPT = -c $(CFLAGS) -UINTEL_MODIFIED -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DFAKE_MPI  -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF)   -I/opt/intel/cc/9.0/include/c++ $(HAM_STORE_DEF)  -Dmem_64BIT
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT = -c -O2 -Kc++ -axP -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DFAKE_MPI  -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF)  -Dmem_64BIT
	endif
	LIB_PARPACK= 

	LD = icc
	    LDOPT  = $(CFLAGS) -Wl,-static -Wall
	    LDOPT  = $(CFLAGS)  -Wall


	    MISCLIB0   = -lpthread -lm /opt/intel/fc/9.0/lib/libifcore.a /opt/intel/fc/9.0/lib/libifport.a -ldl -lsvml 
	    MISCLIB0   = -lpthread -lm /opt/intel/cc/9.0/lib/libcprts.a /opt/intel/cc/9.0/lib/libcxa.a /opt/intel/cc/9.0/lib/libunwind.a /opt/intel/fc/9.0/lib/libifcore.a /opt/intel/fc/9.0/lib/libifport.a -ldl  /opt/intel/cc/9.0/lib/libsvml.a
	MISCINCL0  =

	#XERCES specific options
	XERCES_CCOPT= -c -O1 -Kc++

        #XMLtoC specific options
        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT  -Dmem_64BIT
	    NML_OPTIONS = $(CCOPT)

	FORTRAN =ifort 

	# LAPACK-specific options
	LAPACK_NOOPT  = -O3 -mp -fp_port
	LAPACK_OPTS   = -O3 -mp -fp_port -axP -tpp7 -xW
	LAPACK_LDOPTS = -O3 -mp -fp_port -axP -tpp7 -xW -Vaxlib

	# ARPACK-specific options
	ARPACK_OPTS   = -O3 -mp -axP -tpp7 -xW

endif




# MPICH Intel optimized with SSE3
ifeq ($(BUILD_TYPE),intel9_64_intelfast_91045)
	CCDIR=/apps/rhel5/intel/Compiler/11.1/038
	FCDIR=/apps/rhel5/intel/Compiler/11.1/038
	AR = xiar
	ifeq ($(MPI_HOME),)
		#@ echo "Error - Undefined environment variable, MPI_HOME"
	endif

	CCMKD = icc -MM    -I${CCDIR}/include/c++ 
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC = icc -I${CCDIR}/include/c++

	CFLAGS = -O2 -Kc++ -axP
        ####################################################################
        ##NOTE: -DINTEL_MODIFIED sets the INTEL optimizations ON, 
		##use -UINTEL_MODIFIED in CCOPT (next line) to turn them back OFF!!
        ####################################################################
	CCOPT = -c $(CFLAGS) -DINTEL_MODIFIED -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DFAKE_MPI -UMPI3d -UMPI_TIMING $(MKL_MACRO_DEF)  $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -I${CCDIR}/include/c++ $(HAM_STORE_DEF)  -Dmem_64BIT $(MKL_INCLUDE)

	ifeq ($(POSTPROCESS), Embedded)
		CCOPT = -c -O2 -Kc++ -axP -DINTEL_MODIFIED -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DFAKE_MPI -UMPI3d -UMPI_TIMING $(MKL_MACRO_DEF)  $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess -I${CCDIR}/include/c++ $(HAM_STORE_DEF)  -Dmem_64BIT $(MKL_INCLUDE)
	endif

	LD = $(CC)
	#LDOPT  = $(CFLAGS) -Wl,-static -Wall
	LDOPT  = $(CFLAGS)  -Wall


	MISCLIB0   = -lpthread -lm ${FCDIR}/lib/libifcore.a ${FCDIR}/lib/libifport.a -ldl -lsvml $(MKL_LINK)
	MISCLIB0   = -lpthread -lm ${CCDIR}/lib/libcprts.a ${CCDIR}/lib/libcxa.a ${CCDIR}/lib/libunwind.a ${FCDIR}/lib/libifcore.a ${FCDIR}/lib/libifport.a -ldl  ${CCDIR}/lib/libsvml.a $(MPI_HOME)/lib/libmpichf90.a $(MKL_LINK)
	MISCINCL0  =

	#XERCES specific options
	XERCES_CCOPT= $(CCOPT) -c -O1 -Kc++ -axP
	XERCES_CC= icc

  	#XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)
	XMLtoC_CC = icc

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT  -Dmem_64BIT
	NML_OPTIONS = $(CCOPT)

	FORTRAN   = $(FC)

	# LAPACK-specific options
	LAPACK_NOOPT  = -O3 -mp -fp_port
	LAPACK_OPTS   = -O3 -mp -fp_port -axP -tpp7 
	LAPACK_LDOPTS = -O3 -mp -fp_port -axP -tpp7 -Vaxlib

	# ARPACK-specific options
	ARPACK_OPTS   = -O3 -mp -axP -tpp7 -xW

	##############################################################
	# TRACEMIN-specific options (these are set here only for choice of BUILD_TYPE: intel8,  intel9_64_mpi_intelfast (lear.rcac.purdue.edu))  
	SCTRACEMIN_OPTS = -O3
 	PCTRACEMIN_OPTS = -O3
	SQTRACEMIN_OPTS = -O3
	PQTRACEMIN_OPTS = -O3
	FORTRAN90 = $(F90)
	##############################################################
endif





# MPICH Intel optimized with SSE3
ifeq ($(BUILD_TYPE),intel9_64_mpi_intelfast)
	CCDIR=/opt/intel/cce/9.1.045
	FCDIR=/opt/intel/fce/9.1.045
	AR = xiar
	ifeq ($(MPI_HOME),)
		#@ echo "Error - Undefined environment variable, MPI_HOME"
	endif

	CCMKD = $(MPI_HOME)/bin/mpiCC -MM    -I${CCDIR}/include/c++ 
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC = $(MPI_HOME)/bin/mpiCC -I${CCDIR}/include/c++

	CFLAGS = -O2 -Kc++ -axP
        ####################################################################
        ##NOTE: -DINTEL_MODIFIED sets the INTEL optimizations ON, 
		##use -UINTEL_MODIFIED in CCOPT (next line) to turn them back OFF!!
        ####################################################################
	CCOPT = -c $(CFLAGS) -DINTEL_MODIFIED -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(MKL_MACRO_DEF)  $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -I${CCDIR}/include/c++ $(HAM_STORE_DEF)  -Dmem_64BIT $(MKL_INCLUDE)

	ifeq ($(POSTPROCESS), Embedded)
		CCOPT = -c -O2 -Kc++ -axP -DINTEL_MODIFIED -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(MKL_MACRO_DEF)  $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess -I${CCDIR}/include/c++ $(HAM_STORE_DEF)  -Dmem_64BIT $(MKL_INCLUDE)
	endif

	LD = $(CC)
	LDOPT  = $(CFLAGS) -Wl,-static -Wall
	LDOPT  = $(CFLAGS)  -Wall


	MISCLIB0   = -lpthread -lm ${FCDIR}/lib/libifcore.a ${FCDIR}/lib/libifport.a -ldl -lsvml $(MKL_LINK)
	MISCLIB0   = -lpthread -lm ${CCDIR}/lib/libcprts.a ${CCDIR}/lib/libcxa.a ${CCDIR}/lib/libunwind.a ${FCDIR}/lib/libifcore.a ${FCDIR}/lib/libifport.a -ldl  ${CCDIR}/lib/libsvml.a $(MPI_HOME)/lib/libmpichf90.a $(MKL_LINK)
	MISCINCL0  =

	#XERCES specific options
	XERCES_CCOPT= $(CCOPT) -c -O1 -Kc++ -axP
	XERCES_CC= icc

  	#XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)
	XMLtoC_CC = icc

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT  -Dmem_64BIT
	NML_OPTIONS = $(CCOPT)

	FORTRAN   = $(MPI_HOME)/bin/mpif77

	# LAPACK-specific options
	LAPACK_NOOPT  = -O3 -mp -fp_port
	LAPACK_OPTS   = -O3 -mp -fp_port -axP -tpp7 
	LAPACK_LDOPTS = -O3 -mp -fp_port -axP -tpp7 -Vaxlib

	# ARPACK-specific options
	ARPACK_OPTS   = -O3 -mp -axP -tpp7 -xW

	##############################################################
	# TRACEMIN-specific options (these are set here only for choice of BUILD_TYPE: intel8,  intel9_64_mpi_intelfast (lear.rcac.purdue.edu))  
	SCTRACEMIN_OPTS = -O3
 	PCTRACEMIN_OPTS = -O3
	SQTRACEMIN_OPTS = -O3
	PQTRACEMIN_OPTS = -O3
	FORTRAN90 = $(MPI_HOME)/bin/mpif90  #tracemin fortran compiler
	##############################################################
endif

ifeq ($(BUILD_TYPE),ranger_intel10_mvapich2)

	MPI_HOME=$(MPICH_HOME)
	CCDIR=/opt/apps/intel/10.1/cc
	FCDIR=/opt/apps/intel/10.1/fc
        
	AR = xiar
        
	ifeq ($(MPI_HOME),)
	        #@ echo "Error - Undefined environment variable, MPI_HOME"
	endif
        
	CCMKD = $(MPI_HOME)/bin/mpicxx -MM    -I${CCDIR}/include/c++
        
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
        
	CC = $(MPI_HOME)/bin/mpicxx -I${CCDIR}/include/c++
	CFLAGS = -O2 -Kc++ -axP
        
	####################################################################
	##NOTE: -DINTEL_MODIFIED sets the INTEL optimizations ON,
	        ##use -UINTEL_MODIFIED in CCOPT (next line) to turn them back OFF!!
	####################################################################
        
	CCOPT = -c $(CFLAGS) -DINTEL_MODIFIED -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(MKL_MACRO_DEF)  $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -I${CCDIR}/include/c++ $(HAM_STORE_DEF)  -Dmem_64BIT $(MKL_INCLUDE) -DMPICH_IGNORE_CXX_SEEK
        
	ifeq ($(POSTPROCESS), Embedded)
	        CCOPT = -c -O2 -Kc++ -axP -DINTEL_MODIFIED -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(MKL_MACRO_DEF)  $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess -I${CCDIR}/include/c++ $(HAM_STORE_DEF)  -Dmem_64BIT $(MKL_INCLUDE)
	endif
        
	LD = $(CC)
        
	#LDOPT  = $(CFLAGS) -Wl,-static -Wall
	LDOPT  = $(CFLAGS)  -Wall
        
	#MISCLIB0   = -lpthread -lm ${FCDIR}/lib/libifcore.a ${FCDIR}/lib/libifport.a -ldl -lsvml $(MKL_LINK)
	MISCLIB0   = -lpthread -lm ${CCDIR}/lib/libcprts.a ${CCDIR}/lib/libcxa.a ${CCDIR}/lib/libunwind.a ${FCDIR}/lib/libifcore.a ${FCDIR}/lib/libifport.a -ldl  ${CCDIR}/lib/libsvml.a $(MPI_HOME)/lib/libtvmpich.so $(MKL_LINK)
	MISCINCL0  =

        #XERCES specific options
        
	XERCES_CCOPT= $(CCOPT) -c -O1 -Kc++ -axP
	XERCES_CC= icc

        #XMLtoC specific options
        
	XMLtoC_CCOPT=$(CCOPT)
	XMLtoC_CC = icc

        # NML-specific options
        
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT  -Dmem_64BIT
	NML_OPTIONS = $(CCOPT)
	FORTRAN   = $(MPI_HOME)/bin/mpif77

        # LAPACK-specific options
        
	LAPACK_NOOPT  = -O3 -mp -fp_port
	LAPACK_OPTS   = -O3 -mp -fp_port -axP 
	LAPACK_LDOPTS = -O3 -mp -fp_port -axP -Vaxlib

        # ARPACK-specific options
        
	ARPACK_OPTS   = -O3 -mp -axP -xW
	
        ##############################################################
        # TRACEMIN-specific options (these are set here only for choice of BUILD_TYPE: intel8,  intel9_64_mpi_intelfast (lear.rcac.purdue.edu))
	SCTRACEMIN_OPTS = -O3
	PCTRACEMIN_OPTS = -O3
	SQTRACEMIN_OPTS = -O3
	PQTRACEMIN_OPTS = -O3
	FORTRAN90 = $(MPI_HOME)/bin/mpif90  #tracemin fortran compiler
        ##############################################################

endif

ifeq ($(BUILD_TYPE),ranger_intel10_openmpi)

	MPI_HOME=/opt/apps/intel10_1/openmpi/1.2.4
	CCDIR=/opt/apps/intel/10.1/cc
	FCDIR=/opt/apps/intel/10.1/fc
        
	AR = xiar
        
	ifeq ($(MPI_HOME),)
	        #@ echo "Error - Undefined environment variable, MPI_HOME"
	endif
        
	CCMKD = $(MPI_HOME)/bin/mpiCC -MM    -I${CCDIR}/include/c++
        
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
        
	CC = $(MPI_HOME)/bin/mpiCC -I${CCDIR}/include/c++
	CFLAGS = -O2 -Kc++ -axP
        
	####################################################################
	##NOTE: -DINTEL_MODIFIED sets the INTEL optimizations ON,
	        ##use -UINTEL_MODIFIED in CCOPT (next line) to turn them back OFF!!
	####################################################################
        
	CCOPT = -c $(CFLAGS) -DINTEL_MODIFIED -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(MKL_MACRO_DEF)  $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -I${CCDIR}/include/c++ $(HAM_STORE_DEF)  -Dmem_64BIT $(MKL_INCLUDE)
        
	ifeq ($(POSTPROCESS), Embedded)
	        CCOPT = -c -O2 -Kc++ -axP -DINTEL_MODIFIED -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(MKL_MACRO_DEF)  $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess -I${CCDIR}/include/c++ $(HAM_STORE_DEF)  -Dmem_64BIT $(MKL_INCLUDE)
	endif
        
	LD = $(CC)
        
	#LDOPT  = $(CFLAGS) -Wl,-static -Wall
	LDOPT  = $(CFLAGS)  -Wall
        
	#MISCLIB0   = -lpthread -lm ${FCDIR}/lib/libifcore.a ${FCDIR}/lib/libifport.a -ldl -lsvml $(MKL_LINK)
	MISCLIB0   = -lpthread -lm ${CCDIR}/lib/libcprts.a ${CCDIR}/lib/libcxa.a ${CCDIR}/lib/libunwind.a ${FCDIR}/lib/libifcore.a ${FCDIR}/lib/libifport.a -ldl  ${CCDIR}/lib/libsvml.a ${MPI_HOME}/lib/libmpi_f77.so $(MKL_LINK)
	MISCINCL0  =

        #XERCES specific options
        
	XERCES_CCOPT= $(CCOPT) -c -O1 -Kc++ -axP
	XERCES_CC= icc

        #XMLtoC specific options
        
	XMLtoC_CCOPT=$(CCOPT)
	XMLtoC_CC = icc

        # NML-specific options
        
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT  -Dmem_64BIT
	NML_OPTIONS = $(CCOPT)
	FORTRAN   = $(MPI_HOME)/bin/mpif77

        # LAPACK-specific options
        
	LAPACK_NOOPT  = -O3 -mp -fp_port
	LAPACK_OPTS   = -O3 -mp -fp_port -axP 
	LAPACK_LDOPTS = -O3 -mp -fp_port -axP -Vaxlib

        # ARPACK-specific options
        
	ARPACK_OPTS   = -O3 -mp -axP -xW
	
        ##############################################################
        # TRACEMIN-specific options (these are set here only for choice of BUILD_TYPE: intel8,  intel9_64_mpi_intelfast (lear.rcac.purdue.edu))
	SCTRACEMIN_OPTS = -O3
	PCTRACEMIN_OPTS = -O3
	SQTRACEMIN_OPTS = -O3
	PQTRACEMIN_OPTS = -O3
	FORTRAN90 = $(MPI_HOME)/bin/mpif90  #tracemin fortran compiler
        ##############################################################

endif


# MPICH Intel optimized with SSE3
ifeq ($(BUILD_TYPE),intel10_64_mpi_intelfast)

	CCDIR=/opt/intel/cce/10.0.023
	FCDIR=/opt/intel/fce/10.0.023
	AR = xiar
	ifeq ($(MPI_HOME),)
		#@ echo "Error - Undefined environment variable, MPI_HOME"
	endif

	CCMKD = $(MPI_HOME)/bin/mpiCC -MM    -I${CCDIR}/include/c++ 
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC = $(MPI_HOME)/bin/mpiCC -I${CCDIR}/include/c++

	CFLAGS = -O1 -Kc++ -axP
        ####################################################################
        ##NOTE: -DINTEL_MODIFIED sets the INTEL optimizations ON, 
		##use -UINTEL_MODIFIED in CCOPT (next line) to turn them back OFF!!
        ####################################################################
	CCOPT = -c $(CFLAGS) -DINTEL_MODIFIED -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(MKL_MACRO_DEF)  $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -I${CCDIR}/include/c++ $(HAM_STORE_DEF)  -Dmem_64BIT $(MKL_INCLUDE)

	ifeq ($(POSTPROCESS), Embedded)
		CCOPT = -c -O1 -Kc++ -axP -DINTEL_MODIFIED -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(MKL_MACRO_DEF)  $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess -I${CCDIR}/include/c++ $(HAM_STORE_DEF)  -Dmem_64BIT $(MKL_INCLUDE)
	endif

	LD = $(CC)
	LDOPT  = $(CFLAGS) -Wl,-static -Wall
	LDOPT  = $(CFLAGS)  -Wall


	#MISCLIB0   = -lpthread -lm ${FCDIR}/lib/libifcore.a ${FCDIR}/lib/libifport.a -ldl -lsvml $(MKL_LINK)
	#MISCLIB0   = -lpthread -lm ${CCDIR}/lib/libcprts.a ${CCDIR}/lib/libcxa.a ${CCDIR}/lib/libunwind.a ${FCDIR}/lib/libifcore.a ${FCDIR}/lib/libifport.a -ldl  ${CCDIR}/lib/libsvml.a $(MPI_HOME)/lib/libmpichf90.a $(MKL_LINK)
        MISCLIB0   = -lpthread -lm ${CCDIR}/lib/libcxaguard.a ${FCDIR}/lib/libifcore.a ${FCDIR}/lib/libifport.a -ldl ${CCDIR}/lib/libsvml.a $(MPI_HOME)/lib/libmpichf90.a $(MKL_LINK)
	MISCINCL0  =

	#XERCES specific options
	XERCES_CCOPT= $(CCOPT)
	XERCES_CC= icc

  	#XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)
	XMLtoC_CC = icc

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT  -Dmem_64BIT
	NML_OPTIONS = $(CCOPT)

	FORTRAN   = $(MPI_HOME)/bin/mpif77

	# LAPACK-specific options
	LAPACK_NOOPT  = -O1 -mp -fp_port -mtune=pentium4 -xW -axP
	LAPACK_OPTS   = -O1 -mp -fp_port -mtune=pentium4 -xW -axP
	LAPACK_LDOPTS = -O1 -mp -fp_port -mtune=pentium4 -xW -axP -Vaxlib

	# ARPACK-specific options
	ARPACK_OPTS   = -O1 -mp -mtune=pentium4 -axP -xW

	##############################################################
	# TRACEMIN-specific options (these are set here only for choice of BUILD_TYPE: intel8,  intel9_64_mpi_intelfast (lear.rcac.purdue.edu))  
	SCTRACEMIN_OPTS = -O1
 	PCTRACEMIN_OPTS = -O1
	SQTRACEMIN_OPTS = -O1
	PQTRACEMIN_OPTS = -O1
	FORTRAN90 = $(MPI_HOME)/bin/mpif90  #tracemin fortran compiler
	##############################################################

endif

# MPICH Intel optimized with SSE3
ifeq ($(BUILD_TYPE),intel20_64_openmpi_gadi)
        ICC_HOME=$(INTEL_COMPILER_BASE)
        IFORT_HOME=$(INTEL_COMPILER_BASE)
        MPI_HOME=$(OPENMPI_BASE)
        CCDIR=$(INTEL_COMPILER_BASE)
        FCDIR=$(INTEL_COMPILER_BASE)
        AR = xiar
        ifeq ($(MPI_HOME),)
                #@ echo "Error - Undefined environment variable, MPI_HOME"
        endif

        CCMKD = $(MPI_HOME)/bin/mpicxx -MM    -I${CCDIR}/include/c++
        # Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
        CC = $(MPI_HOME)/bin/mpicxx

        CFLAGS = -O0 -g
        ####################################################################
        ##NOTE: -DINTEL_MODIFIED sets the INTEL optimizations ON,
                ##use -UINTEL_MODIFIED in CCOPT (next line) to turn them back OFF!!
        ####################################################################
        CCOPT = -c $(CFLAGS) -DINTEL_MODIFIED -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(MKL_MACRO_DEF)  $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF)  -Dmem_64BIT $(MKL_INCLUDE)

        ifeq ($(POSTPROCESS), Embedded)
                CCOPT = -c -O0 -g -DINTEL_MODIFIED -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(MKL_MACRO_DEF)  $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF)  -Dmem_64BIT $(MKL_INCLUDE)
        endif

        LD = $(CC)
        LDOPT  = $(CFLAGS) -Wl,-static -Wall
        LDOPT  = $(CFLAGS)  -Wall


        #MISCLIB0   = -lpthread -lm ${FCDIR}/lib/libifcore.a ${FCDIR}/lib/libifport.a -ldl -lsvml $(MKL_LINK)
        #MISCLIB0   = -lpthread -lm ${CCDIR}/lib/libcprts.a ${CCDIR}/lib/libcxa.a ${CCDIR}/lib/libunwind.a ${FCDIR}/lib/libifcore.a ${FCDIR}/lib/libifport.a -ldl  ${CCDIR}/lib/libsvml.a $(MPI_HOME)/lib/libmpichf90.a $(MKL_LINK)
        MISCLIB0   = -lpthread -mkl ${CCDIR}/lib/libifcore.a ${CCDIR}/lib/libifport.a -ldl ${CCDIR}/lib/libsvml.a -L$(MPI_HOME)/lib/ -lmpi -lmpi_mpifh_Intel
        MISCINCL0  =

        #XERCES specific options
        XERCES_CCOPT= $(CCOPT)
        XERCES_CC= icc -Kc++

        #XMLtoC specific options
        XMLtoC_CCOPT=$(CCOPT)
        XMLtoC_CC = icpc

        # NML-specific options
        NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT  -Dmem_64BIT
        NML_OPTIONS = $(CCOPT)

        FORTRAN   = $(MPI_HOME)/bin/mpifort

        # LAPACK-specific options
        LAPACK_NOOPT  = -O2
        LAPACK_OPTS   = -O2
        LAPACK_LDOPTS = -O2

        # ARPACK-specific options
        ARPACK_OPTS   = -O2

        ##############################################################
        # TRACEMIN-specific options (these are set here only for choice of BUILD_TYPE: intel8,  intel9_64_mpi_intelfast (lear.rcac.purdue.edu))
        SCTRACEMIN_OPTS = -O2
        PCTRACEMIN_OPTS = -O2
        SQTRACEMIN_OPTS = -O2
        PQTRACEMIN_OPTS = -O2
        FORTRAN90 = $(MPI_HOME)/bin/mpif90  #tracemin fortran compiler
        ##############################################################

endif

# MPICH Intel optimized with SSE3
ifeq ($(BUILD_TYPE),intel19_64_mpi_intelfast_gadi)
        ICC_HOME=$(INTEL_COMPILER_BASE)
        IFORT_HOME=$(INTEL_COMPILER_BASE)
        MPI_HOME=$(INTEL_MPI_BASE)/intel64
        CCDIR=$(INTEL_COMPILER_BASE)
        FCDIR=$(INTEL_COMPILER_BASE)
        AR = xiar
        ifeq ($(MPI_HOME),)
                #@ echo "Error - Undefined environment variable, MPI_HOME"
        endif

        CCMKD = $(MPI_HOME)/bin/mpiicpc -MM
        # Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
        CC = $(MPI_HOME)/bin/mpiicpc

        CFLAGS = -O0 -g
        ####################################################################
        ##NOTE: -DINTEL_MODIFIED sets the INTEL optimizations ON,
                ##use -UINTEL_MODIFIED in CCOPT (next line) to turn them back OFF!!
        ####################################################################
        CCOPT = -c $(CFLAGS) -DINTEL_MODIFIED -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(MKL_MACRO_DEF)  $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF)  -Dmem_64BIT $(MKL_INCLUDE)

        ifeq ($(POSTPROCESS), Embedded)
                CCOPT = -c -O0 -g -DINTEL_MODIFIED -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(MKL_MACRO_DEF)  $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF)  -Dmem_64BIT $(MKL_INCLUDE)
        endif

        LD = $(CC)
        LDOPT  = $(CFLAGS) -Wl,-static -Wall
        LDOPT  = $(CFLAGS)  -Wall


        #MISCLIB0   = -lpthread -lm ${FCDIR}/lib/libifcore.a ${FCDIR}/lib/libifport.a -ldl -lsvml $(MKL_LINK)
        #MISCLIB0   = -lpthread -lm ${CCDIR}/lib/libcprts.a ${CCDIR}/lib/libcxa.a ${CCDIR}/lib/libunwind.a ${FCDIR}/lib/libifcore.a ${FCDIR}/lib/libifport.a -ldl  ${CCDIR}/lib/libsvml.a $(MPI_HOME)/lib/libmpichf90.a $(MKL_LINK)
        MISCLIB0   = -lpthread -mkl ${CCDIR}/lib/libifcore.a ${CCDIR}/lib/libifport.a -ldl ${CCDIR}/lib/libsvml.a #$(MPI_HOME)/lib/libmpichf90.a
        MISCINCL0  =

        #XERCES specific options
        XERCES_CCOPT= $(CCOPT)
        XERCES_CC= icc -Kc++

        #XMLtoC specific options
        XMLtoC_CCOPT=$(CCOPT)
        XMLtoC_CC = icpc

        # NML-specific options
        NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT  -Dmem_64BIT
        NML_OPTIONS = $(CCOPT)

        FORTRAN   = $(MPI_HOME)/bin/mpiifort

        # LAPACK-specific options
        LAPACK_NOOPT  = -O2
        LAPACK_OPTS   = -O2
        LAPACK_LDOPTS = -O2

        # ARPACK-specific options
        ARPACK_OPTS   = -O2

        ##############################################################
        # TRACEMIN-specific options (these are set here only for choice of BUILD_TYPE: intel8,  intel9_64_mpi_intelfast (lear.rcac.purdue.edu))
        SCTRACEMIN_OPTS = -O2
        PCTRACEMIN_OPTS = -O2
        SQTRACEMIN_OPTS = -O2
        PQTRACEMIN_OPTS = -O2
        FORTRAN90 = $(MPI_HOME)/bin/mpif90  #tracemin fortran compiler
        ##############################################################

endif


# MPICH Intel optimized with SSE3
ifeq ($(BUILD_TYPE),intel15_64_mpi_intelfast)
        ICC_HOME=$(INTEL_COMPILER_BASE)
        IFORT_HOME=$(INTEL_COMPILER_BASE)
        MPI_HOME=$(INTEL_MPI_BASE)/intel64
        CCDIR=$(INTEL_COMPILER_BASE)
        FCDIR=$(INTEL_COMPILER_BASE)
        AR = xiar
        ifeq ($(MPI_HOME),)
                #@ echo "Error - Undefined environment variable, MPI_HOME"
        endif

        CCMKD = $(MPI_HOME)/bin/mpiicpc -MM
        # Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
        CC = $(MPI_HOME)/bin/mpiicpc

        CFLAGS = -O2
        ####################################################################
        ##NOTE: -DINTEL_MODIFIED sets the INTEL optimizations ON,
                ##use -UINTEL_MODIFIED in CCOPT (next line) to turn them back OFF!!
        ####################################################################
        CCOPT = -c $(CFLAGS) -DINTEL_MODIFIED -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(MKL_MACRO_DEF)  $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF)  -Dmem_64BIT $(MKL_INCLUDE)

        ifeq ($(POSTPROCESS), Embedded)
                CCOPT = -c -O1 -DINTEL_MODIFIED -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(MKL_MACRO_DEF)  $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF)  -Dmem_64BIT $(MKL_INCLUDE)
        endif

        LD = $(CC)
        LDOPT  = $(CFLAGS) -Wl,-static -Wall
        LDOPT  = $(CFLAGS)  -Wall


        #MISCLIB0   = -lpthread -lm ${FCDIR}/lib/libifcore.a ${FCDIR}/lib/libifport.a -ldl -lsvml $(MKL_LINK)
        #MISCLIB0   = -lpthread -lm ${CCDIR}/lib/libcprts.a ${CCDIR}/lib/libcxa.a ${CCDIR}/lib/libunwind.a ${FCDIR}/lib/libifcore.a ${FCDIR}/lib/libifport.a -ldl  ${CCDIR}/lib/libsvml.a $(MPI_HOME)/lib/libmpichf90.a $(MKL_LINK)
        MISCLIB0   = -lpthread -mkl ${CCDIR}/lib/libifcore.a ${CCDIR}/lib/libifport.a -ldl ${CCDIR}/lib/libsvml.a #$(MPI_HOME)/lib/libmpichf90.a
        MISCINCL0  =

        #XERCES specific options
        XERCES_CCOPT= $(CCOPT)
        XERCES_CC= icc -Kc++

        #XMLtoC specific options
        XMLtoC_CCOPT=$(CCOPT)
        XMLtoC_CC = icpc

        # NML-specific options
        NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT  -Dmem_64BIT
        NML_OPTIONS = $(CCOPT)

        FORTRAN   = $(MPI_HOME)/bin/mpiifort

        # LAPACK-specific options
        LAPACK_NOOPT  = -O2
        LAPACK_OPTS   = -O2
        LAPACK_LDOPTS = -O2

        # ARPACK-specific options
        ARPACK_OPTS   = -O2

        ##############################################################
        # TRACEMIN-specific options (these are set here only for choice of BUILD_TYPE: intel8,  intel9_64_mpi_intelfast (lear.rcac.purdue.edu))
        SCTRACEMIN_OPTS = -O2
        PCTRACEMIN_OPTS = -O2
        SQTRACEMIN_OPTS = -O2
        PQTRACEMIN_OPTS = -O2
        FORTRAN90 = $(MPI_HOME)/bin/mpif90  #tracemin fortran compiler
        ##############################################################

endif

# MPICH Intel optimized with SSE3
ifeq ($(BUILD_TYPE),intel19_64_mpi_intelfast)
        ICC_HOME=$(INTEL_CC_BASE)/intel64
        IFORT_HOME=$(INTEL_FC_BASE)/intel64
        MPI_HOME=$(INTEL_MPI_BASE)/intel64
        CCDIR=$(INTEL_CC_BASE)
        FCDIR=$(INTEL_FC_BASE)
        AR = xiar
        ifeq ($(MPI_HOME),)
                #@ echo "Error - Undefined environment variable, MPI_HOME"
        endif

        CCMKD = $(MPI_HOME)/bin/mpiicpc -MM
        # Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
        CC = $(MPI_HOME)/bin/mpiicpc

        CFLAGS = -O2 -xCORE-AVX2
        ####################################################################
        ##NOTE: -DINTEL_MODIFIED sets the INTEL optimizations ON,
                ##use -UINTEL_MODIFIED in CCOPT (next line) to turn them back OFF!!
        ####################################################################
        CCOPT = -c $(CFLAGS) -DINTEL_MODIFIED -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(MKL_MACRO_DEF)  $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF)  -Dmem_64BIT $(MKL_INCLUDE)

        ifeq ($(POSTPROCESS), Embedded)
                CCOPT = -c -O1 -DINTEL_MODIFIED -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(MKL_MACRO_DEF)  $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF)  -Dmem_64BIT $(MKL_INCLUDE)
        endif

        LD = $(CC)
        LDOPT  = $(CFLAGS) -Wl,-static -Wall
        LDOPT  = $(CFLAGS)  -Wall


        #MISCLIB0   = -lpthread -lm ${FCDIR}/lib/libifcore.a ${FCDIR}/lib/libifport.a -ldl -lsvml $(MKL_LINK)
        #MISCLIB0   = -lpthread -lm ${CCDIR}/lib/libcprts.a ${CCDIR}/lib/libcxa.a ${CCDIR}/lib/libunwind.a ${FCDIR}/lib/libifcore.a ${FCDIR}/lib/libifport.a -ldl  ${CCDIR}/lib/libsvml.a $(MPI_HOME)/lib/libmpichf90.a $(MKL_LINK)
        MISCLIB0   = -lpthread -mkl ${CCDIR}/../common/lib/libifcore.a ${CCDIR}/../common/lib/libifport.a -ldl ${CCDIR}/../common/lib/libsvml.a #$(MPI_HOME)/lib/libmpichf90.a
        MISCINCL0  =

        #XERCES specific options
        XERCES_CCOPT= $(CCOPT)
        XERCES_CC= icc -Kc++

        #XMLtoC specific options
        XMLtoC_CCOPT=$(CCOPT)
        XMLtoC_CC = icpc

        # NML-specific options
        NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT  -Dmem_64BIT
        NML_OPTIONS = $(CCOPT)

        FORTRAN   = $(MPI_HOME)/bin/mpiifort

        # LAPACK-specific options
        LAPACK_NOOPT  = -O2
        LAPACK_OPTS   = -O2
        LAPACK_LDOPTS = -O2

        # ARPACK-specific options
        ARPACK_OPTS   = -O2

        ##############################################################
        # TRACEMIN-specific options (these are set here only for choice of BUILD_TYPE: intel8,  intel9_64_mpi_intelfast (lear.rcac.purdue.edu))
        SCTRACEMIN_OPTS = -O2
        PCTRACEMIN_OPTS = -O2
        SQTRACEMIN_OPTS = -O2
        PQTRACEMIN_OPTS = -O2
        FORTRAN90 = $(MPI_HOME)/bin/mpif90  #tracemin fortran compiler
        ##############################################################

endif



# MPICH Intel optimized with SSE3
ifeq ($(BUILD_TYPE),steele_intel11_64_mpi_intelfast)

	CCDIR=/apps/rhel6/intel/composer_xe_2013.3.163/compiler/
        FCDIR=/apps/rhel6/intel/composer_xe_2013.3.163/compiler/
	AR = xiar
	ifeq ($(MPI_HOME),)
		#@ echo "Error - Undefined environment variable, MPI_HOME"
	endif

	CCMKD = $(MPI_HOME)/bin/mpiCC -MM    -I${CCDIR}/include/c++ 
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC = $(MPI_HOME)/bin/mpiCC -I${CCDIR}/include/c++

	CFLAGS = -O1 -Kc++ -axP
        ####################################################################
        ##NOTE: -DINTEL_MODIFIED sets the INTEL optimizations ON, 
		##use -UINTEL_MODIFIED in CCOPT (next line) to turn them back OFF!!
        ####################################################################
	CCOPT = -c $(CFLAGS) -DINTEL_MODIFIED -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(MKL_MACRO_DEF)  $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -I${CCDIR}/include/c++ $(HAM_STORE_DEF)  -Dmem_64BIT $(MKL_INCLUDE) -DMPICH_IGNORE_CXX_SEEK

	ifeq ($(POSTPROCESS), Embedded)
		CCOPT = -c -O1 -Kc++ -axP -DINTEL_MODIFIED -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(MKL_MACRO_DEF)  $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess -I${CCDIR}/include/c++ $(HAM_STORE_DEF)  -Dmem_64BIT $(MKL_INCLUDE)
	endif

	LD = $(CC)
	LDOPT  = $(CFLAGS) -Wl,-static -Wall
	LDOPT  = $(CFLAGS)  -Wall


	#MISCLIB0   = -lpthread -lm ${FCDIR}/lib/intel64/libifcore.a ${FCDIR}/lib/intel64/libifport.a -ldl -lsvml $(MKL_LINK)
	#MISCLIB0   = -lpthread -lm ${CCDIR}/lib/intel64/libcprts.a ${CCDIR}/lib/intel64/libcxa.a ${CCDIR}/lib/intel64/libunwind.a ${FCDIR}/lib/intel64/libifcore.a ${FCDIR}/lib/intel64/libifport.a -ldl  ${CCDIR}/lib/intel64/libsvml.a $(MPI_HOME)/lib/intel64/libmpichf90.a $(MKL_LINK)
        MISCLIB0   = -lpthread -lm ${CCDIR}/lib/intel64/libcxaguard.a ${FCDIR}/lib/intel64/libifcore.a ${FCDIR}/lib/intel64/libifport.a -ldl ${CCDIR}/lib/intel64/libsvml.a $(MPI_HOME)/lib/libmpichf90.a $(MKL_LINK)
	MISCINCL0  =

	#XERCES specific options
	XERCES_CCOPT= $(CCOPT)
	XERCES_CC= icc

  	#XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)
	XMLtoC_CC = icc

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT  -Dmem_64BIT
	NML_OPTIONS = $(CCOPT)

	FORTRAN   = $(MPI_HOME)/bin/mpif77

	# LAPACK-specific options
	LAPACK_NOOPT  = -O1 -mp -fp_port -mtune=pentium4 -xW -axP
	LAPACK_OPTS   = -O1 -mp -fp_port -mtune=pentium4 -xW -axP
	LAPACK_LDOPTS = -O1 -mp -fp_port -mtune=pentium4 -xW -axP -Vaxlib

	# ARPACK-specific options
	ARPACK_OPTS   = -O1 -mp -mtune=pentium4 -axP -xW

	##############################################################
	# TRACEMIN-specific options (these are set here only for choice of BUILD_TYPE: intel8,  intel9_64_mpi_intelfast (lear.rcac.purdue.edu))  
	SCTRACEMIN_OPTS = -O1
 	PCTRACEMIN_OPTS = -O1
	SQTRACEMIN_OPTS = -O1
	PQTRACEMIN_OPTS = -O1
	FORTRAN90 = $(MPI_HOME)/bin/mpif90  #tracemin fortran compiler
	##############################################################

endif


# MPICH Intel optimized with SSE3
ifeq ($(BUILD_TYPE),steele_intel10_64_mpi_intelfast)

	CCDIR=/opt/intel/cce/10.1.015
	FCDIR=/opt/intel/fce/10.1.015
	AR = xiar
	ifeq ($(MPI_HOME),)
		#@ echo "Error - Undefined environment variable, MPI_HOME"
	endif

	CCMKD = $(MPI_HOME)/bin/mpiCC -MM    -I${CCDIR}/include/c++ 
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC = $(MPI_HOME)/bin/mpiCC -I${CCDIR}/include/c++

	CFLAGS = -O1 -Kc++ -axP
        ####################################################################
        ##NOTE: -DINTEL_MODIFIED sets the INTEL optimizations ON, 
		##use -UINTEL_MODIFIED in CCOPT (next line) to turn them back OFF!!
        ####################################################################
	CCOPT = -c $(CFLAGS) -DINTEL_MODIFIED -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(MKL_MACRO_DEF)  $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -I${CCDIR}/include/c++ $(HAM_STORE_DEF)  -Dmem_64BIT $(MKL_INCLUDE) -DMPICH_IGNORE_CXX_SEEK

	ifeq ($(POSTPROCESS), Embedded)
		CCOPT = -c -O1 -Kc++ -axP -DINTEL_MODIFIED -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(MKL_MACRO_DEF)  $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess -I${CCDIR}/include/c++ $(HAM_STORE_DEF)  -Dmem_64BIT $(MKL_INCLUDE)
	endif

	LD = $(CC)
	LDOPT  = $(CFLAGS) -Wl,-static -Wall
	LDOPT  = $(CFLAGS)  -Wall


	#MISCLIB0   = -lpthread -lm ${FCDIR}/lib/libifcore.a ${FCDIR}/lib/libifport.a -ldl -lsvml $(MKL_LINK)
	#MISCLIB0   = -lpthread -lm ${CCDIR}/lib/libcprts.a ${CCDIR}/lib/libcxa.a ${CCDIR}/lib/libunwind.a ${FCDIR}/lib/libifcore.a ${FCDIR}/lib/libifport.a -ldl  ${CCDIR}/lib/libsvml.a $(MPI_HOME)/lib/libmpichf90.a $(MKL_LINK)
        MISCLIB0   = -lpthread -lm ${CCDIR}/lib/libcxaguard.a ${FCDIR}/lib/libifcore.a ${FCDIR}/lib/libifport.a -ldl ${CCDIR}/lib/libsvml.a $(MPI_HOME)/lib/libmpichf90.a $(MKL_LINK)
	MISCINCL0  =

	#XERCES specific options
	XERCES_CCOPT= $(CCOPT)
	XERCES_CC= icc

  	#XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)
	XMLtoC_CC = icc

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT  -Dmem_64BIT
	NML_OPTIONS = $(CCOPT)

	FORTRAN   = $(MPI_HOME)/bin/mpif77

	# LAPACK-specific options
	LAPACK_NOOPT  = -O1 -mp -fp_port -mtune=pentium4 -xW -axP
	LAPACK_OPTS   = -O1 -mp -fp_port -mtune=pentium4 -xW -axP
	LAPACK_LDOPTS = -O1 -mp -fp_port -mtune=pentium4 -xW -axP -Vaxlib

	# ARPACK-specific options
	ARPACK_OPTS   = -O1 -mp -mtune=pentium4 -axP -xW

	##############################################################
	# TRACEMIN-specific options (these are set here only for choice of BUILD_TYPE: intel8,  intel9_64_mpi_intelfast (lear.rcac.purdue.edu))  
	SCTRACEMIN_OPTS = -O1
 	PCTRACEMIN_OPTS = -O1
	SQTRACEMIN_OPTS = -O1
	PQTRACEMIN_OPTS = -O1
	FORTRAN90 = $(MPI_HOME)/bin/mpif90  #tracemin fortran compiler
	##############################################################

endif



# MPICH Intel optimized with SSE3
ifeq ($(BUILD_TYPE),intel9_64_mpi_intelfast_instr)
	CCDIR=/opt/intel/cce/9.1.045
	FCDIR=/opt/intel/fce/9.1.040
	AR = xiar
	ifeq ($(MPI_HOME),)
		#@ echo "Error - Undefined environment variable, MPI_HOME"
	endif

	CCMKD = $(MPI_HOME)/bin/mpiCC -MM    -I${CCDIR}/include/c++ 
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC = $(MPI_HOME)/bin/mpiCC -I${CCDIR}/include/c++

	CFLAGS = -O2 -Kc++ -axP
        ####################################################################
        ##NOTE: -DINTEL_MODIFIED sets the INTEL optimizations ON, 
		##use -UINTEL_MODIFIED in CCOPT (next line) to turn them back OFF!!
        ####################################################################
	CCOPT = -c $(CFLAGS) -DINTEL_MODIFIED -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(MKL_MACRO_DEF)  $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -I${CCDIR}/include/c++ $(HAM_STORE_DEF)  -Dmem_64BIT $(MKL_INCLUDE) -DMPI_INSTR -I../../nemo3d/src/base -I../../nemo3d/src/math

	ifeq ($(POSTPROCESS), Embedded)
		CCOPT = -c -O2 -Kc++ -axP -DINTEL_MODIFIED -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(MKL_MACRO_DEF)  $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess -I${CCDIR}/include/c++ $(HAM_STORE_DEF)  -Dmem_64BIT $(MKL_INCLUDE)
	endif

	LD = $(CC)
	LDOPT  = $(CFLAGS) -Wl,-static -Wall
	LDOPT  = $(CFLAGS)  -Wall


	MISCLIB0   = -lpthread -lm ${FCDIR}/lib/libifcore.a ${FCDIR}/lib/libifport.a -ldl -lsvml $(MKL_LINK)
	MISCLIB0   = -lpthread -lm ${CCDIR}/lib/libcprts.a ${CCDIR}/lib/libcxa.a ${CCDIR}/lib/libunwind.a ${FCDIR}/lib/libifcore.a ${FCDIR}/lib/libifport.a -ldl  ${CCDIR}/lib/libsvml.a $(MPI_HOME)/lib/libmpichf90.a $(MKL_LINK)
	MISCINCL0  =

	#XERCES specific options
	XERCES_CCOPT= $(CCOPT) -c -O1 -Kc++ -axP
	XERCES_CC= icc

  	#XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)
	XMLtoC_CC = icc

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT  -Dmem_64BIT
	NML_OPTIONS = $(CCOPT)

	FORTRAN   = $(MPI_HOME)/bin/mpif77

	# LAPACK-specific options
	LAPACK_NOOPT  = -O3 -mp -fp_port
	LAPACK_OPTS   = -O3 -mp -fp_port -axP -tpp7 
	LAPACK_LDOPTS = -O3 -mp -fp_port -axP -tpp7 -Vaxlib

	# ARPACK-specific options
	ARPACK_OPTS   = -O3 -mp -axP -tpp7 -xW

	##############################################################
	# TRACEMIN-specific options (these are set here only for choice of BUILD_TYPE: intel8,  intel9_64_mpi_intelfast (lear.rcac.purdue.edu))  
	SCTRACEMIN_OPTS = -O3
 	PCTRACEMIN_OPTS = -O3
	SQTRACEMIN_OPTS = -O3
	PQTRACEMIN_OPTS = -O3
	FORTRAN90 = $(MPI_HOME)/bin/mpif90  #tracemin fortran compiler
	##############################################################
endif


# MPICH Intel 
ifeq ($(BUILD_TYPE),intel9_64_mpi_inteldebug)
	CCDIR=/opt/intel/cce/9.0.025
	FCDIR=/opt/intel/fce/9.0.025
	AR = xiar
	ifeq ($(MPI_HOME),)
		#@ echo "Error - Undefined environment variable, MPI_HOME"
	endif

	 CCMKD = $(MPI_HOME)/bin/mpiCC -MM    -I${CCDIR}/include/c++ 
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC = $(MPI_HOME)/bin/mpiCC -I${CCDIR}/include/c++ 

	CFLAGS = -g -Kc++ 

	    CCOPT = -c $(CFLAGS) -DINTEL_MODIFIED -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING   $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -I${CCDIR}/include/c++ $(HAM_STORE_DEF)  -Dmem_64BIT
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT = -c -g -Kc++ -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING   $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF)  -Dmem_64BIT
	endif

	LD = $(CC)
	    LDOPT  = $(CFLAGS) -Wl,-static -Wall
	    LDOPT  = $(CFLAGS)  -Wall


	    MISCLIB0   = -lpthread -lm ${FCDIR}/lib/libifcore.a ${FCDIR}/lib/libifport.a -ldl -lsvml 
	    MISCLIB0   = -lpthread -lm ${CCDIR}/lib/libcprts.a ${CCDIR}/lib/libcxa.a ${CCDIR}/lib/libunwind.a ${FCDIR}/lib/libifcore.a ${FCDIR}/lib/libifport.a -ldl  ${CCDIR}/lib/libsvml.a $(MPI_HOME)/lib/libmpichf90.a
	MISCINCL0  =

	#XERCES specific options
	XERCES_CCOPT= -c -g -Kc++ 
	XERCES_CC= icc

        #XMLtoC specific options
        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)
	XMLtoC_CC = icc

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT  -Dmem_64BIT
	    NML_OPTIONS = $(CCOPT)

	FORTRAN = $(MPI_HOME)/bin/mpif77

	# LAPACK-specific options
	LAPACK_NOOPT  = -g -mp -fp_port
	LAPACK_OPTS   = -g -mp -fp_port  
	LAPACK_LDOPTS = -g -mp -fp_port -Vaxlib

	# ARPACK-specific options
	ARPACK_OPTS   = -g -mp 

	##############################################################
	# TRACEMIN-specific options (these are set here only for choice of BUILD_TYPE: intel8,  intel9_64_mpi_intelfast (lear.rcac.purdue.edu))  
	SCTRACEMIN_OPTS = -g
   	PCTRACEMIN_OPTS = -g
	SQTRACEMIN_OPTS = -g
	PQTRACEMIN_OPTS = -g
	FORTRAN90 = $(MPI_HOME)/bin/mpif90  #tracemin fortran compiler
	##############################################################
endif


# MPICH Intel 
ifeq ($(BUILD_TYPE),intel9_32_mpi_intelfast)

	CCDIR=/opt/intel/cc/9.0.032
	FCDIR=/opt/intel/fc/9.0.032
	AR = xiar
	ifeq ($(MPI_HOME),)
		#@ echo "Error - Undefined environment variable, MPI_HOME"
	endif

	 CCMKD = $(MPI_HOME)/bin/mpiCC -MM    -I${CCDIR}/include/c++ 
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC = $(MPI_HOME)/bin/mpiCC -I${CCDIR}/include/c++ 

	CFLAGS = -O2 -Kc++ -axW

	CCOPT = -c $(CFLAGS) -UINTEL_MODIFIED -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING   $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -I${CCDIR}/include/c++ $(HAM_STORE_DEF)  -Umem_64BIT
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT = -c -O2 -Kc++ -axW -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING   $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess -I${CCDIR}/include/c++ $(HAM_STORE_DEF)  -Umem_64BIT
	endif

	LD = $(CC)
	LDOPT  = $(CFLAGS) -Wl,-static -Wall
	LDOPT  = $(CFLAGS)  -Wall


	MISCLIB0   = -lpthread -lm ${FCDIR}/lib/libifcore.a ${FCDIR}/lib/libifport.a -ldl -lsvml 
	MISCLIB0   = -lpthread -lm ${CCDIR}/lib/libcprts.a ${CCDIR}/lib/libcxa.a ${CCDIR}/lib/libunwind.a ${FCDIR}/lib/libifcore.a ${FCDIR}/lib/libifport.a -ldl  ${CCDIR}/lib/libsvml.a
	MISCINCL0  =

	#XERCES specific options
	XERCES_CCOPT= ${CCOPT} -c -O1 -Kc++ -axW
	XERCES_CC= icc

        #XMLtoC specific options
        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)
	XMLtoC_CC = icc

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT  -Umem_64BIT
	NML_OPTIONS = $(CCOPT)

	FORTRAN = $(MPI_HOME)/bin/mpif77

	# LAPACK-specific options
	LAPACK_NOOPT  = -O3 -mp -fp_port
	LAPACK_OPTS   = -O3 -mp -fp_port -xW -tpp7 
	LAPACK_LDOPTS = -O3 -mp -fp_port -xW -tpp7 -Vaxlib

	# ARPACK-specific options
	ARPACK_OPTS   = -O3 -mp -tpp7 -xW

  	##############################################################
  	# TRACEMIN-specific options (these are set here only for choice of BUILD_TYPE equal to intel9_64_mpi_intelfast (lear.rcac.purdue.edu))
	SCTRACEMIN_OPTS = -O3
	PCTRACEMIN_OPTS = -O3
	SQTRACEMIN_OPTS = -O3
	PQTRACEMIN_OPTS = -O3
	FORTRAN90 = $(MPI_HOME)/bin/mpif90  #tracemin fortran compiler
  	##############################################################
endif


# MPICH Intel PentiumIII optimized MPI Static no machine specific options
ifeq ($(BUILD_TYPE),intel_mpi_static)

        AR = xiar
        ifeq ($(MPI_HOME),)
                #@ echo "Error - Undefined environment variable, MPI_HOME"
        endif

        CCMKD = $(MPI_HOME)/bin/mpicxx -MM
        # Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
        CC = $(MPI_HOME)/bin/mpicxx
        CCOPT = -c -O3 -Kc++ -axK -tpp6 -xK -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) 
        ifeq ($(POSTPROCESS), Embedded)
                CCOPT = -c -O3 -Kc++ -axK -tpp6 -xK -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
        endif
        LD = $(MPI_HOME)/bin/mpicxx
        LDOPT  = -O3 -Wl,-static -Wall -axK -tpp6 -xK

        MISCLIB0   = -ldl -lpthread -lm -lPEPCF90 -lIEPCF90 -lF90
        MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

        # NML-specific options
        NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
        NML_OPTIONS = -O3 -Wall -ansi -march=pentiumiii -mcpu=pentiumpro -axK -tpp6 -xK

        FORTRAN = $(MPI_HOME)/bin/mpif77

        # LAPACK-specific options
        LAPACK_NOOPT  = -O3 -mp -fp_port
        LAPACK_OPTS   = -O3 -mp -fp_port -axK -tpp6 -xK
        LAPACK_LDOPTS = -O3 -mp -fp_port -axK -tpp6 -xK -Vaxlib

        # ARPACK-specific options
        ARPACK_OPTS   = -O3 -mp -axK -tpp6 -xK

endif


# MPICH Intel PentiumIII optimized MPI Static
ifeq ($(BUILD_TYPE),intel_piii_mpi_static)

	AR = xiar
	ifeq ($(MPI_HOME),)
		#@ echo "Error - Undefined environment variable, MPI_HOME"
	endif

	CCMKD = $(MPI_HOME)/bin/mpicxx -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC = $(MPI_HOME)/bin/mpicxx
	CCOPT = -c -O3 -Kc++ -march=pentiumiii -mcpu=pentiumpro -axK -tpp6 -xK -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) 
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT = -c -O3 -Kc++ -march=pentiumiii -mcpu=pentiumpro -axK -tpp6 -xK -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	endif
	LD = $(MPI_HOME)/bin/mpicxx
	LDOPT  = -O3 -Wl,-static -Wall -march=pentiumiii -mcpu=pentiumpro -axK -tpp6 -xK

	MISCLIB0   = -ldl -lpthread -lm -lPEPCF90 -lIEPCF90 -lF90
	MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
	NML_OPTIONS = -O3 -Wall -ansi -march=pentiumiii -mcpu=pentiumpro -axK -tpp6 -xK

	FORTRAN = $(MPI_HOME)/bin/mpif77

	# LAPACK-specific options
	LAPACK_NOOPT  = -O3 -mp -fp_port
	LAPACK_OPTS   = -O3 -mp -fp_port -axK -tpp6 -xK
	LAPACK_LDOPTS = -O3 -mp -fp_port -axK -tpp6 -xK -Vaxlib

	# ARPACK-specific options
	ARPACK_OPTS   = -O3 -mp -axK -tpp6 -xK

endif

# MPICH Intel PentiumIII optimized Static
ifeq ($(BUILD_TYPE),intel_piii_static)

	AR = xiar

	CCMKD =  /home/packages/intel/compiler70/ia32/bin/icpc -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC = /home/packages/intel/compiler70/ia32/bin/icpc
	CCOPT = -c -O3 -Kc++ -march=pentiumiii -mcpu=pentiumpro -axK -tpp6 -xK -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DFAKE_MPI -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DNoGnuStatementExpre $(HAM_STORE_DEF) 
	LIB_PARPACK =
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT = -c -O3 -Kc++ -march=pentiumiii -mcpu=pentiumpro -axK -tpp6 -xK -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DFAKE_MPI -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DNoGnuStatementExpre -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	endif
	LD = /home/packages/intel/compiler70/ia32/bin/icpc
	LDOPT  = -O3 -Wl,-static -Wall -march=pentiumiii -mcpu=pentiumpro -axK -tpp6 -xK

	MISCLIB0   = -ldl -lpthread -lm -lPEPCF90 -lIEPCF90 -lF90
	MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
	NML_OPTIONS = -O3 -Wall -ansi -march=pentiumiii -mcpu=pentiumpro -axK -tpp6 -xK

	FORTRAN = /home/packages/intel/compiler70/ia32/bin/ifc

	# LAPACK-specific options
	LAPACK_NOOPT  = -O3 -mp -fp_port
	LAPACK_OPTS   = -O3 -mp -fp_port -axK -tpp6 -xK
	LAPACK_LDOPTS = -O3 -mp -fp_port -axK -tpp6 -xK -Vaxlib

	# ARPACK-specific options
	ARPACK_OPTS   = -O3 -mp -fp_port -axK -tpp6 -xK

endif





# GNU Debug Static
ifeq ($(BUILD_TYPE),itanium2_gnu_debug)

	CCMKD  = g++ -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY -DSANS_PY 
	CC     = g++
	CCOPT  = -c -g -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE  -DFAKE_MPI -DMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF)  $(HAM_STORE_DEF) 
	LIB_PARPACK =
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT  = -c -g -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE  -DFAKE_MPI -DMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	endif

	LD     = g++
	LDOPT  = -g 

	MISCLIB0   = -ldl -lpthread -lm -lg2c
	MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT -Dmem_64BIT
	NML_OPTIONS = -g -Wall -ansi -pedantic

	# LAPACK-specific options
	FORTRAN = g77
	LAPACK_NOOPT =
	LAPACK_OPTS   = -g -funderscoring
	LAPACK_LDOPTS = -g -funderscoring
	ARPACK_OPTS   = -g -funderscoring

	# If -DNO_PARPACK is used in CCOPT, then we need to make sure that the PARPACK library doesn't get linked in
	LIB_PARPACK =

endif




# GNU Debug Static
ifeq ($(BUILD_TYPE),itanium2_gnu_fast)

	CCMKD  = g++ -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY -DSANS_PY 
	CC     = g++
	CCOPT  = -c -O3  -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE  -DFAKE_MPI -DMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF)  $(HAM_STORE_DEF) 
	LIB_PARPACK =
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT  = -c -O3 -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE  -DFAKE_MPI -DMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	endif

	LD     = g++
	LDOPT  = -O3  

	MISCLIB0   = -ldl -lpthread -lm -lg2c
	MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT -Dmem_64BIT
	NML_OPTIONS =  -Wall -ansi -pedantic

	# LAPACK-specific options
	FORTRAN = g77
	LAPACK_NOOPT = -O1
	LAPACK_OPTS   =  -O2 -funderscoring
	LAPACK_LDOPTS =  -O2 -funderscoring
	ARPACK_OPTS   =  -O2 -funderscoring

	# If -DNO_PARPACK is used in CCOPT, then we need to make sure that the PARPACK library doesn't get linked in
	LIB_PARPACK =

endif



# GNU Debug Static
ifeq ($(BUILD_TYPE),opteron_pgi_fast)

	CCMKD  = g++ -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY -DSANS_PY 
	CC     = pgCC
	CCOPT  = -c -fastsse -DNoGnuStatementExpre  -D$(OSTYPE) -DFORTRAN_UNDERSCORE  -DFAKE_MPI -DMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF)  $(HAM_STORE_DEF) 
	LIB_PARPACK =
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT  = -c -fastsse -D$(OSTYPE) -DFORTRAN_UNDERSCORE  -DFAKE_MPI -DMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	endif

	LD     = pgCC
	LDOPT  = -fastsse  

	# MISCLIB0   = -ldl -lpthread -lm -lg2c
	MISCLIB0   = -ldl -lpthread -lm -lpgftnrtl -lstd
	MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT -Dmem_64BIT
	NML_OPTIONS =  -Wall -ansi -pedantic

	# LAPACK-specific options
	FORTRAN = pgf77
	LAPACK_NOOPT = -O1
	LAPACK_OPTS   =  -fastsse 
	LAPACK_LDOPTS =  -fastsse 
	ARPACK_OPTS   =  -fastsse

	# If -DNO_PARPACK is used in CCOPT, then we need to make sure that the PARPACK library doesn't get linked in
	LIB_PARPACK =

endif


ifeq ($(BUILD_TYPE),crayxt_backend_pgi_fast)

	CCMKD  = g++ -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY -DSANS_PY 
	CC     = CC
	CCOPT  = -c -fastsse -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF)  $(HAM_STORE_DEF) -Dmem_64BIT -DDONTBUILDFORMATCONVERTER -DUNSET_SEEK -DMPI_INSTR -I../../nemo3d/src/base -I../../nemo3d/src/math
	LIB_PARPACK =
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT  = -c -fastsse -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	endif

	LD     = CC
	LDOPT  = -fastsse  

	#MISCLIB0   = -ldl -lpthread -lm -lpgftnrtl -lstd
	MISCLIB0   = -ldl -lm -lpgf90 -lpgftnrtl -lpghpf_prof -lpghpf_mpi -lpghpf -lpghpf2 -lstd
	MISCINCL0  =

	XERCES_CCOPT=$(CCOPT) -DAPP_NO_THREADS -DNOLINUXLIMITS
        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT -Dmem_64BIT
	NML_OPTIONS =

	# LAPACK-specific options
	FORTRAN = ftn
	LAPACK_NOOPT = -O1
	LAPACK_OPTS   =  -fastsse 
	LAPACK_LDOPTS =  -fastsse 
	ARPACK_OPTS   =  -fastsse

	# If -DNO_PARPACK is used in CCOPT, then we need to make sure that the PARPACK library doesn't get linked in
	LIB_PARPACK =

endif


ifeq ($(BUILD_TYPE),crayxt_frontend_pgi_fast)

	CCMKD  = g++ -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY -DSANS_PY 
	CC     = pgCC
	CCOPT  = -c -fastsse -D$(OSTYPE) -DFORTRAN_UNDERSCORE  -DFAKE_MPI -DMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF)  $(HAM_STORE_DEF) 
	LIB_PARPACK =
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT  = -c -fastsse -D$(OSTYPE) -DFORTRAN_UNDERSCORE  -DFAKE_MPI -DMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	endif

	LD     = pgCC
	LDOPT  = -fastsse  

	# MISCLIB0   = -ldl -lpthread -lm -lg2c
	MISCLIB0   = -ldl -lpthread -lm -lpgftnrtl -lstd
	MISCINCL0  =

	XERCES_CC=$(CC)
	XERCES_LD=$(LD)

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT -Dmem_64BIT
	NML_OPTIONS =

	# LAPACK-specific options
	FORTRAN = ftn
	LAPACK_NOOPT = -O1
	LAPACK_OPTS   =  -fastsse 
	LAPACK_LDOPTS =  -fastsse 
	ARPACK_OPTS   =  -fastsse

	# If -DNO_PARPACK is used in CCOPT, then we need to make sure that the PARPACK library doesn't get linked in
	LIB_PARPACK =

endif



ifeq ($(BUILD_TYPE),crayxt_backend_pgi_debug)

	CCMKD  = g++ -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY -DSANS_PY 
	CC     = CC
	CCOPT  = -c -g -O0 -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF)  $(HAM_STORE_DEF) -Dmem_64BIT -DDONTBUILDFORMATCONVERTER -DUNSET_SEEK -DMPI_INSTR -I../../nemo3d/src/base -I../../nemo3d/src/math
	LIB_PARPACK =
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT  = -c -g -O0 -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	endif

	LD     = CC
	LDOPT  = -g -O0  

	#MISCLIB0   = -ldl -lpthread -lm -lpgftnrtl -lstd
	MISCLIB0   = -ldl -lm -lpgf90 -lpgftnrtl -lpghpf_prof -lpghpf_mpi -lpghpf -lpghpf2 -lstd
	MISCINCL0  =

	XERCES_CCOPT=$(CCOPT) -DAPP_NO_THREADS -DNOLINUXLIMITS
        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT -Dmem_64BIT
	NML_OPTIONS =

	# LAPACK-specific options
	FORTRAN = ftn
	LAPACK_NOOPT = -g -O0
	LAPACK_OPTS   =  -g -O0 
	LAPACK_LDOPTS =  -g -O0 
	ARPACK_OPTS   =  -g -O0

	# If -DNO_PARPACK is used in CCOPT, then we need to make sure that the PARPACK library doesn't get linked in
	LIB_PARPACK =

endif


ifeq ($(BUILD_TYPE),crayxt_frontend_pgi_debug)

	CCMKD  = g++ -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY -DSANS_PY 
	CC     = pgCC
	CCOPT  = -c -g -O0 -D$(OSTYPE) -DFORTRAN_UNDERSCORE  -DFAKE_MPI -DMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF)  $(HAM_STORE_DEF) 
	LIB_PARPACK =
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT  = -c -g -O0 -D$(OSTYPE) -DFORTRAN_UNDERSCORE  -DFAKE_MPI -DMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	endif

	LD     = pgCC
	LDOPT  = -g -O0  

	# MISCLIB0   = -ldl -lpthread -lm -lg2c
	MISCLIB0   = -ldl -lpthread -lm -lpgftnrtl -lstd
	MISCINCL0  =

	XERCES_CC=$(CC)
	XERCES_LD=$(LD)

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT -Dmem_64BIT
	NML_OPTIONS =

	# LAPACK-specific options
	FORTRAN = ftn
	LAPACK_NOOPT = -g -O0
	LAPACK_OPTS   =  -g -O0 
	LAPACK_LDOPTS =  -g -O0 
	ARPACK_OPTS   =  -g -O0

	# If -DNO_PARPACK is used in CCOPT, then we need to make sure that the PARPACK library doesn't get linked in
	LIB_PARPACK =

endif



ifeq ($(BUILD_TYPE),crayxt_backend_pgi_fast_profile)

	CCMKD  = g++ -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY -DSANS_PY 
	CC     = CC
	CCOPT  = -c -fastsse -Mprof=func -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF)  $(HAM_STORE_DEF) -Dmem_64BIT -DDONTBUILDFORMATCONVERTER -DUNSET_SEEK -DMPI_INSTR -I../../nemo3d/src/base -I../../nemo3d/src/math
	LIB_PARPACK =
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT  = -c -fastsse -Mprof=func -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	endif

	LD     = CC
	LDOPT  = -fastsse -Mprof=func

	#MISCLIB0   = -ldl -lpthread -lm -lpgftnrtl -lstd
#	MISCLIB0   = -ldl -lm -lpgf90 -lpgftnrtl -lpghpf_prof -lpghpf_mpi -lpghpf -lpghpf2 -lstd
	MISCLIB0   = -ldl -lm -lpgf90 -lpgftnrtl              -lpghpf_mpi -lpghpf -lpghpf2 -lstd
	MISCINCL0  =

	XERCES_CCOPT=$(CCOPT) -DAPP_NO_THREADS -DNOLINUXLIMITS
        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT -Dmem_64BIT
	NML_OPTIONS = -Mprof=func

	# LAPACK-specific options
	FORTRAN = ftn
	LAPACK_NOOPT  = -O1 -Mprof=func
	LAPACK_OPTS   = -fastsse -Mprof=func
	LAPACK_LDOPTS = -fastsse -Mprof=func
	ARPACK_OPTS   = -fastsse -Mprof=func

	# If -DNO_PARPACK is used in CCOPT, then we need to make sure that the PARPACK library doesn't get linked in
	LIB_PARPACK =

endif


ifeq ($(BUILD_TYPE),crayxt_frontend_pgi_fast_profile)

	CCMKD  = g++ -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY -DSANS_PY 
	CC     = pgCC
	CCOPT  = -c -fastsse -Mprof=func -D$(OSTYPE) -DFORTRAN_UNDERSCORE  -DFAKE_MPI -DMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF)  $(HAM_STORE_DEF) 
	LIB_PARPACK =
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT  = -c -fastsse -Mprof=func -D$(OSTYPE) -DFORTRAN_UNDERSCORE  -DFAKE_MPI -DMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	endif

	LD     = pgCC
	LDOPT  = -fastsse -Mprof=func  

	# MISCLIB0   = -ldl -lpthread -lm -lg2c
	MISCLIB0   = -ldl -lpthread -lm -lpgftnrtl -lstd
	MISCINCL0  =

	XERCES_CC=$(CC)
	XERCES_LD=$(LD)

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT -Dmem_64BIT
	NML_OPTIONS = -Mprof=func

	# LAPACK-specific options
	FORTRAN = ftn
	LAPACK_NOOPT  = -O1 -Mprof=func
	LAPACK_OPTS   = -fastsse -Mprof=func
	LAPACK_LDOPTS = -fastsse -Mprof=func
	ARPACK_OPTS   = -fastsse -Mprof=func

	# If -DNO_PARPACK is used in CCOPT, then we need to make sure that the PARPACK library doesn't get linked in
	LIB_PARPACK =

endif






# 
ifeq ($(BUILD_TYPE),itanium2_intel8_static)

	AR = xiar

	CCMKD =  icpc -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC = icpc
	#  need to compile the XMLtoC with just -O0
	CCOPT = -c -O3 -Kc++   -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DFAKE_MPI -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DNoGnuStatementExpre -DMPI_TIMING  $(HAM_STORE_DEF) 
	LIB_PARPACK =
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT = -c -O3 -Kc++   -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DFAKE_MPI -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DNoGnuStatementExpre -DEmbedded_Postprocess -DMPI_TIMING  $(HAM_STORE_DEF) 
	endif
	LD = icpc
	LDOPT  = -O3 -Wall   
	LDOPT  = -O3 -Wl,-static -Wall   

	MISCLIB0   = -ldl -lpthread -lm -lPEPCF90 -lIEPCF90 -lF90
	#  libifcore  libifcoremt
	MISCLIB0   = -lpthread -lm /opt/intel_fc_80/lib/libifcore.a /opt/intel_fc_80/lib/libifport.a -ldl
	MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT = -c -O0 -Kc++   -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DFAKE_MPI -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DNoGnuStatementExpre -DMPI_TIMING  $(HAM_STORE_DEF) 

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT 
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
	NML_OPTIONS = -O3 -Kc++ -Wall -ansi  -tpp6 

	FORTRAN = ifort

	# LAPACK-specific options
	LAPACK_NOOPT  = -O3 -mp
	LAPACK_OPTS   = -O3 -mp   
	LAPACK_LDOPTS = -O3 -mp    -Vaxlib

	# ARPACK-specific options
	ARPACK_OPTS   = -O3 -mp  

endif

# MPICH Intel version 7 compilers Itanium2 optimized MPI Static no machine specific options
ifeq ($(BUILD_TYPE),itanium2_intel7_mpi_static)
        AR = xiar
        ifeq ($(MPI_HOME),)
                #@ echo "Error - Undefined environment variable, MPI_HOME"
        endif
        CCMKD = $(MPI_HOME)/bin/mpiCC -MM
        # Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
        CC = $(MPI_HOME)/bin/mpiCC
        CCOPT = -c -O3 -Kc++ -tpp6 -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) 
        ifeq ($(POSTPROCESS), Embedded)
                CCOPT = -c -O3 -Kc++ -tpp6 -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
        endif
        LD = $(MPI_HOME)/bin/mpiCC
        LDOPT  = -O3 -Wl,-static -Wall -tpp6

        MISCLIB0   = -ldl -lpthread -lm -lPEPCF90 -lIEPCF90 -lF90
        MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

        # NML-specific options
        NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
        NML_OPTIONS = -O3 -Kc++ -Wall -ansi -tpp6 -I/users/ux453783/DEVELOP/NEMO3D/nml/include

        FORTRAN = $(MPI_HOME)/bin/mpif90

        # LAPACK-specific options
        LAPACK_NOOPT  = -O3 -mp
        LAPACK_OPTS   = -O3 -mp -tpp6
        LAPACK_LDOPTS = -O3 -mp -tpp6 -Vaxlib

        # ARPACK-specific options
        ARPACK_OPTS   = -O3 -mp -tpp6

endif

# MPICH Intel version 8 compilers Itanium2 optimized MPI Static no machine specific options
ifeq ($(BUILD_TYPE),itanium2_intel8_mpi_static)
        AR = xiar
        ifeq ($(MPI_HOME),)
                #@ echo "Error - Undefined environment variable, MPI_HOME"
        endif
        CCMKD = $(MPI_HOME)/bin/mpiCC -MM -I/opt/intel_cc_80/include/c++
        # Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
        CC = $(MPI_HOME)/bin/mpiCC -I/opt/intel_cc_80/include/c++
        CC = $(MPI_HOME)/bin/mpiCC 
        CCOPT = -c -O3 -Kc++  -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DNoGnuStatementExpre -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) 
        ifeq ($(POSTPROCESS), Embedded)
                CCOPT = -c -O3 -Kc++  -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DNoGnuStatementExpre -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
        endif
        LD = $(MPI_HOME)/bin/mpiCC
        LDOPT  = -O3 -Wl,-static -Wall 

        MISCLIB0   = -ldl -lpthread -lm -lPEPCF90 -lIEPCF90 -lF90
	MISCLIB0   = -lpthread -lm /opt/intel_fc_80/lib/libifcore.a /opt/intel_fc_80/lib/libifport.a -ldl
        MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT= -c -O0 -Kc++  -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DNoGnuStatementExpre -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) 

        # NML-specific options
        NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
        NML_OPTIONS = -O3 -Kc++ -Wall -DNoGnuStatementExpre -ansi 

        FORTRAN = $(MPI_HOME)/bin/mpif77 -intconstant

        # LAPACK-specific options
        LAPACK_NOOPT  = -O3 -mp
        LAPACK_OPTS   = -O3 -mp 
        LAPACK_LDOPTS = -O3 -mp  -Vaxlib

        # ARPACK-specific options
        ARPACK_OPTS   = -O3 -mp 

endif

# MPICH-G2 Intel Itanium2 optimized MPI Static no machine specific options
ifeq ($(BUILD_TYPE),itanium2_mpichg2_static)
        AR = xiar
        ifeq ($(MPI_HOME),)
                #@ echo "Error - Undefined environment variable, MPI_HOME"
        endif
        CCMKD = $(MPI_HOME)/bin/mpiCC -MM
        # Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
        CC = $(MPI_HOME)/bin/mpiCC
        CCOPT = -c -O3 -Kc++ -tpp6 -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) 
        ifeq ($(POSTPROCESS), Embedded)
                CCOPT = -c -O3 -Kc++ -tpp6 -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
        endif
        LD = $(MPI_HOME)/bin/mpiCC
        LDOPT  = -O3 -Wl,-static -Wall -tpp6

        MISCLIB0   = -ldl -lpthread -lm -lPEPCF90 -lIEPCF90 -lF90
        MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

        # NML-specific options
        NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
        NML_OPTIONS = -O3 -Kc++ -Wall -ansi -tpp6 -I/users/ux453783/DEVELOP/NEMO3D/nml/include

        FORTRAN = $(MPI_HOME)/bin/mpif90

        # LAPACK-specific options
        LAPACK_NOOPT  = -O3 -mp
        LAPACK_OPTS   = -O3 -mp -tpp6
        LAPACK_LDOPTS = -O3 -mp -tpp6 -Vaxlib

        # ARPACK-specific options
        ARPACK_OPTS   = -O3 -mp -tpp6

endif


# MACOSX GNU Debug
ifeq ($(BUILD_TYPE),macosx_intel_g95_debug)

	CCMKD  = g++ -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC     = g++
	CCOPT  = -c -g -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE  -DFAKE_MPI -DMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) 
	LIB_PARPACK =
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT  = -c -g -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE  -DFAKE_MPI -DMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	LIB_PARPACK =
	endif
	LD     = g++
	LDOPT  = -g

	MISCLIB0   = /sw/lib/libg95.a  -ldl -lpthread -lm \
		/System/Library/Frameworks/CoreServices.framework/Versions/A/CoreServices \
		/usr/lib/libSystem.B.dylib
	MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
	NML_OPTIONS = -g -Wall -ansi 

	# LAPACK-specific options
    FORTRAN = g95
    LAPACK_NOOPT =
    LAPACK_OPTS   = -g -funderscoring
    LAPACK_LDOPTS = -g -funderscoring
    ARPACK_OPTS   = -g -funderscoring

	# If -DNO_PARPACK is used in CCOPT, then we need to make sure that the PARPACK library doesn't get linked in
	LIB_PARPACK =

endif


# MACOSX GNU Debug
ifeq ($(BUILD_TYPE),macosx_intel_g95_fast)

	CCMKD  = g++ -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC     = g++
	CCOPT  = -c -O3 -Wall  -msse -msse2 -msse3 -D$(OSTYPE) -DFORTRAN_UNDERSCORE  -DFAKE_MPI -UMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) 
	LIB_PARPACK =
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT  = -c -O3  -msse -msse2 -msse3  -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE  -DFAKE_MPI -UMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	LIB_PARPACK =
	endif
	LD     = g++
	LDOPT  = -g

	MISCLIB0   = /sw/lib/libg95.a  -ldl -lpthread -lm \
		/System/Library/Frameworks/CoreServices.framework/Versions/A/CoreServices \
		/usr/lib/libSystem.B.dylib
	MISCINCL0  =

	#XERCES specific options

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT 
	NML_OPTIONS = -O3  -msse -msse2 -msse3  -Wall -ansi 

	# LAPACK-specific options
        FORTRAN = g95 
	LAPACK_NOOPT  = 
	LAPACK_OPTS   = -O2 -msse -msse2 -msse3  -funderscoring -funroll-all-loops 
	LAPACK_LDOPTS = -O2 -msse -msse2 -msse3  -funderscoring -funroll-all-loops
	ARPACK_OPTS   = -O2 -msse -msse2 -msse3  -funderscoring -funroll-all-loops

	# If -DNO_PARPACK is used in CCOPT, then we need to make sure that the PARPACK library doesn't get linked in
	LIB_PARPACK =

endif


# MACOSX GNU Debug
ifeq ($(BUILD_TYPE),macosx_intel_gfortran_fast)

	CCMKD  = g++ -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC     = g++
	CCOPT  = -c -O3 -Wall  -msse -msse2 -msse3 -D$(OSTYPE) -DFORTRAN_UNDERSCORE  -DFAKE_MPI -UMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) 
	LIB_PARPACK =
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT  = -c -O3  -msse -msse2 -msse3  -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE  -DFAKE_MPI -UMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	LIB_PARPACK =
	endif
	LD     = gfortran
	LDOPT  = -O3

	MISCLIB0   =  -lstdc++-static -ldl -lpthread -lm \
		/System/Library/Frameworks/CoreServices.framework/Versions/A/CoreServices \
		/usr/lib/libSystem.B.dylib
	MISCINCL0  =

	#XERCES specific options

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT 
	NML_OPTIONS = -O3  -msse -msse2 -msse3  -Wall -ansi 

	# LAPACK-specific options
        FORTRAN = gfortran 
	LAPACK_NOOPT  = 
	LAPACK_OPTS   = -O2 -msse -msse2 -msse3  -funderscoring -funroll-all-loops 
	LAPACK_LDOPTS = -O2 -msse -msse2 -msse3  -funderscoring -funroll-all-loops
	ARPACK_OPTS   = -O2 -msse -msse2 -msse3  -funderscoring -funroll-all-loops

	# If -DNO_PARPACK is used in CCOPT, then we need to make sure that the PARPACK library doesn't get linked in
	LIB_PARPACK =

endif


ifeq ($(BUILD_TYPE),macosx_intel_gfortran_local)

	CCMKD  = g++ -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC     = g++
	CCOPT  = -c -O3 -Wall  -msse -msse2 -msse3 -D$(OSTYPE) -UFORTRAN_UNDERSCORE  -DFAKE_MPI -UMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) 
	LIB_PARPACK =
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT  = -c -O3  -msse -msse2 -msse3  -Wall -D$(OSTYPE) -UFORTRAN_UNDERSCORE  -DFAKE_MPI -UMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	LIB_PARPACK =
	endif
	LD     = gfortran
	LDOPT  = -O3

	MISCLIB0   =  -lstdc++-static -ldl -lpthread -lm \
		/System/Library/Frameworks/CoreServices.framework/Versions/A/CoreServices \
		/usr/lib/libSystem.B.dylib
	MISCINCL0  =

	#XERCES specific options

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT 
	NML_OPTIONS = -O3  -msse -msse2 -msse3  -Wall -ansi 

	# LAPACK-specific options
        FORTRAN = gfortran 
	LAPACK_NOOPT  = 
	LAPACK_OPTS   = -O2 -msse -msse2 -msse3  -fno-underscoring -funroll-all-loops 
	LAPACK_LDOPTS = -O2 -msse -msse2 -msse3  -fno-underscoring -funroll-all-loops
	ARPACK_OPTS   = -O2 -msse -msse2 -msse3  -fno-underscoring -funroll-all-loops

	# If -DNO_PARPACK is used in CCOPT, then we need to make sure that the PARPACK library doesn't get linked in
	LIB_PARPACK =

endif



# MACOSX GNU Debug

# MACOSX GNU Debug
ifeq ($(BUILD_TYPE),macosx_intel_gfortran_mpi_fast)

	MPI_HOME = /opt/mpich2-1.0.5
	CCMKD  = g++ -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC     = $(MPI_HOME)/bin/mpicxx

	CCOPT  = -c -O3 -Wall -DMPICH_IGNORE_CXX_SEEK -msse -msse2 -msse3 -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d   -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) 
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT  = -c -O3  -DMPICH_IGNORE_CXX_SEEK -msse -msse2 -msse3  -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE  -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
#	LIB_PARPACK =
	endif
	LD     = $(MPI_HOME)/bin/mpif90
	LDOPT  = 

	MISCLIB0   =  $(MPI_HOME)/lib/libmpichcxx.a -lstdc++-static -ldl -lpthread -lm \
		/System/Library/Frameworks/CoreServices.framework/Versions/A/CoreServices \
		/usr/lib/libSystem.B.dylib
	MISCINCL0  =

	#XERCES specific options

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT 
	NML_OPTIONS = -O3  -msse -msse2 -msse3  -Wall -ansi 

	# LAPACK-specific options
        FORTRAN = $(MPI_HOME)/bin/mpif90
 
	LAPACK_NOOPT  = 
	LAPACK_OPTS   = -O2 -msse -msse2 -msse3  -funderscoring -funroll-all-loops 
	LAPACK_LDOPTS = -O2 -msse -msse2 -msse3  -funderscoring -funroll-all-loops

	ARPACK_OPTS   = -O2 -msse -msse2 -msse3  -funderscoring -funroll-all-loops

endif

ifeq ($(BUILD_TYPE),macosx_powerpc_gfortran_mpi_fast)
        
	MPI_HOME = /opt/mpich2-1.0.5
	CCMKD  = g++ -MM
        
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
        
	CC     = $(MPI_HOME)/bin/mpicxx
	CCOPT  = -c -O3 -Wall -DMPICH_IGNORE_CXX_SEEK -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d   -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF)
        
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT  = -c -O3 -Wall -DMPICH_IGNORE_CXX_SEEK -D$(OSTYPE) -DFORTRAN_UNDERSCORE  -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF)
	endif
        
	LD     = $(MPI_HOME)/bin/mpif90
	LDOPT  =
        
	MISCLIB0   =  $(MPI_HOME)/lib/libmpichcxx.a -lstdc++-static -ldl -lpthread -lm \
		/System/Library/Frameworks/CoreServices.framework/Versions/A/CoreServices \
		/usr/lib/libSystem.B.dylib
	MISCINCL0  =
        
	#XERCES specific options
        
	#XMLtoC specific options
        
	XMLtoC_CCOPT=$(CCOPT)
        
	# NML-specific options
        
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
	NML_OPTIONS = -O3 -Wall -ansi
        
	# LAPACK-specific options
        
	FORTRAN = $(MPI_HOME)/bin/mpif90
        
	LAPACK_NOOPT  =
	LAPACK_OPTS   = -O2 -funderscoring -funroll-all-loops
	LAPACK_LDOPTS = -O2 -funderscoring -funroll-all-loops
	ARPACK_OPTS   = -O2 -funderscoring -funroll-all-loops

endif







# MACOSX GNU Debug
ifeq ($(BUILD_TYPE),macosx_gnu_debug)

	CCMKD  = g++ -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC     = g++
	CCOPT  = -c -g -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE  -DFAKE_MPI -DMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) 
	LIB_PARPACK =
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT  = -c -g -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE  -DFAKE_MPI -DMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	LIB_PARPACK =
	endif
	LD     = g++
	LDOPT  = -g

	MISCLIB0   = -ldl -lpthread -lm \
		/System/Library/Frameworks/CoreServices.framework/Versions/A/CoreServices \
		/usr/lib/libSystem.B.dylib
	MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
	NML_OPTIONS = -g -Wall -ansi -pedantic

	# LAPACK-specific options
    FORTRAN = g77
    LAPACK_NOOPT =
    LAPACK_OPTS   = -g -funderscoring
    LAPACK_LDOPTS = -g -funderscoring
    ARPACK_OPTS   = -g -funderscoring

	# If -DNO_PARPACK is used in CCOPT, then we need to make sure that the PARPACK library doesn't get linked in
	LIB_PARPACK =

endif


# MACOSX GNU Debug
ifeq ($(BUILD_TYPE),macosx_gnu_debug_64bit)

	CCMKD  = g++ -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC     = g++
	CCOPT  = -c -g -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE  -DFAKE_MPI -DMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) -Dmem_64BIT
	LIB_PARPACK =
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT  = -c -g -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE  -DFAKE_MPI -DMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) -Dmem_64BIT
	LIB_PARPACK =
	endif
	LD     = g++
	LDOPT  = -g

	MISCLIB0   = -ldl -lpthread -lm \
		/System/Library/Frameworks/CoreServices.framework/Versions/A/CoreServices \
		/usr/lib/libSystem.B.dylib
	MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT -Dmem_64BIT
	NML_OPTIONS = -g -Wall -ansi -pedantic

	# LAPACK-specific options
    FORTRAN = g77
    LAPACK_NOOPT =
    LAPACK_OPTS   = -g -funderscoring
    LAPACK_LDOPTS = -g -funderscoring
    ARPACK_OPTS   = -g -funderscoring

	# If -DNO_PARPACK is used in CCOPT, then we need to make sure that the PARPACK library doesn't get linked in
	LIB_PARPACK =

endif

# MACOSX MPICH Absoft MPI Fast
ifeq ($(BUILD_TYPE),macosx_absoft_mpi_fast)

	ifeq ($(ABSOFT),)
		@ echo "Error - Undefined environment variable, ABSOFT"
	endif

	ifeq ($(MPI_HOME),)
		@ echo "Error - Undefined environment variable, MPI_HOME"
	endif

	CCMKD = g++ -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC = $(MPI_HOME)/bin/mpicxx
	CCOPT = -c -O3 -Wall -D$(OSTYPE) -DFORTRAN_ALLCAPS -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) 
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT = -c -O3 -Wall -D$(OSTYPE) -DFORTRAN_ALLCAPS -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	endif
	LD = $(MPI_HOME)/bin/mpicxx
	LDOPT  =

	MISCLIB0   = $(ABSOFT)/lib/libU77.a $(ABSOFT)/lib/libf77math.a $(ABSOFT)/lib/libf90math.a $(ABSOFT)/lib/libfio.a \
		/System/Library/Frameworks/CoreServices.framework/Versions/A/CoreServices \
		/usr/lib/libSystem.B.dylib \
		-ldl -lpthread -lm
	MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
	NML_OPTIONS = -O3 -Wall -ansi -pedantic

	FORTRAN = $(MPI_HOME)/bin/mpif77

	# LAPACK-specific options
	LAPACK_NOOPT =
    LAPACK_OPTS = -O2 -N11
    LAPACK_LDOPTS = -O2 -N11

    ARPACK_OPTS   = -O2 -N11

endif



# MACOSX MPICH Absoft Fast
ifeq ($(BUILD_TYPE),macosx_absoft_fast)

	ifeq ($(ABSOFT),)
		@ echo "Error - Undefined environment variable, ABSOFT"
	endif

	CCMKD = g++ -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC = g++
	CCOPT = -c -O3 -Wall -D$(OSTYPE) -DFORTRAN_ALLCAPS -DFAKE_MPI  -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF)  $(HAM_STORE_DEF) 
	LIB_PARPACK =
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT = -c -O3 -Wall -D$(OSTYPE) -DFORTRAN_ALLCAPS -DFAKE_MPI  -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	endif
	LD = g++
	LDOPT  =

	MISCLIB0   = $(ABSOFT)/lib/libU77.a $(ABSOFT)/lib/libf77math.a $(ABSOFT)/lib/libf90math.a $(ABSOFT)/lib/libfio.a \
		/System/Library/Frameworks/CoreServices.framework/Versions/A/CoreServices \
		/usr/lib/libSystem.B.dylib \
		-ldl -lpthread -lm
	MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
	NML_OPTIONS = -O3 -Wall -ansi -pedantic

	FORTRAN = $(ABSOFT)/bin/f90

	# LAPACK-specific options
	LAPACK_NOOPT  =
    LAPACK_OPTS   = -O2 -N11
    LAPACK_LDOPTS = -O2 -N11

    ARPACK_OPTS   = -O2 -N11

endif

# MACOSX MPICH Absoft Debug
ifeq ($(BUILD_TYPE),macosx_absoft_debug)

	ifeq ($(ABSOFT),)
		@ echo "Error - Undefined environment variable, ABSOFT"
	endif

	CCMKD = g++ -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC = g++
	CCOPT = -c -g -Wall -D$(OSTYPE) -DFORTRAN_ALLCAPS -DFAKE_MPI  -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) 
	LIB_PARPACK =
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT = -c -g -Wall -D$(OSTYPE) -DFORTRAN_ALLCAPS -DFAKE_MPI  -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	endif
	LD = g++
	LDOPT  =

	MISCLIB0   = $(ABSOFT)/lib/libU77.a $(ABSOFT)/lib/libf77math.a $(ABSOFT)/lib/libf90math.a $(ABSOFT)/lib/libfio.a  \
		/System/Library/Frameworks/CoreServices.framework/Versions/A/CoreServices \
		/usr/lib/libSystem.B.dylib \
		-ldl -lpthread -lm
	MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
	NML_OPTIONS = -g -Wall -ansi -pedantic

	FORTRAN = $(ABSOFT)/bin/f90

	# LAPACK-specific options
	LAPACK_NOOPT =
    LAPACK_OPTS =  -g -N11
    LAPACK_LDOPTS = -g -N11

    ARPACK_OPTS   = -g -N11

endif


# MACOSX Absoft MetroWerksDebug
ifeq ($(BUILD_TYPE),macosx_absoft_mwerks_debug)

	ifeq ($(ABSOFT),)
		@ echo "Error - Undefined environment variable, ABSOFT"
	endif

	CCMKD = g++ -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC = mwcc
	CCOPT = -c -g -Wall -lang cplusplus -D$(OSTYPE) -DFORTRAN_ALLCAPS -DFAKE_MPI  -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) 
	LIB_PARPACK =
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT = -c -g -Wall -lang cplusplus -D$(OSTYPE) -DFORTRAN_ALLCAPS -DFAKE_MPI  -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	endif
	LD = mwld
	LDOPT  =

	MISCLIB0   = $(ABSOFT)/lib/libU77.a $(ABSOFT)/lib/libf77math.a $(ABSOFT)/lib/libf90math.a $(ABSOFT)/lib/libfio.a \
		/System/Library/Frameworks/CoreServices.framework/Versions/A/CoreServices \
		/usr/lib/libSystem.B.dylib \
		-ldl -lpthread -lm
	MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
	NML_OPTIONS = -g -Wall -lang cplusplus -ext o

	FORTRAN = $(ABSOFT)/bin/f90

	# LAPACK-specific options
	LAPACK_NOOPT  =
    LAPACK_OPTS   = -g -N11
    LAPACK_LDOPTS = -g -N11

    ARPACK_OPTS   = -g -N11

endif

# MACOSX IBM Fast
ifeq ($(BUILD_TYPE),macosx_ibm_fast)
	CCMKD  = g++ -MM
	RANLIB = ranlib -s
	#LIB_PATHS = -L$(XERCESCROOT)/lib
	AR = ar

	FORTRAN = xlf_r
	CC = xlC_r
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d,-DFAKE_MPI, -DLIMIT_THEORY
	#CCOPT = -c -O3 -+ -qhot -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DFAKE_MPI -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) 
	CCOPT = -c -O3 -+ -qhot -D$(OSTYPE) -DFAKE_MPI -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) 

	ifeq ($(POSTPROCESS), Embedded)
		CC = g++
		#CCOPT = -c -O3 -+ -D$(OSTYPE) -DFORTRAN_UNDERSCORE $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
		CCOPT = -c -O3 -+ -D$(OSTYPE) $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
		LD = g++
		LDOPT =
	endif

	LD = xlC_r
	LDOPT  = -O3 -w 

##        MISCLIB0   = -ldl -lpthread -lmass -lg2c
	MISCLIB0   = -L/opt/ibmcmp/lib -lxlf90_r.A -lxlfmath.A -ldl -lpthread -framework CoreServices
	MISCINCL0  =

	#XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
	NML_OPTIONS = -+ -O3 -qhot


	# LAPACK-specific options
	LAPACK_NOOPT  = -O3 -mp
	LAPACK_OPTS = -O3 -qcache=auto -qhot -qipa -qtune=auto 
	LAPACK_LDOPTS = -O3 -mp 

	ARPACK_OPTS = -O3 -qcache=auto -qhot -qipa -qtune=auto
	LIB_PARPACK = 
endif

# MACOSX IBM Debug
ifeq ($(BUILD_TYPE),macosx_ibm_debug)
	CCMKD  = g++ -MM
	RANLIB = ranlib -s
	#LIB_PATHS = -L$(XERCESCROOT)/lib
	AR = ar

	FORTRAN = xlf_r
	CC = xlC_r
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d,-DFAKE_MPI, -DLIMIT_THEORY
	#CCOPT = -c -O3 -+ -qhot -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DFAKE_MPI -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) 
	CCOPT = -c -g -+ -D$(OSTYPE) -DFAKE_MPI -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) 

	ifeq ($(POSTPROCESS), Embedded)
		CC = g++
		#CCOPT = -c -O3 -+ -D$(OSTYPE) -DFORTRAN_UNDERSCORE $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
		CCOPT = -c -g -+ -D$(OSTYPE) $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
		LD = g++
		LDOPT =
	endif

	LD = xlC_r
	LDOPT  = -g -w 

##        MISCLIB0   = -ldl -lpthread -lmass -lg2c
	MISCLIB0   = -L/opt/ibmcmp/lib -lxlf90_r.A -lxlfmath.A -ldl -lpthread -framework CoreServices
	MISCINCL0  =

	#XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
	NML_OPTIONS = -+ -g


	# LAPACK-specific options
	LAPACK_NOOPT  = -g
	LAPACK_OPTS = -g -qcache=auto -qtune=auto 
	LAPACK_LDOPTS = -g -mp 

	ARPACK_OPTS = -g -qcache=auto -qtune=auto
	LIB_PARPACK = 
endif

# AIX MPI
ifeq ($(BUILD_TYPE),aix_mpi_static)
	LIB_PATHS = -L$(XERCESCROOT)/lib
	HOSTTYPE=AIX
	MPI_HOME=/usr/lpp/ppe.poe
	MAKE=gmake
	AR = ar
        AR_OPTIONS = X 32 r
	CCMKD = g++ -MM

	FORTRAN = mpxlf_r
	CC = mpCC_r
	LD = mpCC_r

	ifeq ($(POSTPROCESS), Embedded)
		CC = g++
		CCOPT = -c -O3 -+ -D$(OSTYPE) -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	LD = g++
	LDOPT =
	endif

	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d,-DFAKE_MPI, -DLIMIT_THEORY
	CCOPT = -c -O3 -+ -qhot -D$(OSTYPE) -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DAIX_SPECIFIC $(HAM_STORE_DEF) 
	LDOPT  = -O3 -brtl -w -lxlf90 -bmaxdata:0x80000000

##        MISCLIB0   = -ldl -lpthread -lmass -lg2c
	MISCLIB0   = -ldl -lpthread
	MISCINCL0  =

	#XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
	NML_OPTIONS = -+ -O3 -qhot


	# LAPACK-specific options
	LAPACK_NOOPT  = -O3 -mp
	LAPACK_OPTS = -O3 -bmaxdata:0x80000000 -qarch=auto -qcache=auto -qhot -qipa -qtune=auto
	LAPACK_LDOPTS = -O3 -mp -bmaxdata:0x80000000

	ARPACK_OPTS = -O3 -qarch=auto -qcache=auto -qhot -qipa -qtune=auto
endif

# AIX MPI 64-bit
# Must set environment variable OBJECT_MODE=64 for build
ifeq ($(BUILD_TYPE),aix_mpi_static64)
	LIB_PATHS = -L$(XERCESCROOT)/lib
	HOSTTYPE=AIX
	MPI_HOME=/usr/lpp/ppe.poe
	MAKE=gmake
	AR = ar
        AR_OPTIONS = X 64 r
	CCMKD = g++ -MM

	FORTRAN = mpxlf_r
	CC = mpCC_r
	LD = mpCC_r

	ifeq ($(POSTPROCESS), Embedded)
		CC = g++
		CCOPT = -c -O3 -+ -D$(OSTYPE) -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	LD = g++
	LDOPT =
	endif

	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d,-DFAKE_MPI, -DLIMIT_THEORY
	CCOPT = -c -O3 -+ -qhot -D$(OSTYPE) -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DAIX_SPECIFIC $(HAM_STORE_DEF) 
	LDOPT  = -O3 -brtl -w -lxlf90

##        MISCLIB0   = -ldl -lpthread -lmass -lg2c
	MISCLIB0   = -ldl -lpthread
	MISCINCL0  =

	#XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
	NML_OPTIONS = -+ -O3 -qhot


	# LAPACK-specific options
	LAPACK_NOOPT  = -O3 -mp
	LAPACK_OPTS = -O3 -qarch=auto -qcache=auto -qhot -qipa -qtune=auto
	LAPACK_LDOPTS = -O3 -mp 

	ARPACK_OPTS = -O3 -qarch=auto -qcache=auto -qhot -qipa -qtune=auto
endif

# MPICH GNU MPI Fast
ifeq ($(BUILD_TYPE),pathscale_mpi_fast)

	ifeq ($(MPI_HOME),)
		@echo "Error - Undefined environment variable, MPI_HOME"
		exit
	endif

	CCMKD = g++ -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
        CC = $(MPI_HOME)/bin/mpiCC
	CCOPT = -c -O3 -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF)  -m64 -Dmem_64BIT
	ifeq ($(POSTPROCESS), Embedded)
		JUNK = $(CCOPT)
		CCOPT = $(JUNK) -DEmbedded_Postprocess 
	endif
        LD = $(MPI_HOME)/bin/mpiCC
	LDOPT  = 

	MISCLIB0   = -ldl -lpthread -lm -lg2c
	MISCINCL0  =


        FORTRAN = $(MPI_HOME)/bin/mpif77

	#XERCES specific options
	XERCES_CC =pathCC
	XERCES_CCMKD =
	XERCES_CCOPT=
	XERCES_LD=pathCC
	XERCES_LDOPT=

        #XMLtoC specific options
	XMLtoC_CC= pathCC
	XMLtoC_CCOPT=
	XMLtoC_CCMKD=
	XMLtoC_LD=pathCC
	XMLtoC_LDOPT=

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT -Dmem_64BIT
	NML_OPTIONS = -O3 -m64

	FORTRAN = $(MPI_HOME)/bin/mpif77

	# LAPACK-specific options
	LAPACK_NOOPT  = -O1 -m64
	LAPACK_OPTS   = -O3  -m64
	LAPACK_LDOPTS = -O3 -m64

	# ARPACK-specific options
	ARPACK_OPTS   = -O3  -m64

endif

# GNU Debug Static
ifeq ($(BUILD_TYPE),cygwin_gnu_debug)

        CCMKD  = g++ -MM
        # Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
        CC     = g++
        CCOPT  = -c -g -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE  -DFAKE_MPI -UMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) 
	LIB_PARPACK =
        ifeq ($(POSTPROCESS), Embedded)
                CCOPT  = -c -g -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE  -DFAKE_MPI -UMPI_TIMING -DNO_PARPACK $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
        endif
        LD     = g++
        LDOPT  = -g -Wl
        LDOPT  = -g 

        MISCLIB0   = -ldl -lpthread -lm -lg2c
        MISCLIB0   = -lpthread -lm -lg2c
        MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

        # NML-specific options
        NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
        NML_OPTIONS = -g -Wall -ansi 

        # LAPACK-specific options
        ifeq ($(OS),Darwin)
                FORTRAN = f77
                LAPACK_NOOPT =
                LAPACK_OPTS   = -g
                LAPACK_LDOPTS = -g
                ARPACK_OPTS   = -g
        else
                FORTRAN = g77
                LAPACK_NOOPT =
                LAPACK_OPTS   = -g -funderscoring
                LAPACK_LDOPTS = -g -funderscoring
                ARPACK_OPTS   = -g -funderscoring
        endif

        # If -DNO_PARPACK is used in CCOPT, then we need to make sure that the PARPACK library doesn't get linked in
        LIB_PARPACK =

endif


ifeq ($(BUILD_TYPE),ibm64_mpich_instr)
  # Compiler location
  MPI_HOME = $(MPICH_HOME)
  # Compiler/Archiver/Linker definition
  AR     = ar
  # AR_OPTIONS = r
  CC     = mpicxx -+
  CCMKD  = g++ -MM
  CCFAST = -O3 -qhot
  CCOPT  = -c -D$(OSTYPE) -DMPI3d -DMPI_TIMING $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) $(CCFAST) -DMPI_INSTR -I../../nemo3d/src/base -I../../nemo3d/src/math
  ifeq ($(POSTPROCESS), Embedded)
    CCOPT += -DEmbedded_Postprocess
  endif
  LD     = $(CC)
  LDOPT  = $(CCFAST) -w -lxlf90
  LDOPT  = $(CCFAST) -w -L/opt/ibmcmp/xlf/10.1/lib64 -lxlf90_r -lxlfmath
  # Libraries
  MISCLIB0  = -lpthread -ldl -lm
  MISCINCL0 =
  # XERCES specific options
  XERCES_CC    = g++ -m64
  XERCES_CCMKD = $(CCMKD)
  XERCES_CCOPT = $(CCOPT)
  XERCES_LD    = $(XERCES_CC)
  XERCES_LDOPT = $(LDOPT)
  #XMLtoC specific options
  XMLtoC_CC    = $(CC)
  XMLtoC_CCOPT = $(CCOPT)
  XMLtoC_CCMKD = $(CCMKD)
  XMLtoC_LD    = $(XMLtoC_CC)
  XMLtoC_LDOPT = $(LDOPT)
  # NML-specific options
  NML_DEFINES  = -DNML_ERRANT -DNML_MEM_EXT
  NML_OPTIONS  = $(CCOPT)
  # Options for Fortran codes
	# etime with trailing underscore
  FORTRAN       = mpif77 -I$(MPI_HOME)/include -qextname=etime
  LAPACK_NOOPT  = $(CCFAST) -qstrict
  LAPACK_OPTS   = $(CCFAST) -qstrict -qarch=ppc970 -qtune=ppc970 -qipa
  LAPACK_LDOPTS = $(CCFAST) -qstrict
  ARPACK_OPTS   = $(CCFAST) -qstrict -qarch=ppc970 -qtune=ppc970 -qipa
endif

# NANOHUB MPICH Intel 10
ifeq ($(BUILD_TYPE),nanohub_intel10_mpi)

	CCDIR=/apps/intel/cc/10.0.023
	FCDIR=/apps/intel/fc/10.0.023
	AR = xiar

	ifeq ($(MPI_HOME),)
		@ echo "Error - Undefined environment variable, MPI_HOME"
	endif

	CCMKD = $(MPI_HOME)/bin/mpiCC -MM -I/usr/include/c++/4.1 -I/usr/include/c++/4.1/i486-linux-gnu
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
	CC = $(MPI_HOME)/bin/mpiCC -I/usr/include/c++/4.1 -I/usr/include/c++/4.1/i486-linux-gnu -I${CCDIR}/include

	CFLAGS = -O1 -Kc++ -axP
	CCOPT = -c $(CFLAGS) -UINTEL_MODIFIED -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(MKL_MACRO_DEF)  $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -I/usr/include/c++/4.1 -I/usr/include/c++/4.1/i486-linux-gnu -I${CCDIR}/include $(HAM_STORE_DEF)  -Umem_64BIT $(MKL_INCLUDE)

	ifeq ($(POSTPROCESS), Embedded)
		CCOPT = -c -O1 -Kc++ -axP -DINTEL_MODIFIED -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(MKL_MACRO_DEF)  $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess -I${CCDIR}/include/c++ $(HAM_STORE_DEF)  -Umem_64BIT $(MKL_INCLUDE)
	endif

	LD = $(CC)
	LDOPT  = $(CFLAGS) -Wl,-static -Wall
	LDOPT  = $(CFLAGS)  -Wall


	#MISCLIB0   = -lpthread -lm ${FCDIR}/lib/libifcore.a ${FCDIR}/lib/libifport.a -ldl -lsvml $(MKL_LINK)
	#MISCLIB0   = -lpthread -lm ${CCDIR}/lib/libcprts.a ${CCDIR}/lib/libcxa.a ${CCDIR}/lib/libunwind.a ${FCDIR}/lib/libifcore.a ${FCDIR}/lib/libifport.a -ldl  ${CCDIR}/lib/libsvml.a $(MPI_HOME)/lib/libmpichf90.a $(MKL_LINK)
        MISCLIB0   = -lpthread -lm ${CCDIR}/lib/libcxaguard.a ${FCDIR}/lib/libifcore.a ${FCDIR}/lib/libifport.a -ldl ${CCDIR}/lib/libsvml.a $(MPI_HOME)/lib/libmpichf90.a $(MKL_LINK)
	MISCINCL0  =

	#XERCES specific options
	XERCES_CCOPT= $(CCOPT)
	XERCES_CC= icc

  	#XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)
	XMLtoC_CC = icc

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT  -Dmem_64BIT
	NML_OPTIONS = $(CCOPT)

	FORTRAN   = $(MPI_HOME)/bin/mpif77

	# LAPACK-specific options
	LAPACK_NOOPT  = -O1 -mp -fp_port -mtune=pentium4 -xW -axP
	LAPACK_OPTS   = -O1 -mp -fp_port -mtune=pentium4 -xW -axP
	LAPACK_LDOPTS = -O1 -mp -fp_port -mtune=pentium4 -xW -axP -Vaxlib

	# ARPACK-specific options
	ARPACK_OPTS   = -O1 -mp -mtune=pentium4 -axP -xW

	# TRACEMIN-specific options (these are set here only for choice of BUILD_TYPE: intel8,  intel9_64_mpi_intelfast (lear.rcac.purdue.edu))  
	SCTRACEMIN_OPTS = -O1
 	PCTRACEMIN_OPTS = -O1
	SQTRACEMIN_OPTS = -O1
	PQTRACEMIN_OPTS = -O1
	FORTRAN90 = $(MPI_HOME)/bin/mpif90  #tracemin fortran compiler
	##############################################################

endif

# MPICH GNU MPI Fast (nanohub workspace)
ifeq ($(BUILD_TYPE),nanohub_gcc_mpi)

	ifeq ($(MPI_HOME),)
		@echo "Error - Undefined environment variable, MPI_HOME"
		exit
	endif

	CCMKD = g++ -MM
	# Allowed CCOPT flags: -DFORTRAN_UNDERSCORE, -DFORTRAN_ALLCAPS, -DMPI3d, -DFAKE_MPI, -DLIMIT_THEORY
        CC = $(MPI_HOME)/bin/mpiCC
	CCOPT = -c -O3 -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) $(HAM_STORE_DEF) 
	ifeq ($(POSTPROCESS), Embedded)
		CCOPT = -c -O3 -Wall -D$(OSTYPE) -DFORTRAN_UNDERSCORE -DMPI3d -DMPI_TIMING $(TRACEMIN_MACRO_DEF) $(PYTHON_MACRO_DEF) -DEmbedded_Postprocess $(HAM_STORE_DEF) 
	endif
        LD = $(MPI_HOME)/bin/mpiCC
	LDOPT  = 

	MISCLIB0   = -ldl -lpthread -lm -lg2c
	MISCINCL0  =

        #XMLtoC specific options
	XMLtoC_CCOPT=$(CCOPT)

	# NML-specific options
	NML_DEFINES = -DNML_ERRANT  -DNML_MEM_EXT
	NML_OPTIONS = -O3 -Wall -ansi -pedantic
	NML_OPTIONS = -O3 -Wall -ansi 

        FORTRAN = $(MPI_HOME)/bin/mpif77

	# LAPACK-specific options
	LAPACK_NOOPT =
	ifeq ($(OS),Darwin)
		LAPACK_OPTS   = -O2 -N11
		LAPACK_LDOPTS = -O2 -N11
	else
		LAPACK_OPTS   = -O2 -funderscoring -funroll-all-loops
		LAPACK_LDOPTS = -O2 -funderscoring -funroll-all-loops
	endif

        ifeq ($(OS),Darwin)
		ARPACK_OPTS   = -O2 -N11
	else
		ARPACK_OPTS   = -O2 -funderscoring -funroll-all-loops
	endif
endif

# Misc additional library and include flags

ifeq ($(OS),Linux)
	MISCLIB  = $(MISCLIB0) -lutil
ifeq ($(BUILD_TYPE),crayxt_backend_pgi_fast)
	MISCLIB  = $(MISCLIB0)
endif
ifeq ($(BUILD_TYPE),crayxt_backend_pgi_debug)
	MISCLIB  = $(MISCLIB0)
endif
ifeq ($(BUILD_TYPE),crayxt_backend_pgi_fast_profile)
	MISCLIB  = $(MISCLIB0)
endif
else
	MISCLIB  = $(MISCLIB0)
endif

ifeq ($(OS),Darwin)
	MISCINCL = $(MISCINCL0) -I/usr/include/objc
else
	MISCINCL = $(MISCINCL0)
endif

ifeq ($(XERCES_CC),)
	XERCES_CC = $(CC)
endif
ifeq ($(XERCES_CCMKD),)
	XERCES_CCMKD = $(CCMKD)
endif
ifeq ($(XERCES_CCOPT),)
	XERCES_CCOPT = $(CCOPT)
endif
ifeq ($(XERCES_LD),)
	XERCES_LD = $(LD)
endif
ifeq ($(XERCES_LDOPT),)
	XERCES_LDOPT = $(LDOPT)
endif

ifeq ($(XMLtoC_CC),)
	XMLtoC_CC = $(CC)
endif
ifeq ($(XMLtoC_CCOPT),)
	XMLtoC_CCOPT=$(CCOPT)
endif
ifeq ($(XMLtoC_CCMKD),)
	XMLtoC_CCMKD = $(CCMKD)
endif
ifeq ($(XMLtoC_LD),)
	XMLtoC_LD = $(LD)
endif
ifeq ($(XMLtoC_LDOPT),)
	XMLtoC_LDOPT = $(LDOPT)
endif





# ---------------------------------------------------------------
# OVERRIDE default profiles for different systems
# ( You can safely delete this section if you only build on one system )
#----------------------------------------------------------------

ifeq ($(shell hostname),pluto.jpl.nasa.gov)
        ifeq ($(PYTHON_MACRO_DEF),-DUSE_PYTHON)
	PYTHON_ROOT = /usr/lib/python1.5
	PYTHON_EXPORT_INCL = -I/usr/include/python1.5
	LIB_PYTHON  = $(PYTHON_ROOT)/config/libpython1.5.a
	endif
	ifeq ($(BUILD_TYPE),absoft_mpi_fast)
		CC      = /home/packages/mpich-1.2.5/absoft/bin/mpiCC
		LD      = /home/packages/mpich-1.2.5/absoft/bin/mpiCC
		FORTRAN = /home/packages/mpich-1.2.4/absoft/bin/mpif77
	endif
	ifeq ($(BUILD_TYPE),absoft_mpi_fast_static)
		CC      = /home/packages/mpich-1.2.5/absoft/bin/mpiCC
		LD      = /home/packages/mpich-1.2.5/absoft/bin/mpiCC
		FORTRAN = /home/packages/mpich-1.2.4/absoft/bin/mpif77
	endif
	ifeq ($(BUILD_TYPE),absoft_mpi_debug)
		CC      = /home/packages/mpich-1.2.5/absoft/bin/mpiCC
		LD      = /home/packages/mpich-1.2.5/absoft/bin/mpiCC
		FORTRAN = /home/packages/mpich-1.2.4/absoft/bin/mpif77
	endif
	ifeq ($(BUILD_TYPE),intel_piii_mpi_static)
		CC      = /home/packages/mpich-1.2.5/intel/bin/mpicxx
		LD      = /home/packages/mpich-1.2.5/intel/bin/mpicxx
		FORTRAN = /home/packages/mpich-1.2.5/intel/bin/mpif77
	endif
endif

ifeq ($(shell hostname),polo)
	MAKE = /usr/freeware/bin/gmake
	RANLIB = echo
ifeq ($(PYTHON_MACRO_DEF),-DUSE_PYTHON)
	PYTHON_ROOT = /usr/freeware/lib/python2.1
	LIB_PYTHON  = $(PYTHON_ROOT)/config/libpython2.1.a
	PYTHON_EXPORT_INCL = -I/usr/freeware/include/python2.1
endif
	LIB_PATHS = -L$(XERCESCROOT)/lib -Wl,-rpath,$(XERCESCROOT)/lib \
		-L/usr/freeware/lib/gcc-lib/mips-sgi-irix6.5/3.0.4
endif

ifeq ($(shell hostname),newyork)
ifeq ($(PYTHON_MACRO_DEF),-DUSE_PYTHON)
	PYTHON_ROOT = /usr/local/lib/python2.2-Zach
	LIB_PYTHON  = $(PYTHON_ROOT)/config/libpython2.2.a
	PYTHON_EXPORT_INCL = -I/usr/include/python2.2
endif
	LIB_PATHS = -L$(XERCESCROOT)/lib
	MISCLIB  = $(MISCLIB0) \
		/Users/hook/temp/gcc3-1151/obj/ppc/libf2c/.libs/libg2c.a 
	ifeq ($(BUILD_TYPE),gnu_debug_static)
		FORTRAN       = g77
		LAPACK_OPTS   = -g -funderscoring
		LAPACK_LDOPTS = -g -funderscoring
		ARPACK_OPTS   = -g -funderscoring
	endif
	ifeq ($(BUILD_TYPE),gnu_fast)
		FORTRAN = g77
		LAPACK_OPTS   = -g -O2 -funderscoring -funroll-all-loops
		LAPACK_LDOPTS = -g -O2 -funderscoring -funroll-all-loops
		ARPACK_OPTS   = -g -O2 -funderscoring -funroll-all-loops
	endif
	ifeq ($(BUILD_TYPE),absoft_mpi_fast)
		CC      = /usr/local/bin/MPI/mpich-1.2.5/absoft/bin/mpicxx
		LD      = /usr/local/bin/MPI/mpich-1.2.5/absoft/bin/mpicxx
		FORTRAN = /usr/local/bin/MPI/mpich-1.2.5/absoft/bin/mpif77
	endif
	ifeq ($(BUILD_TYPE),absoft_mpi_fast_static)
		CC      = /usr/local/bin/MPI/mpich-1.2.5/absoft/bin/mpicxx
		LD      = /usr/local/bin/MPI/mpich-1.2.5/absoft/bin/mpicxx
		FORTRAN = /usr/local/bin/MPI/mpich-1.2.5/absoft/bin/mpif77
	endif
	ifeq ($(BUILD_TYPE),absoft_mpi_debug)
		CC      = /usr/local/bin/MPI/mpich-1.2.5/absoft/bin/mpicxx
		LD      = /usr/local/bin/MPI/mpich-1.2.5/absoft/bin/mpicxx
		FORTRAN = /usr/local/bin/MPI/mpich-1.2.5/absoft/bin/mpif77
	endif
endif

ifeq ($(shell hostname),nimrod.jpl.nasa.gov)
ifeq ($(PYTHON_MACRO_DEF),-DUSE_PYTHON)
	PYTHON_ROOT        = /home/packages/Python-2.2.2
	LIB_PYTHON         = $(PYTHON_ROOT)/lib/python2.2/config/libpython2.2.a
	PYTHON_EXPORT_INCL = -I$(PYTHON_ROOT)/include/python2.2
endif
	LAPACK_OPTS =  -g -cpu:host
	ARPACK_OPTS =  -g -cpu:host
endif

ifeq ($(shell hostname),akos.jpl.nasa.gov)
ifeq ($(PYTHON_MACRO_DEF),-DUSE_PYTHON)
	PYTHON_ROOT = /sw/lib/python2.2
	LIB_PYTHON  = $(PYTHON_ROOT)/config/libpython2.2.a
	PYTHON_EXPORT_INCL = -I/sw/include/python2.2
endif
	LIB_PATHS = \
		-L$(XERCESCROOT)/lib \
		-L/sw/lib
	MISCLIB0   = -ldl -lpthread -lm
endif


ifeq ($(shell hostname),superman.rcs.purdue.edu)
ifeq ($(PYTHON_MACRO_DEF),-DUSE_PYTHON)
	PYTHON_ROOT        = /usr/lib/python1.5
	LIB_PYTHON         = $(PYTHON_ROOT)/config/libpython1.5.a
	PYTHON_EXPORT_INCL = -I$(PYTHON_ROOT)/include/python1.5
	PYTHON_EXPORT_INCL = -I/usr/include/python1.5
endif
#	LAPACK_OPTS =  -cpu:host
#	ARPACK_OPTS =  -cpu:host
endif

ifeq ($(shell hostname),ecn200-dhcp-7.ecn.purdue.edu)
ifeq ($(PYTHON_MACRO_DEF),-DUSE_PYTHON)
        PYTHON_ROOT = /sw/lib/python2.2
        LIB_PYTHON  = $(PYTHON_ROOT)/config/libpython2.2.a
        LIB_PYTHON  = /Volumes/Data_gekco/DEVELOP/pythonDB/lib/libpython2.2_OSX.a
        PYTHON_EXPORT_INCL = -I/sw/include/python2.2
endif
        LIB_PATHS = \
                -L$(XERCESCROOT)/lib \
                -L/sw/lib
        MISCLIB0   = -ldl -lpthread -lm
endif

ifeq ($(shell hostname),Gerhard-Klimecks-Computer.local)
ifeq ($(PYTHON_MACRO_DEF),-DUSE_PYTHON)
        PYTHON_ROOT = /sw/lib/python2.2
        LIB_PYTHON  = $(PYTHON_ROOT)/config/libpython2.2.a
        PYTHON_EXPORT_INCL = -I/sw/include/python2.2
endif
        LIB_PATHS = \
                -L$(XERCESCROOT)/lib \
                -L/sw/lib
        MISCLIB  = $(MISCLIB0) \
                /System/Library/Frameworks/CoreServices.framework/Versions/A/CoreServices
endif

ifeq ($(shell hostname),belle.ecn.purdue.edu)
ifeq ($(PYTHON_MACRO_DEF),-DUSE_PYTHON)
        PYTHON_ROOT = /sw/lib/python2.2
        LIB_PYTHON  = $(PYTHON_ROOT)/config/libpython2.2.a
        PYTHON_EXPORT_INCL = -I/sw/include/python2.2
endif
        LIB_PATHS = \
                -L$(XERCESCROOT)/lib \
                -L/sw/lib
        MISCLIB  = $(MISCLIB0) \
                /System/Library/Frameworks/CoreServices.framework/Versions/A/CoreServices
endif


ifeq ($(shell hostname),bubba.ecn.purdue.edu)
ifeq ($(PYTHON_MACRO_DEF),-DUSE_PYTHON)
        PYTHON_ROOT = /Volumes/OS9_Apps/sw/lib/python2.2
        LIB_PYTHON  = $(PYTHON_ROOT)/config/libpython2.2.a
        PYTHON_EXPORT_INCL = -I/Volumes/OS9_Apps/sw/include/python2.2
endif
        LIB_PATHS = \
                -L$(XERCESCROOT)/lib \
                -L/sw/lib
        MISCLIB  = $(MISCLIB0) \
                /System/Library/Frameworks/CoreServices.framework/Versions/A/CoreServices
endif

ifeq ($(HOSTNAME),tg-login2.ncsa.teragrid.org)
     CCMKD = $(MPI_HOME)/bin/mpiCC -MM -I/opt/intel/compiler80/include/c++
     MISCLIB0 = -lpthread /usr/local/intel/compiler80/lib/libifcore.a \
     /usr/local/intel/compiler80/lib/libifcoremt.a \
     /usr/local/intel/compiler80/lib/libifport.a -ldl -lm
endif

# Misc additional library and include flags

ifeq ($(OS),Linux)
	MISCLIB  = $(MISCLIB0) -lutil
ifeq ($(BUILD_TYPE),crayxt_backend_pgi_fast)
	MISCLIB  = $(MISCLIB0)
endif
ifeq ($(BUILD_TYPE),crayxt_backend_pgi_debug)
	MISCLIB  = $(MISCLIB0)
endif
ifeq ($(BUILD_TYPE),crayxt_backend_pgi_fast_profile)
	MISCLIB  = $(MISCLIB0)
endif
else
	MISCLIB  = $(MISCLIB0)
endif

ifeq ($(OS),Darwin)
	MISCINCL = $(MISCINCL0) -I/usr/include/objc
else
	MISCINCL = $(MISCINCL0)
endif

ifeq ($(XERCES_CC),)
	XERCES_CC = $(CC)
endif
ifeq ($(XERCES_CCMKD),)
	XERCES_CCMKD = $(CCMKD)
endif
ifeq ($(XERCES_CCOPT),)
	XERCES_CCOPT = $(CCOPT)
endif
ifeq ($(XERCES_LD),)
	XERCES_LD = $(LD)
endif
ifeq ($(XERCES_LDOPT),)
	XERCES_LDOPT = $(LDOPT)
endif

ifeq ($(XMLtoC_CC),)
	XMLtoC_CC = $(CC)
endif
ifeq ($(XMLtoC_CCOPT),)
	XMLtoC_CCOPT=$(CCOPT)
endif
ifeq ($(XMLtoC_CCMKD),)
	XMLtoC_CCMKD = $(CCMKD)
endif
ifeq ($(XMLtoC_LD),)
	XMLtoC_LD = $(LD)
endif
ifeq ($(XMLtoC_LDOPT),)
	XMLtoC_LDOPT = $(LDOPT)
endif
ifeq ($(FORTRAN90),)
	FORTRAN90 = $(FORTRAN)
endif
ifeq ($(SCTRACEMIN_OPTS),)
	SCTRACEMIN_OPTS = $(LAPACK_OPTS)
endif
ifeq ($(PCTRACEMIN_OPTS),)
	PCTRACEMIN_OPTS = $(LAPACK_OPTS)
endif
ifeq ($(SQTRACEMIN_OPTS),)
	SQTRACEMIN_OPTS = $(LAPACK_OPTS)
endif
ifeq ($(PQTRACEMIN_OPTS),)
	PQTRACEMIN_OPTS = $(LAPACK_OPTS)
endif

