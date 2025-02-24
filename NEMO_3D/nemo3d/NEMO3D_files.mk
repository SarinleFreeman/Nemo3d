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
# $Header: /repo/nemo3d/NEMO3D_files.mk,v 1.35 2007/03/08 06:49:38 lee509 Exp $
# *****************************************************************************

#================================================================
# NEMO-3D Specification of files
#================================================================


#----------------------------------------------------------------
# Directory structure
#---------------------------------------------------------------

NEMO3D_SRCD   = $(NEMO3D_ROOT)/src
NEMO3D_OBJD   = $(NEMO3D_ROOT)/obj/$(PLATFORM)
NEMO3D_LIBD   = $(NEMO3D_ROOT)/lib/$(PLATFORM)
NEMO3D_BIND   = $(NEMO3D_ROOT)/bin


#----------------------------------------------------------------
#  Include and library information for compiler
#----------------------------------------------------------------


NEMO3D_INCLD_TOP  = $(NEMO3D_ROOT)/src/top
NEMO3D_INCLD_BASE = $(NEMO3D_ROOT)/src/base
NEMO3D_INCLD_IO   = $(NEMO3D_ROOT)/src/io
NEMO3D_INCLD_MATH = $(NEMO3D_ROOT)/src/math
NEMO3D_INCLD_UTIL = $(NEMO3D_ROOT)/src/util
NEMO3D_INCLD_HAM  = $(NEMO3D_ROOT)/src/ham
NEMO3D_INCLD_POST = $(NEMO3D_ROOT)/src/postprocessing

ifeq ($(POSTPROCESS), Embedded)
	NEMO3D_SRCH_INCLD = \
	-I$(NEMO3D_INCLD_TOP) \
	-I$(NEMO3D_INCLD_BASE) \
	-I$(NEMO3D_INCLD_IO) \
	-I$(NEMO3D_INCLD_MATH) \
	-I$(NEMO3D_INCLD_UTIL) \
	-I$(NEMO3D_INCLD_HAM) \
	-I$(NEMO3D_INCLD_POST) \
	-I$(NML_ROOT)/include \
	-I$(EIGEN_ROOT)/include \
	-I$(XMLtoC_ROOT)/apps/ioNemo3D/src \
	-I$(XMLtoC_ROOT)/src \
	$(XERCESC_EXPORT_INCL) \
	$(PYTHON_EXPORT_INCL) \
	$(MISCINCL)
else
	NEMO3D_SRCH_INCLD = \
	-I$(NEMO3D_INCLD_TOP) \
	-I$(NEMO3D_INCLD_BASE) \
	-I$(NEMO3D_INCLD_IO) \
	-I$(NEMO3D_INCLD_MATH) \
	-I$(NEMO3D_INCLD_UTIL) \
	-I$(NEMO3D_INCLD_HAM) \
	-I$(NML_ROOT)/include \
	-I$(EIGEN_ROOT)/include \
	-I$(XMLtoC_ROOT)/apps/ioNemo3D/src \
	-I$(XMLtoC_ROOT)/src \
	$(XERCESC_EXPORT_INCL) \
	$(PYTHON_EXPORT_INCL) \
	$(MISCINCL)
endif

#----------------------------------------------------------------
#  Library information
#----------------------------------------------------------------

NEMO3D_LIB_TOP  = $(NEMO3D_LIBD)/libnemo3d-top-$(PLATFORM).a
NEMO3D_LIB_BASE = $(NEMO3D_LIBD)/libnemo3d-base-$(PLATFORM).a
NEMO3D_LIB_IO   = $(NEMO3D_LIBD)/libnemo3d-io-$(PLATFORM).a
NEMO3D_LIB_MATH = $(NEMO3D_LIBD)/libnemo3d-math-$(PLATFORM).a
NEMO3D_LIB_UTIL = $(NEMO3D_LIBD)/libnemo3d-util-$(PLATFORM).a
NEMO3D_LIB_HAM  = $(NEMO3D_LIBD)/libnemo3d-ham-$(PLATFORM).a
ifeq ($(POSTPROCESS), Embedded)
	NEMO3D_LIB_POST  = $(NEMO3D_LIBD)/libnemo3d-post-$(PLATFORM).a
else
	NEMO3D_LIB_POST  =
endif

NEMO3D_LIBS = \
	$(NEMO3D_LIB_POST) \
	$(NEMO3D_LIB_IO) \
	$(NEMO3D_LIB_BASE) \
	$(NEMO3D_LIB_HAM) \
	$(NEMO3D_LIB_IO) \
	$(NEMO3D_LIB_MATH) \
	$(NEMO3D_LIB_UTIL)

NEMO3D_LIB_TOP_LNK  = $(NEMO3D_LIBD)/libnemo3d-top.a
NEMO3D_LIB_BASE_LNK = $(NEMO3D_LIBD)/libnemo3d-base.a
NEMO3D_LIB_IO_LNK   = $(NEMO3D_LIBD)/libnemo3d-io.a
NEMO3D_LIB_MATH_LNK = $(NEMO3D_LIBD)/libnemo3d-math.a
NEMO3D_LIB_UTIL_LNK = $(NEMO3D_LIBD)/libnemo3d-util.a
NEMO3D_LIB_HAM_LNK  = $(NEMO3D_LIBD)/libnemo3d-ham.a
NEMO3D_LIB_POST_LNK  = $(NEMO3D_LIBD)/libnemo3d-post.a


NEMO3D_EXTERNLIBS = \
	$(LIB_XMLtoC) \
	$(LIB_XERCESC) \
	$(LIB_EIGEN) \
	$(LIB_NML) \
	$(LIB_PARPACK) \
	$(LIB_ARPACK) \
	$(LIB_TRACEMIN) \
	$(LIB_LAPACK) \
	$(LIB_BLAS) \
	$(LIB_PYTHON) \
#	$(LIB_SVMTL)

LIBS = $(NEMO3D_LIBS) $(NEMO3D_EXTERNLIBS) $(MISCLIB)



#----------------------------------------------------------------
# header information
#----------------------------------------------------------------

NEMO3D_HDRS0_TOP = \
	run3d.h \
	serv3d.h

NEMO3D_HDRS0_BASE = \
	ctracemin_driver.h \
	qtracemin_driver.h \
	arpack_driver.h \
	parpack_driver.h \
	blas_run3d_f.h \
	Boolean.h \
	Cell.h \
	constants_nemo.h \
	eigsys3d.h \
	eigsys3d_par.h \
	geo_struct.h \
	h_cvectr_mult.h \
	h_cmat_mult.h \
	ham_local.h \
	ham3d.h \
	ham3d_final.h \
	HamZincBlende.h \
	HamZincBlendeComplex.h \
	Hmult_stored.h \
	is_def.h \
	mat3d.h \
	mat_def.h \
	MaterialHandle.h \
	matmul.h \
	matsite_init.h \
	MPI_Timing.h \
	nemo3d_entry.h \
	nemo3d_mpi.h \
	opt3d.h \
	qd_struct.h \
	qd_struct_global.h \
	rayleigh_ritz.h \
	nemo3d_mpi.h \
	nemo3d_entry.h \
	run3d_mpi.h \
	strain3d.h \
	StrainParam.h \
	waveFunction.h \
	phon_struct.h \
	phon_output.h \
	phon_Keating.h \
	phon_parpack.h \
	phon_arpack.h \
	dmk_mult.h

NEMO3D_HDRS0_IO = \
	fileformat.h \
	i_3dout_files.h \
	i_mat.h \
	io_utils.h

NEMO3D_HDRS0_MATH = \
	c3tensor.h \
	cmatrix_lib.h \
	cmatrix_nonsym.h \
	cmatrixmatrix.h \
	cvector_lib.h \
	geopoint.h \
	geopointvector.h \
	geopointvectorvector.h \
	i3matrix.h \
	ivectori3matrix.h \
	ivectorvector.h \
	math3d.h \
	math_vdim.h \
	polyshape.h \
	ptinply3.h \
	random.h \
	r3tensor.h \
	realtype.h \
	rmatrix_lib.h \
	rmatrixrvector.h \
	rvector_lib.h \
	tensor.h

NEMO3D_HDRS0_UTIL = \
	mpi_def_local.h \
	util_argv_B.h \
	util_argv_loc_B.h \
	util_command_BX.h \
	util_memory.h \
	util_messages.h \
	util_stringops.h \
	system.h

NEMO3D_HDRS0_HAM = \
	mb_ham_spds.h \
	mb_ham_spds_so.h \
	mb_ham_spds_diag.h \
	mb_ham_spds_E_const.h \
	mb_ham_spds_inclVoffdEdiag.h \
	mb_ham_spds_micro.h \
   mb_ham_spds_VoffdEdiag_constructors.h

ifeq ($(POSTPROCESS), Embedded)
	NEMO3D_HDRS0_POST = \
	AngularMomentum.h \
	Decoherence.h \
	Dipole.h \
	ReadData.h \
	SimpleParser.h \
	Dipole.h \
	Coulomb.h \
	Projection.h \
	TimeReversal.h \
	RandRandom.h \
	Postprocess.h
else 
	NEMO3D_HDRSO_POST =
endif

NEMO3D_HDRS = \
	$(patsubst %,$(NEMO3D_INCLD_TOP)/%, $(NEMO3D_HDRS0_TOP)) \
	$(patsubst %,$(NEMO3D_INCLD_BASE)/%,$(NEMO3D_HDRS0_BASE)) \
	$(patsubst %,$(NEMO3D_INCLD_IO)/%,  $(NEMO3D_HDRS0_IO)) \
	$(patsubst %,$(NEMO3D_INCLD_MATH)/%,$(NEMO3D_HDRS0_MATH)) \
	$(patsubst %,$(NEMO3D_INCLD_UTIL)/%,$(NEMO3D_HDRS0_UTIL)) \
	$(patsubst %,$(NEMO3D_INCLD_HAM)/%, $(NEMO3D_HDRS0_HAM)) \
	$(patsubst %,$(NEMO3D_INCLD_POST)/%, $(NEMO3D_HDRS0_POST))

#----------------------------------------------------------------
# object information
#----------------------------------------------------------------

NEMO3D_OBJS0_TOP = \
	run3d.o

NEMO3D_OBJS0_BASE = \
	ctracemin_driver.o  \
	qtracemin_driver.o  \
	arpack_driver.o \
	parpack_driver.o \
	nemo3d_mpi.o \
	nemo3d_entry.o \
	Cell.o \
	eigsys3d.o \
	eigsys3d_par.o \
	geo_struct.o \
	h_cvectr_mult.o \
	h_cmat_mult.o \
	ham_local.o \
	ham3d.o \
	ham3d_final.o \
	HamZincBlende.o \
	HamZincBlendeComplex.o \
	Hmult_stored.o \
	mat3d.o \
	matmul.o \
	matsite_init.o \
	nemo3d_entry.o \
	nemo3d_mpi.o \
	opt3d.o \
	qd_struct_global.o \
	qd_struct_make.o \
	rayleigh_ritz.o \
	strain3d.o \
	StrainParam.o \
	waveFunction.o \
	phon_struct.o \
	phon_output.o \
	phon_Keating.o \
	phon_parpack.o \
	phon_arpack.o \
	dmk_mult.o

NEMO3D_OBJS0_IO = \
	fileformat.o \
	i_3dout_files.o \
	i_mat.o \
	io_utils.o

NEMO3D_OBJS0_MATH = \
	c3tensor.o \
	cmatrix_lib.o \
	cmatrix_nonsym.o \
	cmatrixmatrix.o \
	cvector_lib.o \
	geopointvector.o \
	geopointvectorvector.o \
	i3matrix.o \
	ivectori3matrix.o \
	ivectorvector.o \
	math3d.o \
	polyshape.o \
	ptinply3.o \
	r3tensor.o \
	random.o \
	rmatrix_lib.o \
	rmatrixrvector.o \
	rvector_lib.o

NEMO3D_OBJS0_UTIL = \
	mpi_def_local.o \
	util_argv_B.o \
	util_command_BX.o \
	util_memory.o \
	util_messages.o \
	util_stringops.o

NEMO3D_OBJS0_HAM = \
	mb_ham_spds.o \
	mb_ham_spds_micro.o

ifeq ($(POSTPROCESS), Embedded)
	NEMO3D_OBJS0_POST = \
	AngularMomentum.o \
	Decoherence.o \
	Dipole.o \
	ReadData.o \
	Coulomb.o \
	Dipole.o \
	Projection.o \
	TimeReversal.o \
	SimpleParser.o \
	RandRandom.o \
	Postprocess.o
else 
	NEMO3D_OBJS0_POST =
endif

NEMO3D_OBJS_TOP  = $(patsubst %,$(NEMO3D_OBJD)/%,$(NEMO3D_OBJS0_TOP))
NEMO3D_OBJS_BASE = $(patsubst %,$(NEMO3D_OBJD)/%,$(NEMO3D_OBJS0_BASE))
NEMO3D_OBJS_IO   = $(patsubst %,$(NEMO3D_OBJD)/%,$(NEMO3D_OBJS0_IO))
NEMO3D_OBJS_MATH = $(patsubst %,$(NEMO3D_OBJD)/%,$(NEMO3D_OBJS0_MATH))
NEMO3D_OBJS_UTIL = $(patsubst %,$(NEMO3D_OBJD)/%,$(NEMO3D_OBJS0_UTIL))
NEMO3D_OBJS_HAM  = $(patsubst %,$(NEMO3D_OBJD)/%,$(NEMO3D_OBJS0_HAM))
NEMO3D_OBJS_POST = $(patsubst %,$(NEMO3D_OBJD)/%,$(NEMO3D_OBJS0_POST))


NEMO3D_OBJS = \
	$(NEMO3D_OBJS_TOP) \
	$(NEMO3D_OBJS_BASE) \
	$(NEMO3D_OBJS_IO) \
	$(NEMO3D_OBJS_MATH) \
	$(NEMO3D_OBJS_UTIL) \
	$(NEMO3D_OBJS_HAM) \
	$(NEMO3D_OBJS_POST)

#----------------------------------------------------------------
# source information
#----------------------------------------------------------------

NEMO3D_SRCD_TOP  = $(NEMO3D_SRCD)/top
NEMO3D_SRCD_BASE = $(NEMO3D_SRCD)/base
NEMO3D_SRCD_IO   = $(NEMO3D_SRCD)/io
NEMO3D_SRCD_MATH = $(NEMO3D_SRCD)/math
NEMO3D_SRCD_UTIL = $(NEMO3D_SRCD)/util
NEMO3D_SRCD_HAM  = $(NEMO3D_SRCD)/ham
ifeq ($(POSTPROCESS), Embedded)
 NEMO3D_SRCD_POST = $(NEMO3D_SRCD)/postprocessing
else
 NEMO3D_SRCD_POST =
endif	

NEMO3D_SRCS0_TOP  = $(patsubst %.o,%.c,$(NEMO3D_OBJS0_TOP))
NEMO3D_SRCS0_BASE = $(patsubst %.o,%.c,$(NEMO3D_OBJS0_BASE)) 
NEMO3D_SRCS0_IO   = $(patsubst %.o,%.c,$(NEMO3D_OBJS0_IO))
NEMO3D_SRCS0_MATH = $(patsubst %.o,%.c,$(NEMO3D_OBJS0_MATH))
NEMO3D_SRCS0_UTIL = $(patsubst %.o,%.c,$(NEMO3D_OBJS0_UTIL))
NEMO3D_SRCS0_HAM  = $(patsubst %.o,%.c,$(NEMO3D_OBJS0_HAM))
NEMO3D_SRCS0_POST = $(patsubst %.o,%.cpp,$(NEMO3D_OBJS0_POST))

NEMO3D_SRCS_TOP  = $(patsubst %,$(NEMO3D_SRCD_TOP)/%,$(NEMO3D_SRCS0_TOP))
NEMO3D_SRCS_BASE = $(patsubst %,$(NEMO3D_SRCD_BASE)/%,$(NEMO3D_SRCS0_BASE))
NEMO3D_SRCS_IO   = $(patsubst %,$(NEMO3D_SRCD_IO)/%,$(NEMO3D_SRCS0_IO))
NEMO3D_SRCS_MATH = $(patsubst %,$(NEMO3D_SRCD_MATH)/%,$(NEMO3D_SRCS0_MATH))
NEMO3D_SRCS_UTIL = $(patsubst %,$(NEMO3D_SRCD_UTIL)/%,$(NEMO3D_SRCS0_UTIL))
NEMO3D_SRCS_HAM  = $(patsubst %,$(NEMO3D_SRCD_HAM)/%,$(NEMO3D_SRCS0_HAM))
NEMO3D_SRCS_POST = $(patsubst %,$(NEMO3D_SRCD_POST)/%,$(NEMO3D_SRCS0_POST))

NEMO3D_SRCS = \
	$(NEMO3D_SRCS_TOP) \
	$(NEMO3D_SRCS_BASE) \
	$(NEMO3D_SRCS_IO) \
	$(NEMO3D_SRCS_MATH) \
	$(NEMO3D_SRCS_UTIL) \
	$(NEMO3D_SRCS_HAM) \
	$(NEMO3D_SRCS_POST)

#----------------------------------------------------------------
#  Executable information
#----------------------------------------------------------------

NEMO3D_EXEC     = $(NEMO3D_BIND)/nemo3d-$(PLATFORM).ex
NEMO3D_EXEC_LNK = $(NEMO3D_BIND)/nemo3d.ex

