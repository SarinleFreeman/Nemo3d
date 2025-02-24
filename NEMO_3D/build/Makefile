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
# NEMO-3D Application Makefile
#================================================================


#----------------------------------------------------------------
# Check for required system environment variables
#----------------------------------------------------------------

ifeq ($(HOSTTYPE),)
#  $(error "Undefined environment variable, HOSTTYPE")
endif

ifeq ($(PYTHONPATH),)
#   $(error "Undefined environment variable, PYTHONPATH")
endif

ifeq ($(XERCESCROOT),)
#   $(error "Undefined environment variable, XERCESCROOT")
endif

ifeq ($(MPI_HOME),)
#   $(@ echo "NOTE - Undefined environment variable, MPI_HOME")
endif



#----------------------------------------------------------------
# System configuration: compiler and libraries
#----------------------------------------------------------------

include ./make.inc


# ---------------------------------------------------------------
# Xerces-C++ runConfigure script
# ---------------------------------------------------------------

XERCES_RUNCONFIGURE = ./configure

ifeq ($(OS),Linux)
	XERCES_RUNCONFIGURE = env CC='$(XERCES_CC)' CFLAGS='$(XERCES_CCOPT)' CXX='$(XERCES_CC)' CXXFLAGS='$(XERCES_CCOPT)' ./configure
endif

ifeq ($(OS),Darwin)
	XERCES_RUNCONFIGURE = ./runConfigure -p macosx -n native -c $(XERCES_CC) -x $(XERCES_CC)
endif

ifeq ($(OS),IRIX64)
	XERCES_RUNCONFIGURE = ./runConfigure -p irix -m inmem -n fileonly -c $(XERCES_CC) -x $(XERCES_CC)
endif

ifeq ($(OS),AIX)
	AR=ar
	XERCES_RUNCONFIGURE = ./runConfigure -d -p aix -m inmem -n fileonly -r pthread -b$${OBJECT_MODE-32} -c mpcc_r -x mpCC_r
endif


ifeq ($(OS),CYGWIN_NT-5.1)
	XERCES_RUNCONFIGURE = ./runConfigure -p cygwin  -c $(XERCES_CC) -x $(XERCES_CC)
endif

#----------------------------------------------------------------
# User rules
#----------------------------------------------------------------

# Default rule
all: nemo3d
#all: nemo3d postprocess

# Libraries --------------------------------

$(LIB_BLAS):
	(cd $(BLAS_ROOT)/SRC; $(MAKE) PLAT='_$(PLATFORM)' MAKE='$(MAKE)' SHELL='$(SHELL)' \
	FORTRAN='$(FORTRAN)' OPTS='$(LAPACK_OPTS)' LOADER='$(FORTRAN)' LOADOPTS='$(LAPACK_LDOPTS)' \
	ARCH='$(AR)' ARCHFLAGS='$(AR_OPTIONS)' RANLIB='$(RANLIB)' \
	BLASLIB='../../$(LIB_BLAS)')
	@ echo "==== BLAS library archive up to date ===="

$(LIB_LAPACK):
	(cd $(LAPACK_ROOT); $(MAKE) install PLAT='_$(PLATFORM)' MAKE='$(MAKE)' SHELL='$(SHELL)' \
	FORTRAN='$(FORTRAN)' OPTS='$(LAPACK_OPTS)' NOOPT='$(LAPACK_NOOPT)' LOADER='$(FORTRAN)' LOADOPTS='$(LAPACK_LDOPTS)')
	(cd $(LAPACK_ROOT)/SRC; $(MAKE) PLAT='_$(PLATFORM)' MAKE='$(MAKE)' SHELL='$(SHELL)' \
	FORTRAN='$(FORTRAN)' OPTS='$(LAPACK_OPTS)' NOOPT='$(LAPACK_NOOPT)' LOADER='$(FORTRAN)' LOADOPTS='$(LAPACK_LDOPTS)' \
	ARCH='$(AR)' ARCHFLAGS='$(AR_OPTIONS)' RANLIB='$(RANLIB)' \
	LAPACKLIB='$(LIB_LAPACK)')
	@ echo "==== LAPACK library archive up to date ===="

##################################################################
$(LIB_SCTRACEMIN):
	(cd $(TRACEMIN_ROOT); $(MAKE) sctracemin_lib STRACEMINLIB='$(LIB_SCTRACEMIN)' DIRS='$(TRACEMIN_ROOT)' FC='$(FORTRAN90)' \
                                           FFLAGS='$(SCTRACEMIN_OPTS)' MAKE='$(MAKE)' RANLIB='$(RANLIB)')
	@ echo "==== SCRACEMIN library archive up to date ===="

$(LIB_PCTRACEMIN):
	(cd $(TRACEMIN_ROOT); $(MAKE) pctracemin_lib PTRACEMINLIB='$(LIB_PCTRACEMIN)' DIRS='$(TRACEMIN_ROOT)' FC='$(FORTRAN90)' \
	                                   FFLAGS='$(PCTRACEMIN_OPTS)' MAKE='$(MAKE)' RANLIB='$(RANLIB)')
	@ echo "==== PCRACEMIN library archive up to date ===="

$(LIB_SQTRACEMIN):
	(cd $(TRACEMIN_ROOT); $(MAKE) sqtracemin_lib STRACEMINLIB='$(LIB_SQTRACEMIN)' DIRS='$(TRACEMIN_ROOT)' FC='$(FORTRAN90)' \
	                                   FFLAGS='$(SQTRACEMIN_OPTS)' MAKE='$(MAKE)' RANLIB='$(RANLIB)')
	@ echo "==== SQRACEMIN library archive up to date ===="

$(LIB_PQTRACEMIN):
	(cd $(TRACEMIN_ROOT); $(MAKE) pqtracemin_lib PTRACEMINLIB='$(LIB_PQTRACEMIN)' DIRS='$(TRACEMIN_ROOT)' FC='$(FORTRAN90)' \
	                                    FFLAGS='$(PQTRACEMIN_OPTS)' MAKE='$(MAKE)' RANLIB='$(RANLIB)')
	@ echo "==== PQRACEMIN library archive up to date ===="
##################################################################

$(LIB_ARPACK):
	(cd $(ARPACK_ROOT); $(MAKE) lib ARPACKLIB='$(LIB_ARPACK)' \
	DIRS='$(ARPACK_ROOT)/UTIL $(ARPACK_ROOT)/SRC' \
	FC='$(FORTRAN)' FFLAGS='$(ARPACK_OPTS)' MAKE='$(MAKE)' RANLIB='$(RANLIB)')
	@ echo "==== ARPACK library archive up to date ===="

$(LIB_PARPACK):
ifeq ($(MPI_HOME),)
else
	(cp $(MPI_HOME)/include/Intel/mpif.h $(PARPACK_ROOT)/SRC/MPI; cp $(MPI_HOME)/include/Intel/mpif.h $(PARPACK_ROOT)/UTIL/MPI)
endif
	(cd $(ARPACK_ROOT); $(MAKE) plib PARPACKLIB='$(LIB_PARPACK)' \
	PUTILdir="$(PARPACK_ROOT)/UTIL/MPI" PSRCdir="$(PARPACK_ROOT)/SRC/MPI" \
	FC='$(FORTRAN)' FFLAGS='$(ARPACK_OPTS)' MAKE='$(MAKE)' RANLIB='$(RANLIB)')
	@ echo "==== PARPACK library archive up to date ===="

$(LIB_NML):
	(cd $(NML_ROOT)/src; $(MAKE) PLATFORM='$(PLATFORM)' \
	 CC='$(CC)' OPTIONS='$(NML_OPTIONS)' DEFINES='$(NML_DEFINES)' \
	 AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)' MAKE='$(MAKE)')
	@ echo "==== NML library archive up to date ===="

$(LIB_EIGEN):
	(cd $(EIGEN_ROOT); $(MAKE) PLATFORM='$(PLATFORM)' \
	CCMKD='$(CCMKD)' CC='$(CC)' CCOPT='$(CCOPT)' LD='$(LD)' AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)' MAKE='$(MAKE)' \
	LIBBLAS='$(LIB_BLAS)' LIBLAPACK='$(LIB_LAPACK)' LIBNML='$(LIB_NML)')
	@ echo "==== EIGEN library archive up to date ===="

$(LIB_XERCESC):
	@ echo $(XERCES_RUNCONFIGURE)
	( cd $(XERCESCROOT)/src/xercesc; $(XERCES_RUNCONFIGURE); $(MAKE) )
	( cd $(XERCESCROOT)/lib/; $(AR) -$(AR_OPTIONS) $(LIB_XERCESC) ../obj/$(OSTYPE)/*.o; $(RANLIB) $(LIB_XERCESC) )
	@ echo "==== XERCES-C++ library archive up to date ===="

$(LIB_XMLtoC):
	( cd $(XMLtoC_ROOT); $(MAKE) applib PLATFORM='$(PLATFORM)' MAKE='$(MAKE)' APP='ioNemo3D' \
	CCMKD='$(XMLtoC_CCMKD)' CC='$(XMLtoC_CC)' CCOPT='$(XMLtoC_CCOPT)' LD='$(XMLtoC_LD)' LD_OPT='$(XMLtoC_LDOPT)' AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)' \
	PYTHON_ROOT='$(PYTHON_ROOT)' PYTHON_EXPORT_INCL='$(PYTHON_EXPORT_INCL)' PYTHON_EXPORT_LIB='$(LIB_PYTHON)' \
	XERCESCROOT='$(XERCESCROOT)' XERCESC_EXPORT_INCL='$(XERCESC_EXPORT_INCL)' XERCESC_EXPORT_LIB='$(LIB_XERCESC)' \
	MISCLIB='$(MISCLIB)' LIB_PATHS='$(LIB_PATHS)' )
	@ echo "==== XMLtoC library archive up to date ===="


# Rules for main components --------------------------------

############################################################
sctracemin: $(LIB_SCTRACEMIN)
pctracemin: $(LIB_PCTRACEMIN)
sqtracemin: $(LIB_SQTRACEMIN)
pqtracemin: $(LIB_PQTRACEMIN)
tracemin:   $(LIB_TRACEMIN)
###########################################################

blas: $(LIB_BLAS)

lapack: $(LIB_LAPACK)

arpack: $(LIB_ARPACK)

parpack: $(LIB_PARPACK)

nml: $(LIB_NML)

eigen: $(LIB_EIGEN)

xercesc: $(LIB_XERCESC)

XMLtoC: $(LIB_XMLtoC)


################################################################# 
# taken out from nemo3d:
nemo3d: blas lapack arpack parpack nml eigen xercesc XMLtoC tracemin 
	(cd $(NEMO3D_ROOT); $(MAKE) nemo3d PLATFORM='$(PLATFORM)' POSTPROCESS='$(POSTPROCESS)' CCMKD='$(CCMKD)' \
	CC='$(CC)' CCOPT='$(CCOPT)' LD='$(LD)' LDOPT='$(LDOPT)' AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)' \
	OSTYPE='$(OSTYPE)')

postprocess:
	(cd $(POSTPROCESS_ROOT); $(MAKE) PLATFORM='$(PLATFORM)' CCMKD='$(CCMKD)' CC='$(CC)' CCOPT='$(CCOPT)' \
        OSTYPE='$(OSTYPE)' LD='$(LD)'  LDOPT='$(LDOPT)' AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)')
################################################################# 
#added test option (runs in ./nemo_3d/example/TestRun)
test:
	(cd  $(TEST_ROOT);  $(MAKE) )


# Cleanup --------------------------------

clean: clean_nemo3d clean_XMLtoC clean_eigen clean_nml clean_arpackparpack clean_lapackblas clean_xercesc clean_test
	@ echo "==== Clean nemo3d, XMLtoC, XERCES-C++, eigen, nml, ARPACK, PARPACK, LAPACK, and BLAS objects & dependencies ===="

################################################################# 
#cleanall:cleanall_nemo3d cleanall_XMLtoC cleanall_eigen cleanall_nml cleanall_arpackparpack cleanall_lapackblas cleanall_xercesc ####################################### 
cleanall: cleanall_nemo3d cleanall_XMLtoC cleanall_eigen cleanall_nml cleanall_arpackparpack cleanall_lapackblas cleanall_xercesc cleanall_tracemin \
	@ echo "==== Clean all nemo3d, XMLtoC, XERCES-C++, eigen, nml, ARPACK, PARPACK, LAPACK, and BLAS executables, library archives, objects, & dependencies ===="


clean_lapack: clean_lapackblas
clean_blas: clean_lapackblas
cleanall_lapack: cleanall_lapackblas
cleanall_blas: cleanall_lapackblas

clean_lapackblas:
	(cd $(LAPACK_ROOT); $(MAKE) cleanlib PLAT='_$(PLATFORM)' MAKE='$(MAKE)' \
	SHELL='$(SHELL)' LAPACKLIB='$(LIB_LAPACK)' BLASLIB='../../$(LIB_BLAS)')
	@ echo "  ...... Finished cleaning LAPACK and BLAS"

cleanall_lapackblas: clean_lapackblas
	rm -f $(LIB_LAPACK)
	rm -f $(LIB_BLAS)
	@ echo "  ...... Finished cleaning all LAPACK and BLAS"

clean_arpack:  clean_arpackparpack 
clean_parpack: clean_arpackparpack
cleanall_arpack:  cleanall_arpackparpack 
cleanall_parpack: cleanall_arpackparpack


clean_arpackparpack:
	(cd $(ARPACK_ROOT); $(MAKE) clean BLASdir="$(ARPACK_ROOT)/BLAS" LAPACKdir="$(ARPACK_ROOT)/LAPACK" \
	UTILdir="$(ARPACK_ROOT)/UTIL" SRCdir="$(ARPACK_ROOT)/SRC" \
	PUTILdir="$(PARPACK_ROOT)/UTIL/MPI" PSRCdir="$(PARPACK_ROOT)/SRC/MPI" \
	MAKE="$(MAKE)")
	@ echo "  ...... Finished cleaning ARPACK and PARPACK"

cleanall_arpackparpack: clean_arpackparpack
	rm -f $(LIB_ARPACK)
	rm -f $(LIB_PARPACK)
	@ echo "  ...... Finished cleaning all ARPACK and PARPACK"

################################################################################
cleanall_tracemin:
	(cd $(TRACEMIN_ROOT); $(MAKE) clean_tracemin)
	@echo "   ...... Finished cleaning TRACEMIN drivers"

cleanall_sctracemin:
	(cd $(TRACEMIN_ROOT); $(MAKE) clean_sctracemin  STRACEMINLIB='$(LIB_SCTRACEMIN)')  
	@ echo "  ...... Finished cleaning all Sequential TRACEMIN"

cleanall_pctracemin:
	(cd $(TRACEMIN_ROOT); $(MAKE) clean_pctracemin  PTRACEMINLIB='$(LIB_PCTRACEMIN)')
	@ echo "  ...... Finished cleaning all Parallel TRACEMIN"

cleanall_sqtracemin:
	(cd $(TRACEMIN_ROOT); $(MAKE) clean_sqtracemin  STRACEMINLIB='$(LIB_SQTRACEMIN)')
	@ echo "  ...... Finished cleaning all Sequential TRACEMIN"

cleanall_pqtracemin:
	(cd $(TRACEMIN_ROOT); $(MAKE) clean_pqtracemin  PTRACEMINLIB='$(LIB_PQTRACEMIN)')
	@ echo "  ...... Finished cleaning all Parallel TRACEMIN"

clean_test:
	(cd  $(TEST_ROOT);  $(MAKE) clean1; $(MAKE) clean2)
	@ echo "  ...... Finished cleaning TestRun"

cleanall_test:
	(cd  $(TEST_ROOT);  $(MAKE) cleanall)
	@ echo "  ...... Finished cleaning all TestRun"
################################################################################

clean_nml:
	(cd $(NML_ROOT)/src; $(MAKE) clean PLATFORM='$(PLATFORM)' MAKE='$(MAKE)')
	@ echo "  ...... Finished cleaning nml"

cleanall_nml:
	(cd $(NML_ROOT)/src; $(MAKE) sterile PLATFORM='$(PLATFORM)' MAKE='$(MAKE)')
	@ echo "  ...... Finished cleaning all nml"

clean_eigen:
	(cd $(EIGEN_ROOT);  $(MAKE) clean  PLATFORM='$(PLATFORM)' MAKE='$(MAKE)')
	@ echo "  ...... Finished cleaning eigen"

cleanall_eigen:
	(cd $(EIGEN_ROOT);  $(MAKE) cleanall PLATFORM='$(PLATFORM)' MAKE='$(MAKE)')
	@ echo "  ...... Finished cleaning all eigen"

clean_xercesc:
	(cd $(XERCESCROOT)/src/xercesc; $(MAKE) clean)
	@ echo "  ...... Finished cleaning XERCES-C++"

cleanall_xercesc:
	(cd $(XERCESCROOT)/src/xercesc; $(MAKE) distclean)
	rm -f $(XERCESCROOT)/obj/libxerces-c.so.*.map
	rm -f $(LIB_XERCESC)
	rm -f $(XERCESCROOT)/lib/libxerces-c.so*
	@ echo "  ...... Finished cleaning all XERCES-C++"

clean_XMLtoC:
	(cd $(XMLtoC_ROOT); $(MAKE) clean PLATFORM='$(PLATFORM)' MAKE='$(MAKE)' APP='ioNemo3D' \
	CCMKD='$(XMLtoC_CCMKD)' CC='$(XMLtoC_CC)' CCOPT='$(XMLtoC_CCOPT)' LD='$(XMLtoC_LD)' AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)' \
	PYTHON_ROOT='$(PYTHON_ROOT)' PYTHON_EXPORT_INCL='$(PYTHON_EXPORT_INCL)' PYTHON_EXPORT_LIB='$(LIB_PYTHON)' \
	XERCESCROOT='$(XERCESCROOT)' XERCESC_EXPORT_INCL='$(XERCESC_EXPORT_INCL)' XERCESC_EXPORT_LIB='$(LIB_XERCESC)' \
	MISCLIB='$(MISCLIB)' LIB_PATHS='$(LIB_PATHS)' )
	@ echo "  ...... Finished cleaning XMLtoC"

cleanall_XMLtoC:
	(cd $(XMLtoC_ROOT); $(MAKE) cleanall PLATFORM='$(PLATFORM)' MAKE='$(MAKE)' APP='ioNemo3D' \
	CCMKD='$(XMLtoC_CCMKD)' CC='$(XMLtoC_CC)' CCOPT='$(XMLtoC_CCOPT)' LD='$(XMLtoC_LD)' AR='$(AR)' AR_OPTIONS='$(AR_OPTIONS)' \
	PYTHON_ROOT='$(PYTHON_ROOT)' PYTHON_EXPORT_INCL='$(PYTHON_EXPORT_INCL)' PYTHON_EXPORT_LIB='$(LIB_PYTHON)' \
	XERCESCROOT='$(XERCESCROOT)' XERCESC_EXPORT_INCL='$(XERCESC_EXPORT_INCL)' XERCESC_EXPORT_LIB='$(LIB_XERCESC)' \
	MISCLIB='$(MISCLIB)' LIB_PATHS='$(LIB_PATHS)' )
	@ echo "  ...... Finished cleaning all XMLtoC"

clean_nemo3d:
	(cd $(NEMO3D_ROOT); $(MAKE) clean PLATFORM='$(PLATFORM)')
	@ echo "  ...... Finished cleaning nemo3d"

cleanall_nemo3d:
	(cd $(NEMO3D_ROOT); $(MAKE) cleanall PLATFORM='$(PLATFORM)')
	@ echo "  ...... Finished cleaning all nemo3d"


# Clear the current BUILD_TYPE's object files of
# blas, lapack, arpack, parpack, and xerces-c
# so it can be safe to build another *coexisting* BUILT_TYPE.
coexisting_build_type:  clean_lapackblas clean_arpackparpack clean_xercesc
	@ echo " ==== DONE! It is safe to build another coexisting BUILT_TYPE now. ===="


# Complete NEMO3D Library --------------------------------

# Creates a complete library of NEMO3D by merging individual sub libraries
# of NEMO3D together.
library: blas lapack arpack parpack nml eigen xercesc XMLtoC nemo3d
	(rm -f $(LIB_NEMO3D); \
	mkdir -p __TEMP; \
	cd __TEMP; \
	if [ "$(LIB_BLAS)" != "" ]; then echo "# Extracting ../$(LIB_BLAS)"; $(AR) -x ../$(LIB_BLAS); fi; \
	if [ "$(LIB_LAPACK)" != "" ]; then echo "# Extracting ../$(LIB_LAPACK)"; $(AR) -x ../$(LIB_LAPACK); fi; \
	if [ "$(LIB_ARPACK)" != "" ]; then echo "# Extracting $(LIB_ARPACK)"; $(AR) -x $(LIB_ARPACK); fi; \
	if [ "$(LIB_PARPACK)" != "" ]; then echo "# Extracting $(LIB_PARPACK)"; $(AR) -x $(LIB_PARPACK); fi; \
	echo "# Extracting $(LIB_NML)"; \
	$(AR) -x $(LIB_NML); \
	echo "# Extracting $(LIB_EIGEN)"; \
	$(AR) -x $(LIB_EIGEN); \
	echo "# Extracting $(LIB_XERCESC)"; \
	$(AR) -x $(LIB_XERCESC); \
	echo "# Extracting $(LIB_XMLtoC)"; \
	$(AR) -x $(LIB_XMLtoC); \
	echo "# Extracting $(NEMO3D_ROOT)/lib/$(PLATFORM)/*.a"; \
	for file in $(NEMO3D_ROOT)/lib/$(PLATFORM)/*.a; do echo "  # $$file"; $(AR) -x $$file; done; \
	echo "# Creating complete static library $(LIB_NEMO3D)"; \
	$(AR) -$(AR_OPTIONS) $(LIB_NEMO3D) *.o; \
	cd ..; \
	echo "# Removing temporary files"; \
	rm -rf __TEMP; \
	echo "# DONE!")


# Archiving --------------------------------

tarball:
	(cd ../..; tar -cvf - NEMO3D | gzip -rv9 > NEMO3D/build/NEMO3D.tar.gz)

dist_tarball:
	(cd ../..; \
	find NEMO3D/* | grep 'CVS' > __OMIT; \
	tar -cvf - NEMO3D --exclude-from __OMIT | gzip -rv9 > NEMO3D.tar.gz; \
	rm -f __OMIT)


