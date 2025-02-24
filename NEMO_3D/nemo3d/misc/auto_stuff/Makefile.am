## Makefile.am -- Process this file with automake to produce Makefile.in



# The name of the main project application
APP_NAME = NEMO3D

SUBDIRS = nml eigen XMLtoC nemo3d

#----------------------------------------------------------------
# System configuration: compiler and libraries
#----------------------------------------------------------------

NML_HOME    = ../nml
EIGEN_HOME  = ../eigen
XMLtoC_HOME = ../XMLtoC
NEMO3D_HOME = ../nemo3d

# -----------------------------------------------------------------------

SRCS = \
        cmacopt.c \
        io.c \
        lanczos.c \
        macopt.c \
        macopt_par.c \
        ritz.c

HDRS = \
        cmacopt.h \
        io.h \
        lanczos.h \
        macopt.h \
        macopt_par.h \
	mpi_local.h \
        ritz.h


# ------------------ Global Assignments ------------------

INCLUDES = \
	-I$(top_srcdir)/include \
	-I$(NML_HOME)/include \
	-I$(MPI_HOME)/include
LDADD = 
LDFLAGS = 
TESTS = $(check_PROGRAMS)


# ------------------ Primitive Assignments ------------------

SUBDIRS = 
lib_LIBRARIES = libeigen.a
check_PROGRAMS = 
include_HEADERS = $(HDRS)


# ------------------ Property Assignments ------------------

libeigen_a_SOURCES = $(SRCS) $(HDRS)

