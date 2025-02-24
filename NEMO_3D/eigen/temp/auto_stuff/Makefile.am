## Makefile.am -- Process this file with automake to produce Makefile.in


# The name of the main project application
PJCT = EIGEN


#----------------------------------------------------------------
#  Default platform specific project build options
#----------------------------------------------------------------

# Default Project Build if none was specified
#if COND_PJCTBLD_null
ifeq ($($(PJCT)_BLD),)
   $(PJCT)_BLD = gnuLinux
endif

# GNU Linux
#if COND_PJCTBLD_gnuLinux
ifeq ($($(PJCT)_BLD),gnuLinux)
   CFLAGS   = -Wall
   CC       = /usr/bin/g++ -c
   CCLD     = /usr/bin/g++
   AR       = /usr/bin/ar -r
   MISCLIB0 = -lm
endif

# Absoft MPI Fast
#if COND_PJCTBLD_absoft_mpi_fast
ifeq ($($(PJCT)_BLD),absoft_mpi_fast)
   CFLAGS   = -Wall -O3 -DFORTRAN_ALLCAPS -DUSEMPI
   CC       = $(MPI_HOME)/bin/mpiCC -c
   CCLD     = $(MPI_HOME)/bin/mpiCC 
   AR       = /usr/bin/ar -r
   MISCLIB0 = -lm
endif

#AM_CONDITIONAL(COND_PJCTBLD_gnuLinux, [test $($(PJCT)_BLD) = gnuLinux])
#AM_CONDITIONAL(COND_PJCTBLD_absoft_mpi_fast, [test $($(PJCT)_BLD) = absoft_mpi_fast])
#AM_CONDITIONAL(COND_PJCTBLD_null, [test $($(PJCT)_BLD) = ])


#----------------------------------------------------------------
# System configuration: compiler and libraries
#----------------------------------------------------------------

NML_HOME = ../nml

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

