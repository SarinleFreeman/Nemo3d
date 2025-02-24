include ../build/make.inc

#Names of the standalone executables
PCTRACEMINNAME = pctracemin
PQTRACEMINNAME = pqtracemin
SCTRACEMINNAME = sctracemin
SQTRACEMINNAME = sqtracemin

#Directories
OBJDIR = Obj
SRCDIR = Src
LIBDIR = Lib
BINDIR = Bin

#Compilers 
#PFC = mpiifort  #Fortran/MPI Intel compiler                                                       (Intel cluster - nocona[1|2|3|4].cs.purdue.edu)        
#SFC = ifort     #Fortran Intel compiler                                                           (Intel cluster - nocona[1|2|3|4].cs.purdue.edu)
#SCC = icc       #C           Intel compiler                                                       (Intel cluster - nocona[1|2|3|4].cs.purdue.edu)
#PFC = mpxlf90   #Fortran/MPI IBM   compiler, should change all "use mpi" into "include "mpif.h""  (IBM SP supercomputer frontend - cloud.cc.purdue.edu)
#PFC = mpxlf90_r #Fortran/MPI IBM   compiler, accepts "use mpi" syntax                             (IBM SP supercomputer frontend - cloud.cc.purdue.edu)
#SFC = xlf90     #Fortran IBM       compiler                                                       (IBM SP supercomputer frontend - cloud.cc.purdue.edu)
#SCC = xlc       #C           IBM   compiler							   (IBM SP supercomputer frontend - cloud.cc.purdue.edu)
PFC = $(FORTRAN90)  #Fortran/MPI Intel compiler                 (NEMO3D - lear.rcac.purdue.edu)               
PCC = $(CC)   #C/MPI       Intel compiler                 (NEMO3D - lear.rcac.purdue.edu)
#SFC = /opt/intel/fce/9.0.025/bin/ifort                 #Fortran Intel compiler (NEMO3D - lear.rcac.purdue.edu)
SFC = $(FORTRAN90)                 #Fortran Intel compiler (NEMO3D - lear.rcac.purdue.edu)
#SCC = /opt/intel/cce/9.0.025/bin/icc                   #C       Intel compiler (NEMO3D - lear.rcac.purdue.edu)
SCC = $(CC)                   #C       Intel compiler (NEMO3D - lear.rcac.purdue.edu)


#Commands used to create a library from the object files
CD            = cd  
AR            = ar
ARFLAGS       = rcvs 
RANLIB        = ranlib
PCTRACEMINLIB = ../$(LIBDIR)/libpctracemin_$(PLATFORM).a
PQTRACEMINLIB = ../$(LIBDIR)/libpqtracemin_$(PLATFORM).a
SCTRACEMINLIB = ../$(LIBDIR)/libsctracemin_$(PLATFORM).a
SQTRACEMINLIB = ../$(LIBDIR)/libsqtracemin_$(PLATFORM).a

#Compiler Flags
FFLAGS_v2 =  $(FFLAGS) -module ../$(OBJDIR) #-nus option for Intel compiler can be used for Fortran-C compatibility
#FFLAGS =  -qsave -bmaxdata:0x80000000 -qarch=pwr3 -qtune=pwr3 \
           -qcache=level=1:type=i:assoc=128:cost=12:line=128:size=32  -qcache=level=1:type=d:assoc=128:cost=12:line=128:size=64     \
           -qcache=level=2:type=i:assoc=1:cost=40:line=2048:size=8192 -qcache=level=2:type=d:assoc=1:cost=40:line=2048:size=8192    \
	   #optimization options for IBM compiler for running on Purdue University IBM SP
CFLAGS_v2 = 

#Include the MKL-Intel libraries (place lapack library first since it depends on blas) 
#INCF = -I/opt/intel/mkl172cluster/include -I/opt/intel_mpi_10/include64                                 #(nocona[1|2|3|4].cs.purdue.edu)
#INCF = -I/opt/intel/mkl/8.0.1/include/ -I/opt/mpich-1.2.6/64/p4-intel-9.0.025/include/f90choice/         #(lear.rcac.purdue.edu)
INCF = /opt/intel/mkl/9.1/include -I/opt/mpich-1.2.6/64/p4-intel-9.0.025/include/f90choice/         #(lear.rcac.purdue.edu)

#Library Path 
#LIBPATH = -L/opt/intel/mkl72cluster/lib/em64t/   #(nocona[1|2|3|4].cs.purdue.edu)      
#LIBPATH = -L/opt/lapack/lib                      #(cloud.cc.purdue.edu) 
#LIBPATH = -L/opt/intel/mkl/8.0.1/lib/em64t/      #(lear.rcac.purdue.edu)
LIBPATH = -L/opt/intel/mkl/9.1/lib/em64t          #(lear.rcac.purdue.edu)

#Libraries: MKL - Lapack, Blas and others
#LIB = -lmkl_lapack -lmkl_em64t -lguide -lpthread #(nocona[1|2|3|4].cs.purdue.edu)                                
#LIB =  -llapack -lessl                           #(cloud.cc.purdue.edu)
LIB = -lmkl_lapack -lmkl_em64t -lguide -lpthread  #(lear.rcac.purdue.edu)

#Compilation Rules
#Standalone
all: pctracemin pqtracemin 
	@echo "==== ALL TRACEMIN STANDALONE archives up to date ===="

pctracemin: 
	($(CD) $(SRCDIR); $(MAKE) pctracemin    CD='$(CD)' SRCDIR='../$(SRCDIR)' OBJDIR='../$(OBJDIR)' LIBDIR='../$(LIBDIR)' BINDIR='../$(BINDIR)' \
						PCTRACEMINNAME='$(PCTRACEMINNAME)'\
						FC='$(PFC)'          \
						CC='$(PCC)'          \
						FFLAGS='$(FFLAGS_v2)'   \
						CFLAGS='$(CFLAGS_v2)'   \
						INCF='$(INCF)'       \
						LIBPATH='$(LIBPATH)' \
						LIB='$(LIB)')
	@echo "==== PCTRACEMIN STANDALONE archives up to date ===="

pqtracemin:
	($(CD) $(SRCDIR); $(MAKE) pqtracemin	CD='$(CD)' SRCDIR='../$(SRCDIR)' OBJDIR='../$(OBJDIR)' LIBDIR='../$(LIBDIR)' BINDIR='../$(BINDIR)' \
						PQTRACEMINNAME='$(PQTRACEMINNAME)'\
						FC='$(PFC)'         \
						CC='$(PCC)'         \
						FFLAGS='$(FFLAGS_v2)'  \
						CFLAGS='$(CFLAGS_v2)'  \
						INCF='$(INCF)'      \
						LIBPATH='$(LIBPATH)'\
						LIB='$(LIB)')
	@echo "==== PQTRACEMIN STANDALONE archives up to date ===="

sctracemin:
	($(CD) $(SRCDIR); $(MAKE) sctracemin	CD='$(CD)' SRCDIR='../$(SRCDIR)' OBJDIR='../$(OBJDIR)' LIBDIR='../$(LIBDIR)' BINDIR='../$(BINDIR)' \
						SCTRACEMINNAME='$(SCTRACEMINNAME)'\
						FC='$(SFC)'         \
						CC='$(SCC)'         \
						FFLAGS='$(FFLAGS_v2)'  \
						CFLAGS='$(CFLAGS_v2)'  \
						INCF='$(INCF)'      \
						LIBPATH='$(LIBPATH)'\
						LIB='$(LIB)')
	@echo "==== SCTRACEMIN STANDALONE archives up to date ===="

sqtracemin: 
	($(CD) $(SRCDIR); $(MAKE) sqtracemin	CD='$(CD)' SRCDIR='../$(SRCDIR)' OBJDIR='../$(OBJDIR)' LIBDIR='../$(LIBDIR)' BINDIR='../$(BINDIR)' \
						SQTRACEMINNAME='$(SQTRACEMINNAME)'\
						FC='$(SFC)'         \
						CC='$(SCC)'         \
						FFLAGS='$(FFLAGS_v2)'  \
						CFLAGS='$(CFLAGS_v2)'  \
						INCF='$(INCF)'      \
						LIBPATH='$(LIBPATH)'\
						LIB='$(LIB)')
	@echo "==== SQTRACEMIN STANDALONE archives up to date ===="


#Libraries
lib: pctracemin_lib pqtracemin_lib sctracemin_lib sqtracemin_lib
	@echo "==== ALL TRACEMIN LIBRARIES archives up to date ===="
 
pctracemin_lib: 
	($(CD) $(SRCDIR); $(MAKE) pctracemin_lib CD='$(CD)' SRCDIR='../$(SRCDIR)' OBJDIR='../$(OBJDIR)' LIBDIR='../$(LIBDIR)' BINDIR='../$(BINDIR)' \
						PCTRACEMINLIB='$(PTRACEMINLIB)'\
						FC='$(PFC)'         \
						CC='$(PCC)'         \
						FFLAGS='$(FFLAGS_v2)'  \
						CFLAGS='$(CFLAGS_v2)'  \
						INCF='$(INCF)'      \
						LIBPATH='$(LIBPATH)'\
						LIB='$(LIB)'        \
						AR='$(AR)'          \
						ARFLAGS='$(ARFLAGS)'\
						RANLIB='$(RANLIB)')
	@echo "==== PCTRACEMIN LIBRARY archives up to date ===="
 
pqtracemin_lib:
	($(CD) $(SRCDIR); $(MAKE) pqtracemin_lib CD='$(CD)' SRCDIR='../$(SRCDIR)' OBJDIR='../$(OBJDIR)' LIBDIR='../$(LIBDIR)' BINDIR='../$(BINDIR)' \
						PQTRACEMINLIB='$(PTRACEMINLIB)'\
						FC='$(PFC)'         \
						CC='$(PCC)'         \
						FFLAGS='$(FFLAGS_v2)'  \
						CFLAGS='$(CFLAGS_v2)'  \
						INCF='$(INCF)'      \
                        LIBPATH='$(LIBPATH)'\
						LIB='$(LIB)'        \
						AR='$(AR)'          \
						ARFLAGS='$(ARFLAGS)'\
						RANLIB='$(RANLIB)')
	@echo "==== PQTRACEMIN LIBRARY archives up to date ===="
 
sctracemin_lib:
	($(CD) $(SRCDIR); $(MAKE) sctracemin_lib CD='$(CD)' SRCDIR='../$(SRCDIR)' OBJDIR='../$(OBJDIR)' LIBDIR='../$(LIBDIR)' BINDIR='../$(BINDIR)' \
						SCTRACEMINLIB='$(STRACEMINLIB)'\
						FC='$(SFC)'         \
						CC='$(SCC)'         \
						FFLAGS='$(FFLAGS_v2)'  \
						CFLAGS='$(CFLAGS_v2)'  \
						INCF='$(INCF)'      \
						LIBPATH='$(LIBPATH)'\
						LIB='$(LIB)'        \
						AR='$(AR)'          \
						ARFLAGS='$(ARFLAGS)'\
						RANLIB='$(RANLIB)')
	@echo "==== SCTRACEMIN LIBRARY archives up to date ===="
 
sqtracemin_lib: 
	($(CD) $(SRCDIR); $(MAKE) sqtracemin_lib CD='$(CD)' SRCDIR='../$(SRCDIR)' OBJDIR='../$(OBJDIR)' LIBDIR='../$(LIBDIR)' BINDIR='../$(BINDIR)' \
						SQTRACEMINLIB='$(STRACEMINLIB)'\
						FC='$(SFC)'         \
						CC='$(SCC)'         \
						FFLAGS='$(FFLAGS_v2)'  \
						CFLAGS='$(CFLAGS_v2)'  \
						INCF='$(INCF)'      \
						LIBPATH='$(LIBPATH)'\
						LIB='$(LIB)'        \
						AR='$(AR)'          \
						ARFLAGS='$(ARFLAGS)'\
						RANLIB='$(RANLIB)')
	@echo "==== SQTRACEMIN LIBRARY archives up to date ===="

#Clean Everything
clean: clean_pctracemin clean_pqtracemin clean_sctracemin clean_sqtracemin 

clean_tracemin:
	($(CD) $(SRCDIR); $(MAKE) clean_tracemin CD='$(CD)' SRCDIR='../$(SRCDIR)' OBJDIR='../$(OBJDIR)' LIBDIR='../$(LIBDIR)' BINDIR='../$(BINDIR)' \
                                                PCTRACEMINLIB='$(PTRACEMINLIB)')
	@echo "==== TRACEMIN cleaned ===="

clean_pctracemin:
	($(CD) $(SRCDIR); $(MAKE) clean_pctracemin CD='$(CD)' SRCDIR='../$(SRCDIR)' OBJDIR='../$(OBJDIR)' LIBDIR='../$(LIBDIR)' BINDIR='../$(BINDIR)' \
                                                PCTRACEMINLIB='$(PTRACEMINLIB)')
	@echo "==== PCTRACEMIN cleaned ===="

clean_pqtracemin:
	($(CD) $(SRCDIR); $(MAKE) clean_pqtracemin CD='$(CD)' SRCDIR='../$(SRCDIR)' OBJDIR='../$(OBJDIR)' LIBDIR='../$(LIBDIR)' BINDIR='../$(BINDIR)' \
                                                PQTRACEMINLIB='$(PTRACEMINLIB)')
	@echo "==== PQTRACEMIN cleaned ===="

clean_sctracemin:
	($(CD) $(SRCDIR); $(MAKE) clean_sctracemin CD='$(CD)' SRCDIR='../$(SRCDIR)' OBJDIR='../$(OBJDIR)' LIBDIR='../$(LIBDIR)' BINDIR='../$(BINDIR)' \
                                                SCTRACEMINLIB='$(STRACEMINLIB)')
	@echo "==== SCTRACEMIN cleaned ===="

clean_sqtracemin:
	($(CD) $(SRCDIR); $(MAKE) clean_sqtracemin  CD='$(CD)' SRCDIR='../$(SRCDIR)' OBJDIR='../$(OBJDIR)' LIBDIR='../$(LIBDIR)' BINDIR='../$(BINDIR)' \
                                                SQTRACEMINLIB='$(STRACEMINLIB)')
	@echo "==== SQTRACEMIN cleaned ===="

clean_manually:
	rm $(OBJDIR)/*.o $(OBJDIR)/NEMO3DWrappers/*.o $(OBJDIR)/*.mod $(LIBDIR)/*.a $(BINDIR)/pctracemin $(BINDIR)/pqtracemin \
    	   $(BINDIR)/sctracemin $(BINDIR)/sqtracemin   

