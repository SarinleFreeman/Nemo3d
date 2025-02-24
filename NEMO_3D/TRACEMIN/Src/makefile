#SL Library
SORTOBJFILES       = sort2.o
PCWrappersOBJFILES = NEMO3DWrappers/p_map_chebyshev_shift_wrapper.o NEMO3DWrappers/p_zgtracemin_wrapper.o NEMO3DWrappers/p_extract_eigenpairs_wrapper.o 
SCWrappersOBJFILES = NEMO3DWrappers/s_map_chebyshev_shift_wrapper.o NEMO3DWrappers/s_zgtracemin_wrapper.o NEMO3DWrappers/s_extract_eigenpairs_wrapper.o
PQWrappersOBJFILES = NEMO3DWrappers/p_modified_conjugate_gradient_wrapper.o NEMO3DWrappers/p_ztracemin_wrapper.o
SQWrappersOBJFILES = NEMO3DWrappers/s_modified_conjugate_gradient_wrapper.o NEMO3DWrappers/s_ztracemin_wrapper.o

#Objects Fortran
TRACEMINTYPES = tracemin_types.o $(SORTOBJFILES)

PCTRACEMINOBJ_core      = pctracemin_functions.o
PCTRACEMINOBJ_example   = pctracemin.o
PCTRACEMINOBJ_library   = $(TRACEMINTYPES) $(PCTRACEMINOBJ_core) $(PCWrappersOBJFILES)
PCTRACEMINOBJ_standalone= $(TRACEMINTYPES) $(PCTRACEMINOBJ_core) $(PCTRACEMINOBJ_example)

PQTRACEMINOBJ_core      = pqtracemin_functions.o 
PQTRACEMINOBJ_example   = pqtracemin.o
PQTRACEMINOBJ_library   = $(TRACEMINTYPES) $(PQTRACEMINOBJ_core) $(PQWrappersOBJFILES)
PQTRACEMINOBJ_standalone= $(TRACEMINTYPES) $(PQTRACEMINOBJ_core) $(PQTRACEMINOBJ_example)

SCTRACEMINOBJ_core      = sctracemin_functions.o 
SCTRACEMINOBJ_example   = sctracemin.o
SCTRACEMINOBJ_library   = $(TRACEMINTYPES) $(SCTRACEMINOBJ_core) $(SCWrappersOBJFILES)
SCTRACEMINOBJ_standalone= $(TRACEMINTYPES) $(SCTRACEMINOBJ_core) $(SCTRACEMINOBJ_example)

SQTRACEMINOBJ_core      = sqtracemin_functions.o 
SQTRACEMINOBJ_example   = sqtracemin.o
SQTRACEMINOBJ_library   = $(TRACEMINTYPES) $(SQTRACEMINOBJ_core) $(SQWrappersOBJFILES)
SQTRACEMINOBJ_standalone= $(TRACEMINTYPES) $(SQTRACEMINOBJ_core) $(SQTRACEMINOBJ_example)

ALLTRACEMINOBJ_library   = $(TRACEMINTYPES) $(PCTRACEMINOBJ_core) $(PCWrappersOBJFILES) $(PQTRACEMINOBJ_core) $(PQWrappersOBJFILES) \
			   $(SCTRACEMINOBJ_core) $(SQTRACEMINOBJ_core)	
ALLTRACEMINOBJ_standalone= $(ALLTRACEMINOBJ_library) $(PCTRACEMINOBJ_example) $(PQTRACEMINOBJ_example) $(SCTRACEMINOBJ_example) $(SQTRACEMINOBJ_example)
   

#Standalone
all: $(ALLTRACEMINOBJ_standalone)
	$(FC) $(FFLAGS) $(ALLTRACEMINOBJ_standalone) $(LIBPATH) $(LIB) -o $(NAME)

pctracemin: $(PCTRACEMINOBJ_standalone)
	($(CD) $(OBJDIR); $(FC) $(FFLAGS) $(PCTRACEMINOBJ_standalone) $(LIBPATH) $(LIB) -o $(BINDIR)/$(PCTRACEMINNAME))

pqtracemin: $(PQTRACEMINOBJ_standalone)
	($(CD) $(OBJDIR); $(FC) $(FFLAGS) $(PQTRACEMINOBJ_standalone) $(LIBPATH) $(LIB) -o $(BINDIR)/$(PQTRACEMINNAME))

sctracemin: $(SCTRACEMINOBJ_standalone)
	($(CD) $(OBJDIR); $(FC) $(FFLAGS) $(SCTRACEMINOBJ_standalone) $(LIBPATH) $(LIB) -o $(BINDIR)/$(SCTRACEMINNAME))

sqtracemin: $(SQTRACEMINOBJ_standalone)
	($(CD) $(OBJDIR); $(FC) $(FFLAGS) $(SQTRACEMINOBJ_standalone) $(LIBPATH) $(LIB) -o $(BINDIR)/$(SQTRACEMINNAME))


#Libraries
lib: $(ALLTRACEMINOBJ_library)
	$(AR) $(ARFLAGS) $(TRACEMINLIB) $(ALLTRACEMINOBJ_library)
	$(RANLIB) $(TRACEMINLIB)
	@echo "Created a library for Trace Minimization Algorithm"

pctracemin_lib: $(PCTRACEMINOBJ_library)
	($(CD) $(OBJDIR); $(AR) $(ARFLAGS) $(PCTRACEMINLIB) $(PCTRACEMINOBJ_library); \
	$(RANLIB) $(PCTRACEMINLIB)) 
	@echo "Created a library for Parallel Trace Minimization Algorithm using Chebyshev Mapping"

pqtracemin_lib: $(PQTRACEMINOBJ_library)
	($(CD) $(OBJDIR); $(AR) $(ARFLAGS) $(PQTRACEMINLIB) $(PQTRACEMINOBJ_library); \
	$(RANLIB) $(PQTRACEMINLIB)) 
	@echo "Created a library for Parallel Trace Minimization Algorithm using Quadratic Mapping"	

sctracemin_lib: $(SCTRACEMINOBJ_library)
	($(CD) $(OBJDIR); $(AR) $(ARFLAGS) $(SCTRACEMINLIB) $(SCTRACEMINOBJ_library); \
	$(RANLIB) $(SCTRACEMINLIB)) 
	@echo "Created a library for Sequential Trace Minimization Algorithm using Chebyshev Mapping"

sqtracemin_lib: $(SQTRACEMINOBJ_library)
	($(CD) $(OBJDIR); $(AR) $(ARFLAGS) $(SQTRACEMINLIB) $(SQTRACEMINOBJ_library); \
	$(RANLIB) $(SQTRACEMINLIB))
	@echo "Created a library for Sequential Trace Minimization Algorithm using Quadratic Mapping"

#Cleaning
clean_tracemin:
	($(CD) $(OBJDIR); rm -f *.o; rm -f *.mod; \
	 $(CD) $(LIBDIR); rm -f *.a; \
	 $(CD) $(OBJDIR)/NEMO3DWrappers; rm -f *.o;)  	

clean_pctracemin:
	($(CD) $(OBJDIR); rm -f $(PCTRACEMINOBJ_library) *.mod; \
	 $(CD) $(LIBDIR); rm -f $(PCTRACEMINLIB))  	

clean_pqtracemin:
	($(CD) $(OBJDIR); rm -f $(PQTRACEMINOBJ_library) *.mod; \
	 $(CD) $(LIBDIR); rm -f $(PQTRACEMINLIB))

clean_sctracemin:
	($(CD) $(OBJDIR); rm -f $(SCTRACEMINOBJ_library) *.mod; \
	 $(CD) $(LIBDIR); rm -f $(SCTRACEMINLIB))

clean_sqtracemin:
	($(CD) $(OBJDIR); rm -f $(SQTRACEMINOBJ_library) *.mod; \
	 $(CD) $(LIBDIR); rm -f $(SQTRACEMINLIB))

#Compilation Rules
.SUFFIXES: .f .f90 
.f90.o:
	$(FC) $(FFLAGS) -c -o $(OBJDIR)/$*.o $*.f90

.f.o:
	$(FC) $(FFLAGS) -c -o $(OBJDIR)/$*.o $*.f

















