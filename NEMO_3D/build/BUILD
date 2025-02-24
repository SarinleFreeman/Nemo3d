
                          Installing NEMO3D: BUILD



Build Instructions
==================


Decompressing the tarball
-------------------------

The NEMO3D application is available as a compressed tape archive NEMO3D.tar.gz.

Copy NEMO3D.tar.gz into your home directory and type

     gzip -dc NEMO3D.tar.gz | tar -xvf -

to decompress and untar the distribution. A "NEMO3D/" directory
containing the NEMO3D application and its dependent packages will be created.



Build Tools
-----------

NEMO3D has been developed using freely available GNU tools to ensure maximum
cross-platform compatibility. Please make sure that you have fairly recent
versions of GNU Make and GNU GCC. NEMO3D has been test built using various
versions of GNU Make including 3.79.1 and 8.0. NEMO3D also compiles with 
GNU GCC versions 2.96, 3.04, 3.1 and 3.2.

NEMO3D has successfully been compiled on RedHat 7.3, RedHat 8.0, IRIX64 6.5,
and MAC OS X 10.2 (Jaguar-Darwin).



Required Environmental Variables for the build
----------------------------------------------

Before building the executable, two things must be done.  One, appropriate
paths must be set in the "NEMO3D/build/make.inc" configuration file.  This
is described later. Two, the following environmental variables must be set:

HOSTTYPE: This variable is usually already set by the system during login 
          and is used in the "NEMO3D/build/make.inc" file to set the variable 
          PLATFORM = $(HOSTTYPE)_$(BUILD_TYPE).  An known issue is that
          different shells may, by default, define HOSTTYPE differently, so
          that a compilation under bash (for example) may not recognize
          a previous compilation under tcsh.

XERCESCROOT: Root path of where xerces-c-src2_1_0 is installed.
             Since xerces-c-src2_1_0 is packaged with NEMO3D, this variable
             should be set to <NEMO3D-root-path>/xerces-c-src2_1_0, unless
             the user desires to use his own version.  However, other versions
             of Xerces (in particular pre-2.0 versions and more recent
             versions that have removed depracated functions) have been known 
             to fail.



Optional Environmental Variables
--------------------------------

For parallel execution, NEMO3D uses coarse-grained MPI message passing
calls.  However, it has been tested only with the MPICH implementation 
of MPI-1.  This freely available software is available at

     http://www-unix.mcs.anl.gov/mpi/mpich

If you wish to compile NEMO3D with MPICH, you need to set the MPI_HOME 
environmental variable.  In addition to gcc NEMO3D has been tested with
two commercial compilers, Absoft and Intel.  If the Absoft compiler is 
installed, then you should also set the ABSOFT variable.  Support for
additional compiler types must be added by the user within the
"NEMO3D/build/make.inc" file (see below).


MPI_HOME: Root of where MPI is installed.
          e.g. "/home/packages/mpich-1.2.5/absoft"

ABSOFT  : Root of where a particular compiler type is installed.
          e.g. "/home/packages/absoft"



Editing "NEMO3D/build/make.inc"
-------------------------------

Before building the executable, you must set certain library paths 
specific to your system.

The make include file "NEMO3D/build/make.inc" contains the master configuration
settings for compiling the NEMO3D application.  This file contains settings for
the subcomponents of NEMO3D including LAPACK/BLAS, NML, eigen, and XMLtoC. 
These settings are recursively passed onto and override the settings of the
corresponding Makefiles of the dependent components.  Thus, one need only 
modify "make.inc".

The dependency on Python for XML data interpretion requires
that the Python libraries and export include paths be set to your
system's Python installation environment.  Since Python is typically installed
by default on gnu/Linux systems, it has not been included in the NEMO3D 
distribution. The following 3 Python variables must be correctly set
in the "NEMO3D/build/make.inc" file.


PYTHON_ROOT        : The root directory where Python is installed.

PYTHON_EXPORT_INCL : The include directory for Python where Python.h and
		     other headers can be found.

LIB_PYTHON         : The path to the Python library archive.


Also, although the environment variable XERCESCROOT must be set in order
to build Xerces as part of the NEMO3D distribution, the user may prefer 
to use a version different from that distributed by NEMO3D.  (CAVEAT:  
because the use of XML-based inputs was implemented during a period while 
Xerces was undergoing major modifications, using a different version of 
Xerces from that which is distriuted with NEMO3D may fail.)  To do this, 
the user must also modify XERCESCROOT within the make.inc file.



Build Types
-----------

The "NEMO3D/build/make.inc" file is organized so as to support multiple
build type configurations. "make.inc" includes several example compiler build 
types designated by the variable BUILD_TYPE, which constists of a head 
string denoting the compiler and/or platform followed by various suffixes:
"_mpi" denotes compilation with MPI support; "_fast"/"_debug" means 
with/without optimizations; and "_static" refers to compilation without
reliance on shared libraries.  Possible head strings include:
gnu, absoft, intel_piii, and macosx.  For example,

BUILD_TYPE=absoft_mpi_fast_static

denotes a static build using the Absoft compiler (with optimization)
linked with an MPI library.  Note that command-line setting of BUILD_TYPE
overrides what is specified in the "make.inc" file (see below).



Compiling the Executable
------------------------

Having set the required environmental variables and correctly configured 
"NEMO3D/build/make.inc", one can then build the executable. From the 
"NEMO3D/" directory, type:

	cd build
	make

This make will run the main Makefile which recurses down into each 
component's Makefile to run.  The "make.inc" is included to load the 
master confugration for the subcomponents.

Optionally, you can pass in the build type when calling make.  For example:

	make BUILD_TYPE='gnu_mpi_fast'

After a successful build of the executable, you should see the following
output:
	...
	...

	===== CREATED nemo3d EXECUTABLE =====
	   NEMO-3D up to date


The executable is located in the nemo3d component's bin directory.

	cd ../nemo3d/bin
	ls -l

You should see 2 new files.

1) "nemo3d.ex" is a symbolic link to the newest build of the executable.

2) The second file is the actual executable named with the following
   convention: "nemo3d-<HOSTTYPE>_<BUILD_TYPE>.ex"
   So if the host computer is a "i386-linux" and the compiler build type
   defined in "NEMO3D/build/make.inc" was called "absoft_mpi_fast",
   then the executable will be "nemo3d-i386-linux_absoft_mpi_fast.ex"



Running NEMO3D
==============

Upon successful compilation of NEMO3D, the executable will be placed
in "NEMO3D/nemo3d/bin/nemo3d-<HOSTTYPE>_<BUILD_TYPE>.ex". A soft link will 
also be created at "NEMO3D/nemo3d/bin/nemo3d.ex" which points to the current 
build's executable.  Prior to running the executable, one must set the 
environment variable PYTHONPATH, the path containing the "funcs.py" file 
(which must exist but may be empty).  The default location is NEMO3D's 
"nemo3d/examples" subdirectory, so that if NEMO3D is installed
in the root of your home directory, PYTHONPATH should point to

            "$HOME/NEMO3D/nemo3d/example"

To run the executable, first create a directory and copy the sample
XML input file to it:

        mkdir $HOME/first_run
        cd $HOME/first_run
        cp $HOME/XMLtoC/xml/ioNemo3D.xml .

If you have compiled a non-MPI version of NEMO3D (e.g. BUILD_TYPE=gnu_fast),
you can easily run the application as follows:

	$HOME/NEMO3D/nemo3d/bin/nemo3d.ex ioNemo3D.xml

If you have compiled the MPICH-enabled version and have 'mpirun' in your path,
you can run the executable on (for example) 2 processors as follows:

	mpirun -np 2 $HOME/NEMO3D/nemo3d/bin/nemo3d.ex ioNemo3D.xml

To specify the subset of nodes that mpirun will use by using a text file 
listing the allowed nodes, run

	mpirun -machinefile mpiMachineFile.txt -np 2 $HOME/NEMO3D/nemo3d/bin/nemo3d.ex ioNemo3D.xml

where mpiMachineFile.txt is located in the current directory and contains the
host names of the allowed nodes to run (at one node name per line).
Please read the mpich documentation for additional options.

For examples on how to run NEMO3D, see the README file in nemo3d/doc.


What does "NEMO3D/build/Makefile" do?
-------------------------------------

Basically, when you run 'make' in the "NEMO3D/build/" directory, it builds
all of the dependent libraries first, then builds the core of NEMO3D itself.
Of course, if the library is already built, it moves on to check and build
the next library. Here is a time-ordered list of what "NEMO3D/build/Makefile" 
does:

	1) build BLAS library
	2) build LAPACK library
	4) build ARPACK library
	5) build PARPACK library
	6) build nml library
	7) build eigen library
	8) build Xerces-C++ library
		a) run Xerces's `./configure'
		b) run Xerces's `make'
	9) build XMLtoC library
		a) build XMLtoC.ex executable
		b) run XMLtoC.ex to generate C++ files from XML file
		c) build the XMLtoC I/O library using the generated files
                   specific to NEMO3D: "ioNemo3D"
	10) build NEMO3D executable

Once the libraries and files for these packages have been built, the Makefile
will begin to build the NEMO3D core and then link to make the executable.




Compiling Different Coexisting BUILT_TYPEs
==========================================

The components developed for NEMO3D can handle simultaneous cross platform 
build types.  Each build type is built in a separate build directory that 
contains the built type and platform name. These components can have 
several BUILD_TYPEs coexisting together without conflicts.

So the following components support *simultaneous* cross platoform 
compilation:
   - nml
   - eigen
   - XMLtoC
   - nemo3d

However, NEMO3D also uses 3rd party components that cannot handle 
simultaneous cross platform builds.  Each time a new platform is built, 
their object files are placed in the same location and will conflict in 
namespace with newer BUILT_TYPEs.  However, the library names do not 
conflict since they contain platform built type info in the name to avoid 
namespace collisions.

The following components do NOT support simultaneous cross platoform 
compilation:
   - BLAS
   - LAPACK
   - ARPACK
   - PARPACK
   - Xerces-C++

For these, the object files must be removed before attempting to build 
another BUILD_TYPE profile.

Thus, before trying to build a different BUILD_TYPE (while maintaining an 
existing BUILD_TYPE), one must remove the objects of the non-cross platform 
compatible object files. Their libraries can be left alone since they don't 
conflict in namespace.

The following must be run from "NEMO3D/build" directory before 
attempting to build a different BUILD_TYPE:

   make coexisting_build_type

This clears all the object files of the components that cannot coexist 
with other BUILD_TYPEs.



TROUBLE SHOOTING
================

NEMO3D has been developed using GNU tools and should be built with GNU tools.
Compiling NEMO3D with non-GNU versions of 'make' and 'g++' is known to cause 
problems.




Setting my own GCC in "NEMO3D/build/make.inc" that is different from my
system default results in Xerces undefined reference link errors
-----------------------------------------------------------------------

The settings for the CC variable in "NEMO3D/build/make.inc" are passed onto
the recursive Makefiles for LAPACK, nml, eigen, XMLtoC, and nemo. However,
for xerces-c-src2_1_0, it uses it's own settings based from running
	"NEMO3D/xerces-c-src2_1_0/src/xercesc/configure"
which probes your system. So no matter what GCC you set in the
"NEMO3D/build/make.inc" file, Xerces's configure script may detect another
GCC installed on your system (most like the one in your path. i.e. `which g++')
and generate all it's Makefiles to use that GCC.

A work-around is to edit "NEMO3D/xerces-c-src2_1_0/src/xercesc/config.cache"
and replace instances of "g++" and "c++" with absolute paths to the version
you want to use.  Then rerun `./configure' and it will generate all the 
Makefiles using the GCC of your choice.  Finally, run `make' manually 
to build the xerces-c shared objct library.  Now when you go back to  
"NEMO3D/build/" to run `make', it will skip building Xerces-C++ since the 
xerces library has just been built.

The bottom line is that the GCC that is used to compile Xerces-C++ should
be the same as the GCC that compiles the rest of the NEMO3D application.
This way, you will have minimized any GCC incompatibilities with 
different GCC versions and formats.




When the Makefile reaches "making the XMLtoC executable", it stops at:
ld: unknown flag: -rpath
----------------------------------------------------------------------

On some older or non-GCC link editors, it may not recognize the flag `-rpath'.
The LD option flag `-Wl,-rpath,/some/shared_lib/path/' passes the `-rpath' 
parameter to the link editor so it compiles the executable with paths to 
where to search for shared objects during run-time.

For compatibility with other ELF linkers, you may also try the `-R'  option
instead of the `-rpath' option. If `-R' is followed  by a directory name,
rather than a file name, the link editor will treated `-R' as the 
`-rpath' option.

If your link editor does not support `-rpath', use the LD_LIBRARY_PATH
environmental variable instead.  Although using `-Wl,-rpath,' is prefered, 
and using the LD_LIBRARY_PATH method is BAD, it is a viable
option if it gets to this point.

Simply remove the `-Wl,-rpath,' flags from the Makefiles, and then
define LD_LIBRARY_PATH environmental variable with the same paths as
the one defined in -rpath. This way, the run-time of the executable
will not have any problems dynamically linking in shared object libraries.

So to get XMLtoC.ex to compile, simply add the path "$(XERCESCROOT)/lib/"
to LD_LIBRARY_PATH environment variable.




When running "NEMO3D/XMLtoC/bin/XMLtoC.ex" to generate C++ code from XML,
it says it cannot find "libcerces-c.so.21"?
-------------------------------------------------------------------------

* This may occur if you have removed the default `-Wl,-rpath,' flags
  that are passed to the link editor in the Makefiles. *

In our testing, the newer GCC 3.2 will compile XMLtoC.ex and link with no
problems.  However during the runtime of the XML-to-C part of Makefile, the
"NEMO3D/XMLtoC/bin/XMLtoC.ex" gets executed and the loader might complain
that it cannot find "libcerces-c.so.21".

This would not have occured if the default `-Wl,-rpath,' flag is used
to passed into the link editor.  However, if for some reason your link
editor does not recognize `-rpath', then the alternative is to simply add
the path "$(XERCESCROOT)/lib/" to runtime library path environment variable.

In Linux, the runtime library path environment variable is "LD_LIBRARY_PATH".
In MACOSX, it is "DYLD_LIBRARY_PATH". For IRIX, it could be LD_LIBRARY_PATH,
LD_LIBRARYN32_PATH, or LD_LIBRARYN64_PATH.

Then during runtime, XMLtoC.ex will look in "$(XERCESCROOT)/lib/" to find
the Xerces-C++ shared object library.

If you are on a Solaris system, yet another alternative is to add the
switch "-R $XERCESCROOT/lib" to the link editor (instead of `-Wl,-rpath,')
for compiling the XMLtoC.ex file. This will create the executable that will
know to look in "$XERCESCROOT/lib" for the "libxerces-c.so.21" DURING 
runtime.  Unfortunately, the -R switch is not recognized by GNU GCC. 
GNU GCC does recognize the `-Wl' flag to pass the following `-rpath' 
to the link editor.



Building Xerces-C++ from the Mac OS X command line gives errors?
----------------------------------------------------------------

(NOTE: The "NEMO3D/build/Makefile" already fixes this case)

If you are trying to build NEMO3D in a Mac OS X command line environment,
please note that Xerces-C++'s make process will stop and complain that
you need a message loader and a native transcoder.

	MacOSPlatformUtils.cpp:1033:6: #error You must provide a message loader
	MacOSPlatformUtils.cpp:1052:6: #error You must provide a transcoding service implementation

The remedy is to run generate the Xerces Makefile with special flags that
compile the code with:

	cd NEMO3D/xerces-c-src2_1_0/src/xercesc
	./runConfigure -p macosx -n native
	make

	( Instructions taken from:
	  http://xml.apache.org/xerces-c/build-other.html#Mac )

The new Makefile will include new defines to force Xerces to use native
transcoders and in-memory message loader. You may notice the following new
defines added:
	-DXML_USE_NATIVE_TRANSCODER -DXML_USE_INMEM_MESSAGELOADER -DXML_USE_PTHREADS -DXML_USE_NETACCESSOR_NATIVE

Now that the Xerces library has been created, go back to "NEMO3D/build"
and run "make" again. Make will check the Xerces library dependency and
will continue onto making NEMO3D since the Xerces library has already been
created.




Compiling Xerces-C++ v2.1.0 on IRIX platform fails at IRIXPlatformUtils.cpp 
with undefined references to openFileToWrite() and writeBufferToFile()?
----------------------------------------------------------------------------

(NOTE: The "NEMO3D/build/Makefile" already fixes this case)

The version 2.1.0 of Xerces-C++ contains bugs and incomplete implementation 
for file read/write operations on IRIX platforms. If you try to build 
Xerces-C++ v2.1.0 on IRIX, you may see error messages about incomplete 
functions openFileToWrite() and writeBufferToFile().

This is a known bug acknowledged in the bug list for Xerces:
http://nagoya.apache.org/bugzilla/show_bug.cgi?id=14598

Until a newer version of Xerces-C++ is release (newer than version 2.1.0), 
please replace the original file:
"NEMO3D/xerces-c-src2_1_0/src/xercesc/util/Platforms/IRIX/IRIXPlatformUtils.cpp"
with a patched version distributed along with this NEMO3D application.
The patch file is located at:
"NEMO3D/nemo3d/misc/IRIXPlatformUtils.cpp".




Building nemo3d.ex with build type macosx_absoft_mpi_fast fails.  The
libraries were created successfully but the linking stage failed.
The error message is
   ld: Undefined symbols: _PZNAPPSZSORTC
------------------------------------------------------------------------

This is an issue between the Fortran compiler on the MACOSX and the source
code of PARPACK. PARPACK's "PARPACK/SRC/MPI/pznaup2.f" has unusually long
line spacing on lines 253-254. 
         
         external   zcopy,pzgetv0,pznaitr,pzneigh,pzngets,pznapps,
         &          zsortc,zswap,pzmout,pzvout,pivout,second

There is a bug with the Fortran compiler handling the number of total chars 
per line.  The compiler will think think "pznapps" and "zsortc" are really 
one variable "pznappszsortc". We have provided a patch in build/patches 
that shortens the line spacing.

===========================================================
