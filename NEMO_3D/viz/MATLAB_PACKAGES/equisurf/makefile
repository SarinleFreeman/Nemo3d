# Makefile for SUN Solaris 2.x

# You probably have to modify these four lines
CC = gcc
CCOPT = -O3
INCLUDES = /usr/local/matlab/extern/include
SUFFIX = mexsol 

CFLAGS = $(CCOPT) -I$(INCLUDES) -DMATLAB_MEX_FILE

OBJECTS1 = eyta.o konf.o eyta_matlab.o

all:	polygons.$(SUFFIX) interp.$(SUFFIX)

polygons.$(SUFFIX):	$(OBJECTS1)
	mex $(OBJECTS1)
	mv eyta.$(SUFFIX) polygons.$(SUFFIX)

eyta.o:		eyta.h
eyta_matlab.o:	eyta.h

OBJECTS2 = interp.o

interp.$(SUFFIX):	$(OBJECTS2)
	mex $(OBJECTS2)

interp.o:	eyta.h

clean:
	rm -f *.o *~
