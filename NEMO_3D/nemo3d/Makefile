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
# $Header: /repo/nemo3d/Makefile,v 1.14 2006/10/17 14:10:14 hoonryu Exp $
# *****************************************************************************

#================================================================
# NEMO-3D Makefile
#================================================================


#----------------------------------------------------------------
# System configuration: compiler and libraries
#----------------------------------------------------------------

include ../build/make.inc

#----------------------------------------------------------------
# Makefile includes
#----------------------------------------------------------------

NEMO3D_TOPMK   = $(NEMO3D_ROOT)/Makefile
NEMO3D_FILMK   = $(NEMO3D_ROOT)/NEMO3D_files.mk
NEMO3D_DEPMK   = $(NEMO3D_ROOT)/NEMO3D_dependencies.mk

# sets variables containing file information
include $(NEMO3D_FILMK)

#----------------------------------------------------------------
# user rules
#----------------------------------------------------------------
nemo3d : $(NEMO3D_EXEC)
	@ echo "   NEMO-3D up to date"

depend : $(NEMO3D_DEPMK)
	@ echo "   " $(NEMO3D_DEPMK) "is now up-to-date"


# Nemo3D library components

base : $(NEMO3D_LIB_BASE)
	@ echo "   NEMO3D::base library up to date"

ham : $(NEMO3D_LIB_HAM)
	@ echo "   NEMO3D::ham library up to date"

io : $(NEMO3D_LIB_IO)
	@ echo "   NEMO3D::io library up to date"

math : $(NEMO3D_LIB_MATH)
	@ echo "   NEMO3D::math library up to date"

top : $(NEMO3D_LIB_TOP)
	@ echo "   NEMO3D::top library up to date"

util : $(NEMO3D_LIB_UTIL)
	@ echo "   NEMO3D::util library up to date"



# Cleanup

# cleanup executables
cleane :
	rm -rf $(NEMO3D_EXEC)
	rm -rf $(NEMO3D_EXEC_LNK)
	@ echo "   All executable files deleted"

# cleanup libraries
cleanl :
	rm -rf $(NEMO3D_LIBD)
	@ echo "   All library archive files deleted"

# cleanup object files
cleano :
	rm -rf $(NEMO3D_OBJD)
	@ echo "   All object files deleted"

# cleanup dependencies
cleand :
	rm -f $(NEMO3D_DEPMK)
	@ echo "   NEMO-3D dependency file deleted"

clean : cleano cleand
cleanall : cleane cleano cleanl cleand

#----------------------------------------------------------------
# local rules
#----------------------------------------------------------------


$(NEMO3D_DEPMK) : $(NEMO3D_HDRS)
	@ echo "===== MAKING DEPENDENCIES ====="
	rm -f $(NEMO3D_DEPMK)
	$(CCMKD) $(NEMO3D_SRCH_INCLD) $(NEMO3D_SRCS) > $(NEMO3D_DEPMK)-temp
	sed "s?^.*\.o:?$(NEMO3D_OBJD)/&?" $(NEMO3D_DEPMK)-temp > $(NEMO3D_DEPMK)
	rm $(NEMO3D_DEPMK)-temp

$(NEMO3D_EXEC) : $(NEMO3D_OBJS_TOP) $(NEMO3D_LIBS) $(NEMO3D_EXTERNLIBS)
	$(LD) $(LDOPT)  -o $(@) $(NEMO3D_OBJS_TOP) $(LIB_PATHS) $(LIBS)
	ln -s -f $(NEMO3D_EXEC) $(NEMO3D_EXEC_LNK)
	@ echo "===== CREATED nemo3d EXECUTABLE ====="

# Nemo3D library components

$(NEMO3D_LIB_BASE) : $(NEMO3D_OBJS_BASE)
	@ mkdir -p $(NEMO3D_LIBD)
	$(AR) -$(AR_OPTIONS) $(@) $(NEMO3D_OBJS_BASE)
	@ echo "===== CREATED nemo3d::base LIBRARY ====="

$(NEMO3D_LIB_HAM) : $(NEMO3D_OBJS_HAM)
	@ mkdir -p $(NEMO3D_LIBD)
	$(AR) -$(AR_OPTIONS) $(@) $(NEMO3D_OBJS_HAM)
	@ echo "===== CREATED nemo3d::ham LIBRARY ====="

$(NEMO3D_LIB_POST) : $(NEMO3D_OBJS_POST)
	@ mkdir -p $(NEMO3D_LIBD)
	$(AR) -$(AR_OPTIONS) $(@) $(NEMO3D_OBJS_POST)
	@ echo "===== CREATED nemo3d::postprocessing LIBRARY ====="

$(NEMO3D_LIB_IO) : $(NEMO3D_OBJS_IO)
	@ mkdir -p $(NEMO3D_LIBD)
	$(AR) -$(AR_OPTIONS) $(@) $(NEMO3D_OBJS_IO)
	@ echo "===== CREATED nemo3d::io LIBRARY ====="

$(NEMO3D_LIB_MATH) : $(NEMO3D_OBJS_MATH)
	@ mkdir -p $(NEMO3D_LIBD)
	$(AR) -$(AR_OPTIONS) $(@) $(NEMO3D_OBJS_MATH)
	@ echo "===== CREATED nemo3d::math LIBRARY ====="

$(NEMO3D_LIB_TOP) : $(NEMO3D_OBJS_TOP)
	@ mkdir -p $(NEMO3D_LIBD)
	$(AR) -$(AR_OPTIONS) $(@) $(NEMO3D_OBJS_TOP)
	@ echo "===== CREATED nemo3d::top LIBRARY ====="

$(NEMO3D_LIB_UTIL) : $(NEMO3D_OBJS_UTIL)
	@ mkdir -p $(NEMO3D_LIBD)
	$(AR) -$(AR_OPTIONS) $(@) $(NEMO3D_OBJS_UTIL)
	@ echo "===== CREATED nemo3d::util LIBRARY ====="


# Nemo3D object componenets

$(NEMO3D_OBJD)/%.o: $(NEMO3D_SRCD_TOP)/%.c
	@ mkdir -p $(NEMO3D_OBJD)
	$(CC) $(CCOPT) $(NEMO3D_SRCH_INCLD) $(<) -o $(@)

$(NEMO3D_OBJD)/%.o: $(NEMO3D_SRCD_BASE)/%.c
	@ mkdir -p $(NEMO3D_OBJD)
	$(CC) $(CCOPT) $(NEMO3D_SRCH_INCLD) $(<) -o $(@)

$(NEMO3D_OBJD)/%.o: $(NEMO3D_SRCD_IO)/%.c
	@ mkdir -p $(NEMO3D_OBJD)
	$(CC) $(CCOPT) $(NEMO3D_SRCH_INCLD) $(<) -o $(@)

$(NEMO3D_OBJD)/%.o: $(NEMO3D_SRCD_MATH)/%.c
	@ mkdir -p $(NEMO3D_OBJD)
	$(CC) $(CCOPT) $(NEMO3D_SRCH_INCLD) $(<) -o $(@)

$(NEMO3D_OBJD)/%.o: $(NEMO3D_SRCD_UTIL)/%.c
	@ mkdir -p $(NEMO3D_OBJD)
	$(CC) $(CCOPT) $(NEMO3D_SRCH_INCLD) $(<) -o $(@)

$(NEMO3D_OBJD)/%.o: $(NEMO3D_SRCD_HAM)/%.c
	@ mkdir -p $(NEMO3D_OBJD)
	$(CC) $(CCOPT) $(NEMO3D_SRCH_INCLD) $(<) -o $(@)

$(NEMO3D_OBJD)/%.o: $(NEMO3D_SRCD_POST)/%.cpp
	@ mkdir -p $(NEMO3D_OBJD)
	$(CC) $(CCOPT) $(NEMO3D_SRCH_INCLD) $(<) -o $(@)


#----------------------------------------------------------------
# dependency rules
#----------------------------------------------------------------

include $(NEMO3D_DEPMK)

