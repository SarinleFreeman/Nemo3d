# *****************************************************************************
# The Jet Propulsion Laboratory (JPL) XML-to-C++ package.
# Copyright (C) 2002 California Institute of Technology (Caltech)
#
# This library is free software, which you can redistribute and/or modify
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
#   Fabiano Oyafuso (fabiano@jpl.nasa.gov)
#
# Written by:  Fabiano Oyafuso
#              Hook Hua
#              Ed Vinyard
#
# This product includes software developed by the Apache Software Foundation
# (http://www.apache.org/).
#
# *****************************************************************************
# $Header: /repo/XMLtoC/application.mk,v 1.4 2003/10/08 16:15:18 hook Exp $
# *****************************************************************************


#================================================================
# Makefile to generate the XMLtoC application library
#================================================================

ifeq ($(PJCT),)
	$(error "Undefined environment variable, $(PJCT)")
endif

ifeq ($($(PJCT)_ROOT),)
	$(error "Undefined environment variable, $(PJCT_ROOT)")
endif

ifeq ($($(PJCT)_BLD),)
	$(error "Undefined environment variable, $(PJCT_BLD)")
endif

# sets directory structure and default compiler options
include ./global.mk


$(PJCT)_HDRS = $(wildcard $($(PJCT)_SRCD)/*.H) $(XMLC_ROOT)/src/wiglafNode.H
$(PJCT)_SRCS = $(wildcard $($(PJCT)_SRCD)/*.C)
$(PJCT)_SRCS4LIB = $(patsubst $($(PJCT)_SRCD)/main.C,,$($(PJCT)_SRCS))

$(PJCT)_OBJS = \
	$(patsubst $($(PJCT)_SRCD)/%.C,$($(PJCT)_OBJD)/%.o,$($(PJCT)_SRCS))
$(PJCT)_OBJS4LIB = $(patsubst $($(PJCT)_OBJD)/main.o,,$($(PJCT)_OBJS))

$(PJCT)_LIB = $($(PJCT)_LIBD)/lib$(PJCT)-$(PLATFORM).a
$(PJCT)_LIB_LNK = $($(PJCT)_LIBD)/lib$(PJCT).a

$(PJCT)_EXEC      = $($(PJCT)_BIND)/$(PJCT)-driver-$(PLATFORM).ex
$(PJCT)_EXEC_LNK  = $($(PJCT)_BIND)/$(PJCT)-driver.ex

#----------------------------------------------------------------
#  Platform specific options
#----------------------------------------------------------------


# Search paths to shared objects for compiler and link editor
LIB_PATHS = -L$(XERCESCROOT)/lib -Wl,-rpath,$(XERCESCROOT)/lib

$(PJCT)_SRCH_INCLD = \
	-I $(XMLC_ROOT)/src \
	$(XERCESC_EXPORT_INCL) \
	$(PYTHON_EXPORT_INCL)

$(PJCT)_EXTERNLIBS_NOSHARE = \
	$(XMLC_ROOT)/obj/$(PLATFORM)/wiglafNode.o \
	$(XMLC_ROOT)/obj/$(PLATFORM)/DOMTreeErrorReporter.o \

$(PJCT)_EXTERNLIBS = \
	$($(PJCT)_EXTERNLIBS_NOSHARE) \
	$(XERCESC_EXPORT_LIB) \
	$(PYTHON_EXPORT_LIB)

LIBS   = $(LIB_PATHS) $(XMLC_LIBS) $($(PJCT)_EXTERNLIBS) $(MISCLIB)


#----------------------------------------------------------------
# User rules
#----------------------------------------------------------------

exec : $($(PJCT)_EXEC)
	@ echo "\nProgram '$($(PJCT)_EXEC)' is up-to-date\n"

depend : $($(PJCT)_DEPMK)
	@ echo "   " $($(PJCT)_DEPMK) "is now up-to-date"

lib: $($(PJCT)_LIB)
	@ echo "\nLibrary '$($(PJCT)_LIB)' is up-to-date\n"


# cleanup executables
cleane:
	rm -f $(PJCT)_EXEC
	rm -f $(PJCT)_EXEC_LNK
	@ echo "   All executable files deleted"

# cleanup libraries
cleanl : 
	rm -rf $($(PJCT)_LIBD)
	@ echo "   All library archive files deleted"

# cleanup object files
cleano : 
	rm -rf $($(PJCT)_OBJD)
	@ echo "   All object files deleted"

# cleanup dependencies
cleand :
	rm -f $($(PJCT)_DEPMK)
	@ echo "   $(PJCT) dependency file deleted"

clean : cleano cleand
cleanall : cleane cleanl cleano cleand

#----------------------------------------------------------------
# Local rules
#----------------------------------------------------------------

$($(PJCT)_DEPMK) : $($(PJCT)_HDRS)
	@ echo "===== MAKING DEPENDENCIES ====="
	@ echo $($(PJCT)_OBJS)
	@ /bin/rm -f $($(PJCT)_DEPMK)
	@ $(CCMKD) $($(PJCT)_SRCH_INCLD) $($(PJCT)_SRCS) > $($(PJCT)_DEPMK)-temp
	@ sed "s?^.*\.o:?$($(PJCT)_OBJD)/&?" $($(PJCT)_DEPMK)-temp > $($(PJCT)_DEPMK)
	@ /bin/rm $($(PJCT)_DEPMK)-temp

$($(PJCT)_LIB) : $($(PJCT)_OBJS4LIB)
	@ echo "===== MAKING $(PJCT) LIBRARY ====="
	$(AR) -$(AR_OPTIONS) $(@) $($(PJCT)_OBJS4LIB) $($(PJCT)_EXTERNLIBS_NOSHARE)
	ln -s -f $(@) $($(PJCT)_LIB_LNK)

$($(PJCT)_EXEC) : $($(PJCT)_OBJS)
	@ echo "===== MAKING $(PJCT) EXECUTABLE ====="
	$(LD) $(LD_OPT) -o $(@) $($(PJCT)_OBJS) $(LIBS)
	ln -s -f ../$($(PJCT)_EXEC) $($(PJCT)_EXEC_LNK)

$($(PJCT)_OBJD)/%.o : $($(PJCT)_SRCD)/%.C
	@ mkdir -p $($(PJCT)_OBJD)
	$(CC) $(CCOPT) $($(PJCT)_SRCH_INCLD) $(<) -o $(@)

#----------------------------------------------------------------
# dependencies rules
#----------------------------------------------------------------

-include $($(PJCT)_DEPMK)

