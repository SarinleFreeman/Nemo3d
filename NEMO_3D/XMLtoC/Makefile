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
# $Header: /repo/XMLtoC/Makefile,v 1.4 2006/01/17 02:59:39 gekco Exp $
# *****************************************************************************


#================================================================
# Makefile to generate XMLtoC
#================================================================

# some boolean options
_USE_PYTHON = 1

# very top level stuff 
PJCT = XMLC
XMLC_ROOT = $(shell pwd)

ifeq ($(XMLC_BLD),)
	XMLC_BLD = absoft_mpi_fast
endif

# sets directory structure and default compiler options
include $(XMLC_ROOT)/global.mk

# sets variables containing file information
include $(XMLC_FILMK)

# sets compiler options
include $(XMLC_CMPMK)

#----------------------------------------------------------------
# User rules
#----------------------------------------------------------------

#default application
ifeq ($(APP),)
	APP = ioNemo3D
endif

$(APP)_ROOT = $(XMLC_ROOT)/apps/$(APP)
$(APP)_XMLD = $(XMLC_ROOT)/xml
$(APP)_XML0 = $(APP).xml

xmlc : dir $(XMLC_BIND)/$(XMLC_EXEC)
	@ echo "   XMLtoC up to date"

depend : $(XMLC_DEPMK)
	@ echo "   " $(XMLC_DEPMK) "is now up-to-date"

dir:
	@ mkdir -p $(XMLC_OBJD)
	@ mkdir -p $(XMLC_LIBD)



bldapp: xmlc
	mkdir -p $($(APP)_ROOT) $($(APP)_ROOT)/src $($(APP)_ROOT)/obj $($(APP)_ROOT)/lib $($(APP)_ROOT)/bin \
		$($(APP)_ROOT)/lib/$(PLATFORM); \
	cp $($(APP)_XMLD)/$($(APP)_XML0) $($(APP)_ROOT)/$($(APP)_XML0); \
	cp application.mk $($(APP)_ROOT)/Makefile; \
	cp $(XMLC_SRCD)/app_drivers/$(APP)_driver.C $($(APP)_ROOT)/src; \
	cp $(XMLC_SRCD)/app_drivers/main.C $($(APP)_ROOT)/src; \
	cp $(XMLC_SRCD)/extension.[CH] $($(APP)_ROOT)/src; \
	cp $(XMLC_SRCD)/app_drivers/$(APP)_funcs.py $($(APP)_ROOT)/bin/funcs.py; \
	ln -f $(XMLC_ROOT)/global.mk $($(APP)_ROOT); \
	cd $(XMLC_BIND); \
	cp $($(APP)_ROOT)/$($(APP)_XML0) .; \
	./$(XMLC_EXEC) $($(APP)_XML0); \
	mv *.[CH] $($(APP)_ROOT)/src;
#	
#	modification for Cygwin compilation Jan 16, 2006
#	compile line 
#	has difficulty in path on cygwin
#	can copy the XML file to the bin directory and 
#	compile in place instead of having the long pathname in front of the 
#       xml input deck
#	./$(XMLC_EXEC) $($(APP)_ROOT)/$($(APP)_XML0); \

applib: bldapp
	cd $($(APP)_ROOT); \
	$(MAKE) lib "PLATFORM=$(PLATFORM)" \
	     "PJCT=$(APP)" \
	     "$(APP)_ROOT=$($(APP)_ROOT)" \
	     "$(APP)_BLD=$(XMLC_BLD)" \
	     "XMLC_ROOT=$(XMLC_ROOT)" \
	     "XMLC_LIBS=$(LIBS)" \
	     "PYTHON_ROOT=$(PYTHON_ROOT)" "PYTHON_EXPORT_INCL=$(PYTHON_EXPORT_INCL)" "PYTHON_EXPORT_LIB=$(PYTHON_EXPORT_LIB)" \
	     "XERCESCROOT=$(XERCESCROOT)" "XERCESC_EXPORT_INCL=$(XERCESC_EXPORT_INCL)" "XERCESC_EXPORT_LIB=$(XERCESC_EXPORT_LIB)"

appexec: bldapp
	cd $($(APP)_ROOT); \
	$(MAKE) exec "PLATFORM=$(PLATFORM)" \
	     "PJCT=$(APP)" \
	     "$(APP)_ROOT=$($(APP)_ROOT)" \
	     "$(APP)_BLD=$(XMLC_BLD)" \
	     "XMLC_ROOT=$(XMLC_ROOT)" \
	     "XMLC_LIBS=$(LIBS)" \
	     "PYTHON_ROOT=$(PYTHON_ROOT)" "PYTHON_EXPORT_INCL=$(PYTHON_EXPORT_INCL)" "PYTHON_EXPORT_LIB=$(PYTHON_EXPORT_LIB)" \
	     "XERCESCROOT=$(XERCESCROOT)" "XERCESC_EXPORT_INCL=$(XERCESC_EXPORT_INCL)" "XERCESC_EXPORT_LIB=$(XERCESC_EXPORT_LIB)"

applib0:
	cd $($(APP)_ROOT); \
	$(MAKE) lib "PLATFORM=$(PLATFORM)" \
	     "PJCT=$(APP)" \
	     "$(APP)_ROOT=$($(APP)_ROOT)" \
	     "$(APP)_BLD=$(XMLC_BLD)" \
	     "XMLC_ROOT=$(XMLC_ROOT)" \
	     "XMLC_LIBS=$(LIBS)" \
	     "PYTHON_ROOT=$(PYTHON_ROOT)" "PYTHON_EXPORT_INCL=$(PYTHON_EXPORT_INCL)" "PYTHON_EXPORT_LIB=$(PYTHON_EXPORT_LIB)" \
	     "XERCESCROOT=$(XERCESCROOT)" "XERCESC_EXPORT_INCL=$(XERCESC_EXPORT_INCL)" "XERCESC_EXPORT_LIB=$(XERCESC_EXPORT_LIB)"

appexec0:
	cd $($(APP)_ROOT); \
	$(MAKE) exec "PLATFORM=$(PLATFORM)" \
	     "PJCT=$(APP)" \
	     "$(APP)_ROOT=$($(APP)_ROOT)" \
	     "$(APP)_BLD=$(XMLC_BLD)" \
	     "XMLC_ROOT=$(XMLC_ROOT)" \
	     "XMLC_LIBS=$(LIBS)" \
	     "PYTHON_ROOT=$(PYTHON_ROOT)" "PYTHON_EXPORT_INCL=$(PYTHON_EXPORT_INCL)" "PYTHON_EXPORT_LIB=$(PYTHON_EXPORT_LIB)" \
	     "XERCESCROOT=$(XERCESCROOT)" "XERCESC_EXPORT_INCL=$(XERCESC_EXPORT_INCL)" "XERCESC_EXPORT_LIB=$(XERCESC_EXPORT_LIB)"

appclean:
	cd $($(APP)_ROOT); \
	$(MAKE) clean "PLATFORM=$(PLATFORM)" \
	     "PJCT=$(APP)" \
	     "$(APP)_ROOT=$($(APP)_ROOT)" \
	     "$(APP)_BLD=$(XMLC_BLD)" \
	     "XMLC_ROOT=$(XMLC_ROOT)" \
	     "XMLC_LIBS=$(LIBS)" \
	     "PYTHON_ROOT=$(PYTHON_ROOT)" "PYTHON_EXPORT_INCL=$(PYTHON_EXPORT_INCL)" "PYTHON_EXPORT_LIB=$(PYTHON_EXPORT_LIB)" \
	     "XERCESCROOT=$(XERCESCROOT)" "XERCESC_EXPORT_INCL=$(XERCESC_EXPORT_INCL)" "XERCESC_EXPORT_LIB=$(XERCESC_EXPORT_LIB)"

appcleanall:
	rm -rf $($(APP)_ROOT)




# cleanup executables
cleane:
	rm -f $(XMLC_BIND)/$(XMLC_EXEC)
	rm -f $(XMLC_BIND)/$(XMLC_EXEC_LNK)
	@ echo "   All executable files deleted"

# cleanup libraries
cleanl : 
	rm -rf $(XMLC_LIBD)
	@ echo "   All library archive files deleted"

# cleanup object files
cleano : 
	rm -rf $(XMLC_OBJD)
	@ echo "   All object files deleted"

# cleanup dependencies
cleand :
	rm -f $(XMLC_DEPMK)
	@ echo "   XMLtoC dependency file deleted"

clean : cleano cleand appclean
cleanall : cleane cleano cleanl cleand appcleanall

#----------------------------------------------------------------
# Local rules
#----------------------------------------------------------------

$(XMLC_DEPMK) : $(XMLC_HDRS)
	@ echo "===== MAKING DEPENDENCIES ====="
	rm -f $(XMLC_DEPMK)
	$(CCMKD) $(XMLC_SRCH_INCLD) $(XMLC_SRCS) > $(XMLC_DEPMK)-temp
	sed "s?^.*\.o:?$(XMLC_OBJD)/&?" $(XMLC_DEPMK)-temp > $(XMLC_DEPMK)
	@ /bin/rm $(XMLC_DEPMK)-temp

$(XMLC_BIND)/$(XMLC_EXEC) : $(XMLC_OBJS)
	@ echo "===== MAKING XMLtoC EXECUTABLE ====="
	$(LD) $(LD_OPT) -o $(@) $(XMLC_OBJS) $(LIBS)
	ln -s -f $(XMLC_BIND)/$(XMLC_EXEC) $(XMLC_BIND)/$(XMLC_EXEC_LNK)

$(XMLC_OBJD)/%.o : $(XMLC_SRCD)/%.C
	@ mkdir -p $(XMLC_OBJD)
	$(CC) $(CCOPT) $(XMLC_SRCH_INCLD) $(<) -o $(@)

#----------------------------------------------------------------
# dependency rules
#----------------------------------------------------------------

include $(XMLC_DEPMK)

