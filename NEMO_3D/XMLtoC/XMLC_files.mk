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
# $Header: /repo/XMLtoC/XMLC_files.mk,v 1.2 2003/10/08 16:15:17 hook Exp $
# *****************************************************************************


#================================================================
# XMLtoC Specification of files
#================================================================



XMLC_HDRS0 = \
        DOMTreeErrorReporter.H \
        wiglafNode.H \
        domStuff.H

XMLC_HDRS = $(patsubst %,$(XMLC_INCLD)/%,$(XMLC_HDRS0_TOP))


XMLC_OBJS0 = \
        DOMTreeErrorReporter.o \
        XMLtoC.o \
        wiglafNode.o

XMLC_OBJS = $(patsubst %,$(XMLC_OBJD)/%,$(XMLC_OBJS0))


XMLC_SRCS0 = $(patsubst %.o,%.C,$(XMLC_OBJS0))
XMLC_SRCS  = $(patsubst %,$(XMLC_SRCD)/%,$(XMLC_SRCS0))


XMLC_EXEC      = XMLtoC-$(PLATFORM).ex
XMLC_EXEC_LNK  = XMLtoC.ex
