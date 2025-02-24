#!/usr/bin/env python


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
# *****************************************************************************


two = 2.0
pi = 3.14159265

def multiply(a,b):
    print "Python:: multiply()"
    return a * b

def addPi(a):
    print "Python:: addPi()"
    return a + pi

def factorial(n):
    print "Python:: factorial()"
    if n < 2:
        return 1
    else:
        return n*factorial(n-1)

def squareEXT(x):
    print "Python:: squareEXT()"
    import extension
    return extension.square(x)

def extSum(x):
    print "Python:: extSum()"
    import extension
    return extension.sum(x)


def twice(f):
    return lambda x: f(f(x))

class C:
    def __init__(self, i=None):
        if i is None:
            print "No arguments"
        else:
            print "Argument is", i

# a higher order function; now no need
class linear:
    def __init__(self, a, b):
        self.a, self.b = a,b
        def __call__(self, x):
            return self.a * x + self.b

def linear2(a,b):
    return lambda x: a*x + b

linear2_ex = linear2(2,4)
