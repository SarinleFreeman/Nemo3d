/*****************************************************************************
The Jet Propulsion Laboratory (JPL) XML-to-C++ package.
Copyright (C) 2002 California Institute of Technology (Caltech)

This library is free software, which you can redistribute and/or modify
under the terms of the GNU Lesser General Public License as published by the
Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this library; see the file COPYING. If not, write to the
Free Software Foundation, Inc.,
59 Temple Place, Suite 330,
Boston, MA  02111-1307  USA

For additional information, please contact
  Fabiano Oyafuso (fabiano@jpl.nasa.gov)

Written by:  Fabiano Oyafuso
             Hook Hua
             Ed Vinyard

This product includes software developed by the Apache Software Foundation
(http://www.apache.org/).
*****************************************************************************/

#include "extension.H"

#ifdef USE_PYTHON

extern "C" {

    PyObject* extension_system(PyObject* self, PyObject* args) {
        char *command;
        int sts;

        if (!PyArg_ParseTuple(args, "s", &command))
            return NULL;
        sts = system(command);

        return Py_BuildValue("i", sts);
    }

    PyObject* extension_square(PyObject* self, PyObject* args) {
        double x;

        if (!PyArg_ParseTuple(args, "d", &x))
            return NULL;

        cout << "extension_square():: x = " << x << endl;

        double xx = x*x;

        return Py_BuildValue("d", xx);
    }

    PyObject* extension_sum(PyObject* self, PyObject* args) {
        double sum=0;
        double x;

        if (!PyArg_ParseTuple(args, "d", &x))
            return NULL;

        cout << "extension_sum():: x = " << x << endl;

        for (int i=1; i<=x; i++) {
            sum += i;
        }

        return Py_BuildValue("d", sum);
    }

    // Extension method definition mappings
    void initExtension(void) {

        cout << "initExtension():: calling Py_InitModule()... " << endl;

        (void) Py_InitModule("extension", ExtensionMethods);
    }

}
#endif /* USE_PYTHON */
