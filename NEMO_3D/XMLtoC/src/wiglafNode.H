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

#ifndef _WIGLAFNODE_H
#define _WIGLAFNODE_H

#include "domStuff.H"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <fstream>

#ifdef USE_PYTHON
#include <Python.h>
#endif /* USE_PYTHON */

/* GNU GCC 3.0 and later requires explicit specification that <string> is from std */
using namespace std;


// ---------------------------------------------------------------------------

// A wrapper around DOM Nodes that provides the kind of convenience
// you'd get from a JDOM node (without JDOM).

#define INDENT "    "
#define TYPEDEF_SUFFIX "_struct"

const bool DEBUG = false; // Enable to turn on debug statements to trace sequence of instantiations, tree traversals, flow, etc.
const int MaxExpressionLength = 100; // Maximum length of the expression in the XML

class wiglafNode {
protected:
   DOM_Node node;

public:

   wiglafNode() {;}
   wiglafNode(const DOM_Node&);
   wiglafNode(const DOM_Node&, string);
   wiglafNode(const wiglafNode&);
   ~wiglafNode() {;}

   wiglafNode& operator =(const wiglafNode&);
   const DOM_Node& getDOMNode() const;
   DOMString getElementName() const;

   bool isNull() const;
   void testNullException() const;
   bool isObj() const;
   bool isTopLevel() const;
   wiglafNode getDocument() const;
   wiglafNode getParent() const;
   DOMString getAttribute(const DOMString&) const;
   DOMString getName() const;
   DOMString getText() const;
   bool isNamed() const;


   friend bool isOption(const wiglafNode&);
   friend bool isParam(const wiglafNode&);
   friend bool isParamButNotDeleteBranch(const wiglafNode&);
   friend bool isGroup(const wiglafNode&);
   friend bool isGroupOrParam(const wiglafNode&);


   wiglafNode getChild(string) const;
   vector<wiglafNode> getChildren( bool (*)(const wiglafNode&) ) const;
   vector<wiglafNode> getChildren() const;
   DOM_Node getChildDOMNode(const DOM_Node&, string);
   DOM_Node getInputElement(const DOM_Node&);
   wiglafNode getPath(string) const;
   DOMString getValue() const;
   void setValue(string);
   string getStringValue() const;
   int getIntValue() const;
   double getDoubleValue() const;
   bool getBooleanValue() const;
   vector<double> getVectorDoubleValue() const;
   vector<int> getVectorIntegerValue() const;
   vector<string> getVectorStringValue() const;
#ifdef USE_PYTHON
   double getExpressionDoubleValue(PyObject *) const;
   int getExpressionIntegerValue(PyObject *) const;
   void addStatement(PyObject *, string, string) const;
#endif/* USE_PYTHON */
   DOMString getDescription() const;
   DOMString get_cTag() const;
   DOMString getEnumType() const;
   DOMString getType() const;
   DOMString getEnumTypedef() const;
   DOMString meth_stringToEnum(string) const;
   DOMString meth_enumToString(string) const;
   wiglafNode getGroupParamAncestor() const;
   string getPathToRoot() const;
   string genMethod_constructor() const;
   string genMethod_set() const;
   string genMethod_setFromArrayElement() const;
   string genMethod_toString() const;
   string genClassDefinition(vector<string>&) const;
   void traverse() const;
   string topOfFile() const;
   void XMLtoC() const;

}; // End wiglaf class definition

// make these friends accessible outside class:
bool isParamButNotDeleteBranch(const wiglafNode&);
bool isGroup(const wiglafNode&);

#endif
