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

#include "wiglafNode.H"

// ---------------------------------------------------------------------------

// A wrapper around DOM Nodes that provides the kind of convenience
// you'd get from a JDOM node (without JDOM).


   // Create a Node associated with the specified DOM node. 
   wiglafNode::wiglafNode(const DOM_Node& dn) { 
       if (!dn.isNull()) {
           node = dn;
       }
       else cerr << "wiglafNode(dn):: dn is null!\n";
   }


   wiglafNode::wiglafNode(const DOM_Node& dn, string s)  {
       // wiglafNode w should be at the <input> level, where the children are <group> and/or <param> elements
       wiglafNode w( getInputElement(dn) );

       wiglafNode n = w.getPath(s);
       if (!n.isNull()) {
           if (DEBUG) cout << "wiglafNode(dn,s):: Created wiglafNode n from path= " << s << endl;
           if (DEBUG) cout << "        n.getDOMNode().getNodeName().transcode()= " << n.getDOMNode().getNodeName().transcode() << endl;
           node = n.getDOMNode();
       }
       else {
           cerr << "wiglafNode(dn,s):: w.getPath(s) returned a NULL wiglafNode n!\n";
       }
   }

   wiglafNode::wiglafNode(const wiglafNode& wn) {
      node = wn.node;
   }
   
   wiglafNode& wiglafNode::operator =(const wiglafNode& wn) {
      node = wn.node;
      return *this;
   }
   
   
   const DOM_Node& wiglafNode::getDOMNode() const { return node; }

   DOMString wiglafNode::getElementName() const { return node.getNodeName(); }
   
   bool wiglafNode::isNull() const { return node.isNull(); }

   void wiglafNode::testNullException() const {
      if (isNull()) cerr << "node must not be null.";;
   }

   bool wiglafNode::isObj() const {
      return getElementName().equals("group") && 
             getAttribute("type").equals("obj");
   }

   bool wiglafNode::isTopLevel() const {
      return node.getNodeType() == DOM_Node::DOCUMENT_NODE;
   }

   wiglafNode wiglafNode::getDocument() const {
      return wiglafNode(node.getOwnerDocument());
   }

   wiglafNode wiglafNode::getParent() const { 
      DOM_Node parent=getDOMNode().getParentNode();
      return ( parent.isNull() ? wiglafNode() : wiglafNode(parent) );
   }

   DOMString wiglafNode::getAttribute(const DOMString& name) const {
      return node.getAttributes().getNamedItem(name).getNodeValue();
   }

   DOMString wiglafNode::getName() const {
      if (!isNamed()) cerr << "getName():: node must be a named node.\n";
      return getChild("name").getText();
   }

   DOMString wiglafNode::getText() const {
      if ( getChild("#text").isNull() )  
         cerr << "getText():: wiglafNode must have a #text child.\n";
      return getChild("#text").getDOMNode().getNodeValue();
   }

   // ------------------------------------------------
   bool wiglafNode::isNamed() const {
      testNullException();
      return !getChild("name").isNull();
   }

   bool isOption(const wiglafNode& wn) {
      return wn.getElementName().equals("option");
   }

   bool isParam(const wiglafNode& wn) {
      return wn.getElementName().equals("param");
   }

   bool isParamButNotDeleteBranch(const wiglafNode& wn) {
      bool is_param = wn.getElementName().equals("param");
      bool is_deleteBranch = string(wn.getType().transcode())=="deletebranch";
      return is_param && !is_deleteBranch;
   }

   bool isGroup(const wiglafNode& wn) {
      return wn.getElementName().equals("group");
   }

   bool isGroupOrParam(const wiglafNode& wn) {
      return isGroup(wn) || isParam(wn);
   }
   // ------------------------------------------------

   wiglafNode wiglafNode::getChild(string name) const {
      if (name.empty()) cerr << "getChild(name):: name is empty!\n";

      DOMString name_ds = DOMString(name.c_str());
      const DOM_NodeList children = node.getChildNodes();

      for (unsigned int i = 0; i < children.getLength(); i++) {
         wiglafNode child = children.item(i);
         
         if (child.isNull()) continue;

         if ( child.isNamed() ) {
             if ( child.getName().equals(name_ds) ) return child;
         }

         // used when name="value" for going from value -> group
         if ( child.getElementName().equals(name_ds) ) {
            return child;
         }
      }
      return wiglafNode();
   }

   vector<wiglafNode> wiglafNode::getChildren( bool (*test)(const wiglafNode&) ) const {
      DOM_NodeList children = node.getChildNodes();
      int sz = children.getLength();

      vector<wiglafNode> result;
      //result.reserve(sz);
      
      for (int i = 0; i < sz; i++) {
         const wiglafNode& child = children.item(i);
         if( !test || test(child) )
            result.push_back(wiglafNode(child));
      }
      
      return result;
    }

   vector<wiglafNode> wiglafNode::getChildren() const {
      DOM_NodeList children = node.getChildNodes();
      int sz = children.getLength();
      
      vector<wiglafNode> result;
      result.reserve(sz);
      
      for (int i = 0; i < sz; i++)
         result.push_back(wiglafNode(children.item(i)));
	
      return result;
   }


   // Returns the child DOM_Node of dn whose getNodeName() equals DOMNodeName
   DOM_Node wiglafNode::getChildDOMNode(const DOM_Node& dn, string DOMNodeName) {

       DOM_Node child;
       DOM_NodeList children = dn.getChildNodes();

       for (unsigned int i = 0; i < children.getLength(); i++) {
           child = children.item(i);
           if (DEBUG) cout << " " << i << " child.getNodeName().transcode(): " << child.getNodeName().transcode() << endl;
           if ( child.getNodeName().equals( DOMNodeName.c_str() ) ) {
	       if (DEBUG) cout << " Found DOM_Node: " << DOMNodeName << endl;
               return child;
           }
       }

       cerr << "getChildDOMNode(dn,DOMNodeName):: DOMNodeName= " << DOMNodeName << " not found! returning dn." << endl;
       return dn;
   }

   // It traverses down to the DOM_Node: document root -> application -> input
   // "dn" must be the root node.
   // Returns the <input> node so that all the children are the input <groups> and/or <params> elements
   DOM_Node wiglafNode::getInputElement(const DOM_Node& dn) {

       DOM_Node child;

       // find "application" DOM_Node
       child = getChildDOMNode(dn, "application");

       // find "input" DOM_Node
       child = getChildDOMNode(child, "input");

       return child;
   }


   wiglafNode wiglafNode::getPath(string path) const {

      // remove possible leading "/"
      if (path.substr(0,1) == "/") { path.erase(0,1); }

      if (path.empty()) cerr << "getPath(path):: path is empty!\n";

      int i=0;
      int pos;
      string rempath = path;
      wiglafNode current = *this;
      wiglafNode child;

      while (!rempath.empty()) {

         // find position of first "/" in path
         pos = rempath.find("/", i);

         if ( current.isNamed() ) {
             if (DEBUG) cout << "getPath(path):: current.getName().transcode()= " << current.getName().transcode() << endl;
         }
         else {
             if (DEBUG) cout << "getPath(path):: current has no name!\n";
         }

         // if no more "/" in path
         if (pos==-1) {
            if (DEBUG) cout << "getPath(path):: last downpath, rempath = \"" << rempath << "\"\n";

            //return current.getParent();
	    
            child = current.getChild(rempath);
            if ( child.isNull() ) {

                // now check for case when we are at the <value> node
   	        child = current.getChild("value");
                if ( !child.getChild(rempath).isNull() ) {
                    // ...then the <value> element contains the child <group> with name=rempath 
              	    if (DEBUG) cout << "getPath(path):: reached <value> element where it's child is group with name= " << rempath << endl;
                    return child.getChild(rempath);
                }
                else {
                    cerr << "getPath(path):: LAST child is null! returning current.\n";
                    if (DEBUG) cout << "getPath(path):: current.getDOMNode().getNodeName().transcode()= " << current.getDOMNode().getNodeName().transcode() << endl;
                    return current;
                }

            }
            else {
                if (DEBUG) cout << "getPath(path):: LAST child exists. Returning child, where child.getDOMNode().getNodeName().transcode()= " << child.getDOMNode().getNodeName().transcode() << endl;
                return child;
            }
	    

         }

         // consume first token
         else {
            string downpath = rempath.substr(0,pos);
            rempath.erase(0,pos+1);
            if (DEBUG) cout << "getPath(path):: downpath = " << downpath << endl;

            // Just in case path starts from root. If at "application" or "input", do not have getNames() for XML element children
            if ( (downpath == "application") || (downpath == "input") ) {
                vector<wiglafNode> g = current.getChildren();
                for (vector<wiglafNode>::iterator iChild=g.begin(); iChild != g.end(); iChild++) {
		    if (DEBUG) cout << " searching children DOM_Node names: " << iChild->getDOMNode().getNodeName().transcode() << endl;
		    if ( iChild->getDOMNode().getNodeName().transcode() == downpath ) {
		        if (DEBUG) cout << " => Found DOM_Node: " << iChild->getDOMNode().getNodeName().transcode() << endl;
                        wiglafNode& found = *iChild;
                        current = found;
                        break;
                    }
                }
	    }

	    // else, we are in Group or Params, so we do have getNames() available
            else {
                child = current.getChild(downpath);
                if ( child.isNull() ) {
                   cerr << "getPath(path):: child is null, returning current.\n";
                   return current;
                }
                else {
                current = child;
                }
            }
         }
      }
      // shouldn't get here!
      return current;

   }

   DOMString wiglafNode::getValue() const {
      if ( getChild("value").isNull() ) 
         cerr << "getValue():: child <value> node is Null.";
      return getChild("value").getText();
   }

   void wiglafNode::setValue(string value) {
      DOMString value_ds = DOMString(value.c_str());
      getChild("value").getChild("#text").node.setNodeValue(value_ds);      
   }

   string wiglafNode::getStringValue() const {
      return string(getValue().transcode());
   }

   int wiglafNode::getIntValue() const {
       int val;
       sscanf(getStringValue().c_str(), "%d", &val);
       return val;
   }

   double wiglafNode::getDoubleValue() const {
      double val;
      sscanf(getStringValue().c_str(), "%lf", &val);
      return val;
   }

   bool wiglafNode::getBooleanValue() const {
      string s = getStringValue();
      if      (s=="true")
         return true;
      else if (s=="false")
         return false;
      else
         cerr << "invalid boolean value";

      return false;
   }

   vector<double> wiglafNode::getVectorDoubleValue() const {

      int indexStart;
      int indexEnd;
      string value = getStringValue();

      // To be formal, only look at everything between the { }
      indexStart = value.find("{", 0) + 1;
      indexEnd   = value.find("}", indexStart) - 1;
      value = value.substr(indexStart,indexEnd);

      char* s = new char[value.length()+1];
      char *token;
      char *delimiters = "\t ";
      double tokenDouble;
      vector<double> tokenList;

      strcpy( s, value.c_str() );
      token = strtok(s,delimiters);
      while (token) {
          if (DEBUG) cout << "getVectorDoubleValue():: token= " << token << endl;
          tokenDouble = strtod(token,NULL);
          //sscanf(token, "%lf", &tokenDouble);
          tokenList.push_back(tokenDouble);
          token = strtok(NULL,delimiters); // Next token
      }
      delete [] s;
      return tokenList;

   }

   vector<int> wiglafNode::getVectorIntegerValue() const {

      int indexStart;
      int indexEnd;
      string value = getStringValue();

      // To be formal, only look at everything between the { }
      indexStart = value.find("{", 0) + 1;
      indexEnd   = value.find("}", indexStart) - 1;
      value = value.substr(indexStart,indexEnd);
      if (DEBUG) cout << "getVectorIntegerValue():: value.length()= " << value.length() << endl;

      char* s = new char[value.length()+1];
      char *token;
      char *delimiters = "\t ";
      int tokenInteger;
      vector<int> tokenList;

      strcpy( s, value.c_str() );
      token = strtok(s,delimiters);
      while (token) {
          if (DEBUG) cout << "getVectorIntegerValue():: token= " << token << endl;
          tokenInteger = (int)strtol(token,NULL,10); // Base10
          //sscanf(token, "%ld", &tokenInteger);
          tokenList.push_back(tokenInteger);
          token = strtok(NULL,delimiters); // Next token
      }
      delete [] s;
      return tokenList;
   }


   vector<string> wiglafNode::getVectorStringValue() const {

      int indexStart;
      int indexEnd;
      string value = getStringValue();

      // To be formal, only look at everything between the { }
      indexStart = value.find("{", 0) + 1;
      indexEnd   = value.find("}", indexStart) - 1;
      value = value.substr(indexStart,indexEnd);
      if (DEBUG) cout << "getVectorStringValue():: value.length()= " << value.length() << endl;

      char* s = new char[value.length()+1];
      char *token;
      char *delimiters = "\t ";
      string tokenString;
      vector<string> tokenList;

      strcpy( s, value.c_str() );
      token = strtok(s,delimiters);
      while (token) {
          if (DEBUG) cout << "getVectorStringValue():: token= " << token << endl;
          tokenString = string(token); // Base10
          //sscanf(token, "%ld", &tokenString);
          tokenList.push_back(tokenString);
          token = strtok(NULL,delimiters); // Next token
      }
      delete [] s;
      return tokenList;
   }


#ifdef USE_PYTHON
   double wiglafNode::getExpressionDoubleValue(PyObject *pDict) const {

       string sExpression = getStringValue();
       //sExpression.insert(0,"renv.r_eval(\"");
       //sExpression.append("\")");

       char cExpression[sExpression.length()+1];
       strcpy( cExpression, sExpression.c_str() );

       // Return the value of the evaluated identifier EXPRESSION from the dictionary.
       PyObject *v = PyRun_String(cExpression, Py_eval_input, pDict, pDict);
       char *s = PyString_AsString( PyObject_Str( v ) );

       double result = 0.0;
       sscanf(s,"%lf", &result);

       if (DEBUG) cout << "getExpressionDoubleValue():: sExpression: " << sExpression << ",  result: " << result << endl;

       return result;

   }
#endif  /* USE_PYTHON */

#ifdef USE_PYTHON
   int wiglafNode::getExpressionIntegerValue(PyObject *pDict) const {

       string sExpression = getStringValue();
       //sExpression.insert(0,"renv.r_eval(\"");
       //sExpression.append("\")");

       char cExpression[sExpression.length()+1];
       strcpy( cExpression, sExpression.c_str() );

       // Return the value of the evaluated identifier EXPRESSION from the dictionary.
       PyObject *v = PyRun_String(cExpression, Py_eval_input, pDict, pDict);
       char *s = PyString_AsString( PyObject_Str( v ) );

       int result = 0;
       sscanf(s,"%d", &result);

       if (DEBUG) cout << "getExpressionIntegerValue():: sExpression: " << sExpression << ",  result: " << result << endl;

       return result;

   }
#endif  /* USE_PYTHON */


   // Adds a Python statement into the Dictionary. This allows relative nodes to access the values
   // stored in other relative nodes.
   // e.g.: a=1, b=2, c=3, sum=a+b+c
#ifdef USE_PYTHON
   void wiglafNode::addStatement(PyObject *pDict, string key, string value) const {

       string sExpression;
       sExpression = key + string(" = ") + value;
       //sExpression.insert(0,"renv.r_exec(\"");
       //sExpression.append("\")");

       char cExpression[sExpression.length()+1];
       strcpy( cExpression, sExpression.c_str() );
       PyObject *v = PyRun_String(cExpression, Py_file_input, pDict, pDict);

       if (DEBUG) {
           char *s = PyString_AsString( PyObject_Str( v ) );
           double result = 0.0;
           sscanf(s,"%lf", &result);
           cout << "addStatement():: sExpression Statement: " << sExpression << ",  Run Result: " << result << endl;
       }
   }
#endif /* USE_PYTHON */

   DOMString wiglafNode::getDescription() const {
      if ( getChild("desc").isNull() ) 
         cerr << "getDescription():: child <desc> node is Null.";
      return getChild("desc").getText();
   }

   DOMString wiglafNode::get_cTag() const {
      if (!isNamed()) cerr << "get_cTag():: node does NOT have a name.";
      return getChild("cTag").getText();
   }

   DOMString wiglafNode::getEnumType() const {
      wiglafNode parent = getParent();

      if ( !(getElementName().equals("param") &&
             getAttribute("type").equals("option") &&
             !parent.isNull() && parent.getElementName().equals("group")) )
         cerr << "this must be a <param type=\"option\"> AND parent must exist AND be of type group";

      return get_cTag() + "Type";
   }

   DOMString wiglafNode::getType() const {
      DOMString elementName = getElementName();

      // Param elements

      if ( elementName.equals("param") ) {
         DOMString t = getAttribute("type");

         if (      t.equals("integer") )
            return DOMString("int");
         else if ( t.equals("real") )
            return DOMString("double");
         else if ( t.equals("boolean") )
            return DOMString("short");
         else if ( t.equals("option") )
            return getEnumType();
         else if ( t.equals("string") )
            return DOMString("string");
         else if ( t.equals("array_real") )
            return DOMString("vector<double>");
         else if ( t.equals("array_integer") )
            return DOMString("vector<int>");
         else if ( t.equals("array_string") )
            return DOMString("vector<string>");
         else if ( t.equals("expr_real") )
            return DOMString("double");
         else if ( t.equals("expr_integer") )
            return DOMString("int");
         else
            return t;
      } 


      // Group elements

      else if ( elementName.equals("group") ) {
           // NB:  We've now required groups to have the attribute "type",
           //      although it would be better to allow for a
           //      default type (ie a struct) for groups; for now we retain
           //      the exception raising of a typeless group
         DOMString t = getAttribute("type");

         if      (t.equals("struct"))
            return get_cTag() + TYPEDEF_SUFFIX;

         else if (t.equals("obj"))
            return DOMString("NULL");

         else if (t.equals("array")) {
            // Get the array element type from the "newbranch" group
            vector<wiglafNode> paramChildren = getChildren(isParam);
            for (vector<wiglafNode>::iterator paramChild=paramChildren.begin(); paramChild != paramChildren.end(); paramChild++) {
                if ( paramChild->getAttribute("type").equals("newbranch") ) {
                    wiglafNode groupInNewBranch(paramChild->getPath("value/group").getDOMNode());
                    if ( groupInNewBranch.getAttribute("type").equals("struct") ) {
                        return groupInNewBranch.get_cTag() + TYPEDEF_SUFFIX;
                    }
                    else { cerr << "getType():: Group in <param type=\"newbranch\"> is NOT of type=struct!\n"; }
                }
            }
            return DOMString("void*"); 
         }

         else
            return get_cTag() + TYPEDEF_SUFFIX;
      }

      else { // Not param nor group
         return DOMString();
      }
   }

   // Generates the enumeration type definition. i.e. "enum {elem1, elem2, elem3} "
   DOMString wiglafNode::getEnumTypedef() const {
      DOMString result = DOMString("enum ");
      result += getEnumType();
      result += " {";

      vector<wiglafNode> e = getChildren(isOption);

      vector<wiglafNode>::iterator it;
      for (it=e.begin(); it != e.end(); it++) {
         wiglafNode& opt = *it;
         result.appendData(opt.getText());
         if (it+1 != e.end())
            result.appendData(", ");
      }

      result += "} ";
      
      return result;
   }

   DOMString wiglafNode::meth_stringToEnum(string ind) const {
      DOMString result = DOMString("get");
      result.appendData(getEnumType().transcode());
      result.appendData("(const string& s) {\n");
      
      if (DEBUG) {
          result.appendData( "        cout << \"getEnumType___():: s = \" << s << endl;\n");
      }
      
      vector<wiglafNode> e = getChildren(isOption);

      vector<wiglafNode>::iterator it;
      for (it=e.begin(); it != e.end(); it++) {

         wiglafNode& opt = *it;
         result += ind.c_str();
         result += INDENT;
         result += "if (s==\"";
         result += opt.getText();
         result += "\")  return ";
         result += opt.getText();
         result += ";\n";
      }
      result += INDENT;
      result += INDENT;
      result += "cerr << \"Invalid string!  ( returning arbitrary value! )\" << endl;\n";
      result += INDENT;
      result += INDENT;
      result += "return ";
      result += e.begin()->getText();
      result += ";\n";
      result += INDENT;
      result += "}\n";
      
      return result;
   }


   DOMString wiglafNode::meth_enumToString(string ind) const {
      DOMString result = DOMString("getEnumString");
      result.appendData(getEnumType().transcode());
      result.appendData("(const ");
      result.appendData(getEnumType().transcode());
      result.appendData("& e) {\n");

      if (DEBUG) {
          result.appendData( "        cout << \"getEnumString_____():: (int)e = \" << (int)e << endl;\n");
      }

      vector<wiglafNode> e = getChildren(isOption);

      vector<wiglafNode>::iterator it;
      for (it=e.begin(); it != e.end(); it++) {

         wiglafNode& opt = *it;
         result += ind.c_str();
         result += INDENT;
         result += "if (e==";
         result += opt.getText();
         result += ")  return string(\"";
         result += opt.getText();
         result += "\");\n";
      }
      result += INDENT;
      result += INDENT;
      result += "return string(\"\");\n";
      result += INDENT;
      result += "}\n";
      
      return result;
   }


   wiglafNode wiglafNode::getGroupParamAncestor() const {
      for (wiglafNode w=getParent(); !w.isNull(); w=w.getParent()) {
         if ( isGroupOrParam(w) )
            return w;
      }
      return wiglafNode();
   }


   // Returns a "/" delimited path of Group/Param child element <name>s 
   string wiglafNode::getPathToRoot() const {
      string path;
      for (wiglafNode wn=*this; !wn.isNull() && wn.isNamed(); wn=wn.getGroupParamAncestor()) {
         path.insert(0,wn.getName().transcode()).insert(0,"/");
      }

      return path;
   }

   
  // Returns Implementation output code for the constructor
   string wiglafNode::genMethod_constructor() const {
      if (getAttribute("type").equals("array")) return string("");

      string thisType = getType().transcode();
#ifdef USE_PYTHON
      string top = thisType + "::" + thisType + "(const DOM_Node& dn, PyObject *pDict) {\n";
#else
      string top = thisType + "::" + thisType + "(const DOM_Node& dn) {\n";
#endif
      string bottom = "}\n\n";
      string body;

#ifdef USE_PYTHON
      body.append( string(INDENT) + "set(dn,pDict,\"" + getPathToRoot() + "\");\n");
#else
      body.append( string(INDENT) + "set(dn,\"" + getPathToRoot() + "\");\n");
#endif

      return top + body + bottom;
   }


  // Returns Implementation output code for method set()
   string wiglafNode::genMethod_set() const {
      if (getAttribute("type").equals("array")) return string("");

      string thisType = getType().transcode();
#ifdef USE_PYTHON
      string top = "void " + thisType + "::set(const DOM_Node& dn, PyObject *pDict, string s) {\n";
#else
      string top = "void " + thisType + "::set(const DOM_Node& dn, string s) {\n";
#endif
      string bottom = "}\n\n";
      string body;

      if (DEBUG) body = string(INDENT) + "cout << \"set():: setting " + thisType + " from path: \" << s << endl;\n";

      string path = getPathToRoot();
      if (DEBUG) cout << "genMethod_set():: getPathToRoot(): " << getPathToRoot() << endl;


      // Params
#ifdef USE_PYTHON
      body.append( INDENT + string("char buffer[20];\n") );
#endif
      body.append( INDENT + string("wiglafNode w(dn,s);\n") );
      body.append( INDENT + string("vector<wiglafNode> g = w.getChildren(isParamButNotDeleteBranch);\n") );
      body.append( INDENT + string("vector<wiglafNode>::iterator child=g.begin();\n\n") );

      vector<wiglafNode> p = getChildren(isParamButNotDeleteBranch);
      string strType;
      string strIdent;
      for (vector<wiglafNode>::iterator wnp=p.begin(); wnp != p.end(); wnp++) {
             strType = string(wnp->getType().transcode());
             strIdent = string(wnp->get_cTag().transcode());

             // Check for expression first since getType() returns doubles the C/C++ code type.
             // Double from reals or expression reals cannot be distinguisted in this version with getType().
             // Expression
#ifdef USE_PYTHON
             if ( wnp->getAttribute("type").equals("expr_real") ) {
                body.append(string(INDENT) + strIdent + string(" = (child)->getExpressionDoubleValue(pDict);\n"));
                body.append(string(INDENT) + string("sprintf(buffer, \"%f\", " + strIdent + ");\n"));
                body.append(string(INDENT) + string("(child++)->addStatement(pDict, \"" + strIdent + "\", string(buffer) );\n"));
             }
             else if ( wnp->getAttribute("type").equals("expr_integer") ) {
                body.append(INDENT + strIdent + string(" = (child)->getExpressionIntegerValue(pDict);\n"));
                body.append(string(INDENT) + string("sprintf(buffer, \"%d\", " + strIdent + ");\n"));
                body.append(string(INDENT) + string("(child++)->addStatement(pDict, \"" + strIdent + "\", string(buffer) );\n"));
             }
             // Primitive
             else 
#endif  /* USE_PYTHON */
	     if (strType=="int") {
#ifdef USE_PYTHON
                body.append(INDENT + strIdent + string(" = (child)->getIntValue();\n"));
                body.append(string(INDENT) + string("sprintf(buffer, \"%d\", " + strIdent + ");\n"));
                body.append(string(INDENT) + string("(child++)->addStatement(pDict, \"" + strIdent + "\", string(buffer) );\n"));
#else
                body.append(INDENT + strIdent + string(" = (child++)->getIntValue();\n"));
#endif
             }
             else if (strType=="double") {
#ifdef USE_PYTHON
                body.append(INDENT + strIdent + string(" = (child)->getDoubleValue();\n"));
                body.append(string(INDENT) + string("sprintf(buffer, \"%f\", " + strIdent + ");\n"));
                body.append(string(INDENT) + string("(child++)->addStatement(pDict, \"" + strIdent + "\", string(buffer) );\n"));
#else
                body.append(INDENT + strIdent + string(" = (child++)->getDoubleValue();\n"));
#endif
             }
             else if (strType=="short") {
                body.append(INDENT + strIdent + string(" = (child++)->getBooleanValue();\n"));
             }
             else if (strType=="string") {
                body.append(INDENT + strIdent + string(" = (child++)->getStringValue();\n"));
             }

             // Vector
             else if (strType=="vector<double>") {
                body.append(INDENT + strIdent + string(" = (child++)->getVectorDoubleValue();\n"));
             }
             else if (strType=="vector<int>") {
                body.append(INDENT + strIdent + string(" = (child++)->getVectorIntegerValue();\n"));
             }
             else if (strType=="vector<string>") {
                body.append(INDENT + strIdent + string(" = (child++)->getVectorStringValue();\n"));
             }

             // Option
             else if (wnp->getAttribute("type").equals("option")) {
                body.append(INDENT + strIdent + string(" = ") +
                            string("get") + string(wnp->getEnumType().transcode()) + 
                            string("((child++)->getStringValue());\n"));
                if (DEBUG) body.append( string(INDENT) + string(INDENT) + string("cout << \"Just set " + strIdent + " = \" << (int)" + strIdent + " << endl;\n") );
             }

             // Default
             else {
                body.append(INDENT + string("child++; // ") + strIdent + " " + string(wnp->getElementName().transcode()) + ";\n" );
             }
      }
      body.append("\n");


      // Groups

      vector<wiglafNode> g = getChildren(isGroup);
      for (vector<wiglafNode>::iterator child=g.begin(); child != g.end(); child++) {

          // Struct type, One instance
          if ( child->getAttribute("type").equals("struct") ) {
#ifdef USE_PYTHON
              body.append(INDENT + string(child->get_cTag().transcode()) + string(".set(dn,pDict,\"" + child->getPathToRoot() + "\");\n") );
#else
              body.append(INDENT + string(child->get_cTag().transcode()) + string(".set(dn,\"" + child->getPathToRoot() + "\");\n") );
#endif
          }

          // Array type, multiple instances
          else if ( child->getAttribute("type").equals("array") ) {

              string elementType = child->getType().transcode();

              body.append( string(INDENT) + string("vector<wiglafNode> allGroupsInTypeStruct = w.getChildren(isGroup);\n") );
              body.append( string(INDENT) + elementType + string(" item;\n") );
	      body.append( string(INDENT) + string("for (vector<wiglafNode>::iterator groupInTypeStruct = allGroupsInTypeStruct.begin(); groupInTypeStruct != allGroupsInTypeStruct.end(); groupInTypeStruct++ ) {\n") );
	      body.append( string(INDENT) + string(INDENT) + string("if ( groupInTypeStruct->getAttribute(\"type\").equals(\"array\") ) {\n") );
	      body.append( string(INDENT) + string(INDENT) + string(INDENT) + string("int offset = 0;\n") );
	      body.append( string(INDENT) + string(INDENT) + string(INDENT) + string("vector<wiglafNode> allGroupsInTypeArray = groupInTypeStruct->getChildren(isGroup);\n") );
	      body.append( string(INDENT) + string(INDENT) + string(INDENT) + string("for (vector<wiglafNode>::iterator groupInTypeArray = allGroupsInTypeArray.begin(); groupInTypeArray != allGroupsInTypeArray.end(); groupInTypeArray++ ) {\n") );
	      body.append( string(INDENT) + string(INDENT) + string(INDENT) + string(INDENT) + string("if ( groupInTypeArray->getAttribute(\"type\").equals(\"obj\") ) {\n") );
	      if (DEBUG) body.append( string(INDENT) + string(INDENT) + string(INDENT) + string(INDENT) + string(INDENT) + string("cout << \"" + thisType + "():: setting " + child->get_cTag().transcode() + " vector item, getPathToRoot()= \" << groupInTypeArray->getPathToRoot() << endl;\n") );
#ifdef USE_PYTHON
	      body.append( string(INDENT) + string(INDENT) + string(INDENT) + string(INDENT) + string(INDENT) + string("item.setFromArrayElement(dn,pDict,groupInTypeStruct->getPathToRoot(),offset++);\n") );
#else
	      body.append( string(INDENT) + string(INDENT) + string(INDENT) + string(INDENT) + string(INDENT) + string("item.setFromArrayElement(dn,groupInTypeStruct->getPathToRoot(),offset++);\n") );
#endif
	      body.append( string(INDENT) + string(INDENT) + string(INDENT) + string(INDENT) + string(INDENT) + child->get_cTag().transcode() + string(".push_back(item);\n") );
	      body.append( string(INDENT) + string(INDENT) + string(INDENT) + string(INDENT) + string("}\n") );
	      body.append( string(INDENT) + string(INDENT) + string(INDENT) + string("}\n") );
	      body.append( string(INDENT) + string(INDENT) + string("}\n") );
	      body.append( string(INDENT) + string("}\n") );

          }
      }

      return top + body + bottom;
   }


  // Returns Implementation output code for method setFromArrayElement()
   string wiglafNode::genMethod_setFromArrayElement() const {
      if (getAttribute("type").equals("array")) return string("");

      string thisType = getType().transcode();
#ifdef USE_PYTHON
      string top = "void " + thisType + "::setFromArrayElement(const DOM_Node& dn, PyObject *pDict, string pathToArray, int offset) {\n";
#else
      string top = "void " + thisType + "::setFromArrayElement(const DOM_Node& dn,  string pathToArray, int offset) {\n";
#endif
      string bottom = "}\n\n";
      string body;

      if (DEBUG) body = string(INDENT) + "cout << \"setFromArrayElement():: setting " + thisType + " from element \" << offset << \" of path: \" << pathToArray << endl;\n";

      string path = getPathToRoot();
      if (DEBUG) cout << "genMethod_setFromArrayElement_():: getPathToRoot(): " << getPathToRoot() << endl;


      // Params
#ifdef USE_PYTHON
      body.append( INDENT + string("char buffer[20];\n\n") );
#endif
      // string pathToArray is the Path to the <group type="array"> node. It will get the index offset child element
      body.append( INDENT + string("wiglafNode w(dn,pathToArray);\n\n") );
      body.append( INDENT + string("if (!( w.getElementName().equals(\"group\") || w.getAttribute(\"type\").equals(\"array\") )) {\n") );
      body.append( INDENT + string(INDENT) + string("cerr << \"" + thisType + "::setFromArrayElement(): wiglafNode is not <group type=array> !\\n\";\n") );
      body.append( INDENT + string(INDENT) + string("return;\n") );
      body.append( INDENT + string("}\n\n") );
      body.append( INDENT + string("int index = 0;\n") );
      body.append( INDENT + string("vector<wiglafNode> allGroupsInTypeArray = w.getChildren(isGroup);\n") );

      body.append( string(INDENT) + string("for (vector<wiglafNode>::iterator groupInTypeArray = allGroupsInTypeArray.begin(); groupInTypeArray != allGroupsInTypeArray.end(); groupInTypeArray++ ) {\n") );
      body.append( string(INDENT) + string(INDENT) + string("if ( groupInTypeArray->getAttribute(\"type\").equals(\"obj\") && (index == offset) ) {\n\n") );
      if (DEBUG) body.append( string(INDENT) + string(INDENT) + string(INDENT) + string("cout << \"Material_struct.setArrayElement():: setting element offset: \" << offset << endl;\n") );

      body.append( string(INDENT) + string(INDENT) + string(INDENT) + string("vector<wiglafNode> g = groupInTypeArray->getChildren(isParamButNotDeleteBranch);\n") );
      body.append( string(INDENT) + string(INDENT) + string(INDENT) + string("vector<wiglafNode>::iterator child=g.begin();\n") );

      vector<wiglafNode> p = getChildren(isParamButNotDeleteBranch);
      for (vector<wiglafNode>::iterator wnp=p.begin(); wnp != p.end(); wnp++) {
             string strType = string(wnp->getType().transcode());
             string strIdent = string(wnp->get_cTag().transcode());

#ifdef USE_PYTHON
             // Expression
             if ( wnp->getAttribute("type").equals("expr_double") ) {
                body.append( string(INDENT) + string(INDENT) + string(INDENT) + strIdent + string(" = (child)->getExpressionDoubleValue(pDict);\n"));
                body.append( string(INDENT) + string(INDENT) + string(INDENT) + string("sprintf(buffer, \"%f\", " + strIdent + ");\n"));
                body.append( string(INDENT) + string(INDENT) + string(INDENT) + string("(child++)->addStatement(pDict, \"" + strIdent + "\", string(buffer) );\n"));
             }
             else if ( wnp->getAttribute("type").equals("expr_integer") ) {
                body.append( string(INDENT) + string(INDENT) + string(INDENT) + strIdent + string(" = (child++)->getExpressionIntegerValue(pDict);\n"));
                body.append( string(INDENT) + string(INDENT) + string(INDENT) + string("sprintf(buffer, \"%d\", " + strIdent + ");\n"));
                body.append( string(INDENT) + string(INDENT) + string(INDENT) + string("(child++)->addStatement(pDict, \"" + strIdent + "\", string(buffer) );\n"));
             }
             else 
#endif /*  USE_PYTHON */
	     if (strType=="int") {
#ifdef USE_PYTHON
                body.append( string(INDENT) + string(INDENT) + string(INDENT) + strIdent + string(" = (child)->getIntValue();\n"));
                body.append( string(INDENT) + string(INDENT) + string(INDENT) + string("sprintf(buffer, \"%d\", " + strIdent + ");\n"));
                body.append( string(INDENT) + string(INDENT) + string(INDENT) + string("(child++)->addStatement(pDict, \"" + strIdent + "\", string(buffer) );\n"));
#else
                body.append( string(INDENT) + string(INDENT) + string(INDENT) + strIdent + string(" = (child++)->getIntValue();\n"));
#endif
             }
             else if (strType=="double") {
#ifdef USE_PYTHON
                body.append( string(INDENT) + string(INDENT) + string(INDENT) + strIdent + string(" = (child)->getDoubleValue();\n"));
                body.append( string(INDENT) + string(INDENT) + string(INDENT) + string("sprintf(buffer, \"%f\", " + strIdent + ");\n"));
                body.append( string(INDENT) + string(INDENT) + string(INDENT) + string("(child++)->addStatement(pDict, \"" + strIdent + "\", string(buffer) );\n"));
#else
                body.append( string(INDENT) + string(INDENT) + string(INDENT) + strIdent + string(" = (child++)->getDoubleValue();\n"));
#endif
             }
             else if (strType=="short") {
                body.append( string(INDENT) + string(INDENT) + string(INDENT) + strIdent + string(" = (child++)->getBooleanValue();\n"));
             }
             else if (strType=="string") {
                body.append( string(INDENT) + string(INDENT) + string(INDENT) + strIdent + string(" = (child++)->getStringValue();\n"));
             }
             else if (strType=="vector<double>") {
                body.append( string(INDENT) + string(INDENT) + string(INDENT) + strIdent + string(" = (child++)->getVectorDoubleValue();\n"));
             }
	     else if (strType=="vector<int>") {
                body.append( string(INDENT) + string(INDENT) + string(INDENT) + strIdent + string(" = (child++)->getVectorIntegerValue();\n"));
             }
	     else if (strType=="vector<string>") {
                body.append( string(INDENT) + string(INDENT) + string(INDENT) + strIdent + string(" = (child++)->getVectorStringValue();\n"));
             }
             else if (wnp->getAttribute("type").equals("option")) {
                body.append( string(INDENT) + string(INDENT) + string(INDENT) + strIdent + string(" = ") +
                             string("get") + string(wnp->getEnumType().transcode()) + 
                             string("((child++)->getStringValue());\n"));
                if (DEBUG) body.append( string(INDENT) + string(INDENT) + string("cout << \"Just set " + strIdent + " = \" << (int)" + strIdent + " << endl;\n") );
             }
             else
                body.append( string(INDENT) + string(INDENT) + string(INDENT) + string("child++; // ") + strIdent + " " + string(wnp->getElementName().transcode()) + ";\n" );
      }
      body.append( string(INDENT) + string(INDENT) + string(INDENT) + string("break;\n") );
      body.append("\n");
      body.append( string(INDENT) + string(INDENT) + string("}\n") );
      body.append( string(INDENT) + string(INDENT) + string("index++;\n") );
      body.append( string(INDENT) + string("}\n") );


      // Groups

      return top + body + bottom;
   }


  // Returns Implementation output code for method toString()
   string wiglafNode::genMethod_toString() const {
      if (getAttribute("type").equals("array")) return string("");

      string thisType = getType().transcode();
      string top = "string " + thisType + "::toString(int indentLength) {\n";
      string bottom = "}\n\n";
      string body;

      if (DEBUG) body = string(INDENT) + "cout << \"" + thisType + "::toString()\"<< endl;\n";

      string path = getPathToRoot();
      if (DEBUG) cout << "genMethod_set():: getPathToRoot(): " << path << endl;

      //body.append( INDENT + string("wiglafNode w(dn,s);\n") );

      body.append( INDENT + string("int i;\n") );
      body.append( INDENT + string("char buffer[20];\n") );
      body.append( INDENT + string("const string Indent(\"    \");\n\n") );

      body.append( INDENT + string("string indentSpaces;\n") );
      body.append( INDENT + string("for (i=0; i<indentLength; i++) {;\n") );
      body.append( INDENT + string(INDENT) + string("indentSpaces.append(Indent);\n") );
      body.append( INDENT + string("}\n") );
      body.append( INDENT + string("string indentSpacesMore = string(indentSpaces+Indent);\n") );
      body.append( INDENT + string("indentLength++;\n\n") );


      body.append( INDENT + string("string output( string(indentSpaces) + \"---- Begin " + thisType + " -----\\n\");\n") );

      // Params

      vector<wiglafNode> p = getChildren(isParamButNotDeleteBranch);
      for (vector<wiglafNode>::iterator wnp=p.begin(); wnp != p.end(); wnp++) {
             string strType = string(wnp->getType().transcode());
             string strIdent = string(wnp->get_cTag().transcode());


             if ( (strType == "int") || (strType == "short")  ) {
                 body.append( string(INDENT) + "sprintf(buffer, \"%d\", " +  strIdent + ");\n");
                 body.append( string(INDENT) + "output.append( string(indentSpacesMore) + string(\"" + strIdent + ": \") + string(buffer) + string(\"\\n\"));\n");
             }
             else if (strType == "double") {
                 body.append( string(INDENT) + "sprintf(buffer, \"%f\", " +  strIdent + ");\n");
                 body.append( string(INDENT) + "output.append( string(indentSpacesMore) + string(\"" + strIdent + ": \") + string(buffer) + string(\"\\n\"));\n");
             }
             else if (strType == "string") {
                 body.append( string(INDENT) + "output.append( string(indentSpacesMore) + string(\"" + strIdent + ": \\\"\") + " + strIdent + " + string(\"\\\"\\n\"));\n");
             }
             else if (strType == "vector<double>") {
                     body.append( string(INDENT) + "i = 0;\n");
                     body.append( string(INDENT) + "for (vector<double>::iterator item=" + strIdent + ".begin(); item != " + strIdent + ".end(); item++) {\n");
                     body.append( string(INDENT) + string(INDENT) + "sprintf(buffer, \"%d\", i++);\n");
                     body.append( string(INDENT) + string(INDENT) + "output.append( string(indentSpacesMore) + \"" + strIdent + "[\" + string(buffer) + \"]: \");\n");
                     body.append( string(INDENT) + string(INDENT) + "sprintf(buffer, \"%f\", (double)*item);\n");
                     body.append( string(INDENT) + string(INDENT) + "output.append( string(buffer) + \"\\n\");\n");
                     body.append( string(INDENT) + "}\n");
             }
             else if (strType == "vector<int>") {
                     body.append( string(INDENT) + "i = 0;\n");
                     body.append( string(INDENT) + "for (vector<int>::iterator item=" + strIdent + ".begin(); item != " + strIdent + ".end(); item++) {\n");
                     body.append( string(INDENT) + string(INDENT) + "sprintf(buffer, \"%d\", i++);\n");

                     body.append( string(INDENT) + string(INDENT) + "output.append( string(indentSpacesMore) + \"" + strIdent + "[\" + string(buffer) + \"]: \");\n");
                     body.append( string(INDENT) + string(INDENT) + "sprintf(buffer, \"%d\", (int)*item);\n");
                     body.append( string(INDENT) + string(INDENT) + "output.append( string(indentSpacesMore) + string(buffer) + \"\\n\");\n");

                     body.append( string(INDENT) + "}\n");
             }
             else if (strType == "vector<string>") {
                     body.append( string(INDENT) + "i = 0;\n");
                     body.append( string(INDENT) + "for (vector<string>::iterator item=" + strIdent + ".begin(); item != " + strIdent + ".end(); item++) {\n");
                     body.append( string(INDENT) + string(INDENT) + "sprintf(buffer, \"%d\", i++);\n");
                     body.append( string(INDENT) + string(INDENT) + "output.append( string(indentSpacesMore) + \"" + strIdent + "[\" + string(buffer) + \"]: \");\n");
                     body.append( string(INDENT) + string(INDENT) + "output.append(*item + \"\\n\");\n");
                     body.append( string(INDENT) + "}\n");
             }
             else if (wnp->getAttribute("type").equals("option")) {
                     body.append( string(INDENT) + "output.append( string(indentSpacesMore) + string(\"" + strIdent + ": enum \") + string(getEnumString" + strType + "(" + strIdent + ")) + string(\"\\n\") );\n");
             }

      }
      body.append("\n");


      // Groups

      vector<wiglafNode> g = getChildren(isGroup);
      for (vector<wiglafNode>::iterator child=g.begin(); child != g.end(); child++) {

          // Struct type, One instance
          if ( child->getAttribute("type").equals("struct") ) {
              body.append( string(INDENT) + "output.append( " + string(child->get_cTag().transcode()) + string(".toString(indentLength) );\n") );
          }

          // Array type, multiple instances
          else if ( child->getAttribute("type").equals("array") ) {

              string elementType = child->getType().transcode();

	      body.append( string(INDENT) + string("for (vector<" + elementType + ">::iterator item = " + child->get_cTag().transcode() + ".begin(); item != " + child->get_cTag().transcode() + ".end(); item++ ) {\n") );
	      if (DEBUG) body.append( string(INDENT) + string(INDENT) + string("cout << \"" + thisType + "():: toString(\" << indentLength << \") for " + child->get_cTag().transcode() + " vector item\" << endl;\n") );
	      body.append( string(INDENT) + string(INDENT) + string("output.append( item->toString(indentLength) );\n") );
	      body.append( string(INDENT) + string("}\n") );

          }
      }

      body.append( string(INDENT) + "output.append( string(indentSpaces) + \"---- End " + thisType + " -----\\n\");\n");
      body.append( string(INDENT) + "output.append(\"\\n\");\n");
      body.append( string(INDENT) + "return output;\n");

      return top + body + bottom;
   }



  // Returns Header output code
   string wiglafNode::genClassDefinition(vector<string> &dependencies) const {

      string class_comment = string("/* ") + string(getDescription().transcode()) + string("*/\n");
      string class_top = string("class ") + getType().transcode() + string(" {\n");
      class_top += "public:\n";
      string class_end = "\n};\n";
      string class_body;
      string methods_pub;


      // Generate class methods
#ifdef USE_PYTHON
      methods_pub.append( INDENT + string(getType().transcode()) + "(const DOM_Node&, PyObject *);\n" );   // Constructor
#else
      methods_pub.append( INDENT + string(getType().transcode()) + "(const DOM_Node&);\n" );   // Constructor
#endif  /* USE_PYTHON */
      methods_pub.append( INDENT + string(getType().transcode()) + "()  {;}\n" );                 // Constructor-default
      methods_pub.append( INDENT + string("~") + string(getType().transcode()) + "() {;}\n\n" );  // Destructor-default

#ifdef USE_PYTHON
      methods_pub.append( INDENT + string("void set(const DOM_Node&, PyObject *, string);\n") );
      methods_pub.append( INDENT + string("void setFromArrayElement(const DOM_Node&, PyObject *, string, int);\n") );
#else
      methods_pub.append( INDENT + string("void set(const DOM_Node&, string);\n") );
      methods_pub.append( INDENT + string("void setFromArrayElement(const DOM_Node&,  string, int);\n") );
#endif
      methods_pub.append( INDENT + string("string toString(int);\n") );


      // Generate class body
      vector<wiglafNode> e = getChildren(isGroupOrParam);
      for (vector<wiglafNode>::iterator wnp=e.begin(); wnp != e.end(); wnp++) {
          string comment(INDENT);
          comment.append( string("/* ") + string(wnp->getDescription().transcode()) + string(" */") );

          if (!wnp->getType().equals("deletebranch")) {
              class_body.append(string("\n") + INDENT);

              // <param>
              if ( wnp->getElementName().equals("param") ) {

                  // type=option
                  if ( wnp->getAttribute("type").equals("option") ) {

                      class_body.append(string(wnp->getEnumTypedef().transcode()));

                      methods_pub.append( string("\n") + string(INDENT) + 
                         wnp->getEnumType().transcode() + string(" ") + 
                         wnp->meth_stringToEnum(string(INDENT)).transcode() );

                      methods_pub.append( string("\n") + string(INDENT) + 
                         string("string ") + wnp->meth_enumToString(string(INDENT)).transcode() );

                  }

                  // other "primitive" param types
                  else {
                      class_body.append(string(wnp->getType().transcode()) + " ");
                  }
              }

              // <group>
              else if ( wnp->getElementName().equals("group") ) {

                  dependencies.push_back( wnp->getType().transcode() );

                  // type=struct
                  if ( wnp->getAttribute("type").equals("struct") ) {
                      class_body.append(string(wnp->getType().transcode()) + " ");
                  }

                  // type=array
                  else if ( wnp->getAttribute("type").equals("array") ) {
                      class_body.append("vector<" + string(wnp->getType().transcode()) + "> " );
                  }
              }

              class_body.append( string(wnp->get_cTag().transcode()) + ";   " + comment );
          }
      }

      return class_comment + class_top + class_body + "\n\n" + methods_pub + class_end;
   }


   void wiglafNode::traverse() const {
      vector<wiglafNode>::iterator child;

      // Recursively translate children
      // First params ...
      vector<wiglafNode> p = getChildren(isParam);
      for (child=p.begin(); child != p.end(); child++) {
         if( child->getAttribute("type").equals("newbranch") ) {
            wiglafNode childp(child->getPath("value/group").getDOMNode());
            childp.traverse();
         }
      }

      // ... now groups
      vector<wiglafNode> g = getChildren(isGroup);
      for (child=g.begin(); child != g.end(); child++) {
         if (child->getAttribute("type").equals("struct") ||
             child->getAttribute("type").equals("array") )
            child->traverse();
      }


      if (!getAttribute("type").equals("array")) {

          string thisType = getType().transcode();
          string fileNameHeader = string(thisType+".H");
          string fileNameImplementation = string(thisType+".C");
          ofstream ofstreamHeader(fileNameHeader.c_str());
          ofstream ofstreamImplmentation(fileNameImplementation.c_str());


          ofstreamImplmentation << topOfFile();
          ofstreamImplmentation << "#include \"domStuff.H\"\n";
          ofstreamImplmentation << "#include \"" << fileNameHeader << "\"\n";
          ofstreamImplmentation << "#include \"wiglafNode.H\"\n\n";

#ifdef AIX_SPECIFIC
#define   EXPLICIT_FRIEND_FUNCTION
#ifdef    EXPLICIT_FRIEND_FUNCTION
          ofstreamImplmentation << "bool isOption(const wiglafNode&);\nbool isParam(const wiglafNode&);\nbool isParamButNotDeleteBranch(const wiglafNode&);\nbool isGroup(const wiglafNode&);\nbool isGroupOrParam(const wiglafNode&);";
#endif
#endif

          ofstreamImplmentation << genMethod_constructor() << genMethod_set() << genMethod_setFromArrayElement() << genMethod_toString() << endl;

          vector<string> dependencies;
          string classDefinition = genClassDefinition(dependencies);
          ofstreamHeader << topOfFile();
          ofstreamHeader << "#ifndef _" + thisType + "_H\n";
          ofstreamHeader << "#define _" + thisType + "_H\n\n";
          ofstreamHeader << "#include \"domStuff.H\"\n";
          ofstreamHeader << "#include <vector>\n";
          ofstreamHeader << "#include <string>\n";
#ifdef USE_PYTHON
          ofstreamHeader << "#include <Python.h>\n\n";
#endif /* USE_PYTHON */
          // include all user-defined object dependencies for this object
          for (vector<string>::iterator item=dependencies.begin(); item!=dependencies.end(); item++) {
              ofstreamHeader << "#include \"" + ((string)*item)  + ".H\"\n";
          }
          ofstreamHeader << "\n" << classDefinition << "\n\n";
          ofstreamHeader << "#endif\n";


          ofstreamImplmentation.close();
          ofstreamHeader.close();
      }


   }

   string wiglafNode::topOfFile() const {
      return string("/* This file is automatically generated.  DO NOT EDIT! */\n\n");
   }


   // This method is called from XML root,
   // so need to traverse down to 1st group
   void wiglafNode::XMLtoC() const {


      // -------------------------
      DOM_Node child;
      DOM_NodeList children;

       // find "application" DOM_Node
       children = node.getChildNodes();
       for (unsigned int i = 0; i < children.getLength(); i++) {
           child = children.item(i);
           if (DEBUG) cout << " " << i << " child.getNodeName().transcode(): " << child.getNodeName().transcode() << endl;
           if ( child.getNodeName().equals("application") ) {
               if (DEBUG) cout << " Found DOM_Node: application" << endl;
               break;
           }
       }

       // find "input" DOM_Node
       children = child.getChildNodes();
       for (unsigned int i = 0; i < children.getLength(); i++) {
           child = children.item(i);
           if (DEBUG) cout << " " << i << " child.getNodeName().transcode(): " << child.getNodeName().transcode() << endl;
           if ( child.getNodeName().equals("input") ) {
               if (DEBUG) cout << " Found DOM_Node: input" << endl;
               break;
           }
       }
       //DOM_Node nodeAtInput = getInputElement(node.clone());
       // -------------------------


       // Now wiglafNode(nodeAtInput) is at the <input> element, whose children are <groups> and/or <params>
       (wiglafNode( child )).getPath("/top").traverse();


   }

