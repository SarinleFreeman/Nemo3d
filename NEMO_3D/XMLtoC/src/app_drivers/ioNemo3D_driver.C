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

#include "top_struct.H"
#include "Dev_struct.H"
#include "ExecParam_struct.H"
#include "Keating_struct.H"
#include "Lanczos_struct.H"
#include "Material_struct.H"
#include "Output_struct.H"
#include "PhonCalc_struct.H"
#include "RayleighRitz_struct.H"
#include "Shape_struct.H"


static char*                    gXmlFile               = 0;
static bool                     gDoNamespaces          = false;
static bool                     gDoSchema              = false;
static bool                     gSchemaFullChecking    = false;
static bool                     gDoCreate              = false;
static XMLCh*                   gEncodingName          = 0;
static XMLFormatter::UnRepFlags gUnRepFlags            = XMLFormatter::UnRep_CharRef;
static DOMParser::ValSchemes    gValScheme             = DOMParser::Val_Auto;


// ---------------------------------------------------------------------------
//  ostream << DOMString
//
//  Stream out a DOM string. Doing this requires that we first transcode
//  to char * form in the default code page for the system
// ---------------------------------------------------------------------------
ostream& operator<< (ostream& target, const DOMString& s)
{
    char *p = s.transcode();
    target << p;
    delete [] p;
    return target;
}



void usage()
{
    cout << "\nUsage:\n"
            "    DOMPrint [options] <XML file>\n\n"
            "This program invokes the DOM parser, and builds the DOM tree.\n"
            "It then traverses the DOM tree and prints the contents of the\n"
            "tree for the specified XML file.\n\n"
            "Options:\n"
            "    -e          create entity reference nodes. Default is no expansion.\n"
            "    -u=xxx      Handle unrepresentable chars [fail | rep | ref*].\n"
            "    -v=xxx      Validation scheme [always | never | auto*].\n"
            "    -n          Enable namespace processing. Default is off.\n"
            "    -s          Enable schema processing. Default is off.\n"
            "    -f          Enable full schema constraint checking. Defaults to off.\n"
            "    -x=XXX      Use a particular encoding for output. Default is\n"
            "                the same encoding as the input XML file. UTF-8 if\n"
            "                input XML file has not XML declaration.\n"
		      "    -?          Show this help.\n\n"
            "  * = Default if not provided explicitly.\n\n"
            "The parser has intrinsic support for the following encodings:\n"
            "    UTF-8, USASCII, ISO8859-1, UTF-16[BL]E, UCS-4[BL]E,\n"
            "    WINDOWS-1252, IBM1140, IBM037.\n"
          <<  endl;
}


int app_input_driver(int argC, char* argV[])
{
    int retval = 0;

    // Initialize the XML4C2 system
    try {
       XMLPlatformUtils::Initialize();
    }
    catch(const XMLException& toCatch) {
       cerr << "Error during Xerces-c Initialization.\n"
            << "  Exception message:"
            << DOMString(toCatch.getMessage()) << endl;
       return 1;
    }

    // Check command line and extract arguments.
    if (argC < 2)
    {
        usage();
        XMLPlatformUtils::Terminate();
        return 1;
    }

    // See if non validating dom parser configuration is requested.
    int parmInd;
    for (parmInd = 1; parmInd < argC; parmInd++)
    {
        // Break out on first parm not starting with a dash
        if (argV[parmInd][0] != '-')
            break;

        // Watch for special case help request
        if (!strcmp(argV[parmInd], "-?"))
        {
            usage();
            XMLPlatformUtils::Terminate();
            return 2;
        }
         else if (!strncmp(argV[parmInd], "-v=", 3)
              ||  !strncmp(argV[parmInd], "-V=", 3))
        {
            const char* const parm = &argV[parmInd][3];

            if (!strcmp(parm, "never"))
                gValScheme = DOMParser::Val_Never;
            else if (!strcmp(parm, "auto"))
                gValScheme = DOMParser::Val_Auto;
            else if (!strcmp(parm, "always"))
                gValScheme = DOMParser::Val_Always;
            else
            {
                cerr << "Unknown -v= value: " << parm << endl;
                XMLPlatformUtils::Terminate();
                return 2;
            }
        }
         else if (!strcmp(argV[parmInd], "-n")
              ||  !strcmp(argV[parmInd], "-N"))
        {
            gDoNamespaces = true;
        }
         else if (!strcmp(argV[parmInd], "-s")
              ||  !strcmp(argV[parmInd], "-S"))
        {
            gDoSchema = true;
        }
         else if (!strcmp(argV[parmInd], "-f")
              ||  !strcmp(argV[parmInd], "-F"))
        {
            gSchemaFullChecking = true;
        }
         else if (!strcmp(argV[parmInd], "-e")
              ||  !strcmp(argV[parmInd], "-E"))
        {
            gDoCreate = true;
        }
         else if (!strncmp(argV[parmInd], "-x=", 3)
              ||  !strncmp(argV[parmInd], "-X=", 3))
        {
             // Get out the encoding name
             gEncodingName = XMLString::transcode( &(argV[parmInd][3]) );
        }
         else if (!strncmp(argV[parmInd], "-u=", 3)
              ||  !strncmp(argV[parmInd], "-U=", 3))
        {
            const char* const parm = &argV[parmInd][3];

            if (!strcmp(parm, "fail"))
                gUnRepFlags = XMLFormatter::UnRep_Fail;
            else if (!strcmp(parm, "rep"))
                gUnRepFlags = XMLFormatter::UnRep_Replace;
            else if (!strcmp(parm, "ref"))
                gUnRepFlags = XMLFormatter::UnRep_CharRef;
            else
            {
                cerr << "Unknown -u= value: " << parm << endl;
                XMLPlatformUtils::Terminate();
                return 2;
            }
        }
        //  else if (!strcmp(argV[parmInd], "-NoEscape"))
        // {
        //     gDoEscapes = false;
        // }
         else
        {
            cerr << "Unknown option '" << argV[parmInd]
                 << "', ignoring it.\n" << endl;
        }
    }

    //
    //  And now we have to have only one parameter left and it must be
    //  the file name.
    //
    if (parmInd + 1 != argC)
    {
        usage();
        XMLPlatformUtils::Terminate();
        return 1;
    }
    gXmlFile = argV[parmInd];

    //
    //  Create our parser, then attach an error handler to the parser.
    //  The parser will call back to methods of the ErrorHandler if it
    //  discovers errors during the course of parsing the XML document.
    //
    DOMParser *parser = new DOMParser;
    parser->setValidationScheme(gValScheme);
    parser->setDoNamespaces(gDoNamespaces);
    parser->setDoSchema(gDoSchema);
    parser->setValidationSchemaFullChecking(gSchemaFullChecking);
    DOMTreeErrorReporter *errReporter = new DOMTreeErrorReporter();
    parser->setErrorHandler(errReporter);
    parser->setCreateEntityReferenceNodes(gDoCreate);
    parser->setToCreateXMLDeclTypeNode(true);

    //
    //  Parse the XML file, catching any XML exceptions that might propogate
    //  out of it.
    //
    bool errorsOccured = false;
    try
    {
        parser->parse(gXmlFile);
        int errorCount = parser->getErrorCount();
        if (errorCount > 0)
            errorsOccured = true;
    }
    catch (const XMLException& e)
    {
        cerr << "An error occured during parsing\n   Message: "
             << DOMString(e.getMessage()) << endl;
        errorsOccured = true;
    }
    catch (const DOM_DOMException& e)
    {
       cerr << "A DOM error occured during parsing\n   DOMException code: "
             << int(e.code) << endl;
        errorsOccured = true;
    }
    catch (...)
    {
        cerr << "An error occured during parsing\n " << endl;
        errorsOccured = true;
    }

    // If the parse was successful, output the document data from the DOM tree
    if (!errorsOccured && !errReporter->getSawErrors())
    {
        char *filePy = "funcs";

        cout << " Parsing XML and creating DOM_Node..." << endl;
        DOM_Node doc = parser->getDocument();

        cout << " ____________________________________________________ " << endl;
        cout << " creating top_struct instance ts..." << endl;

#ifdef USE_PYTHON
        // Initialize embedded Python
        Py_Initialize();
        PyObject* pFile = PyString_FromString(filePy);
        PyObject *pModule = PyImport_Import(pFile);
        if (!pModule)  cerr << "Cannot find file " << filePy << ".py" << endl;
        PyObject *pDict = PyModule_GetDict(pModule);

        top_struct ts(doc,pDict);
        cout << " --------------------- " << endl;
        cout << " ts.toString():\n\n" << ts.toString(0) << endl;

        // Finalize embedded Python
        Py_DECREF(pModule);
        Py_DECREF(pFile);
        Py_Finalize();
#else
        top_struct ts(doc);
        cout << " --------------------- " << endl;
        cout << " ts.toString():\n\n" << ts.toString(0) << endl;
#endif /* USE_PYTHON */

	/*
        cout << " --------------------- " << endl;
        cout << "ts.ExecParam.StrainStruct: "    << ts.ExecParam.StrainStruct << endl;
        cout << "ts.ExecParam.ResFind: "         << ts.ExecParam.ResFind << endl;
        cout << "ts.ExecParam.BoundaryCond: "    << ts.ExecParam.BoundaryCond << endl;
        cout << "ts.ExecParam.Hamiltonian: "     << ts.ExecParam.Hamiltonian << endl;
        cout << "ts.ExecParam.RandomSeed: "      << ts.ExecParam.RandomSeed << endl;
        cout << "ts.ExecParam.BondShift_s: "     << ts.ExecParam.BondShift_s << endl;
        cout << "ts.ExecParam.BondShift_p: "     << ts.ExecParam.BondShift_p << endl;
        cout << "ts.ExecParam.BondShift_d: "     << ts.ExecParam.BondShift_d << endl;
        cout << "ts.ExecParam.BondShift_sstar: " << ts.ExecParam.BondShift_sstar << endl;
        cout << endl;

        cout << " --------------------- " << endl;
        cout << "ts.ExecParam.keating.BoundCond: "     << ts.ExecParam.keating.BoundCond << endl;
        cout << "ts.ExecParam.keating.MinWrtLatt: "    << ts.ExecParam.keating.MinWrtLatt << endl;
        cout << "ts.ExecParam.keating.tol "            << ts.ExecParam.keating.tol << endl;
        cout << "ts.ExecParam.keating.ConvCrit: "      << ts.ExecParam.keating.ConvCrit << endl;
        cout << "ts.ExecParam.keating.StartFileRead: " << ts.ExecParam.keating.StartFileRead << endl;
        cout << endl;


        cout << " ____________________________________________________ " << endl;
        cout << " creating  Shape_struct instance s..." << endl;
        Shape_struct s(doc);
        cout << "s.xorigin: " << s.xorigin << endl;
        cout << "s.yorigin: " << s.yorigin << endl;
        cout << "s.zorigin: " << s.zorigin << endl; 
        cout << "s.dx: " << s.dx << endl;
        cout << "s.mat: " << s.mat << endl;
        cout << "s.x: " << s.x << endl;
        cout << "s.mol: " << s.mol << endl;
        cout << "s.mol: " << s.mol << endl;
        cout << "s.pol: " << s.pol << endl;
        cout << "s.AlloyTreatment: " << s.AlloyTreatment << endl;
        cout << "s.AlloyClusterRadius: " << s.AlloyClusterRadius << endl;
        cout << endl;

        cout << " ____________________________________________________ " << endl;
        cout << " creating  Dev_struct instance d..." << endl;
        Dev_struct d(doc);
        cout << "d.band_model: "    << d.band_model << endl;
        cout << "d.substrate: "     << d.substrate  << endl;
        cout << "d.a_unstrained: "  << d.a_unstrained  << endl;
        cout << "d.CrystalStruct: " << d.CrystalStruct  << endl;
        cout << "d.ElStrain: "      << d.ElStrain  << endl;
        cout << "d.a_lattice_x: "   << d.a_lattice_x << endl;
        cout << "d.a_lattice_y: "   << d.a_lattice_y << endl;
        cout << "d.a_lattice_z: "   << d.a_lattice_z << endl;
        cout << endl;

        cout << "d.ShapeList is an array of " << d.ShapeList.size() << " Shape_structs: " << endl;
        for (vector<Shape_struct>::iterator item=d.ShapeList.begin(); item != d.ShapeList.end(); item++) {
	    
            cout << "      ---------- Shape_struct ---------- " << endl;
            cout << "      xorigin: " << item->xorigin << endl;
            cout << "      yorigin: " << item->yorigin << endl;
            cout << "      zorigin: " << item->zorigin << endl; 
            cout << "      dx: "  << item->dx << endl;
            cout << "      mat: " << item->mat << endl;
            cout << "      x: "   << item->x << endl;
            cout << "      mol: " << item->mol << endl;
            cout << "      mol: " << item->mol << endl;
            cout << "      pol: " << item->pol << endl;
            cout << "      AlloyTreatment: "     << item->AlloyTreatment << endl;
            cout << "      AlloyClusterRadius: " << item->AlloyClusterRadius << endl;
            cout << endl;
        }



        cout << " --------------------- " << endl;
        cout << "ts.Dev.band_model: "        << ts.Dev.band_model     << endl;
        cout << "ts.Dev.substrate: "         << ts.Dev.substrate     << endl;
        cout << "ts.Dev.a_unstrained: "      << ts.Dev.a_unstrained     << endl;
        cout << "ts.Dev.a_lattice_x: "       << ts.Dev.a_lattice_x     << endl;
        cout << "ts.Dev.a_lattice_y: "       << ts.Dev.a_lattice_y     << endl;
        cout << "ts.Dev.a_lattice_z: "       << ts.Dev.a_lattice_z     << endl;
        cout << "ts.Dev.CrystalStruct: " << ts.Dev.CrystalStruct << endl;
        cout << endl;

        cout << "ts.Dev.ShapeList is an array of " << ts.Dev.ShapeList.size() << " Shape_structs: " << endl;
        for (vector<Shape_struct>::iterator item=ts.Dev.ShapeList.begin(); item != ts.Dev.ShapeList.end(); item++) {
	    
            cout << "      ---------- Shape_struct ---------- " << endl;
            cout << "      xorigin: " << item->xorigin << endl;
            cout << "      yorigin: " << item->yorigin << endl;
            cout << "      zorigin: " << item->zorigin << endl; 
            cout << "      dx: "  << item->dx << endl;
            cout << "      mat: " << item->mat << endl;
            cout << "      x: "   << item->x << endl;
            cout << "      mol: " << item->mol << endl;
            cout << "      mol: " << item->mol << endl;
            cout << "      pol: " << item->pol << endl;
            cout << "      AlloyTreatment: "     << item->AlloyTreatment << endl;
            cout << "      AlloyClusterRadius: " << item->AlloyClusterRadius << endl;
            cout << endl;
        }

        cout << "ts.material_list an array of " << ts.material_list.size() << " Material_structs: " << endl;
        for (vector<Material_struct>::iterator item=ts.material_list.begin(); item != ts.material_list.end(); item++) {

            cout << "      ---------- Material_struct ---------- " << endl;
            cout << "      material: " << item->material << endl;
            cout << "      xval: " << item->xval << endl;
            cout << "      yval: " << item->yval << endl;
            cout << "      tkel: " << item->tkel << endl; 
            cout << "      substrate: "  << item->substrate << endl;
            cout << "      band_model: " << item->band_model << endl;
            cout << "      Ec: " << item->Ec << endl;
            cout << "      Ev: " << item->Ev << endl;
            cout << "      epsilon: " << item->epsilon << endl;
            cout << "      strain_alpha: " << item->strain_alpha << endl;
            cout << "      strain_beta: "  << item->strain_beta << endl;
            cout << "      Ev_offset: "    << item->Ev_offset << endl;

            cout << "      mb_pv is an array of " << item->mb_pv.size() << " doubles: " << endl;
            double mb_pv;
            for (vector<double>::iterator mb_pv_item=item->mb_pv.begin(); mb_pv_item != item->mb_pv.end(); mb_pv_item++) {
                cout << "            mb_pv: " << (double)*mb_pv_item << endl;
            }

            cout << "      mb_strain is an array of " << item->mb_strain.size() << " integers: " << endl;
            for (vector<int>::iterator mb_strain_item=item->mb_strain.begin(); mb_strain_item != item->mb_strain.end(); mb_strain_item++) {
	        cout << "            mb_strain : "    << (int)*mb_strain_item << endl;
            }

            cout << endl;
        }

	*/

        cout << endl;
    }
    else
        retval = 4;

    delete [] gEncodingName;
    //
    //  Clean up the error handler. The parser does not adopt handlers
    //  since they could be many objects or one object installed for multiple
    //  handlers.
    //
    delete errReporter;

    //
    //  Delete the parser itself.  Must be done prior to calling Terminate, below.
    //
    delete parser;

    // And call the termination method
    XMLPlatformUtils::Terminate();
    // DomMemDebug().print();

    //
    //  The DOM document and its contents are reference counted, and need
    //  no explicit deletion.
    //

    return 1;
    return retval;
}
