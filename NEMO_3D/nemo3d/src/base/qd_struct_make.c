/*****************************************************************************
The Jet Propulsion Laboratory (JPL) NanoElectronicMOdeling-3D package.
Copyright (C) 2002 California Institute of Technology (Caltech)

This application is free software, which you can redistribute and/or modify
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
  Gerhard Klimeck (gekco@jpl.nasa.gov)
  Fabiano Oyafuso (fabiano@jpl.nasa.gov)

Written by:  Chris Bowen
             Gerhard Klimeck
             Fabiano Oyafuso
             Seungwon Lee
             Olga Lazarenkova
             Hook Hua

This product includes software developed by the Apache Software Foundation
(http://www.apache.org/).

*****************************************************************************
$Header: /repo/nemo3d/src/base/qd_struct_make.c,v 1.39 2007/04/23 23:30:35
musman Exp $
*****************************************************************************/

#include <fstream>

#include "cmatrix.h"
#include "cvector.h"
#include "ivector.h"
#include "mb_ham_spds_micro.h"
#include "qd_struct.h"
#include "util_stringops.h"

#ifdef MPI3d
#include "run3d_mpi.h"
#endif

#define VERSION_ID "0.19"

ostream &operator<<(ostream &target, const DOMString &s);

QD_struct::QD_struct() {
  this->inputfile = NULL;
  this->Atom_to_lm = NULL;
  this->surfaceAtoms = NULL;
  this->mstar = NULL;
  this->te = NULL;
  this->tindx = NULL;
  this->wf = NULL;
  this->opmat = NULL;
  this->opmat_pxyz = NULL;
  this->wfc = NULL;
  this->eigv = NULL;
  this->mat = NULL;
  this->Basis_Spin = NULL;
  this->BasisAnion = NULL;
  this->BasisCation = NULL;
  this->param = NULL;
  this->hd = NULL;
  this->ho = NULL;
  this->parmat = NULL;
  this->phi = NULL;
  this->Ham = NULL;
  this->Ham_neighbor = NULL;
  this->Ham_offset_x = NULL;
  this->Ham_offset_y = NULL;
  this->Ham_yc = NULL;
  this->Ham_transp = NULL;
  this->Ham_atom_nbr_2_m = NULL;
  this->Ham_atom_nbr_2_l = NULL;
  this->l_m_2_atom = NULL;
  this->strnmat = NULL;
  this->cvec_aux = NULL;
  this->rc = NULL;
  this->r0 = NULL;
  this->q0 = NULL;
  this->q_m10 = NULL;
  this->r0_start = NULL;
  this->cell_s = NULL;
  this->cell_ln = NULL;
  this->seg_s = NULL;
  this->seg_ln = NULL;
  this->nvmap = NULL;
  this->commsize = NULL;
  this->commorder = NULL;
  this->vecmap = NULL;
  this->drq = NULL;

  // sunnylee
  // this is for sparse storage system for cubic 1s (used for both sparse and
  // full modes
  this->CubicSparse1s = NULL;
  this->is_CubicSparse1sSet = NULL;
  this->CubicMaxNbr = NULL;

  // sunnylee
  // this is used for matrix-matrix multiplication
  this->mm = NULL;
}

void QD_struct::cleanup() {

  nml_memory_report("Before QD_structu::cleanup");
  this->Hzb.Deallocate();
  this->HzbComplex.Deallocate();

  if (this->Ham_offset_x) {
    rm_imatrix(&this->Ham_offset_x);
    this->Ham_offset_x = NULL;
  }
  if (this->Ham_offset_y) {
    rm_imatrix(&this->Ham_offset_y);
    this->Ham_offset_y = NULL;
  }
  if (this->Ham_yc) {
    rm_imatrix(&this->Ham_yc);
    this->Ham_yc = NULL;
  }
  if (this->Ham_neighbor) {
    rm_ivectr(&this->Ham_neighbor);
    this->Ham_neighbor = NULL;
  }
  if (this->Ham_transp) {
    rm_imatrix(&this->Ham_transp);
    this->Ham_transp = NULL;
  }
  if (this->Ham_atom_nbr_2_l) {
    rm_imatrix(&this->Ham_atom_nbr_2_l);
    this->Ham_atom_nbr_2_l = NULL;
  }
  if (this->Ham_atom_nbr_2_m) {
    rm_imatrix(&this->Ham_atom_nbr_2_m);
    this->Ham_atom_nbr_2_m = NULL;
  }

  if (this->l_m_2_atom) {
    rm_imatrix(&this->l_m_2_atom);
    this->l_m_2_atom = NULL;
  }
  if (this->vecmap) {
    rm_ivectrvectr(&this->vecmap);
    this->vecmap = NULL;
  }

  if (this->Atom_to_lm) {
    rm_imatrix(&this->Atom_to_lm);
    this->Atom_to_lm = NULL;
  }

  if (this->r0) {
    rm_cvectr(&this->r0);
    this->r0 = NULL;
  }
  if (this->cvec_aux) {
    rm_cvectr(&this->cvec_aux);
    this->cvec_aux = NULL;
  }
  if (this->rc) {
    rm_cvectr(&this->rc);
    this->rc = NULL;
  }
  if (this->q0) {
    rm_cvectr(&this->q0);
    this->q0 = NULL;
  }
  if (this->q_m10) {
    rm_cvectr(&this->q_m10);
    this->q_m10 = NULL;
  }

  if (this->ho) {
    rm_cmatrix(&this->ho);
    this->ho = NULL;
  }
  if (this->hd) {
    rm_cmatrix(&this->hd);
    this->hd = NULL;
  }

  /* mpi_segment_init */
  if (this->commorder) {
    rm_imatrix(&this->commorder);
    this->commorder = NULL;
  }
  if (this->commsize) {
    rm_imatrix(&this->commsize);
    this->commsize = NULL;
  }

  if (this->seg_s) {
    rm_ivectr(&this->seg_s);
    this->seg_s = NULL;
  }
  if (this->seg_ln) {
    rm_ivectr(&this->seg_ln);
    this->seg_ln = NULL;
  }
  if (this->cell_s) {
    rm_ivectr(&this->cell_s);
    this->cell_s = NULL;
  }
  if (this->cell_ln) {
    rm_ivectr(&this->cell_ln);
    this->cell_ln = NULL;
  }
  if (this->nvmap) {
    rm_ivectr(&this->nvmap);
    this->nvmap = NULL;
  }

  /* param_init variables */
  if (this->param) {
    rm_rvectr(&this->param);
    this->param = NULL;
  }
  if (this->strnmat) {
    rm_rmatrixrvectr(&this->strnmat);
    this->strnmat = NULL;
  }

  /* mat_param_init */
  if (this->parmat) {
    rm_r3tensor(&this->parmat, 0);
    this->parmat = NULL;
  }

  /* basis_init variables */
  if (this->Basis_Spin) {
    rm_ivectr(&this->Basis_Spin);
    this->Basis_Spin = NULL;
  }
  if (this->BasisAnion) {
    rm_ivectr(&this->BasisAnion);
    this->BasisAnion = NULL;
  }
  if (this->BasisCation) {
    rm_ivectr(&this->BasisCation);
    this->BasisCation = NULL;
  }

  this->geo.Deallocate();

  if (this->mat) {
    rm_i3tensor(&this->mat, 0);
    this->mat = NULL;
  }

  /* trivial_mat_init */
  if (this->mstar) {
    rm_rvectr(&this->mstar);
    this->mstar = NULL;
  }

  /* lanc_driver_c_par */
  if (this->eigv) {
    rm_cvectr(&this->eigv);
    this->eigv = NULL;
  }

  /*Atomic masses array*/
  if (this->Amass) {
    rm_rvectr(&this->Amass);
    this->Amass = NULL;
  }

  // sunnylee
  if (this->CubicSparse1s) {
    rm_rvectr(&this->CubicSparse1s);
    this->CubicSparse1s = NULL;
    this->is_CubicSparse1sSet = NULL;
    this->CubicMaxNbr = NULL;
  }

  // sunnylee
  if (this->mm) {
    this->mm = NULL;
  }

  nml_memory_report("After QD_structu::cleanup");
  return;
}

/*
This function sets up QD_struct.opt, an object that contains the entire
tree of input parameters.  Two arguments are taken.  The first is the
name of the XML file from which data are read.  Typically this is a fifo.
The second is the name of the original XML file used to create the
first argument.  It is used to establish a naming convention for all
the files generated by this code.
*/
void QD_struct::setTopStruct(const char *f_xml, const char *f_in) {
  //******* Added to check VERSION_ID *******
  // check if the version_ID is the good one!
  cout_master << "\nInitial check: reading version # from file " << f_in
              << endl;
  cout_master << "\nCurrently using version ID '" << VERSION_ID
              << "' (line 24 of" << f_in << "file)" << endl;
  FILE *fp;
  char first23lines[200];
  char restOF24line[7];
  char VersionIDread[4];
  int ii;
  fp = fopen(f_in, "r");
  if (fp == NULL)
    cout_master << "\n No file " << f_in << " in current directory" << endl;
  for (ii = 0; ii < 23;) {
    fgets(first23lines, 200, fp);
    ii++;
  }
  fgets(restOF24line, 24, fp);
  fgets(VersionIDread, 5, fp);
  fclose(fp);
  bool IsGoodVersion = true;
  for (ii = 0; ii < 4;) {
    if (VERSION_ID[ii] != VersionIDread[ii])
      IsGoodVersion = false;
    ii++;
  }
  if (!IsGoodVersion) {
    cout << "\n********************************************************** \n"
            "****************!!! INVALID INPUT !!!********************* \n"
            "********************************************************** \n"
            "---------------------------------------------------------- \n"
            "********************************************************** \n"
            "* ======================================================== \n"
            "* The version of the input file you are using is invalid \n"
            "* (%s) and probably outdated. You should update your  \n"
            "* input to version %s.  \n"
            "*  \n"
            "* One can update the grammar to the latest version by \n"
            "* using the script: \n"
            " \n$HOME_NEMO3D/nemo3d/misc/update-nemo3d-input.pl %s \n\n"
            "* with your input deck (%s) as an argument; then \n"
            "* it will be overwritten and updated to version %s. \n"
            "*  \n"
            "* Please see the commented lines at the top for usage \n"
            "* instructions.  After applying the script one should \n"
            "* carefully examine the newly generated input to make sure \n"
            "* there are no unintended alterations to the physics or to \n"
            "* the execution flow.\n"
            "* ======================================================== \n"
            "********************************************************** \n"
         << VersionIDread << VERSION_ID << f_in << f_in << VERSION_ID << endl;
  }
  //******* End check VERSION_ID *******
  if (!mpi_n3d_id)
    cout << "\nParsing input file:  '" << f_xml << "' ('" << f_in << "')"
         << endl;

  static string pyFuncFile = "funcs";

#ifdef USE_PYTHON
  const char *pyEnv_ptr = getenv("PYTHONPATH");
  if (!pyEnv_ptr) {
    cerr << "ERROR:  PYTHONPATH undefined on proc " << mpi_n3d_id << endl;
    MPI_Finalize();
    exit(1);
  }

  string pyEnv(pyEnv_ptr);
#endif

  // Initialize the XML4C2 system
  try {
    XMLPlatformUtils::Initialize();
  } catch (const XMLException &toCatch) {
    cerr << "Error during Xerces-c Initialization.\n"
         << "  Exception message:"
         << DOMString(toCatch.getMessage()).transcode() << endl;
    return;
  }

  DOMParser *parser = new DOMParser;
  parser->setValidationScheme(DOMParser::Val_Auto);
  parser->setDoNamespaces(false);
  parser->setDoSchema(false);
  parser->setValidationSchemaFullChecking(false);

  DOMTreeErrorReporter *errReporter = new DOMTreeErrorReporter();
  parser->setErrorHandler(errReporter);
  parser->setCreateEntityReferenceNodes(false);
  parser->setToCreateXMLDeclTypeNode(true);

  //  Parse the XML file and catch any XML exceptions
  bool errorsOccured = false;
  try {
    parser->parse(f_xml);
    if (parser->getErrorCount() > 0)
      errorsOccured = true;
  } catch (const XMLException &e) {
    cerr << "An error occured during parsing\n   Message: "
         << DOMString(e.getMessage()).transcode() << endl;
    errorsOccured = true;
  } catch (const DOM_DOMException &e) {
    cerr << "A DOM error occured during parsing\n   DOMException code: "
         << int(e.code) << endl;
    errorsOccured = true;
  } catch (...) {
    cerr << "An error occured during parsing\n " << endl;
    errorsOccured = true;
  }
  if (!errorsOccured && !errReporter->getSawErrors()) {
    DOM_Node doc = parser->getDocument();

#ifdef USE_PYTHON
    // Initialize embedded Python
    Py_Initialize();
    PyObject *pFile = PyString_FromString(pyFuncFile.c_str());
    if (!pFile)
      cerr << "UNABLE TO LOCATE PYTHON FILE:  " << pyFuncFile << endl;
    PyObject *pModule = PyImport_Import(pFile);
    if (!pModule)
      cerr << "pModule is NULL -- verify that PYTHONPATH is set correctly"
           << endl;
    PyObject *pDict = PyModule_GetDict(pModule);
    if (!pDict)
      cerr << "pDict is NULL!!!" << endl;
    this->opt = top_struct(doc, pDict);
    // Finalize embedded Python
    Py_DECREF(pModule);
    Py_DECREF(pFile);
    Py_Finalize();
#else
    this->opt = top_struct(doc);
    /*if (this->opt.Version != VERSION_ID){
       die("\n****************************\n"\
           "****!!! INVALID INPUT !!!***\n"\
           "****************************\n"\
       "********************************************************** \n"\
       "* ======================================================== \n"\
       "* The version of the input file you are using is invalid (%s)\n"\
       "* and probably outdated. You should update your input to \n"\
       "* version %s.  \n"\
       "* One can update the grammar to the latest version by \n"\
       "* using the script: \n"\
       "* $HOME_NEMO3D/nemo3d/misc/update-nemo3d-input.pl %s \n"\
       "* with your input deck as an argument; then your .xml file \n"
       "* will be overwritten and updated to version %s. \n"
       "* Please see the commented lines at the top for usage \n"\
       "* instructions.  After applying the script one should \n"\
       "* carefully examine the newly generated input to make sure \n"\
       "* there are no unintended alterations to the physics or to \n"\
       "* the execution flow.\n"\
       "* ======================================================== \n"\
       "********************************************************** \n",
       this->opt.Version,VERSION_ID, VERSION_ID);
    }*/

#endif /* NOPYTHON */
  }
  // end of the errrorr loacation

  //  Clean up the error handler. The parser does not adopt handlers
  //  since they could be many objects or one object installed for multiple
  //  handlers.
  //
  delete errReporter;

  //  first delete the parser ; then call the termination method
  delete parser;
  XMLPlatformUtils::Terminate();

  if (errorsOccured) {
    if (!mpi_n3d_id)
      cout << "\nQuitting..." << endl;
    die("");
  } else if (!mpi_n3d_id) {
    cout << "   * parsing successful" << endl;
    if (this->opt.ExecParam.Output.RunFile) {
      ofstream fildb;
      char *fil_debug = n3d_strdup_n(f_in);
      n3d_FileTypeSet(&fil_debug, "nd_xmlparse", true);
      cout << "   * writing to file '" << fil_debug << "'" << endl;

      fildb.open(fil_debug);
#ifdef USE_PYTHON
      fildb << "PYTHONPATH='" << pyEnv_ptr << "'" << endl;
#endif
      fildb << this->opt.toString(5) << endl;
      fildb.close();
      str_free(fil_debug);
    }
  }

  /*
        {
        die("!!! INVALID INPUT !!!\n"\
  "The version of the input file you are using (%s) is invalid and probably
  outdated."\
  "You should update your input to version %s.  One can update the grammar to
  the latest "\
  "version by using the script $HOME_NEMO3D/nemo3d/misc/update-nemo3d-input.pl.
  "\
  "Please see the commented lines at the top for usage instructions.  After
  applying the "\
  "script one should carefully examine the newly generated input to make sure
  there are "\ "no unintended alterations to the physics or to the execution
  flow.\n", this->opt.Version.c_str(),VERSION_ID); cout_master  <<
  "====================================================== \n"
        << "****************************************************** \n"
        << "The version of the input file you are using
  (this->opt.Version.c_str()) \n"
        << "is invalid and probably outdated. You should update \n"
        << "your input to version VERSION_ID. One can update the grammar \n"
        << "to the latest version by using the script: \n"
        << "$HOME_NEMO3D/nemo3d/misc/update-nemo3d-input.pl .xml\n"
        << "with your input deck as an argument, then your old\n"
        << "xml file will be overwritten and updated. Please, \n"
        << "see the commented lines at the top for usage instructions. \n"
        << " After applying the script one should carefully examine"
        << " the newly generated input to make sure there are no\n"
        << " unintended alterations to the physics or to the execution flow.\n"
        << "****************************************************** \n"
        << "=================================================\n";
        }
  */
}

// do preliminary checks on consistency of inputs
void QD_struct::checkInput() {
  const char bullet[7] = "   * ";
  cout_master << "\nPerforming consistency check on input file...\n";

  // Geometric info
  if (this->opt.Dev.band_model == Dev_struct::Bands_1_s_nospin)
    cout_master
        << "WARNING:  check MPI communication for single orbital model.  Also "
           "look at d->n_atom_tot, d->geo.AtomsPerCellMax()";

  // shape list
  if (this->opt.Dev.ShapeList.size() == 0)
    masterPrint("%sNo shapes were found!  exiting...", bullet), die("");

  // material list
  for (unsigned int i = 0; i < this->opt.material_list.size(); i++) {
    const Material_struct &m = this->opt.material_list[i];

    if (m.mb_strain.size() != sp3d5s_nstrain) {
      die("size of (%dth) material %s strain parameter list (%d) != %d", i,
          m.material.c_str(), m.mb_strain.size(), sp3d5s_nstrain);
    }

    if (m.band_model == "Bands_20_sp3d5ss_spin" &&
        m.mb_pv.size() != sp3d5s_nparam)
      die("Incorrect length of TB parameter array %d != %d\n", m.mb_pv.size(),
          sp3d5s_nparam);
    else if (m.band_model == "Bands_10_sp3ss_spin" &&
             m.mb_pv.size() != sp3s_nparam)
      die("Incorrect length of TB parameter array %d != %d\n", m.mb_pv.size(),
          sp3s_nparam);
    else if (m.band_model == "Bands_10_sp3d5ss_nospin" &&
             m.mb_pv.size() != sp3d5s_nparam)
      die("Incorrect length of TB parameter array %d != %d\n", m.mb_pv.size(),
          sp3d5s_nparam);
    else if (m.band_model == "Bands_1_s_nospin" && m.mb_pv.size() != 5)
      die("Incorrect length of TB parameter array %d != %d\n", m.mb_pv.size(),
          5);

    /* only the sophisticated sp3d5s* model does automatically contain
       the unstrained lattice constant in its parameter array.  The single
       band models and the sp3s* model do not contain that parameter directly
       and we must not check for consistency there.  */
    if (m.band_model == "Bands_20_sp3d5ss_spin") {
      if (m.mb_pv[pV_unstr_latt] != m.unstrnd_cubic_cell_length) {
        die("p[pV_unstr_latt]=%g and m.unstrnd_cubic_cell_length=%g differ.",
            m.mb_pv[pV_unstr_latt], m.unstrnd_cubic_cell_length);
      }
    }
    if (m.band_model == "Bands_10_sp3d5ss_nospin") {
      if (m.mb_pv[pV_unstr_latt] != m.unstrnd_cubic_cell_length) {
        die("p[pV_unstr_latt]=%g and m.unstrnd_cubic_cell_length=%g differ.",
            m.mb_pv[pV_unstr_latt], m.unstrnd_cubic_cell_length);
      }
    }
  }

  // execution parameters
  if (this->opt.ExecParam.Phys.MagneticFieldOn)
    masterPrint("%sMagnetic field is on  B=(%f,%f,%f) T", bullet,
                this->opt.ExecParam.Phys.Bx, this->opt.ExecParam.Phys.By,
                this->opt.ExecParam.Phys.Bz);

  // outputs
}
