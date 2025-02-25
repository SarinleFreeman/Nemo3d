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
$Header: /repo/nemo3d/src/postprocessing/ReadData.cpp,v 1.2 2007/02/23 18:59:38
hoonryu Exp $
*****************************************************************************/

#include "ReadData.h"
#define MAX_VARS 10

void read_natoms_norbitals(const char *base, // name of files to read
                           int *Natoms,      // number of atoms
                           int *Norbitals    // number of orbitals
) {
  FILE *fp;
  char label[100], *dim_lbl[MAX_VARS], *usr_comment[MAX_VARS];
  int Ndim = 0, Nuser = 0, dim[MAX_VARS];
  char filename[200];
  int List = 1;
  for (int i = 0; i < MAX_VARS; i++) {
    dim_lbl[i] = new char[100];
    usr_comment[i] = new char[100];
  }

  // determine file sizes based on header values found in first file
  sprintf(filename, "%s.nd_evec_%d", base, List);
  fp = readHeader(filename, label, &Ndim, dim, dim_lbl, &Nuser, usr_comment);
  fclose(fp);
  *Norbitals = int(dim[0] / 2);
  *Natoms = dim[1];

  for (int i = 0; i < MAX_VARS; i++) {
    delete[] dim_lbl[i];
    delete[] usr_comment[i];
  }
}

void read_wavefunction_only(const char *filename, int Natoms, int Norbitals,
                            complex **psi) {
  FILE *fp;
  char label[100], *dim_lbl[MAX_VARS], *usr_comment[MAX_VARS];
  int Nrow = 1, Ndim = 0, Nuser = 0, dim[MAX_VARS];
  for (int i = 0; i < MAX_VARS; i++) {
    dim_lbl[i] = new char[100];
    usr_comment[i] = new char[100];
  }
  fp = readHeader(filename, label, &Ndim, dim, dim_lbl, &Nuser, usr_comment);
  for (int i = 0; i < Ndim; i++)
    Nrow *= dim[i];
  if (Nrow != Natoms * Norbitals * 2) {
    cout_master << "Nrows from XML and Evec file " << Nrow << " "
                << Natoms * Norbitals * 2 << endl;
    die("XML file and Evec file do not have the same dimension of the wave "
        "function.");
  }
  *psi = Cvectr(Nrow);
  fread(*psi, 2 * sizeof(double), Nrow, fp);
  fclose(fp);
}

void read_data_for_postprocessing(const char *base, // 'name' of files to read
                                  int NeList,       // size of index list
                                  const int *eList, // e- eigenindex list
                                  int NhList,       // size of index list
                                  const int *hList, // h  eigenindex list
                                  complex **psi_e,  // wave func
                                  complex **psi_h,  // wave func
                                  int **idAtom,     // list of atom id's
                                  int **idShape,    // list of atom id's
                                  int **nbr,        // neighbor index
                                  double **lattice  // atomic positions
) {
  FILE *fp;
  char label[100], *dim_lbl[MAX_VARS], *usr_comment[MAX_VARS];
  int Nrow = 1, Ndim = 0, Nuser = 0, dim[MAX_VARS];
  int Nint, Nreal, Nfloat, Ncol_var[MAX_VARS], sz_var[MAX_VARS], Nvar;
  char fieldType[21];
  char filename[200];

  for (int i = 0; i < MAX_VARS; i++) {
    dim_lbl[i] = new char[100];
    usr_comment[i] = new char[100];
  }

  // determine file sizes based on header values found in first file
  if (NeList)
    sprintf(filename, "%s.nd_evec_%d", base, eList[0]);
  else if (NhList)
    sprintf(filename, "%s.nd_evec_%d", base, hList[0]);
  fp = readHeader(filename, label, &Ndim, dim, dim_lbl, &Nuser, usr_comment);
  fclose(fp);
  for (int i = 0; i < Ndim; i++)
    Nrow *= dim[i];
  int Natoms = dim[1];
  parseRecordFormat(label, &Nint, &Nreal, &Nfloat, Ncol_var, sz_var, &Nvar,
                    fieldType);
  cout_master << "size of wavefunction=" << Nrow << " * ( " << Nreal
              << "*sizeof(real) + " << Nint << "*sizeof(int) )" << endl;

  if (Nint)
    die("Can't handle integers yet");

  // temporary kludge: allocate space for output arrays here
  *psi_e = Cvectr(Nrow * NeList);
  *psi_h = Cvectr(Nrow * NhList);
  *idAtom = new int[Natoms];
  *idShape = new int[Natoms];
  *nbr = new int[4 * Natoms]; // assuming zinc-blende here
  *lattice = new double[3 * Natoms];

  // read in atom species data
  sprintf(filename, "%s.nd_aType", base);
  cout_master << "Reading atomic species:  " << filename << endl;
  fp = readHeader(filename, label, &Ndim, dim, dim_lbl, &Nuser, usr_comment);
  fread(*idAtom, sizeof(int), Natoms, fp);
  fclose(fp);
  for (int ia = 0; ia < Natoms; ia++) {
    if ((*idAtom)[ia] == 12)
      (*idAtom)[ia] = 0; // As atom
    else if ((*idAtom)[ia] == 10)
      (*idAtom)[ia] = 1; // Ga atom
    else if ((*idAtom)[ia] == 14)
      (*idAtom)[ia] = 2; // In atom
    else if ((*idAtom)[ia] == 7)
      (*idAtom)[ia] = 3; // Si atom
    else if ((*idAtom)[ia] == 11)
      (*idAtom)[ia] = 4; // Ge atom
  }

  // read in shape rank data
  sprintf(filename, "%s.nd_sType", base);
  FILE *fp_tmp = fopen(filename, "rb");
  if (!fp_tmp) {
    sprintf(filename, "%s.nd_aType", base);
    cout_master << "Reading fake shape rank from atom_type file:  " << filename
                << endl;
  } else
    cout_master << "Reading shape rank:  " << filename << endl;
  fp = readHeader(filename, label, &Ndim, dim, dim_lbl, &Nuser, usr_comment);
  fread(*idShape, sizeof(int), Natoms, fp);
  fclose(fp);

  // read in atom neighbor info
  sprintf(filename, "%s.nd_nbrIndx", base);
  FILE *fp_tmp2 = fopen(filename, "rb");
  if (!fp_tmp2) {
    cout_master << "No neighborhood file available. Neighborhood info will be "
                   "not read!!!!!"
                << filename << endl;
  } else {
    cout_master << "Reading atomic neighborhood info:  " << filename << endl;
    fp = readHeader(filename, label, &Ndim, dim, dim_lbl, &Nuser, usr_comment);
    fread(*nbr, sizeof(int), 4 * Natoms, fp);
    fclose(fp);
  }

  // read in atomic positions
  sprintf(filename, "%s.nd_rAtom", base);
  cout_master << "Reading atomic positions:  " << filename << endl;
  fp = readHeader(filename, label, &Ndim, dim, dim_lbl, &Nuser, usr_comment);
  fread(*lattice, sizeof(double), 3 * Natoms, fp);
  fclose(fp);

  // read in electron wavevector data
  cout_master << "Reading electron eigenvectors:  " << filename << endl;
  for (int f = 0; f < NeList; f++) {
    sprintf(filename, "%s.nd_evec_%d", base, eList[f]);
    fp = readHeader(filename, label, &Ndim, dim, dim_lbl, &Nuser, usr_comment);
    fread(*psi_e + f * Nrow, 2 * sizeof(double), Nrow, fp);
    fclose(fp);
  }

  cout_master << "Reading electron eigenvectors:  " << filename << endl;
  // read in hole wavevector data
  for (int f = 0; f < NhList; f++) {
    sprintf(filename, "%s.nd_evec_%d", base, hList[f]);
    fp = readHeader(filename, label, &Ndim, dim, dim_lbl, &Nuser, usr_comment);
    fread(*psi_h + f * Nrow, 2 * sizeof(double), Nrow, fp);
    fclose(fp);
  }

  for (int i = 0; i < MAX_VARS; i++) {
    delete[] dim_lbl[i];
    delete[] usr_comment[i];
  }
}
