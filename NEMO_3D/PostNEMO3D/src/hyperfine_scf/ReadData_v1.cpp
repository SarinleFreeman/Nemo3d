/*****************************************************************************
The Jet Propulsion Laboratory (JPL) NanoElectronicMOdeling-3D
PostProcessing package.
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
  Seungwon Lee (Seungwon.Lee@jpl.nasa.gov)

Written by:  Seungwon Lee
*****************************************************************************/


#include "ReadData.h"

#define MAX_VARS 10

void die(char* s, ...)
{
   fprintf(stderr,"ERROR:  %s\n",s);
   fflush(stderr);
   exit(1);
}


void parseRecordFormat(char* label, int* Nint, int* Nreal, int* Ncol_var, 
                       int* sz_var, int* Nvar, char* fieldType) 
{
   const char strR[] = "real ";
   const char strI[] = "int ";
   int markR=0, markI=0;

   *Nint = *Nreal = *Nvar = 0;
   
   for (char* c=label; (*c)!='\0'; c++) {
      if (*c=='{') {
	 Ncol_var[*Nvar] = 1;
	 sz_var[*Nvar] = 0;
	 (*Nvar)++;
      }
      else if (*c==',') {
	 Ncol_var[*Nvar-1]++;
      }
      /* check for match with real */
      if (*c==' ' && markR==4) {
	 markR=0;
	 fieldType[*Nint + *Nreal] = 'R';
	 sz_var[*Nvar-1] += sizeof(double);
	 (*Nreal)++;
      }
      else if (*c==strR[markR]) {
	 markR++;
      }
      else
	 markR=0;
      /* check for match with int */
      if (*c==' ' && markI==3) {
	 markI=0;
	 fieldType[*Nint + *Nreal] = 'Z';
	 sz_var[*Nvar-1] += sizeof(int);
	 (*Nint)++;
      }
      else if (*c==strI[markI]) {
	 markI++;
      }
      else
	 markI=0;
   }

   if (*Nvar > MAX_VARS)
      die("Must specify less than %d variable arguments to 'writeSiteInfo'\n", 
          MAX_VARS);
}


FILE* readHeader(const char* filename, char* recFmt, int* Ndim, int* dim, 
                 char* dim_lbl[], int* Nuser, char* usr_comment[]) 
{
   const char strB[]="<HDR>";
   const char strE[]="<\\HDR>";
   const char strRD[]="record dimension";
   const char strUD[]="user data";
   char c;
   int markB=0, markE=0, markRD=0, markUD=0;
   
   FILE* fp = fopen(filename, "rb");
   if (!fp) die("Could not open file %s for reading\n", filename);

   *Ndim = *Nuser = 0;

   /* first get to beginning of header info */
   do {
      c=getc(fp); 

      if (c=='>' && markB==4) {
	 break;
      }
      else if (c==strB[markB]) {
	 markB++;
      }
      else
	 markB=0;
   } 
   while(c!=EOF);

   /* now look for '=' and copy what follows into recFmt */
   do { c=getc(fp); }  while (c!=EOF && c!='=');
   fscanf(fp,"%[^\n]s", recFmt);

   do {
      c=getc(fp); 
      /* match record dimension */
      if (c=='n' && markRD==15) {
	 do { c=getc(fp); }  while (c!=EOF && c!='=');
	 fscanf(fp,"%d %[^\n]s", dim + (*Ndim), dim_lbl[*Ndim]);
	 (*Ndim)++;
      }
      else if (c==strRD[markRD]) {
	 markRD++;
      }
      else
	 markRD=0;
      /* match user data */
      if (c=='a' && markUD==8) {
	 do { c=getc(fp); }  while (c!=EOF && c!='=');
	 fscanf(fp,"%[^\n]s", usr_comment[*Nuser]);
	 (*Nuser)++;
      }
      else if (c==strUD[markUD]) {
	 markUD++;
      }
      else
	 markUD=0;
      /* match end tag */
      if (c=='>' && markE==5) {
	 break;
      }
      else if (c==strE[markE]) {
	 markE++;
      }
      else
	 markE=0;

   } 
   while(c!=EOF);

   dim_lbl[*Ndim] = NULL;
   dim[*Ndim] = 0;
   usr_comment[*Nuser] = NULL;
   
   return fp;
}

void read_natoms_norbitals(const char* base, // name of files to read
	                   int* Natoms, // number of atoms
			   int* Norbitals //number of orbitals
			   )
{ 
   FILE* fp;
   char label[100], *dim_lbl[MAX_VARS], *usr_comment[MAX_VARS];
   int Ndim=0, Nuser=0, dim[MAX_VARS];
   char filename[200];
   int List = 1;
   for (int i=0; i<MAX_VARS; i++) {
      dim_lbl[i] = new char[100];
      usr_comment[i] = new char[100];
   }

   // determine file sizes based on header values found in first file 
   sprintf(filename, "%s.nd_evec_%d", base, List);
   fp = readHeader(filename, label, &Ndim, dim, dim_lbl, &Nuser, usr_comment);
   if (!fp) die("Could not find file %s\n", filename);
   fclose(fp);
   *Norbitals = int(dim[0]/2);
   *Natoms = dim[1];

   for (int i=0; i<MAX_VARS; i++) {
      delete dim_lbl[i];
      delete usr_comment[i];
   }
}
    
void read_data_for_postprocessing(const char* base, // 'name' of files to read
                                  int NeList, // size of index list
                                  const int* eList, // e- eigenindex list
                                  int NhList, // size of index list
                                  const int* hList, // h  eigenindex list
                                  complex<double>** psi_e, // wave func
                                  complex<double>** psi_h,  // wave func
                                  int** idAtom, // list of atom id's
                                  int** nbr, // neighbor index
                                  double** lattice, // atomic positions
                                  double** Eel, // electron eigenenergies
                                  double** Ehl  // hole eigenenergies
                                  ) 
{
   FILE* fp;
   char label[100], *dim_lbl[MAX_VARS], *usr_comment[MAX_VARS];
   int Nrow=1, Ndim=0, Nuser=0, dim[MAX_VARS];
   int Nint, Nreal, Ncol_var[MAX_VARS], sz_var[MAX_VARS], Nvar;
   char fieldType[21];
   char filename[200];

   for (int i=0; i<MAX_VARS; i++) {
      dim_lbl[i] = new char[100];
      usr_comment[i] = new char[100];
   }

   // determine file sizes based on header values found in first file 
   if(NeList) sprintf(filename, "%s.nd_evec_%d", base, eList[0]);
   else if(NhList) sprintf(filename, "%s.nd_evec_%d", base, hList[0]);
   fp = readHeader(filename, label, &Ndim, dim, dim_lbl, &Nuser, usr_comment);
   if (!fp) die("Could not find file %s\n", filename);
   fclose(fp);
   for (int i=0; i<Ndim; i++) Nrow *= dim[i];
   int Natoms = dim[1];
   parseRecordFormat(label, &Nint, &Nreal, Ncol_var, sz_var, &Nvar, fieldType);
   cout << "size of wavefunction=" << Nrow << " * ( " << Nreal 
        << "*sizeof(real) + " << Nint << "*sizeof(int) )" << endl;

   if (Nint) die("Can't handle integers yet");

   // temporary kludge: allocate space for output arrays here
   *psi_e = new complex<double>[Nrow*NeList];
   //*psi_h = new complex<double>[Nrow*NhList];
   *idAtom = new int[Natoms];
   *nbr = new int[4*Natoms];  // assuming zinc-blende here
   *lattice = new double[3*Natoms];
   *Eel = new double[NeList];
   *Ehl = new double[NhList];

   // read in atom species data
   sprintf(filename, "%s.nd_aType", base);
   cout << "Reading atomic species:  " << filename << endl;
   fp=readHeader(filename, label, &Ndim, dim, dim_lbl, &Nuser, usr_comment);
   if (!fp) die("Could not find file %s\n", filename);
   fread(*idAtom, sizeof(int), Natoms, fp);
   fclose(fp);

   // read in atom neighbor info
   sprintf(filename, "%s.nd_nbrIndx", base);
   cout << "Reading atomic neighborhood info:  " << filename << endl;
   fp=readHeader(filename, label, &Ndim, dim, dim_lbl, &Nuser, usr_comment);
   if (!fp) die("Could not find file %s\n", filename);
   fread(*nbr, sizeof(int), 4*Natoms, fp);
   fclose(fp);

   // read in atomic positions
   sprintf(filename, "%s.nd_rAtom", base);
   cout << "Reading atomic positions:  " << filename << endl;
   fp=readHeader(filename, label, &Ndim, dim, dim_lbl, &Nuser, usr_comment);
   if (!fp) die("Could not find file %s\n", filename);
   fread(*lattice, sizeof(double), 3*Natoms, fp);
   fclose(fp);

   // read in eigenvalues
   sprintf(filename, "%s.nd_Ek", base);
   cout << "Reading eigenvalues:  " << filename << endl;
   fp=readHeader(filename, label, &Ndim, dim, dim_lbl, &Nuser, usr_comment);
   if (!fp) die("Could not find file %s\n", filename);
   double* Eread = new double[ 4*dim[0] ];
   fread(Eread, sizeof(double), 4*dim[0], fp);
   fclose(fp);
   for (int i=0; i<NeList; i++) (*Eel)[i] = Eread[ 4*eList[i]-1 ];
   for (int i=0; i<NhList; i++) (*Ehl)[i] = Eread[ 4*hList[i]-1 ];
   delete Eread;

   // read in electron wavevector data
   for (int f=0; f<NeList; f++) {
      sprintf(filename, "%s.nd_evec_%d", base, eList[f]);
      cout << "Reading electron eigenvectors:  " << filename << endl;
      fp=readHeader(filename, label, &Ndim, dim, dim_lbl, &Nuser, usr_comment);
      if (!fp) die("Could not find file %s\n", filename);
      fread(*psi_e + f*Nrow, 2*sizeof(double), Nrow, fp);
      fclose(fp);
   }

   // read in hole wavevector data
   /*for (int f=0; f<NhList; f++) {
      sprintf(filename, "%s.nd_evec_%d", base, hList[f]);
      cout << "Reading hole eigenvectors:  " << filename << endl;
      fp=readHeader(filename, label, &Ndim, dim, dim_lbl, &Nuser, usr_comment);
      if (!fp) die("Could not find file %s\n", filename);
      fread(*psi_h + f*Nrow, 2*sizeof(double), Nrow, fp);
      fclose(fp);
   }*/

   for (int i=0; i<MAX_VARS; i++) {
      delete dim_lbl[i];
      delete usr_comment[i];
   }
}

